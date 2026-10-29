# Preparation is separate from spectral estimation. No input data.table is mutated.
.scalar_number <- function(x, name, minimum, maximum = Inf, integer = FALSE) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) ||
      x < minimum || x > maximum || (integer && x != floor(x)))
    stop(sprintf("'%s' must be one finite %s between %s and %s.", name,
                 if (integer) "integer" else "number", minimum, maximum), call. = FALSE)
  x
}

.unused_name <- function(preferred, existing) {
  candidate <- preferred
  i <- 1L
  while (candidate %in% existing) {
    i <- i + 1L
    candidate <- paste0(preferred, "_", i)
  }
  candidate
}

.prepare_scree_data <- function(data, num_bins, bin_multiplier, low_entropy_thresh,
                               purge_constants, check_collinearity,
                               categorical_missing, numeric_missing,
                               return_analysis_data = FALSE) {
  if (nrow(data) < 2L || ncol(data) < 1L)
    stop("Provide at least two records and one measurement; at least two measurements must remain after preparation.")
  if (anyDuplicated(names(data)) || anyNA(names(data)) || any(!nzchar(names(data))))
    stop("Measurement names must be nonempty and unique.")
  .scalar_number(bin_multiplier, "bin_multiplier", .Machine$double.eps)
  .scalar_number(low_entropy_thresh, "low_entropy_thresh", 0)
  if (is.null(num_bins)) num_bins <- max(2, ceiling(bin_multiplier * nrow(data)^(1/3)))
  .scalar_number(num_bins, "num_bins", 2, 46340, integer = TRUE)
  dt <- data.table::copy(data)
  categorical_info <- list()
  numeric_info <- list()
  flags <- list()
  parents <- list()
  for (nm in names(dt)) {
    x <- dt[[nm]]
    categorical <- is.factor(x) || is.character(x) || is.logical(x)
    if (!categorical && !is.numeric(x))
      stop(sprintf("Unsupported type for '%s'. Use numeric, integer, factor, character, or logical measurements.", nm))
    unavailable <- is.na(x)
    if (is.numeric(x) && any(is.infinite(x)))
      stop(sprintf("Numeric measurement '%s' contains infinite values. Resolve them explicitly before analysis.", nm))
    if (!any(unavailable)) next
    if (categorical) {
      if (categorical_missing == "error")
        stop(sprintf("Categorical measurement '%s' contains missing values. Set categorical_missing = 'category' to represent them as an explicit state.", nm))
      label <- .unused_name("<missing>", c(as.character(x[!unavailable]), levels(x)))
      values <- as.character(x)
      values[unavailable] <- label
      if (is.factor(x)) values <- factor(values, levels = c(levels(x), label))
      data.table::set(dt, j = nm, value = values)
      categorical_info[[nm]] <- list(missing_count = sum(unavailable), label = label)
    } else {
      if (numeric_missing == "error")
        stop(sprintf("Numeric measurement '%s' contains missing values. Set numeric_missing = 'median_indicator' to fill with its observed median and add a missingness indicator.", nm))
      if (all(unavailable))
        stop(sprintf("Numeric measurement '%s' has no observed values; an observed median cannot be calculated.", nm))
      fill <- stats::median(x[!unavailable])
      values <- as.numeric(x)
      values[unavailable] <- fill
      data.table::set(dt, j = nm, value = values)
      pattern <- as.integer(unavailable)
      matching <- names(flags)[vapply(flags, identical, logical(1), pattern)]
      if (length(matching)) {
        id <- matching[[1L]]
        parents[[id]] <- c(parents[[id]], nm)
      } else {
        id <- .unused_name(paste0(nm, "__missing"), c(names(data), names(flags)))
        flags[[id]] <- pattern
        parents[[id]] <- nm
      }
      numeric_info[[nm]] <- list(missing_count = sum(unavailable),
                                median_fill = fill, indicator = id)
    }
  }
  # Original measurements alone determine parent retention. Companions cannot
  # displace their parents in constant, duplicate, QR, or entropy filtering.
  prepared_originals <- if (length(flags) || return_analysis_data) data.table::copy(dt) else NULL
  exclusions <- list(constants = character(), identical = character(),
                     collinear = character(), discretized_duplicates = character(),
                     low_entropy = character())
  message("[1/10] Checking original measurements for constants and duplicates...")
  if (purge_constants) {
    constants <- names(dt)[vapply(dt, data.table::uniqueN, integer(1)) <= 1L]
    exclusions$constants <- constants
    if (length(constants)) dt[, (constants) := NULL]
    duplicates <- duplicated(as.list(dt))
    exclusions$identical <- names(dt)[duplicates]
    if (any(duplicates)) dt <- dt[, !duplicates, with = FALSE]
  }
  message("[2/10] Checking original numeric measurements for linear dependencies...")
  if (check_collinearity && ncol(dt)) {
    num_cols <- names(dt)[vapply(dt, is.numeric, logical(1))]
    if (length(num_cols) > 1L && nrow(dt) > length(num_cols)) {
      target_rows <- min(nrow(dt), length(num_cols) + 500)
      if (as.numeric(target_rows) * length(num_cols) < 2147000000) {
        sample_idx <- sample(seq_len(nrow(dt)), target_rows)
        qr_mat <- cbind(Intercept = 1, as.matrix(dt[sample_idx, num_cols, with = FALSE]))
        qr_decomp <- qr(qr_mat, tol = 1e-7)
        if (qr_decomp$rank < ncol(qr_mat)) {
          drop_indices <- qr_decomp$pivot[(qr_decomp$rank + 1L):ncol(qr_mat)]
          combos <- setdiff(colnames(qr_mat)[drop_indices], "Intercept")
          exclusions$collinear <- combos
          if (length(combos)) dt[, (combos) := NULL]
        }
      } else message("      -> Matrix exceeds LINPACK bounds; skipping the QR check.")
    }
  }
  message("[3/10] Discretizing numeric measurements and coding categorical states...")
  if (ncol(dt)) {
    cols <- names(dt)
    dt[, (cols) := lapply(.SD, function(x) {
      if (is.numeric(x) && data.table::uniqueN(x) > num_bins)
        infotheo::discretize(x, disc = "equalfreq", nbins = num_bins)[[1L]]
      else data.table::frank(x, ties.method = "dense")
    })]
    message("[4/10] Removing duplicate encoded original measurements...")
    duplicates <- duplicated(as.list(dt))
    exclusions$discretized_duplicates <- names(dt)[duplicates]
    if (any(duplicates)) dt <- dt[, !duplicates, with = FALSE]
    cols <- names(dt)
    dt[, (cols) := lapply(.SD, function(x) as.integer(as.factor(x)))]
    message("[5/10] Filtering original measurements by marginal entropy...")
    entropy <- vapply(dt, infotheo::entropy, numeric(1))
    keep <- names(entropy)[is.finite(entropy) & entropy > 0 & entropy >= low_entropy_thresh]
    exclusions$low_entropy <- setdiff(names(dt), keep)
    dt <- dt[, keep, with = FALSE]
  }
  retained_originals <- data.table::copy(names(dt))
  raw_kept <- if (!is.null(prepared_originals))
    as.list(prepared_originals[, retained_originals, with = FALSE]) else list()
  protected <- names(flags)[vapply(parents, function(ids) any(ids %in% retained_originals), logical(1))]
  indicator_info <- list()
  # Process indicators with at least one retained numeric parent first.
  # Retain them despite low entropy or duplicate values; a constant required
  # indicator is invalid. With no retained parent, apply ordinary binary
  # constant, duplicate, and entropy filters against the retained columns.
  for (id in c(protected, setdiff(names(flags), protected))) {
    values <- flags[[id]]
    codes <- values + 1L
    entropy <- infotheo::entropy(codes)
    constant <- data.table::uniqueN(values) <= 1L
    identical_raw <- any(vapply(raw_kept, identical, logical(1), values))
    duplicate_codes <- any(vapply(as.list(dt), identical, logical(1), codes))
    reasons <- c(if (constant) "constants", if (identical_raw) "identical",
                 if (duplicate_codes) "discretized_duplicates",
                 if (entropy < low_entropy_thresh) "low_entropy")
    mandatory <- id %in% protected
    if (mandatory && constant) stop("A required missingness indicator is constant.")
    retained <- mandatory || !length(reasons)
    if (retained) {
      data.table::set(dt, j = id, value = codes)
      raw_kept[[id]] <- values
    } else exclusions[[reasons[[1L]]]] <- c(exclusions[[reasons[[1L]]]], id)
    indicator_info[[id]] <- list(parent_names = parents[[id]],
      retained_parent_names = intersect(parents[[id]], retained_originals),
      retained = retained, protected = mandatory, ordinary_filter_findings = reasons,
      entropy_nats = entropy, raw_observed_state = 0L, raw_missing_state = 1L,
      encoded_observed_state = 1L, encoded_missing_state = 2L)
  }
  final_order <- c(names(data), names(flags))
  dt <- dt[, final_order[final_order %in% names(dt)], with = FALSE]
  for (nm in names(categorical_info)) {
    categorical_info[[nm]]$retained <- nm %in% names(dt)
    categorical_info[[nm]]$encoded_missing_state <- if (nm %in% names(dt))
      dt[[nm]][which(is.na(data[[nm]]))[1L]] else NA_integer_
  }
  for (nm in names(numeric_info)) {
    numeric_info[[nm]]$parent_retained <- nm %in% retained_originals
    numeric_info[[nm]]$indicator_retained <- numeric_info[[nm]]$indicator %in% names(dt)
  }
  if (ncol(dt) < 2L) stop("Fewer than two informative measurements remain after preparation.")
  state_counts <- vapply(dt, data.table::uniqueN, integer(1))
  too_large <- state_counts > num_bins
  if (any(too_large)) stop(sprintf(
    "Encoded category capacity exceeds num_bins=%d: %s. Set num_bins to at least %d; doing so also changes the bin target for numeric measurements.",
    num_bins, paste(sprintf("%s (%d states)", names(state_counts)[too_large], state_counts[too_large]), collapse = ", "), max(state_counts)))
  info <- list(categorical_missing = categorical_missing, numeric_missing = numeric_missing,
    input_records = nrow(data), prepared_records = nrow(data), num_bins = num_bins,
    categorical_changes = categorical_info, numeric_changes = numeric_info,
    missingness_indicators = indicator_info, exclusions = exclusions,
    missingness_policy = "protected-numeric-missingness-1")
  analysis_data <- NULL
  if (return_analysis_data) {
    analysis_data <- data.table::copy(prepared_originals[, retained_originals, with = FALSE])
    for (id in names(flags)[names(flags) %in% names(dt)])
      data.table::set(analysis_data, j = id, value = flags[[id]])
    analysis_data <- analysis_data[, names(dt), with = FALSE]
  }
  list(data = dt, analysis_data = analysis_data,
       entropy = vapply(dt, infotheo::entropy, numeric(1)),
       num_bins = num_bins, information = info,
       raw_indicators = if (any(names(flags) %in% names(dt)))
         data.table::as.data.table(flags[names(flags) %in% names(dt)]) else NULL)
}

.analysis_data_export <- function(original, prepared, requested) {
  if (!requested) return(list(data = NULL, status = list(requested = FALSE,
    changed = NA, returned = FALSE, message = "Analysis-data return was not requested.")))
  shared <- intersect(names(original), names(prepared))
  equal_column <- function(nm) {
    a <- original[[nm]]; b <- prepared[[nm]]
    # Integer/double storage alone is not a change to numeric observations.
    if (is.numeric(a) && is.numeric(b)) identical(as.numeric(a), as.numeric(b))
    else identical(a, b)
  }
  changed_columns <- shared[!vapply(shared, equal_column, logical(1))]
  dropped <- setdiff(names(original), names(prepared))
  added <- setdiff(names(prepared), names(original))
  changed <- length(changed_columns) > 0L || length(dropped) > 0L || length(added) > 0L ||
    nrow(original) != nrow(prepared) || !identical(names(original), names(prepared))
  text <- if (changed) "Prepared analysis data are included because preprocessing changed columns or values." else
    "Input data were unchanged during preprocessing. No duplicate analysis table was returned; use the original data."
  if (!changed) message(text)
  list(data = if (changed) data.table::copy(prepared) else NULL,
       status = list(requested = TRUE, changed = changed, returned = changed,
         changed_columns = changed_columns, dropped_columns = dropped,
         added_columns = added, message = text))
}
