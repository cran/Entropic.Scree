library(Entropic.Scree)
grDevices::pdf(file = NULL)
set.seed(20261029)
ns <- asNamespace("Entropic.Scree")
must_error <- function(expr, pattern) {
  err <- tryCatch({force(expr); NULL}, error = conditionMessage)
  stopifnot(is.character(err), grepl(pattern, err, fixed = TRUE))
}
fit <- function(d, ...) suppressMessages(Entropic.Scree(d, interactive_mode = FALSE,
  num_cores = 1, return_processed_data = TRUE, ...))
d <- data.table::as.data.table(matrix(rnorm(300 * 14), 300))
base <- fit(d, check_collinearity = FALSE)

# public and native calls reject missing/out-of-capacity encoded states.
cat_data <- data.table::data.table(a = rep(sprintf("state%02d", 1:12), 25),
  b = sample(letters[1:3], 300, TRUE), c = sample(letters[1:4], 300, TRUE))
must_error(fit(cat_data), "category capacity")
valid <- fit(cat_data, num_bins = 12)
native <- get("fast_parallel_MI", ns)
must_error(native(matrix(c(1L, 8L, 2L, 1L), 2), 7L, 1L), "Encoded state")
must_error(native(matrix(c(1L, NA_integer_, 2L, 1L), 2), 7L, 1L), "Encoded state")

# explicit categorical opt-in, collision-free state, and no input mutation.
cm <- data.table::copy(cat_data)
cm$a <- factor(sample(c("yes", "<missing>", NA), 300, TRUE),
               levels = c("yes", "<missing>", "<missing>_2"))
saved <- data.table::copy(cm)
must_error(fit(cm), "categorical_missing")
cc <- fit(cm, categorical_missing = "category")
stopifnot(identical(cm, saved), cc$preparation_info$categorical_changes$a$label == "<missing>_3",
  !anyNA(cc$processed_data), nrow(cc$processed_data) == 300,
  is.null(cc$missingness_indicators))

# Observed median plus associated binary indicator; missingness is not inferred
# from the filled values (observed values equal to the median remain observed).
nm <- data.table::copy(d)
nm$V1[1:20] <- NA_real_
nm$V2[21] <- NaN
saved <- data.table::copy(nm)
must_error(fit(nm), "numeric_missing")
nn <- fit(nm, numeric_missing = "median_indicator", check_collinearity = FALSE)
info <- nn$preparation_info$numeric_changes$V1
stopifnot(identical(nm, saved), info$median_fill == median(nm$V1, na.rm = TRUE),
  identical(nn$missingness_indicators[[info$indicator]], as.integer(is.na(nm$V1))),
  identical(nn$processed_data[[info$indicator]], as.integer(is.na(nm$V1))),
  all(nn$processed_data$V1[1:20] == median(nm$V1, na.rm=TRUE)),
  !anyNA(nn$processed_data), nrow(nn$processed_data) == nrow(nm))
bad <- data.table::copy(nm); bad$V1 <- NA_real_
must_error(fit(bad, numeric_missing = "median_indicator"), "no observed values")
bad <- data.table::copy(d); bad$V1[1] <- Inf
must_error(fit(bad, numeric_missing = "median_indicator"), "infinite")

# Identical missingness patterns share one indicator. Filter parents first;
# indicators bypass QR. If any parent remains, retain its indicator despite
# low entropy or duplicate values; otherwise apply ordinary binary filters.
parent <- rnorm(2000); parent[1] <- NA_real_
protect <- data.table::data.table(a = parent, b = 2 * parent + 3,
  unrelated = rnorm(2000), other = rnorm(2000), a__missing = rnorm(2000))
pp <- fit(protect, numeric_missing = "median_indicator")
maps <- pp$preparation_info$missingness_indicators
stopifnot(length(maps) == 1L)
flag <- names(maps)[1]
stopifnot(flag == "a__missing_2", identical(maps[[flag]]$parent_names, c("a", "b")),
  maps[[flag]]$protected, maps[[flag]]$retained, maps[[flag]]$entropy_nats < .05,
  "low_entropy" %in% maps[[flag]]$ordinary_filter_findings,
  "b" %in% pp$preparation_info$exclusions$collinear,
  "a" %in% pp$retained_features, flag %in% pp$retained_features)
dup <- data.table::data.table(parent = c(NA, rnorm(299)),
  observed_flag = as.character(c(1L, rep(0L, 299))), other = rnorm(300), fourth = rnorm(300))
dp <- fit(dup, numeric_missing = "median_indicator", low_entropy_thresh = 0)
fi <- dp$preparation_info$missingness_indicators[[1]]
stopifnot(fi$protected, fi$retained, "discretized_duplicates" %in% fi$ordinary_filter_findings)

# Orphaned flags face normal binary filtering, including low entropy. A common
# orphaned flag may survive even though its constant parent was removed.
orphan <- data.table::data.table(a = c(NA_real_, rep(4, 1999)),
  b = c(rep(NA_real_,1000), rep(8,1000)), c = rnorm(2000), e = rnorm(2000))
oo <- fit(orphan, numeric_missing = "median_indicator")
mi <- oo$preparation_info$missingness_indicators
stopifnot(!mi$a__missing$protected, !mi$a__missing$retained,
  "low_entropy" %in% mi$a__missing$ordinary_filter_findings,
  !mi$b__missing$protected, mi$b__missing$retained)
all_orphans <- data.table::data.table(a = c(rep(NA_real_,100), rep(4,200)),
  b = rep(c(NA_real_,8,8),100))
ao <- fit(all_orphans, numeric_missing = "median_indicator")
stopifnot(ncol(ao$processed_data) == 2L, nrow(ao$processed_data) == 300L)

# primary and all auxiliary allocations use current ranks and recorded
# inclusive cutoffs, including the no-macro mean-reference fallback.
u <- Update.Entropic.Scree(base, 1, 14, interactive_mode = FALSE)
stopifnot(length(u$FSIG_final) == 1L, length(u$FSIG_auxiliary_bulk) == 1L,
  length(u$FSIG_auxiliary_mean_reference) == 1L,
  abs(sum(u$FSIG_auxiliary_bulk) - 14) < 1e-12,
  abs(sum(u$FSIG_auxiliary_mean_reference) - 14) < 1e-12,
  !any(c("SCDR", "FSIG_extended_bulk", "FSIG_extended_kaiser", "top_of_bulk") %in% names(u)),
  identical(names(u$bipolar_modules), "Axis_1"))
aux <- get(".auxiliary_gravity", ns)(c(rep(1,11),1e-9),11,2,2,11/12,NA_integer_)
stopifnot(aux$diagnostics$macro_endpoint == 11,
  aux$diagnostics$bulk_allocation_cutoff == 11,
  identical(aux$bulk, aux$mean_reference))
legacy <- base; legacy$auxiliary_diagnostics <- NULL; legacy$raw_eigenvalues <- NULL
lu <- Update.Entropic.Scree(legacy,1,14,interactive_mode=FALSE)
stopifnot(identical(lu$FSIG_auxiliary_bulk,u$FSIG_auxiliary_bulk),
  lu$positive_axis_diagnostics$status == "unavailable_in_saved_result")

# Return the documented field names without duplicating quantities.
old_names <- base
old_names$SCDR <- old_names$NSMR_percent
old_names$NSMR_percent <- NULL
old_names$NSMR <- NULL
old_names$FSIG_extended_bulk <- old_names$FSIG_auxiliary_bulk
old_names$FSIG_extended_kaiser <- old_names$FSIG_auxiliary_mean_reference
old_names$top_of_bulk <- 2L
names(old_names$bipolar_modules) <- paste0("Factor_", seq_along(old_names$bipolar_modules))
old_snapshot <- old_names
current <- Update.Entropic.Scree(old_names, 1, 14, interactive_mode=FALSE)
stopifnot(identical(old_names, old_snapshot),
  !any(c("SCDR", "FSIG_extended_bulk", "FSIG_extended_kaiser", "top_of_bulk") %in% names(current)),
  identical(current$NSMR_percent, 100 * current$NSMR),
  identical(current$NSMR_percent, old_names$SCDR),
  identical(names(current$bipolar_modules), "Axis_1"))

# reject invalid/inverted ranks consistently, including the interactive path.
for (x in list(0, 15, 2.5, NA_real_, Inf, c(1,2)))
  must_error(Update.Entropic.Scree(base, x, 14, interactive_mode=FALSE), "K_roots")
must_error(Update.Entropic.Scree(base,4,1,interactive_mode=FALSE), "K_extended must")
env <- new.env(parent=ns); env$interactive <- function() TRUE
answers <- c("4","1","1","2","yes"); ai <- 0L
env$readline <- function(prompt="") { ai <<- ai+1L; answers[[ai]] }
fun <- get("Entropic.Scree",ns); environment(fun) <- env
messages <- character()
ir <- withCallingHandlers(fun(d,check_collinearity=FALSE,num_cores=1),
  message=function(c) {messages <<- c(messages,conditionMessage(c)); invokeRestart("muffleMessage")})
stopifnot(ai == 5L, ir$K_roots == 1, ir$K_extended == 2,
  ir$rank_diagnostics$selection == "user_confirmed",
  any(grepl("Please choose valid ranks",messages,fixed=TRUE)))

# a floor-only zero geometry cannot become an asserted resolved result.
g <- rep(c("a","b"),120)
must_error(fit(data.table::data.table(a=g,b=ifelse(g=="a","b","a")),num_bins=2), "Degenerate centered spectrum")
geom <- get(".positive_axis_diagnostics",ns)(c(2,0,-.1),2)
stopifnot(identical(geom$positive_axes,1L), !geom$selected_axes_positive)
warnings <- character()
full <- withCallingHandlers(Update.Entropic.Scree(base,14,14,interactive_mode=FALSE),
  warning=function(c) { warnings <<- c(warnings,conditionMessage(c)); invokeRestart("muffleWarning")})
stopifnot(ncol(full$measurement_loadings)==14,
  !full$positive_axis_diagnostics$selected_axes_positive,
  any(grepl("algebraic coefficient export",warnings,fixed=TRUE)))

# honest selection status and exact raw-trace threshold, preserved across updates.
messages <- capture.output(suppressMessages(NULL))
messages <- capture.output(Entropic.Scree(d,interactive_mode=TRUE,
  num_cores=1,check_collinearity=FALSE),type="message")
stopifnot(any(grepl("Automatic Candidate",messages,fixed=TRUE)),
  !any(grepl("User-Confirmed",messages,fixed=TRUE)),
  any(grepl("Primary-Boundary Eigenvalue Ratio",messages,fixed=TRUE)),
  base$auxiliary_diagnostics$mean_raw_eigenvalue == base$centered_trace/14,
  base$auxiliary_diagnostics$mean_raw_eigenvalue_count == sum(base$eigenvalues>base$centered_trace/14))

# Returning analysis data means original-unit changes, never temporary MI codes.
unchanged <- fit(d,check_collinearity=FALSE)
stopifnot(is.null(unchanged$processed_data),
  unchanged$processed_data_status$requested,
  identical(unchanged$processed_data_status$changed,FALSE),
  grepl("use the original data",unchanged$processed_data_status$message,fixed=TRUE))
chars <- fit(cat_data,num_bins=12)
stopifnot(is.null(chars$processed_data), !chars$processed_data_status$changed)
drop <- data.table::copy(d); drop$constant <- 1
changed <- fit(drop,check_collinearity=FALSE)
stopifnot(changed$processed_data_status$changed,
  "constant" %in% changed$processed_data_status$dropped_columns,
  identical(changed$processed_data$V1,d$V1))
stopifnot(nn$processed_data_status$changed,
  nn$preparation_info$numeric_changes$V1$indicator %in% nn$processed_data_status$added_columns)
grDevices::dev.off()
cat("All metric safeguards, missing-data policies, protected indicators, and data-export tests passed.\n")
