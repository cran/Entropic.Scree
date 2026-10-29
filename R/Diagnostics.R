.validate_ranks <- function(k, e, m) {
  .scalar_number(k, "K_roots", 1, m, integer = TRUE)
  .scalar_number(e, "K_extended", 1, m, integer = TRUE)
  if (e < k) stop("K_extended must be at least K_roots. Choose an inclusive extended endpoint at or above the primary rank.", call. = FALSE)
  invisible(TRUE)
}

.positive_axis_diagnostics <- function(raw, k = NULL) {
  if (is.null(raw)) return(list(status = "unavailable_in_saved_result",
    tolerance = NA_real_, positive_axes = integer(), selected_axes_positive = NA,
    working_floor_only = NA, floor_weight_fraction = NA_real_))
  tolerance <- length(raw) * .Machine$double.eps * max(1, max(abs(raw)))
  positive <- which(raw > tolerance)
  working <- pmax(raw, 1e-9)
  list(status = if (length(positive)) "positive_geometry_available" else "degenerate",
    tolerance = tolerance, positive_axes = positive,
    selected_axes_positive = if (is.null(k)) NA else all(seq_len(k) %in% positive),
    working_floor_only = all(raw <= 1e-9),
    floor_weight_fraction = sum(working - pmax(raw, 0))/sum(working))
}

.check_selected_geometry <- function(raw, k) {
  info <- .positive_axis_diagnostics(raw, k)
  if (identical(info$status, "degenerate"))
    stop("Degenerate centered spectrum: no positive eigenvalue exceeds the documented numerical tolerance. Flooring cannot establish resolved informational axes.", call. = FALSE)
  if (identical(info$selected_axes_positive, FALSE))
    warning("The selected loading table includes numerically null or negative raw directions. It is an algebraic coefficient export, not a wholly positive-spectrum coordinate map; see positive_axis_diagnostics.", call. = FALSE)
  if (isTRUE(info$working_floor_only))
    warning("Every working eigenvalue equals the numerical floor. Volume and rank candidates do not establish resolved informational structure.", call. = FALSE)
  info
}

# Derive auxiliary macro-gap endpoint metadata from the saved working
# spectrum using the Appendix A rule.
.macro_endpoint <- function(eta, mean_raw) {
  b <- sum(eta > mean_raw)
  active <- eta[eta > 1e-8]
  bulk <- NA_integer_
  if (length(active) > 10L) {
    q <- length(active)
    gaps <- abs(diff(active))
    s <- min(b + max(3, floor(q * 0.10)), q - 5L)
    tail <- gaps[s:(q - 1L)]
    quartiles <- stats::quantile(tail, c(.25, .75), names = FALSE)
    typical <- tail[tail <= quartiles[2L] + 3 * (quartiles[2L] - quartiles[1L])]
    baseline <- if (length(typical)) max(typical) else stats::median(tail)
    hits <- which(gaps[seq_len(s - 1L)] > max(1e-6, 1.5 * baseline))
    if (length(hits)) bulk <- max(hits) + 1L
  }
  list(macro_boundary_index = bulk, mean_raw_eigenvalue = mean_raw,
       mean_raw_eigenvalue_count = b,
       macro_endpoint = if (is.na(bulk)) b else bulk - 1L,
       macro_endpoint_source = if (is.na(bulk)) "mean_raw_eigenvalue_count" else "detected_macro_gap")
}

.auxiliary_gravity <- function(eta, effective_rank, k, e, mean_raw, macro_index) {
  b <- sum(eta > mean_raw)
  endpoint <- if (is.na(macro_index)) b else macro_index - 1L
  bulk_cutoff <- max(e, endpoint)
  mean_cutoff <- max(e, b)
  proportions <- eta[seq_len(k)] / sum(eta[seq_len(k)])
  allocate <- function(cutoff) proportions *
    (length(eta) - effective_rank + effective_rank * sum(eta[seq_len(cutoff)])/sum(eta))
  list(bulk = allocate(bulk_cutoff), mean_reference = allocate(mean_cutoff),
    diagnostics = list(macro_boundary_index = macro_index,
      macro_endpoint = endpoint,
      macro_endpoint_source = if (is.na(macro_index)) "mean_raw_eigenvalue_count" else "detected_macro_gap",
      mean_raw_eigenvalue = mean_raw, mean_raw_eigenvalue_count = b,
      bulk_allocation_cutoff = bulk_cutoff, mean_reference_allocation_cutoff = mean_cutoff))
}
