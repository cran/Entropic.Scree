library(Entropic.Scree)

grDevices::pdf(file = NULL)
set.seed(20261028)
dat <- data.table::as.data.table(matrix(rnorm(240 * 12), nrow = 240))
names(dat) <- paste0("Measure_", seq_len(ncol(dat)))
dat[, Measure_12 := sample(c("absent", "present", "unknown"), .N, TRUE)]

fit <- function(...) suppressMessages(Entropic.Scree(
  dat, interactive_mode = FALSE, num_cores = 1,
  check_collinearity = FALSE, ...))
res <- fit()
m <- length(res$retained_features)
stopifnot(
  !any(c("SCDR", "FSIG_extended_bulk", "FSIG_extended_kaiser", "top_of_bulk") %in% names(res)),
  identical(names(res$bipolar_modules), paste0("Axis_", seq_len(res$K_roots))),
  is.matrix(res$measurement_loadings),
  identical(dim(res$measurement_loadings), c(m, as.integer(res$K_roots))),
  identical(rownames(res$measurement_loadings), res$retained_features),
  identical(colnames(res$measurement_loadings), paste0("Axis_", seq_len(res$K_roots))),
  identical(res$measurement_loadings, res$eigenvectors[, seq_len(res$K_roots), drop = FALSE]),
  identical(dim(res$eigenvectors), c(m, m)),
  length(res$raw_eigenvalues) == m,
  isTRUE(all.equal(res$eigenvalues, pmax(1e-9, res$raw_eigenvalues))),
  isTRUE(all.equal(res$centered_trace, sum(diag(res$similarity_matrix)))),
  isTRUE(all.equal(res$NSMR_percent, 100 * res$NSMR)),
  isTRUE(all.equal(res$NSMR,
    sum(abs(res$raw_eigenvalues[res$raw_eigenvalues < 0])) /
    sum(res$raw_eigenvalues[res$raw_eigenvalues > 0]))),
  max(abs(res$similarity_matrix %*% res$eigenvectors -
    sweep(res$eigenvectors, 2L, res$raw_eigenvalues, "*"))) < 1e-10
)

check_poles <- function(x, limit = Inf) {
  stopifnot(length(x$bipolar_modules) == x$K_roots)
  for (k in seq_len(x$K_roots)) {
    loadings <- x$measurement_loadings[, k]
    positive <- sort(loadings[loadings > 0], decreasing = TRUE)
    negative <- sort(loadings[loadings < 0], decreasing = FALSE)
    if (is.finite(limit)) {
      positive <- head(positive, limit)
      negative <- head(negative, limit)
    }
    stopifnot(
      identical(x$bipolar_modules[[paste0("Axis_", k)]]$Positive_Anchor, positive),
      identical(x$bipolar_modules[[paste0("Axis_", k)]]$Negative_Anchor, negative)
    )
  }
}
check_poles(res)
limited <- fit(bipolar_top_n = 0.20)
stopifnot(identical(res$measurement_loadings, limited$measurement_loadings))
check_poles(limited, max(1, floor(0.20 * m)))
one_per_pole <- fit(bipolar_top_n = 1)
check_poles(one_per_pole, 1)
stopifnot(identical(res$measurement_loadings, one_per_pole$measurement_loadings))

one_axis <- Update.Entropic.Scree(res, new_K_roots = 1, interactive_mode = FALSE)
stopifnot(identical(dim(one_axis$measurement_loadings), c(m, 1L)),
          identical(one_axis$measurement_loadings, res$eigenvectors[, 1L, drop = FALSE]))
check_poles(one_axis)
more_axes <- Update.Entropic.Scree(one_axis, new_K_roots = 3,
  new_K_extended = max(3, res$K_extended), interactive_mode = FALSE)
stopifnot(identical(more_axes$measurement_loadings, res$eigenvectors[, 1:3, drop = FALSE]),
          identical(more_axes$eigenvectors, res$eigenvectors))
check_poles(more_axes)
unchanged <- Update.Entropic.Scree(res, interactive_mode = FALSE)
stopifnot(identical(unchanged$measurement_loadings, res$measurement_loadings))

no_vectors <- fit(extract_eigenvectors = FALSE, extract_bipolar_modules = FALSE)
stopifnot(is.null(no_vectors$eigenvectors), is.null(no_vectors$measurement_loadings),
          is.null(no_vectors$bipolar_modules))
old_without_vectors <- Update.Entropic.Scree(no_vectors, interactive_mode = FALSE)
stopifnot(is.null(old_without_vectors$measurement_loadings))
vectors_only <- fit(extract_bipolar_modules = FALSE)
stopifnot(is.matrix(vectors_only$measurement_loadings), is.null(vectors_only$bipolar_modules))

# Exactly zero coefficients belong in the matrix, but neither sign-specific pole.
zero_obj <- res
zero_obj$eigenvectors[1L, 1L] <- 0
zero_obj <- Update.Entropic.Scree(zero_obj, new_K_roots = 1, interactive_mode = FALSE)
stopifnot(zero_obj$measurement_loadings[1L, 1L] == 0,
          nrow(zero_obj$measurement_loadings) == m,
          !(res$retained_features[1L] %in% names(zero_obj$bipolar_modules[[1L]]$Positive_Anchor)),
          !(res$retained_features[1L] %in% names(zero_obj$bipolar_modules[[1L]]$Negative_Anchor)))

for (bad in list(0, -1, NA_real_, Inf, 1.5, c(0.2, 0.3), "all")) {
  stopifnot(inherits(try(fit(bipolar_top_n = bad), silent = TRUE), "try-error"))
}

console <- capture.output(suppressWarnings(Entropic.Scree(dat, num_cores = 1,
  check_collinearity = FALSE, interactive_mode = TRUE)), type = "message")
stopifnot(any(grepl("Primary Informational Rank", console, fixed = TRUE)),
          any(grepl("Extended Signal Tail boundary", console, fixed = TRUE)),
          any(grepl("Negative Spectral Mass Ratio", console, fixed = TRUE)),
          any(grepl("https://doi.org/10.5281/zenodo.22028086", console, fixed = TRUE)),
          !any(grepl("Observed Generative Rank|Synergistic Curvature|Idiosyncratic Informational Variance", console)))
grDevices::dev.off()
