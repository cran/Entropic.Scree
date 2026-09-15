utils::globalVariables(c(".SD", "Eigenvalue", "Rank", ":="))

#' Entropic Scree Dimensionality Estimation (v1.0.1)
#'
#' @description An information-theoretic diagnostic technique for estimating the
#' intrinsic dimensionality of tabular datasets. Evaluates shared probability mass
#' via a transformed mutual information metric. Aims to extract the Intrinsic
#' Generative Rank (r) and structural topology. For full methodological details,
#' see the preprint at <https://zenodo.org/records/22028087>.
#'
#' @param data A data.table containing the dataset to evaluate. Base data.frames are not supported.
#' @param low_entropy_thresh Numeric. Minimum marginal entropy threshold. Default is 0.05.
#' @param num_bins Integer. Number of bins for discretization. Default is NULL (auto-calculated).
#' @param bin_multiplier Numeric. Multiplier for the automatic bin calculation. Default is 1.0.
#' @param num_cores Integer. Number of CPU cores for OpenMP parallelization.
#' @param interactive_mode Logical. If TRUE, displays plot and prompts user to confirm ranks.
#' @param purge_constants Logical. If TRUE, automatically removes zero-variance variables.
#' @param check_collinearity Logical. If TRUE, purges perfectly collinear variables.
#' @param triple_tap_window Integer. Window size for the triple-tap heuristic. Default is 20.
#' @param fwer_alpha Numeric. Family-wise error rate alpha threshold. Default is 0.01.
#' @param extract_eigenvectors Logical. If TRUE, extracts and returns eigenvectors.
#' @param extract_bipolar_modules Logical. If TRUE, forces eigenvector extraction to build topological poles.
#' @param bipolar_top_n Numeric. Proportion of top variables to include in poles.
#' @param return_processed_data Logical. If TRUE, includes the processed data.table in the output.
#'
#' @return A list of class \code{"entropic_scree"} containing the following components:
#' \item{eigenvalues}{A numeric vector of extracted eigenvalues.}
#' \item{similarity_matrix}{The double-centered mutual information matrix.}
#' \item{retained_features}{A character vector of variables that passed entropy thresholds.}
#' \item{bin_distributions}{A table of bin counts for the discretized data.}
#' \item{R_eff}{Total Unique Probabilistic Volume scalar.}
#' \item{K_log_gap}{Integer estimate for the Observed Generative Rank.}
#' \item{K_triple_tap}{Integer estimate for the Extended Signal Tail.}
#' \item{triple_tap_multiplier}{The local dynamic t-multiplier used in Engine B.}
#' \item{K_roots}{The finalized integer representing the Intrinsic Generative Rank.}
#' \item{K_extended}{The finalized integer representing the Extended Signal Tail.}
#' \item{top_of_bulk}{Index marking the boundary of idiosyncratic variance.}
#' \item{total_signal_volume}{Total shared signal volume scalar.}
#' \item{unique_signal_volume}{Unique signal volume scalar.}
#' \item{redundant_signal_volume}{Redundant signal volume scalar.}
#' \item{idiosyncratic_volume}{Idiosyncratic informational variance scalar.}
#' \item{AIG}{Average Informational Gravity variable equivalent.}
#' \item{FSIG_final}{A numeric vector of Factor-Specific Informational Gravity.}
#' \item{structural_topology_profile}{A numeric vector of relative factor gravities.}
#' \item{FSIG_extended_bulk}{Extended FSIG metric utilizing the macro bulk bound.}
#' \item{FSIG_extended_kaiser}{Extended FSIG metric utilizing the Kaiser bound.}
#' \item{eigenvectors}{A matrix of extracted eigenvectors, if requested.}
#' \item{bipolar_modules}{A structured list of topological extraction poles, if requested.}
#' \item{processed_data}{The pre-processed data.table, if requested.}
#'
#' @examples
#' # Generate a small random dataset
#' dummy_data <- data.table::data.table(
#'   V1 = rnorm(50),
#'   V2 = rnorm(50),
#'   V3 = rnorm(50)
#' )
#'
#' # Run the estimation (interactive_mode MUST be FALSE for automated checks)
#' res <- Entropic.Scree(dummy_data, interactive_mode = FALSE, num_cores = 1)
#'
#' @useDynLib Entropic.Scree, .registration = TRUE
#' @import data.table
#' @importFrom Rcpp sourceCpp
#' @importFrom stats lm.fit na.omit qt sd median quantile setNames
#' @importFrom utils capture.output head
#' @export
Entropic.Scree <- function(data
                           , low_entropy_thresh = 0.05
                           , num_bins = NULL
                           , bin_multiplier = 1.0
                           , num_cores = parallel::detectCores() - 2
                           , interactive_mode = TRUE
                           , purge_constants = TRUE
                           , check_collinearity = TRUE
                           , triple_tap_window = 20
                           , fwer_alpha = 0.01
                           , extract_eigenvectors = FALSE
                           , extract_bipolar_modules = FALSE
                           , bipolar_top_n = 0.20
                           , return_processed_data = FALSE) {

  # Enforce data.table requirement
  if (!data.table::is.data.table(data)) {
    stop("Error: 'data' must be a data.table. Please convert your dataset using data.table::as.data.table() before running Entropic.Scree.")
  }

  # Enforce eigenvector extraction if bipolar modules are requested
  if (extract_bipolar_modules) {
    extract_eigenvectors <- TRUE
  }

  start_time <- Sys.time()
  dt <- data.table::copy(data)

  # ----------------------------------------------------------------------------
  # [0/9] INITIAL DIMENSION CHECK
  # ----------------------------------------------------------------------------
  if (ncol(dt) < 2) {
    stop("Execution Halted: The input dataset must contain at least 2 columns to calculate mutual information.")
  }

  # ----------------------------------------------------------------------------
  # [1/9] PURGE CONSTANTS & DUPLICATES
  # ----------------------------------------------------------------------------
  if (purge_constants) {
    message("[1/9] Purging constants and identical duplicates...")
    const_cols <- names(dt)[sapply(dt, function(x)
      data.table::uniqueN(x, na.rm = TRUE) <= 1)]
    if (length(const_cols) > 0) dt[, (const_cols) := NULL]

    dup_cols <- duplicated(as.list(dt))
    if (any(dup_cols)) dt <- dt[, !dup_cols, with = FALSE]
  } else {
    message("[1/10] Skipping constant and duplicate purge (user requested)...")
  }

  # ----------------------------------------------------------------------------
  # [2/9] MULTIVARIATE COLLINEARITY CHECK
  # ----------------------------------------------------------------------------
  if (check_collinearity) {
    message("[2/10] Checking for perfect multivariate linear combinations...")
    num_cols <- names(dt)[sapply(dt, is.numeric)]
    if (length(num_cols) > 1 && nrow(dt) > length(num_cols)) {

      p_cols <- length(num_cols)
      target_rows <- min(nrow(dt), p_cols + 500)

      if ((as.numeric(target_rows) * as.numeric(p_cols)) < 2147000000) {
        sample_idx <- sample(seq_len(nrow(dt)), target_rows)
        num_mat <- as.matrix(na.omit(dt[sample_idx, num_cols, with = FALSE]))

        if (nrow(num_mat) > 0) {
          qr_mat <- cbind(Intercept = 1, num_mat)
          qr_decomp <- qr(qr_mat, tol = 1e-7)

          if (qr_decomp$rank < ncol(qr_mat)) {
            drop_indices <- qr_decomp$pivot[(qr_decomp$rank + 1):ncol(qr_mat)]
            lin_combos <- setdiff(colnames(qr_mat)[drop_indices], "Intercept")
            if (length(lin_combos) > 0) {
              message(sprintf("      -> Purged %d perfectly collinear variables.", length(lin_combos)))
              dt[, (lin_combos) := NULL]
            }
          }
        }
      } else {
        message("      -> Matrix exceeds absolute LINPACK bounds even when sub-sampled. Bypassing...")
      }
    }
  } else {
    message("[2/10] Skipping collinearity check (user requested)...")
  }

  message("[3/10] Discretizing continuous data and dense-ranking categoricals...")
  if (is.null(num_bins)) {
    num_bins <- max(2, ceiling(bin_multiplier * nrow(dt)^(1/3)))
  }

  cols <- names(dt)
  dt[, (cols) := lapply(.SD, function(x) {
    if (is.numeric(x) && data.table::uniqueN(x) > num_bins) {
      infotheo::discretize(x, disc = "equalfreq", nbins = num_bins)[[1]]
    } else {
      data.table::frank(x, ties.method = "dense")
    }
  })]

  message("[4/10] Purging non-linear monotonic duplicates and locking types...")
  dup_cols_post <- duplicated(as.list(dt))
  if (any(dup_cols_post)) dt <- dt[, !dup_cols_post, with = FALSE]

  valid_cols <- names(dt)
  dt[, (valid_cols) := lapply(.SD, function(x) as.integer(as.factor(x)))]

  message("[5/10] Calculating marginal entropies and purging near-constants...")
  H_vec <- sapply(dt, infotheo::entropy)
  valid_vars <- names(H_vec)[H_vec >= low_entropy_thresh]
  if (length(H_vec) > length(valid_vars)) {
    dt <- dt[, valid_vars, with = FALSE]
    H_vec <- H_vec[valid_vars]
  }

  p <- ncol(dt)
  if (p < 2) stop("Execution Halted: Less than 2 valid variables remain.")

  message(sprintf("[6/10] Computing %d x %d Mutual Information Matrix (C++ OpenMP)...", p, p))

  if (return_processed_data) {
    processed_data_out <- data.table::copy(dt)
  } else {
    processed_data_out <- NULL
  }

  mat_data <- as.matrix(dt)
  bin_sample_sizes <- apply(mat_data, 2, tabulate)
  rm(dt)
  gc(verbose = FALSE)

  safe_target_cores <- max(1, num_cores)
  MI_mat <- fast_parallel_MI(mat_data, num_bins = num_bins, cores = safe_target_cores)
  rownames(MI_mat) <- colnames(MI_mat) <- valid_vars

  message("[7/10] Applying Joint Entropy (Jaccard) Normalization...")
  sum_H_mat <- outer(H_vec, H_vec, FUN = "+")
  joint_H_mat <- sum_H_mat - MI_mat
  joint_H_mat[joint_H_mat < 1e-9] <- 1e-9

  NMI_mat <- MI_mat / joint_H_mat
  diag(NMI_mat) <- 1.0

  message("[8/10] Applying Double-Centering...")
  m_valid <- length(valid_vars)

  # Vectorized Double-Centering
  row_means <- rowMeans(NMI_mat)
  grand_mean <- mean(row_means)

  NMI_mat_c <- NMI_mat - outer(row_means, row_means, FUN = "+") + grand_mean

  # Calculate trace of the centered matrix
  Tr_Mc <- sum(diag(NMI_mat_c))
  mean_trace <- Tr_Mc / m_valid

  message("[9/10] Extracting Entropic Latent Factors (Eigen Decomposition)...")
  if (extract_eigenvectors) {
    eigen_res <- eigen(NMI_mat_c, symmetric = TRUE)
    raw_eig_vals <- eigen_res$values
    vectors_out <- eigen_res$vectors
  } else {
    eigen_res <- eigen(NMI_mat_c, symmetric = TRUE, only.values = TRUE)
    raw_eig_vals <- eigen_res$values
    vectors_out <- NULL
  }

  # Calculate SCDR before clipping
  m_plus <- sum(raw_eig_vals[raw_eig_vals > 0])
  m_minus <- sum(abs(raw_eig_vals[raw_eig_vals < 0]))
  SCDR <- (m_minus / m_plus) * 100

  eig_vals <- pmax(raw_eig_vals, 1e-9)

  # Count true positive dimensions before the geometric zero to bound plots
  max_plot_idx <- sum(raw_eig_vals > 1e-12)

  # Constructive Spectral Mass (sum of positive clipped eigenvalues)
  m_plus <- sum(eig_vals)

  message("[10/10] Calculating Total Unique Probabilistic Volume (R_eff) and Estimating Phase Transitions...")
  sig_vals <- eig_vals[eig_vals > 0]
  if (length(sig_vals) > 0) {
    p_vals <- sig_vals / sum(sig_vals)
    H_spec <- -sum(p_vals * log(p_vals))
    R_eff <- exp(H_spec)
  } else {
    R_eff <- 1
  }

  # ==========================================================================
  # STEP 1: MACRO GAP BOUNDARY (IDIOSYNCRATIC INFORMATIONAL VARIANCE CLIFF DETECTION)
  # ==========================================================================
  n_total <- length(eig_vals)
  valid_k <- sum(eig_vals > mean_trace)

  macro_max_bulk_gap <- NA_real_
  macro_actual_gap <- NA_real_
  macro_gap_ratio <- NA_real_
  top_of_bulk_idx <- NA_integer_

  valid_search_space <- eig_vals[eig_vals > 1e-8]

  if (length(valid_search_space) > 10) {
    all_gaps_diag <- abs(diff(valid_search_space))
    n_active <- length(valid_search_space)
    noise_start_idx <- min(valid_k + max(3, floor(n_active * 0.10)), n_active - 5)
    noise_tail_idx <- noise_start_idx:(n_active - 1)

    if(length(noise_tail_idx) > 0) {
      noise_gaps <- all_gaps_diag[noise_tail_idx]
      # Calculate the baseline volatility of the noise tail, ignoring extreme secondary cliffs
      q_noise <- quantile(noise_gaps, c(0.25, 0.75), names = FALSE)
      iqr_noise <- q_noise[2] - q_noise[1]
      upper_fence <- q_noise[2] + 3.0 * iqr_noise

      # Filter out outlier gaps (the secondary cliffs) to find the true noise floor
      typical_noise_gaps <- noise_gaps[noise_gaps <= upper_fence]

      if (length(typical_noise_gaps) > 0) {
        max_noise_gap <- max(typical_noise_gaps)
      } else {
        max_noise_gap <- median(noise_gaps)
      }

      macro_multiplier <- 1.5
      gap_threshold <- max(1e-6, max_noise_gap * macro_multiplier)

      # Restrict the search space to prevent triggering on outliers INSIDE the noise tail
      valid_macro_search_idx <- 1:(noise_start_idx - 1)

      macroscopic_gap_indices <- which(all_gaps_diag[valid_macro_search_idx] > gap_threshold)
      if (length(macroscopic_gap_indices) > 0) {
        top_of_bulk_idx <- max(macroscopic_gap_indices) + 1
        macro_max_bulk_gap <- max_noise_gap
        macro_actual_gap <- all_gaps_diag[top_of_bulk_idx - 1]
        if (max_noise_gap > 1e-9) macro_gap_ratio <- macro_actual_gap / max_noise_gap
      }
    }
  }

  search_start_idx <- if (!is.na(top_of_bulk_idx)) max(2, top_of_bulk_idx - 1) else valid_k

  # ==========================================================================
  # STEP 2: ENGINE A - MAXIMUM SECONDARY EIGENVALUE RATIO (LOG-GAP)
  # ==========================================================================
  # This engine isolates the Observed Generative Rank (K_roots)
  K_log_gap <- NA_integer_
  log_gap_pct_drop <- NA_real_
  log_gap_ratio <- NA_real_
  log_gap_fallback <- FALSE
  n_valid_search <- length(valid_search_space)

  if (n_valid_search >= 3) {
    all_log_gaps <- abs(diff(log(valid_search_space)))

    if (!is.na(top_of_bulk_idx) && top_of_bulk_idx < 4) {
      # Ultra-low rank edge case
      ordered_gaps <- order(all_log_gaps, decreasing = TRUE)
      second_largest_gap_idx <- ordered_gaps[2]
      if (second_largest_gap_idx >= top_of_bulk_idx) {
        K_log_gap <- max(1, top_of_bulk_idx - 1)
        log_gap_fallback <- TRUE
      } else {
        K_log_gap <- second_largest_gap_idx
      }
    } else {
      # Standard Evaluation
      search_limit <- if (!is.na(top_of_bulk_idx)) max(2, top_of_bulk_idx - 1) else valid_k
      if (search_limit >= 2) {
        secondary_gaps <- all_log_gaps[2:search_limit]
        K_log_gap <- which.max(secondary_gaps) + 1
      } else {
        K_log_gap <- 1
        log_gap_fallback <- TRUE
      }
    }
  } else {
    K_log_gap <- max(1, valid_k)
  }

  if (K_log_gap < length(eig_vals)) {
    val_current <- eig_vals[K_log_gap]
    val_next <- eig_vals[K_log_gap + 1]
    if (val_next > 1e-12) {
      log_gap_ratio <- val_current / val_next
      log_gap_pct_drop <- (1 - (val_next / val_current)) * 100
    }
  }

  # ==========================================================================
  # STEP 3: ENGINE B - TRIPLE-TAP (DYNAMIC LINEAR WITH MACRO-STITCH)
  # ==========================================================================

  # Dynamically scale window: 10% of valid variables, bounded between 5 and user-input (default 20)
  dynamic_window <- floor(0.10 * m_valid)
  triple_tap_window <- max(5, min(as.integer(triple_tap_window), dynamic_window))

  K_triple_tap <- NA_integer_
  triple_tap_multiplier <- NA_real_
  triple_tap_actual_sigma <- NA_real_
  triple_tap_expected_val <- NA_real_
  triple_tap_actual_val <- NA_real_
  triple_tap_local_sd <- NA_real_
  prob_target <- NA_real_
  stitch_applied <- FALSE
  triple_tap_fallback <- FALSE

  if (search_start_idx >= 3) {
    # 1. Apply the Topological Stitch directly to the eigenvalues
    stitched_eig_vals <- eig_vals
    if (!is.na(top_of_bulk_idx) && !is.na(macro_actual_gap)) {
      stitch_constant <- 0.5 * macro_actual_gap
      # Add constant to all eigenvalues at or inside the noise bulk
      stitched_eig_vals[top_of_bulk_idx:length(stitched_eig_vals)] <- stitched_eig_vals[top_of_bulk_idx:length(stitched_eig_vals)] + stitch_constant
      stitch_applied <- TRUE
    }

    # 2. Dynamically scale sigma using t-dist to target a family-wise false positive rate of 1/100
    # Bonferroni correction: alpha = fwer_alpha / number of tests (search_start_idx)
    prob_target <- 1 - (fwer_alpha / search_start_idx)

    # Pre-compute design matrix for the dynamic reference window - LINEAR ONLY
    x_ref <- 0:(triple_tap_window - 1)
    X_mat <- cbind(1, x_ref)

    # Scan backwards from the top of the bounded space down to K_log_gap
    for (i in seq(search_start_idx, K_log_gap, by = -1)) {
      target_val <- stitched_eig_vals[i]

      # Reference window: next N eigenvalues to calculate linear expectation
      ref_start <- i + 1
      ref_end <- min(length(stitched_eig_vals), i + triple_tap_window)

      if (ref_start > length(stitched_eig_vals)) next # Safety check at the absolute tail

      ref_vals <- stitched_eig_vals[ref_start:ref_end]

      if (length(ref_vals) == triple_tap_window) {
        # Fit 1st degree polynomial (2 params). df = window - 2
        df_local <- triple_tap_window - 2
        dynamic_t <- qt(prob_target, df = df_local)
        current_multiplier <- max(2.0, dynamic_t)

        fit <- lm.fit(x = X_mat, y = ref_vals)

        # Intercept is the prediction at the first reference point (x = 0).
        # We predict the target eigenvalue one step backward (x = -1).
        expected_val <- fit$coefficients[1] - fit$coefficients[2]

        # Protect against linear overshoot pulling the expected value negative
        expected_val <- max(1e-12, expected_val)

        # Calculate local standard deviation of residuals (RSS / df)
        rss <- sum(fit$residuals^2)
        local_sd <- max(1e-12, sqrt(rss / df_local))

        if (target_val > (expected_val + current_multiplier * local_sd)) {
          K_triple_tap <- i
          triple_tap_multiplier <- current_multiplier
          triple_tap_expected_val <- expected_val
          triple_tap_actual_val <- target_val
          triple_tap_actual_sigma <- (target_val - expected_val) / local_sd
          triple_tap_local_sd <- local_sd
          break
        }
      } else if (length(ref_vals) >= 2) {
        # Fallback for the absolute tail where full window points are not available. df = n - 1
        df_local <- length(ref_vals) - 1
        dynamic_t <- qt(prob_target, df = df_local)
        current_multiplier <- max(2.0, dynamic_t)

        local_mean <- mean(ref_vals)
        local_sd <- max(1e-12, sd(ref_vals))

        if (target_val > (local_mean + current_multiplier * local_sd)) {
          K_triple_tap <- i
          triple_tap_multiplier <- current_multiplier
          triple_tap_expected_val <- local_mean
          triple_tap_actual_val <- target_val
          triple_tap_actual_sigma <- (target_val - local_mean) / local_sd
          triple_tap_local_sd <- local_sd
          break
        }
      }
    }
  }

  # Fallback if Triple-Tap never triggers (or if search_start_idx < 3)
  if (is.na(K_triple_tap)) {
    K_triple_tap <- K_log_gap
    triple_tap_fallback <- TRUE
  }

  # Safely capture the multiplier for the printout if it never broke but we had space
  if (is.na(triple_tap_multiplier) && search_start_idx >= 3) {
    triple_tap_multiplier <- max(2.0, qt(1 - (fwer_alpha / search_start_idx), df = triple_tap_window - 2))
  }

  # ============================================================================
  # WAVE 1: INITIAL OUTPUT & METRICS
  # ============================================================================
  pct_prob_volume <- (R_eff / m_valid) * 100
  pct_redundant_signal <- (1 - (R_eff / m_valid)) * 100
  redundant_signal_volume <- m_valid - R_eff

  n_eigen_gt_mean <- valid_k
  n_eigen_le_mean <- n_total - valid_k

  if (requireNamespace("ggplot2", quietly = TRUE)) {

    # Define dual-lines for GGplot dynamically based on convergence
    if (K_log_gap == K_triple_tap) {
      lines_geom <- ggplot2::geom_vline(xintercept = K_log_gap, color = "forestgreen", linetype = "dashed", linewidth = 1.2)
      sub_title_text <- sprintf("Converged: Observed Generative Rank & Extended Signal Tail = %d", K_log_gap)
    } else {
      lines_geom <- list(
        ggplot2::geom_vline(xintercept = K_log_gap, color = "#D55E00", linetype = "dashed", linewidth = 1.2),
        ggplot2::geom_vline(xintercept = K_triple_tap, color = "purple", linetype = "twodash", linewidth = 1.2)
      )
      sub_title_text <- sprintf("Observed Generative Rank (Orange) = %d | Extended Signal Tail (Purple) = %d", K_log_gap, K_triple_tap)
    }

    # ZOOMED VIEW
    max_k <- max(K_log_gap, K_triple_tap)
    zoom_start <- max(1, min(K_log_gap, K_triple_tap) - 5)
    zoom_end <- min(max_plot_idx, max_k + 15)
    plot_df_zoom <- data.frame(Rank = zoom_start:zoom_end, Eigenvalue = eig_vals[zoom_start:zoom_end])

    p_scree_zoom <- ggplot2::ggplot(plot_df_zoom, ggplot2::aes(x = Rank, y = Eigenvalue)) +
      ggplot2::geom_line(color = "dodgerblue", linewidth = 1) +
      ggplot2::geom_point(color = "dodgerblue", size = 2) +
      lines_geom +
      ggplot2::scale_y_continuous(trans = 'log10') +
      ggplot2::scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(min(x), max(x)))))) +
      ggplot2::labs(title = "Zoomed View", x = "Eigenvalue Index", y = "Log(Eigenvalue)") +
      ggplot2::theme_minimal(base_size = 14) +
      ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 12))

    # MACRO VIEW (Dynamically extended)
    macro_base <- max(50, max_k * 10)
    if (!is.na(top_of_bulk_idx)) {
      macro_base <- max(macro_base, top_of_bulk_idx + 25)
    }
    macro_end <- min(max_plot_idx, macro_base)

    plot_df_macro <- data.frame(Rank = 1:macro_end, Eigenvalue = eig_vals[1:macro_end])
    macro_y_max <- if(length(eig_vals) >= 2) eig_vals[2] * 1.1 else max(eig_vals)

    p_scree_macro <- ggplot2::ggplot(plot_df_macro, ggplot2::aes(x = Rank, y = Eigenvalue)) +
      ggplot2::geom_line(color = "dodgerblue", linewidth = 1) +
      ggplot2::geom_point(color = "dodgerblue", size = 2) +
      lines_geom +
      ggplot2::scale_y_continuous(trans = 'log10') +
      ggplot2::coord_cartesian(ylim = c(NA, macro_y_max)) +
      ggplot2::labs(title = "Macro View", x = "Eigenvalue Index", y = "Log(Eigenvalue)") +
      ggplot2::theme_minimal(base_size = 14) +
      ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 12))

    if (requireNamespace("patchwork", quietly = TRUE)) {
      combined_plot <- (p_scree_macro + p_scree_zoom) +
        patchwork::plot_annotation(
          title = "Entropic Scree Results",
          subtitle = sub_title_text,
          theme = ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 16, hjust = 0.5),
                                 plot.subtitle = ggplot2::element_text(size = 14, hjust = 0.5))
        )
      print(combined_plot)
    } else {
      print(p_scree_macro)
      print(p_scree_zoom)
    }
  }

  if (interactive_mode) {
    message("\n===================================================================================")
    message(" STRUCTURAL COMPOSITION")
    message("===================================================================================")
    message(sprintf(" -> %-50s : %d", "Valid Variables (m)", m_valid))
    message(sprintf(" -> %-50s : %.3f", "Centered Trace (Tr_Mc)", Tr_Mc))
    message(sprintf(" -> %-50s : %.3f%%", "Synergistic Curvature Deficit Ratio (SCDR)", SCDR))
    message("-----------------------------------------------------------------")
    message(sprintf(" -> %-50s : %.3f", "Total Unique Probabilistic Volume (R_eff)", R_eff))
    message(sprintf(" -> %-50s : %.3f%%", "%", pct_prob_volume))
    message("      (Unique Signal Volume + Structural Uncertainty")
    message("      + Independent Measurement Error + Unshared Signal Geometry)")
    message("-----------------------------------------------------------------")
    message(sprintf(" -> %-50s : %.3f", "Redundant Signal Volume (m - R_eff)", redundant_signal_volume))
    message(sprintf(" -> %-50s : %.3f%%", "%", pct_redundant_signal))
    message("-----------------------------------------------------------------")
    message(sprintf(" -> %-50s : %d", "Eigenvalues > Mean Trace", n_eigen_gt_mean))
    message(sprintf(" -> %-50s : %d", "Eigenvalues <= Mean Trace", n_eigen_le_mean))
    message("===================================================================================\n")

    message("===================================================================================")
    message(" AUTOMATED ELBOW DETECTION (DUAL-DIAGNOSTIC ENSEMBLE)")
    message("===================================================================================")
    message(" [Diagnostic: Continuous Bulk Boundary (Macro Cliff)]")
    if (!is.na(top_of_bulk_idx)) {
      message(sprintf(" -> %-43s : Index %d", "Identified Top of Idiosyncratic", top_of_bulk_idx))
      message("    Informational Variance Bulk")
      message(sprintf(" -> %-43s : %.6f", "Macro Gap Baseline (Max Bulk Gap)", macro_max_bulk_gap))
      message(sprintf(" -> %-43s : %.2fx Baseline", "Actual Macro Gap Magnitude", macro_gap_ratio))
    } else {
      message(sprintf(" -> %-43s : %s", "Identified Top of Idiosyncratic Informational Variance Bulk", "Failed (Defaulted to Kaiser)"))
    }
    message("-----------------------------------------------------------------")
    message(" [Engine A: Log-Gap]")
    message(sprintf(" -> %-43s : %d%s", "Observed Generative Rank (K_roots)", K_log_gap, ifelse(log_gap_fallback, " [Forced by Macro Boundary]", "")))
    if (!is.na(log_gap_ratio)) {
      message(sprintf(" -> %-43s : %.2fx (%.1f%% Drop)", "Log-Gap Magnitude", log_gap_ratio, log_gap_pct_drop))
    }
    message("-----------------------------------------------------------------")
    message(sprintf(" [Engine B: Triple-Tap (%d-Point Linear)]", triple_tap_window))
    message(sprintf(" -> %-43s : %d%s", "Extended Signal Tail Rank (K_extended)", K_triple_tap, ifelse(triple_tap_fallback, " [Failed to trigger; defaulted to Log-Gap]", "")))
    message(sprintf(" -> %-43s : %s", "Macro-Stitch Applied", ifelse(stitch_applied, "Yes (+0.5x Macro Gap to Bulk)", "No")))
    if (!is.na(triple_tap_expected_val)) {
      message(sprintf(" -> %-43s : %.5f", "Expected Local Eigenvalue", triple_tap_expected_val))
      message(sprintf(" -> %-43s : %.5f", "Actual Local Eigenvalue", triple_tap_actual_val))
      message(sprintf(" -> %-43s : %.6f", "Calculated Local Baseline Sigma (SD)", triple_tap_local_sd))
    }
    if (!is.na(triple_tap_multiplier)) {
      message(sprintf(" -> %-43s : %.2f (1-alpha = %.6f)", "Required Local t-Multiplier", triple_tap_multiplier, prob_target))
    }
    if (!is.na(triple_tap_actual_sigma)) {
      message(sprintf(" -> %-43s : %.2f Sigma Breakout", "Actual vs Expected", triple_tap_actual_sigma))
    }
    message("===================================================================================")
  }

  # ============================================================================
  # WAVE 2: INTERACTIVE USER OVERRIDE (OR DEFAULTING)
  # ============================================================================
  # Default Assignments for non-interactive execution
  K_roots <- K_log_gap
  K_extended <- max(K_log_gap, K_triple_tap)

  if (interactive_mode) {
    if (!interactive()) {
      # CRAN SAFETY NET: Force bypass if run in a headless test environment
      message("[-] Non-interactive session detected. Bypassing manual rank confirmation.")
    } else {
      message("\n[WARNING]: The automated extractors rely on statistical heuristics and")
      message("may not perfectly align with the true structural elbow of your specific dataset.")
      message("Please visually examine the generated scree plot.\n")

      while (TRUE) {
        # 1. Prompt for Observed Generative Rank
        ans_r <- trimws(readline(prompt = sprintf("Enter the Observed Generative Rank (K_roots) [Press Enter to keep %d]: ", K_roots)))
        if (ans_r != "") {
          parsed_r <- suppressWarnings(as.integer(ans_r))
          if (!is.na(parsed_r) && parsed_r > 0 && parsed_r <= m_valid) K_roots <- parsed_r
        }

        # 2. Prompt for Extended Signal Tail Rank
        ans_ext <- trimws(readline(prompt = sprintf("Enter the rank of the Extended Signal Tail [Press Enter to keep %d]: ", K_extended)))
        if (ans_ext != "") {
          parsed_ext <- suppressWarnings(as.integer(ans_ext))
          if (!is.na(parsed_ext) && parsed_ext > 0 && parsed_ext <= m_valid) K_extended <- parsed_ext
        }

        if (K_roots > K_extended) {
          message("[-] Warning: The Extended Signal Tail Rank should be >= the Observed Generative Rank.")
        }

        # ====================================================================
        # INTERACTIVE GRAPH PREVIEW & METRICS
        # ====================================================================
        if (requireNamespace("ggplot2", quietly = TRUE)) {

          # Setup previous engine lines (faded)
          if (K_log_gap == K_triple_tap) {
            prev_lines <- ggplot2::geom_vline(xintercept = K_log_gap, color = "gray60", linetype = "dashed", linewidth = 1)
          } else {
            prev_lines <- list(
              ggplot2::geom_vline(xintercept = K_log_gap, color = "gray60", linetype = "dashed", linewidth = 1),
              ggplot2::geom_vline(xintercept = K_triple_tap, color = "gray70", linetype = "twodash", linewidth = 1)
            )
          }

          # Draw user lines: Solid for roots, dashed for extended tail
          root_line <- ggplot2::geom_vline(xintercept = K_roots, color = "forestgreen", linetype = "solid", linewidth = 1.2)
          if (K_roots != K_extended) {
            ext_line <- ggplot2::geom_vline(xintercept = K_extended, color = "forestgreen", linetype = "dashed", linewidth = 1.2)
          } else {
            ext_line <- NULL
          }

          final_lines <- list(root_line, ext_line)

          # 1. ZOOMED VIEW
          zoom_start_upd <- max(1, min(K_log_gap, K_triple_tap, K_roots) - 5)
          zoom_end_upd <- min(max_plot_idx, max(K_log_gap, K_triple_tap, K_extended) + 15)

          plot_df_zoom_upd <- data.frame(Rank = zoom_start_upd:zoom_end_upd, Eigenvalue = eig_vals[zoom_start_upd:zoom_end_upd])

          p_scree_zoom_upd <- ggplot2::ggplot(plot_df_zoom_upd, ggplot2::aes(x = Rank, y = Eigenvalue)) +
            ggplot2::geom_line(color = "dodgerblue", linewidth = 1) +
            ggplot2::geom_point(color = "dodgerblue", size = 2) +
            prev_lines + final_lines +
            ggplot2::scale_y_continuous(trans = 'log10') +
            ggplot2::scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(min(x), max(x)))))) +
            ggplot2::labs(title = "Zoomed View", x = "Eigenvalue Index (m)", y = "Log(Eigenvalue)") +
            ggplot2::theme_minimal(base_size = 14) +
            ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 12))

          # 2. MACRO VIEW
          macro_base_upd <- max(50, K_extended * 10, K_log_gap * 10)
          if (!is.na(top_of_bulk_idx)) {
            macro_base_upd <- max(macro_base_upd, top_of_bulk_idx + 25)
          }
          macro_end_upd <- min(max_plot_idx, macro_base_upd)

          plot_df_macro_upd <- data.frame(Rank = 1:macro_end_upd, Eigenvalue = eig_vals[1:macro_end_upd])
          macro_y_max <- if(length(eig_vals) >= 2) eig_vals[2] * 1.1 else max(eig_vals)

          p_scree_macro_upd <- ggplot2::ggplot(plot_df_macro_upd, ggplot2::aes(x = Rank, y = Eigenvalue)) +
            ggplot2::geom_line(color = "dodgerblue", linewidth = 1) +
            ggplot2::geom_point(color = "dodgerblue", size = 2) +
            prev_lines + final_lines +
            ggplot2::scale_y_continuous(trans = 'log10') +
            ggplot2::coord_cartesian(ylim = c(NA, macro_y_max)) +
            ggplot2::labs(title = "Macro View", x = "Eigenvalue Index (m)", y = "Log(Eigenvalue)") +
            ggplot2::theme_minimal(base_size = 14) +
            ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 12))

          # Render Side-by-Side
          if (requireNamespace("patchwork", quietly = TRUE)) {
            combined_plot_upd <- (p_scree_macro_upd + p_scree_zoom_upd) +
              patchwork::plot_annotation(
                title = "Entropic Scree Results",
                subtitle = sprintf("User Confirmed: Observed Generative Rank (Solid) = %d | Extended Signal Tail (Dashed) = %d", K_roots, K_extended),
                theme = ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 16, hjust = 0.5),
                                       plot.subtitle = ggplot2::element_text(size = 14, hjust = 0.5))
              )
            print(combined_plot_upd)
          } else {
            print(p_scree_macro_upd)
            print(p_scree_zoom_upd)
          }
        }

        # --- PREVIEW GRAVITY CALCULATIONS ---
        # 1. Total Volume is aggregated up to K_extended
        signal_variance_prev <- sum(eig_vals[1:K_extended])
        signal_weight_prev <- signal_variance_prev / m_plus
        unique_signal_volume_prev <- R_eff * signal_weight_prev
        total_signal_volume_prev <- unique_signal_volume_prev + redundant_signal_volume

        # 2. Rebundle Gravity purely into K_roots
        AIG_prev <- total_signal_volume_prev / K_roots
        core_eigenvals_prev <- eig_vals[1:K_roots]
        p_core_prev <- core_eigenvals_prev / sum(core_eigenvals_prev)
        FSIG_prev <- p_core_prev * total_signal_volume_prev

        message("\n===================================================================================")
        message(sprintf(" (PREVIEW) ELBOW METRICS (K_roots = %d, K_extended = %d) [K_extended Volume Rebundled to K_roots]", K_roots, K_extended))
        message("===================================================================================")
        STP_prev <- FSIG_prev / FSIG_prev[1]
        message(sprintf(" -> Average Informational Gravity (AIG): %.3f Variable Equivalents", AIG_prev))
        message(" -> Factor-Specific Informational Gravity (FSIG):")
        message(paste(capture.output(print(round(FSIG_prev, 3))), collapse = "\n"))
        message(" -> Structural Topology Profile (Relative to FSIG_1):")
        message(paste(capture.output(print(round(STP_prev, 3))), collapse = "\n"))
        message("===================================================================================\n")

        ans_confirm <- trimws(readline(prompt = "Do you want to finalize these ranks? (Y/N): "))
        if (tolower(ans_confirm) %in% c("y", "yes")) {
          message("\n[+] Finalizing rank selections.")
          break
        }
      }
    }
  }

  # --- FINAL GRAVITY CALCULATIONS ---
  # 1. Capture Total Shared Signal Volume strictly up to the Extended Signal Tail
  signal_variance <- sum(eig_vals[1:K_extended])
  signal_weight <- signal_variance / m_plus
  unique_signal_volume <- R_eff * signal_weight
  total_signal_volume <- unique_signal_volume + redundant_signal_volume

  # 2. Rebundle AIG and FSIG strictly into the Observed Generative Rank (K_roots)
  AIG <- total_signal_volume / K_roots
  core_eigenvals <- eig_vals[1:K_roots]
  p_core <- core_eigenvals / sum(core_eigenvals)
  FSIG_final <- p_core * total_signal_volume

  elbow_label <- ifelse(interactive_mode, "User-Confirmed", "Automated Default")
  STP_final <- FSIG_final / FSIG_final[1]

  idiosyncratic_variance <- if (K_extended < m_valid) sum(eig_vals[(K_extended + 1):m_valid]) else 0
  idiosyncratic_weight <- idiosyncratic_variance / m_plus
  idiosyncratic_volume <- R_eff * idiosyncratic_weight

  if (interactive_mode) {
    message("\n===================================================================================")
    message(sprintf(" (FINAL) ENTROPIC SCREE METRICS (based on %s)", elbow_label))
    message("===================================================================================")
    message(" [Rank Estimates]")
    message(sprintf(" -> Observed Generative Rank (K_roots)     : %d", K_roots))
    message(sprintf(" -> Extended Signal Tail Rank (K_extended) : %d", K_extended))
    message("-----------------------------------------------------------------")
    message(" [Informational Gravity Metrics Estimates] (K_extended Volume Rebundled to K_roots)")
    message(sprintf(" -> Average Informational Gravity (AIG): %.3f Variable Equivalents", AIG))
    message(" -> Factor-Specific Informational Gravity (FSIG):")
    message(paste(capture.output(print(round(FSIG_final, 3))), collapse = "\n"))
    message(" -> Structural Topology Profile (Relative to FSIG_1):")
    message(paste(capture.output(print(round(STP_final, 3))), collapse = "\n"))
    message("===================================================================================\n")

    # ============================================================================
    # WAVE 3: FINAL TRIPARTITE STRUCTURAL COMPOSITION
    # ============================================================================
    message("===================================================================================")
    message(sprintf(" (FINAL) TRIPARTITE STRUCTURAL COMPOSITION (Rebundling Based on K_extended = %d)", K_extended))
    message("===================================================================================")
    message(sprintf(" -> %-50s : %d", "Valid Variables (m)", m_valid))
    message(sprintf(" -> %-50s : %.2f", "Centered Trace (Tr_Mc)", Tr_Mc))
    message("-----------------------------------------------------------------")
    message(sprintf(" -> %-50s : %.3f", "Total Shared Signal Volume (Unique + Redundant)", total_signal_volume))
    message(sprintf(" -> %-40s : %.3f", "   (Unique Signal Volume)", unique_signal_volume))
    message(sprintf(" -> %-40s : %.3f", "   (Redundant Signal Volume)", redundant_signal_volume))
    message("-----------------------------")
    message(sprintf(" -> %-50s : %.3f", "Idiosyncratic Informational Variance", idiosyncratic_volume))
    message("      (Structural Uncertainty + Independent Measurement Error")
    message("      + Unshared Signal Geometry)")
    message("-----------------------------------------------------------------")
    message(sprintf(" -> %-50s : %.3f%%", "% Total Shared Signal", (total_signal_volume / m_valid) * 100))
    message(sprintf(" -> %-40s : %.3f%%", "   (% Unique Signal)", (unique_signal_volume / m_valid) * 100))
    message(sprintf(" -> %-40s : %.3f%%", "   (% Redundant Signal)", (redundant_signal_volume / m_valid) * 100))
    message("-----------------------------")
    message(sprintf(" -> %-50s : %.3f%%", "% Idiosyncratic Informational Variance", (idiosyncratic_volume / m_valid) * 100))
    message("      (Structural Uncertainty + Independent Measurement Error")
    message("      + Unshared Signal Geometry)")
    message("===================================================================================\n")

    # ============================================================================
    # WAVE 4: METHODOLOGICAL REFERENCE
    # ============================================================================
    message("===================================================================================")
    message(" ENTROPIC SCREE (v1.0.1) - METHODOLOGICAL REFERENCE & LICENSE")
    message("===================================================================================")
    message(" -> Framework developed by Terrence J. Lee-St. John (Enli)")
    message(" -> Released under the Apache License 2.0 (Open Source)")
    message(" -> For full methods and metric definitions, see:")
    message("    The Entropic Scree (2026) - https://doi.org/10.5281/zenodo.22028087")
    message("===================================================================================\n")
  }

  # ============================================================================
  # EXTENDED Factor-Specific Informational Gravity (FSIG) CALCULATIONS
  # ============================================================================

  # --- EXTENDED MODEL A: MACRO BULK BOUNDARY ---
  top_bulk_safe <- if(!is.na(top_of_bulk_idx)) top_of_bulk_idx else valid_k
  extended_bulk_k <- max(K_extended, top_bulk_safe - 1)
  sig_var_bulk <- sum(eig_vals[1:extended_bulk_k])
  total_sig_vol_bulk <- (R_eff * (sig_var_bulk / m_plus)) + redundant_signal_volume
  # Rebundle gravity strictly into the Observed Generative Rank (K_roots):
  FSIG_extended_bulk <- p_core * total_sig_vol_bulk

  # --- EXTENDED MODEL B: KAISER RULE BOUNDARY ---
  extended_kaiser_k <- max(K_extended, valid_k)
  sig_var_kaiser <- sum(eig_vals[1:extended_kaiser_k])
  total_sig_vol_kaiser <- (R_eff * (sig_var_kaiser / m_plus)) + redundant_signal_volume
  # Rebundle gravity strictly into the Observed Generative Rank (K_roots):
  FSIG_extended_kaiser <- p_core * total_sig_vol_kaiser

  # ============================================================================
  # EXTRACT BIPOLAR MODULES (TOPOLOGICAL POLES)
  # ============================================================================
  bipolar_modules_out <- NULL
  if (extract_bipolar_modules && !is.null(vectors_out) && K_roots > 0) {
    bipolar_modules_out <- list()

    # Dynamically resolve top N (proportion vs absolute count)
    resolved_top_n <- if (bipolar_top_n > 0 && bipolar_top_n < 1) {
      max(1, floor(bipolar_top_n * m_valid))
    } else {
      max(1, as.integer(bipolar_top_n))
    }

    for (k in 1:K_roots) {
      # Bind variable names to the specific eigenvector
      loadings <- setNames(vectors_out[, k], valid_vars)

      # Extract the extreme positive anchor
      pos_pole <- sort(loadings[loadings > 0], decreasing = TRUE)
      pos_pole <- head(pos_pole, resolved_top_n)

      # Extract the extreme negative anchor (structural estrangement)
      neg_pole <- sort(loadings[loadings < 0], decreasing = FALSE)
      neg_pole <- head(neg_pole, resolved_top_n)

      bipolar_modules_out[[paste0("Factor_", k)]] <- list(
        Positive_Anchor = pos_pole,
        Negative_Anchor = neg_pole
      )
    }
  }

  return(list(
    eigenvalues = eig_vals,
    similarity_matrix = NMI_mat_c,
    retained_features = valid_vars,
    bin_distributions = bin_sample_sizes,
    R_eff = R_eff,
    K_log_gap = K_log_gap,
    K_triple_tap = K_triple_tap,
    triple_tap_multiplier = triple_tap_multiplier,
    K_roots = K_roots,
    K_extended = K_extended,
    top_of_bulk = top_bulk_safe,
    total_signal_volume = total_signal_volume,
    unique_signal_volume = unique_signal_volume,
    redundant_signal_volume = redundant_signal_volume,
    idiosyncratic_volume = idiosyncratic_volume,
    AIG = AIG,
    FSIG_final = FSIG_final,
    structural_topology_profile = STP_final,
    FSIG_extended_bulk = FSIG_extended_bulk,
    FSIG_extended_kaiser = FSIG_extended_kaiser,
    eigenvectors = vectors_out,
    bipolar_modules = bipolar_modules_out,
    processed_data = processed_data_out
  ))
}

# ==============================================================================
# Update.Entropic.Scree (Post-Hoc Rank Modification)
# ==============================================================================

#' Update Entropic Scree Rank Estimates
#'
#' @description A companion function to modify the rank estimates of a previously
#' evaluated Entropic Scree object. This allows users to recalculate informational
#' gravity and structural composition without re-running the computationally
#' expensive C++ mutual information engine.
#'
#' @param scree_obj A valid output list from \code{Entropic.Scree()}.
#' @param new_K_roots Integer. The new Observed Generative Rank. If NULL, retains original.
#' @param new_K_extended Integer. The new Extended Signal Tail Rank. If NULL, retains original.
#' @param bipolar_top_n Numeric. Proportion of top variables to include in newly extracted poles.
#' @param interactive_mode Logical. If TRUE, prints updated metric dashboards to the console.
#'
#' @return A mutated list of class \code{"entropic_scree"} containing the updated
#' structural topology profiles, signal volumes, and informational gravity metrics
#' based on the new rank thresholds.
#'
#' @export
Update.Entropic.Scree <- function(scree_obj, new_K_roots = NULL, new_K_extended = NULL, bipolar_top_n = 0.20, interactive_mode = TRUE) {

  # 1. Input Validation
  if (!is.list(scree_obj) || is.null(scree_obj$eigenvalues)) {
    stop("Error: 'scree_obj' must be a valid output list from Entropic.Scree().")
  }

  eig_vals <- scree_obj$eigenvalues
  m_valid <- length(scree_obj$retained_features)

  # Capture original state for the comparison printout
  orig_K_roots <- scree_obj$K_roots
  orig_K_extended <- scree_obj$K_extended
  orig_AIG <- scree_obj$AIG

  orig_total_signal <- scree_obj$total_signal_volume
  orig_unique_signal <- scree_obj$unique_signal_volume
  orig_redundant_signal <- scree_obj$redundant_signal_volume
  orig_idio_vol <- scree_obj$idiosyncratic_volume

  orig_pct_signal <- (orig_total_signal / m_valid) * 100
  orig_pct_unique <- (orig_unique_signal / m_valid) * 100
  orig_pct_redundant <- (orig_redundant_signal / m_valid) * 100
  orig_pct_idio <- (orig_idio_vol / m_valid) * 100

  # Resolve new ranks (default to the object's current ranks if not provided)
  K_roots <- if (!is.null(new_K_roots)) as.integer(new_K_roots) else scree_obj$K_roots
  K_extended <- if (!is.null(new_K_extended)) as.integer(new_K_extended) else scree_obj$K_extended

  # Safety bounds
  K_roots <- max(1, min(K_roots, m_valid))
  K_extended <- max(K_roots, min(K_extended, m_valid))

  if (new_K_roots > new_K_extended && !is.null(new_K_roots) && !is.null(new_K_extended)) {
    warning("K_roots is greater than K_extended. K_extended has been bounded to K_roots.")
  }

  # 2. Extract Static Core Physics from Object
  R_eff <- scree_obj$R_eff
  m_plus <- sum(eig_vals)
  redundant_signal_volume <- scree_obj$redundant_signal_volume
  vectors_out <- scree_obj$eigenvectors
  valid_vars <- scree_obj$retained_features
  Tr_Mc <- sum(diag(scree_obj$similarity_matrix))

  # ============================================================================
  # 3. RECALCULATE GRAVITY & STRUCTURAL COMPOSITION
  # ============================================================================

  # Capture Total Shared Signal Volume strictly up to the NEW Extended Signal Tail
  signal_variance <- sum(eig_vals[1:K_extended])
  signal_weight <- signal_variance / m_plus
  unique_signal_volume <- R_eff * signal_weight
  total_signal_volume <- unique_signal_volume + redundant_signal_volume
  new_pct_signal <- (total_signal_volume / m_valid) * 100
  new_pct_unique <- (unique_signal_volume / m_valid) * 100
  new_pct_redundant <- (redundant_signal_volume / m_valid) * 100

  # Rebundle AIG and FSIG strictly into the NEW Observed Generative Rank
  AIG <- total_signal_volume / K_roots
  core_eigenvals <- eig_vals[1:K_roots]
  p_core <- core_eigenvals / sum(core_eigenvals)
  FSIG_final <- p_core * total_signal_volume
  STP_final <- FSIG_final / FSIG_final[1]

  # Idiosyncratic Variance
  idiosyncratic_variance <- if (K_extended < m_valid) sum(eig_vals[(K_extended + 1):m_valid]) else 0
  idiosyncratic_weight <- idiosyncratic_variance / m_plus
  idiosyncratic_volume <- R_eff * idiosyncratic_weight
  new_pct_idio <- (idiosyncratic_volume / m_valid) * 100

  # ============================================================================
  # 4. RE-EXTRACT BIPOLAR MODULES
  # ============================================================================
  bipolar_modules_out <- NULL
  if (!is.null(vectors_out) && K_roots > 0) {
    bipolar_modules_out <- list()

    resolved_top_n <- if (bipolar_top_n > 0 && bipolar_top_n < 1) {
      max(1, floor(bipolar_top_n * m_valid))
    } else {
      max(1, as.integer(bipolar_top_n))
    }

    for (k in 1:K_roots) {
      loadings <- setNames(vectors_out[, k], valid_vars)

      pos_pole <- sort(loadings[loadings > 0], decreasing = TRUE)
      pos_pole <- head(pos_pole, resolved_top_n)

      neg_pole <- sort(loadings[loadings < 0], decreasing = FALSE)
      neg_pole <- head(neg_pole, resolved_top_n)

      bipolar_modules_out[[paste0("Factor_", k)]] <- list(
        Positive_Anchor = pos_pole,
        Negative_Anchor = neg_pole
      )
    }
  }

  # ============================================================================
  # 5. PRINT UPDATED DASHBOARDS
  # ============================================================================
  if (interactive_mode) {
    message("\n===================================================================================")
    message(" (UPDATED) ENTROPIC SCREE METRICS (User Override)")
    message("===================================================================================")
    message(" [Rank Estimates]")
    message(sprintf(" -> Observed Generative Rank (K_roots)     : %d", K_roots))
    message(sprintf(" -> Extended Signal Tail Rank (K_extended) : %d", K_extended))
    message("-----------------------------------------------------------------")
    message(" [Informational Gravity Metrics Estimates] (K_extended Volume Rebundled to K_roots)")
    message(sprintf(" -> Average Informational Gravity (AIG): %.3f Variable Equivalents", AIG))
    message(" -> Factor-Specific Informational Gravity (FSIG):")
    message(paste(capture.output(print(round(FSIG_final, 3))), collapse = "\n"))
    message(" -> Structural Topology Profile (Relative to FSIG_1):")
    message(paste(capture.output(print(round(STP_final, 3))), collapse = "\n"))
    message("===================================================================================\n")

    message("===================================================================================")
    message(sprintf(" (UPDATED) TRIPARTITE STRUCTURAL COMPOSITION (Rebundling Based on K_extended = %d)", K_extended))
    message("===================================================================================")
    message(sprintf(" -> %-50s : %d", "Valid Variables (m)", m_valid))
    message(sprintf(" -> %-50s : %.2f", "Centered Trace (Tr_Mc)", Tr_Mc))
    message("-----------------------------------------------------------------")
    message(sprintf(" -> %-50s : %.3f", "Total Shared Signal Volume (Unique + Redundant)", total_signal_volume))
    message(sprintf(" -> %-40s : %.3f", "   (Unique Signal Volume)", unique_signal_volume))
    message(sprintf(" -> %-40s : %.3f", "   (Redundant Signal Volume)", redundant_signal_volume))
    message("-----------------------------")
    message(sprintf(" -> %-50s : %.3f", "Idiosyncratic Informational Variance", idiosyncratic_volume))
    message("      (Structural Uncertainty + Independent Measurement Error")
    message("      + Unshared Signal Geometry)")
    message("-----------------------------------------------------------------")
    message(sprintf(" -> %-50s : %.3f%%", "% Total Shared Signal", new_pct_signal))
    message(sprintf(" -> %-40s : %.3f%%", "   (% Unique Signal)", new_pct_unique))
    message(sprintf(" -> %-40s : %.3f%%", "   (% Redundant Signal)", new_pct_redundant))
    message("-----------------------------")
    message(sprintf(" -> %-50s : %.3f%%", "% Idiosyncratic Informational Variance", new_pct_idio))
    message("      (Structural Uncertainty + Independent Measurement Error")
    message("      + Unshared Signal Geometry)")
    message("===================================================================================\n")

    message("===================================================================================")
    message(" (DELTA) METRIC SHIFT COMPARISON")
    message("===================================================================================")
    message(sprintf(" -> %-37s : %d -> %d", "Observed Generative Rank (K_roots)", orig_K_roots, K_roots))
    message(sprintf(" -> %-37s : %d -> %d", "Extended Signal Tail (K_extended)", orig_K_extended, K_extended))
    message("-----------------------------------------------------------------")
    message(sprintf(" -> %-37s : %.3f -> %.3f", "Average Info Gravity (AIG)", orig_AIG, AIG))
    message("-----------------------------------------------------------------")
    message(sprintf(" -> %-37s : %.3f -> %.3f", "Total Shared Signal Volume", orig_total_signal, total_signal_volume))
    message(sprintf(" -> %-37s : %.3f -> %.3f", "   (Unique Signal Volume)", orig_unique_signal, unique_signal_volume))
    message(sprintf(" -> %-37s : %.3f -> %.3f", "   (Redundant Signal Volume)", orig_redundant_signal, redundant_signal_volume))
    message("-----------------------------")
    message(sprintf(" -> %-37s : %.3f -> %.3f", "Idiosyncratic Variance Volume", orig_idio_vol, idiosyncratic_volume))
    message("-----------------------------------------------------------------")
    message(sprintf(" -> %-37s : %.3f%% -> %.3f%%", "% Total Shared Signal", orig_pct_signal, new_pct_signal))
    message(sprintf(" -> %-37s : %.3f%% -> %.3f%%", "   (% Unique Signal)", orig_pct_unique, new_pct_unique))
    message(sprintf(" -> %-37s : %.3f%% -> %.3f%%", "   (% Redundant Signal)", orig_pct_redundant, new_pct_redundant))
    message("-----------------------------")
    message(sprintf(" -> %-37s : %.3f%% -> %.3f%%", "% Idiosyncratic Variance", orig_pct_idio, new_pct_idio))
    message("===================================================================================\n")

    message("===================================================================================")
    message(" ENTROPIC SCREE (v1.0.1) - METHODOLOGICAL REFERENCE & LICENSE")
    message("===================================================================================")
    message(" -> Framework developed by Terrence J. Lee-St. John (Enli)")
    message(" -> Released under the Apache License 2.0 (Open Source)")
    message(" -> For full methods and metric definitions, see:")
    message("    The Entropic Scree (2026) - https://doi.org/10.5281/zenodo.22028087")
    message("===================================================================================\n")
  }

  # 6. Mutate and Return Object
  scree_obj$K_roots <- K_roots
  scree_obj$K_extended <- K_extended
  scree_obj$total_signal_volume <- total_signal_volume
  scree_obj$unique_signal_volume <- unique_signal_volume
  scree_obj$idiosyncratic_volume <- idiosyncratic_volume
  scree_obj$AIG <- AIG
  scree_obj$FSIG_final <- FSIG_final
  scree_obj$structural_topology_profile <- STP_final
  scree_obj$bipolar_modules <- bipolar_modules_out

  return(invisible(scree_obj))
}
