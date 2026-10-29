utils::globalVariables(c(".SD", "Eigenvalue", "Rank", ":="))

#' Entropic Scree Measurement Mapping and Rank Estimation
#'
#' @description Maps relationships among measurements using a double-centered,
#' joint-entropy-normalized mutual information matrix. Provides informational-axis
#' loadings, primary and extended spectral-boundary diagnostics, and spectral-volume
#' and informational-gravity summaries for discrete and discretized tabular data.
#' The selected Primary Informational Rank is a candidate boundary in the
#' information spectrum; a returned candidate does not establish a resolved elbow
#' or correspondence with generating dimension.
#' For the methods paper, see \doi{10.5281/zenodo.22028086}.
#'
#' @param data A data.table with numeric/integer, factor, character, or logical columns and unique names. Numeric category labels should be supplied as factors if they represent nominal states. Base data.frames are not supported.
#' @param low_entropy_thresh Numeric. Minimum marginal entropy in nats for original measurements and unprotected indicators. Default 0.05. Zero-entropy original columns are always excluded; an indicator with a retained parent is protected.
#' @param num_bins Integer from 2 through 46340. Numeric discretization target and native counting capacity. NULL uses max(2, ceiling(bin_multiplier * N^(1/3))). It must cover every retained categorical vocabulary, including an opted-in missing state; violations are rejected. Increasing it also changes the target for numeric measurements.
#' @param bin_multiplier Numeric. Multiplier for the automatic bin calculation. Default is 1.0.
#' @param num_cores Integer. Number of CPU cores for OpenMP parallelization.
#' @param interactive_mode Logical. If TRUE, displays diagnostics and, in an interactive R session, prompts for the primary rank and Extended Signal Tail boundary.
#' @param purge_constants Logical. If TRUE, removes constant measurements and identical duplicates before discretization.
#' @param check_collinearity Logical. If TRUE, checks numeric measurements for exact linear dependencies and removes dependent measurements when the check can be performed.
#' @param triple_tap_window Positive integer requested window w0 (default 20). The actual window is max(5, min(w0, floor(0.1 * m))).
#' @param fwer_alpha Numeric. Alpha used in the multiple-comparison-adjusted local t threshold. Default is 0.01; this does not establish whole-pipeline error control.
#' @param extract_eigenvectors Logical. Defaults to TRUE. Returns all retained measurements' loadings on every selected primary informational axis, plus the full eigenvector matrix. To omit both, set this and extract_bipolar_modules to FALSE.
#' @param extract_bipolar_modules Logical. Defaults to TRUE. Returns positive and negative loading lists for each primary informational axis; TRUE also enables eigenvector extraction.
#' @param bipolar_top_n NULL or numeric. NULL (default) keeps every nonzero loading in its corresponding pole. A number strictly between zero and one selects up to that proportion of the retained measurement count per pole; a positive integer selects up to that many per pole. This option never truncates measurement_loadings or eigenvectors.
#' @param return_processed_data Logical. If TRUE, returns the retained analysis table in original units only if measurements or values changed through filtering or opted-in missing-data handling. Internal MI binning/coding alone does not count as a data change. If unchanged, processed_data is NULL and processed_data_status explicitly says to use the original data. No file is written automatically.
#' @param categorical_missing Character. "error" (default) rejects categorical missing values. "category" explicitly represents missing factor, character, or logical values as a distinct, collision-free missing state; this is not imputation of an observed category.
#' @param numeric_missing Character. "error" (default) rejects numeric NA/NaN. "median_indicator" fills them with the column's observed median and adds a binary original-missingness companion. A wholly missing numeric column, or an infinite numeric value, is rejected under either policy.
#'
#' @details
#' Rows of \code{measurement_loadings} identify retained measurements and columns
#' \code{Axis_1}, ..., \code{Axis_K} identify the selected primary informational
#' axes, where \code{K = K_roots}. Every loading is retained, including zeros.
#' These are raw signed eigenvector coefficients (equal-axis coordinates), not
#' eigenvalue-scaled coordinates, measurement-to-score correlations, or native
#' observation-level factor scores. For eigenvalue-scaled coordinates on a chosen
#' positive-axis set, multiply each eigenvector column by the square root of its
#' corresponding \code{raw_eigenvalues} value. Do not use floored working weights
#' to treat nonpositive raw directions as positive geometry.
#'
#' The full \code{eigenvectors} matrix includes all spectral directions for later
#' inspection or rank updates; its presence does not mean all directions were
#' selected. Eigenvector signs are arbitrary, and repeated eigenvalues can admit
#' different bases for the same subspace. Opposite poles express informational
#' contrasts, not positive versus negative correlation between measurements.
#'
#' The Extended Signal Tail lies after \code{K_roots} and ends inclusively at
#' \code{K_extended}; it is empty when the boundaries coincide. Spectral-volume
#' allocations are distinct from coordinate geometry and are not measured
#' proportions of causal signal or measurement error. Factor-Specific
#' Informational Gravity refers to an informational axis, not an identified
#' latent generating factor.
#'
#' \strong{Missing-data preparation.} Both missing-data options are off by
#' default. Numeric medians use observed finite values only. Each generated
#' companion is 0 for an originally observed value and 1 for an originally
#' missing value. Identical row-by-row missingness patterns share one companion
#' with all associated parent names, avoiding duplicate indicator weighting.
#' Names and missing-category labels never overwrite or merge existing values.
#' The original input object is not modified. Imputation and explicit missing
#' states change the represented data; results are not missingness-conditioned MI.
#'
#' Original measurements undergo the usual constant, duplicate, numeric QR, and
#' entropy checks first. Indicators do not enter the numeric QR check. An
#' indicator is mandatory when any associated numeric parent survives, even
#' below the entropy threshold or duplicating an encoded original measurement.
#' Mandatory indicators precede orphaned indicators in filtering. When no parent
#' survives, a companion faces normal binary constant, duplicate and entropy
#' filtering against retained originals and mandatory companions. Constant
#' required indicators are invalid. Complete columns generate no indicator.
#' Protection depends on retaining at least one associated numeric parent,
#' not on the indicator's marginal entropy or uniqueness. Numeric filling
#' uses the ordinary sample median of observed values.
#'
#' \strong{Metric definitions and paper notation.} On common prepared records,
#' marginal entropy is \eqn{H(X)=-\sum p(x)\ln p(x)} in nats, and
#' \eqn{I(X;Y)=H(X)+H(Y)-H(X,Y)}. Joint-entropy normalization gives
#' \eqn{\mathcal{M}_{ij}=I(X_i;X_j)/H(X_i,X_j)}. Double-centering produces
#' \eqn{\mathcal{M}_c=\mathbf{H}\mathcal{M}\mathbf{H}}, with
#' \eqn{\mathbf{H}=\mathbf{I}-\mathbf{1}\mathbf{1}^{T}/m}; its trace is the sum
#' of the raw eigenvalues \eqn{\lambda_i}. The detector's mean reference is
#' \eqn{\mu=\operatorname{Tr}(\mathcal{M}_c)/m}; its baseline count is
#' \eqn{b=\#\{i:\eta_i>\mu\}}, not a comparison with the mean working weight.
#'
#' Working weights are \eqn{\eta_i=\max(\varepsilon,\lambda_i)}, with
#' \eqn{\varepsilon=10^{-9}}, \eqn{S_{\varepsilon}=\sum_{i=1}^{m}\eta_i},
#' and \eqn{p_i=\eta_i/S_{\varepsilon}} over all m directions.
#' \eqn{R_{eff}=\exp(-\sum_i p_i\ln p_i)}. With inclusive boundaries
#' \eqn{1\le k=K_{roots}\le e=K_{extended}\le m}, define
#' \eqn{P_{sig}=\sum_{i=1}^{K_{extended}}\eta_i/S_{\varepsilon}}.
#' Redundant Signal Volume is \eqn{m-R_{eff}}, Unique Signal Volume is
#' \eqn{P_{sig}R_{eff}}, and Idiosyncratic Informational Volume is
#' \eqn{(1-P_{sig})R_{eff}}. Their sum is m. Total Shared Signal Volume is
#' \eqn{TSV=(m-R_{eff})+P_{sig}R_{eff}}. Average Informational Gravity is
#' \eqn{\mathrm{AIG}=TSV/K_{roots}}. Factor-Specific Informational Gravity is
#' \eqn{FSIG_i=TSV\eta_i/\sum_{j=1}^{K_{roots}}\eta_j}, and the Structural
#' Topology Profile is \eqn{T_i=FSIG_i/FSIG_1=\eta_i/\eta_1} over primary axes.
#' Printed allocation percentages use m as denominator; the spectral fraction
#' \eqn{P_{sig}} is different from \eqn{TSV/m}.
#'
#' Raw positive mass is \eqn{m_+=\sum_i\max(0,\lambda_i)} and raw negative
#' mass is \eqn{m_-=\sum_{\lambda_i<0}|\lambda_i|}. The Negative Spectral
#' Mass Ratio is \eqn{\mathrm{NSMR}=m_-/m_+} for \eqn{m_+>0}, before flooring;
#' \code{NSMR_percent} is exactly \eqn{100\,\mathrm{NSMR}}. Raw mass \eqn{m_+} is distinct
#' from working-weight sum \eqn{S_{\varepsilon}}. Internally, \code{m_plus}
#' denotes only raw positive mass and \code{S_epsilon} denotes the working sum.
#'
#' \strong{R identifiers.} \code{raw_eigenvalues} corresponds to \eqn{\lambda_i},
#' \code{eigenvalues} to \eqn{\eta_i}, and \code{similarity_matrix} to
#' \eqn{\mathcal{M}_c}. \code{K_log_gap} and \code{K_triple_tap} correspond
#' to the paper's automatic recommendations \eqn{K_{log}} and \eqn{K_{tap}};
#' \code{K_roots} and \code{K_extended} are the selected boundaries.
#' \code{total_signal_volume} corresponds to \eqn{TSV}, \code{FSIG_final}
#' to \eqn{FSIG_i}, and \code{structural_topology_profile} to \eqn{T_i}.
#' \code{measurement_loadings} contains \eqn{v_{ik}} on the selected primary
#' axes, with every retained measurement included. \code{Axis_1},
#' \code{Axis_2}, etc. name the corresponding informational-axis pole lists.
#' Auxiliary allocations use \code{FSIG_auxiliary_bulk} and
#' \code{FSIG_auxiliary_mean_reference}; neither is a unit-eigenvalue Kaiser rule.
#'
#' The auxiliary bulk endpoint is the last axis before a detected macro gap;
#' without such a gap it is b. Its allocation cutoff is max(e, endpoint).
#' The mean-reference auxiliary cutoff is max(e, b). Each auxiliary applies the
#' same TSV formula at its own cutoff and distributes that volume across the
#' current k primary proportions. These are implementation diagnostics, not
#' separately defined primary metrics in the paper. All cutoffs and their sources
#' are returned and all auxiliary vectors refresh on a rank update.
#'
#' \strong{Numerical geometry and interpretation.} A raw axis is numerically
#' positive when its eigenvalue exceeds
#' \eqn{m\epsilon_{machine}\max(1,\max_i|\lambda_i|)}. An entirely degenerate
#' spectrum is rejected before floor-based ranks or gravity can be presented.
#' The tolerance is a numerical guard, not a spectral boundary estimator.
#' Selected coefficients on null or negative raw directions may still be exported
#' algebraically, with a warning and explicit status; they do not constitute a
#' wholly positive-spectrum coordinate map. A working spectrum entirely at the
#' floor also produces a warning. The working-weight and allocation formulas are
#' not altered by this diagnostic. Local t thresholds are nominal heuristic
#' thresholds, not established family-wise error control.
#'
#' \strong{Data return.} processed_data contains retained original-unit values,
#' opted-in median fills or explicit categorical missing labels, and retained
#' 0/1 companions.
#' Internal discretization and category codes are solely MI-estimation inputs
#' and do not trigger export. Integer versus double storage alone is not a
#' change in numeric observations. The status records changed, added and dropped
#' columns, or explains that the original input can be reused.
#'
#' @section Output structure:
#' Axis-specific pole lists use \code{Axis_1}, \code{Axis_2}, and so on,
#' matching the columns of \code{measurement_loadings}. Each pole list contains
#' \code{Positive_Anchor} and \code{Negative_Anchor} coefficients.
#' \code{NSMR} is a ratio; \code{NSMR_percent} is 100 times that ratio.
#'
#' Within \code{auxiliary_diagnostics}, \code{macro_boundary_index} identifies
#' the first axis after a detected macro gap (NA when no gap is detected).
#' \code{macro_endpoint} is the inclusive endpoint, and \code{macro_endpoint_source}
#' identifies whether it comes from a detected gap or the mean-reference count.
#' \code{bulk_allocation_cutoff} is the inclusive auxiliary allocation endpoint.
#' \code{Update.Entropic.Scree()} returns a result without modifying its input
#' object in place.
#'
#' @return A list containing the following components:
#' \item{eigenvalues}{Working spectral weights: raw eigenvalues floored at 1e-9, in descending order.}
#' \item{raw_eigenvalues}{Eigenvalues of the centered information matrix before flooring.}
#' \item{similarity_matrix}{The double-centered, joint-entropy-normalized mutual information matrix.}
#' \item{retained_features}{Names of measurements retained after preparation and entropy filtering.}
#' \item{bin_distributions}{Per-measurement lists or arrays of category/bin counts.}
#' \item{centered_trace}{Sum of the diagonal of the centered information matrix.}
#' \item{NSMR}{Negative Spectral Mass Ratio: absolute negative raw eigenvalue mass divided by positive raw eigenvalue mass, as a ratio.}
#' \item{NSMR_percent}{Negative Spectral Mass Ratio in percentage units: exactly 100 * NSMR.}
#' \item{R_eff}{Total Unique Probabilistic Volume: entropy effective rank of the entire working spectrum.}
#' \item{K_log_gap}{Automatic Log-Gap recommendation for the Primary Informational Rank.}
#' \item{K_triple_tap}{Automatic Triple-Tap recommendation for the inclusive Extended Signal Tail boundary.}
#' \item{triple_tap_multiplier}{Nominal local t-multiplier at an accepted departure, or the full-window reference multiplier when the scan does not trigger. NA when the search is unavailable; rank_diagnostics records fallback status.}
#' \item{K_roots}{Selected Primary Informational Rank; not by definition the Intrinsic Generative Rank r.}
#' \item{K_extended}{Selected inclusive Extended Signal Tail boundary used for volume allocation.}
#' \item{total_signal_volume}{Total Shared Signal Volume, in variable equivalents.}
#' \item{unique_signal_volume}{Unique Signal Volume allocated to the selected signal region.}
#' \item{redundant_signal_volume}{Redundant Signal Volume: retained measurement count minus R_eff.}
#' \item{idiosyncratic_volume}{Idiosyncratic Informational Volume: remaining effective-volume allocation outside the selected signal region.}
#' \item{AIG}{Average Informational Gravity, in variable equivalents per primary informational axis.}
#' \item{FSIG_final}{Factor-Specific Informational Gravity allocated across the primary informational axes.}
#' \item{structural_topology_profile}{Structural Topology Profile: primary-axis FSIG divided by the leading-axis FSIG.}
#' \item{measurement_loadings}{Named matrix of every retained measurement's loading on every selected primary informational axis, returned by default. NULL only when both extraction options are FALSE.}
#' \item{eigenvectors}{Named full eigenvector matrix, including unselected directions, returned by default. NULL only when both extraction options are FALSE.}
#' \item{bipolar_modules}{Positive and negative loading lists for each selected primary axis, containing every nonzero loading per pole by default. Zero coefficients remain in measurement_loadings but belong to neither pole. List keys Axis_1, Axis_2, etc. match the measurement_loadings columns. Each contains Positive_Anchor and Negative_Anchor loading lists for that informational axis. NULL when extract_bipolar_modules is FALSE.}
#' \item{processed_data}{Retained original-unit analysis data.table only when requested and changed; otherwise NULL. Internal MI codes are not returned here.}
#' \item{processed_data_status}{List with requested, changed (NA if not checked), returned, message, and, when requested, changed_columns, dropped_columns and added_columns. An unchanged input produces an explicit message and no duplicate table.}
#' \item{preparation_info}{Policies, input/prepared row counts, bin target, exclusions, categorical missing labels/counts/codes, observed numeric medians and parent links. Each indicator receipt lists its parents, surviving parents, retained/protected status, ordinary filter findings, entropy in nats, and 0/1 meanings (internally coded 1/2).}
#' \item{missingness_indicators}{Retained raw 0/1 companion columns when analysis-data return is requested; NULL when none or not requested. Identical missingness patterns share a column.}
#' \item{positive_axis_diagnostics}{Numerical tolerance, positive raw-axis indices, whether all selected axes are positive, whether all working weights equal the floor, and floor_weight_fraction = sum(eta - pmax(lambda,0))/sum(eta). These diagnostics are unavailable when the input result has no raw eigenvalues.}
#' \item{rank_diagnostics}{Selection status, primary detector path, automatic boundaries, extended fallback status, search endpoint, actual local window, nominal alpha, boundary eigenvalue ratio and percent drop, and macro-stitch status. These do not assert a resolved elbow.}
#' \item{auxiliary_diagnostics}{Macro boundary index (first after the gap, NA if absent), inclusive macro endpoint and source, unfloored mean-eigenvalue reference and its working-weight count, and both current inclusive auxiliary allocation cutoffs.}
#' \item{FSIG_auxiliary_bulk}{Current primary-axis allocation using the broader of K_extended and the explicit macro endpoint, or mean-reference count when no macro gap exists.}
#' \item{FSIG_auxiliary_mean_reference}{Current primary-axis allocation using the broader of K_extended and the count of working weights above the mean raw eigenvalue.}
#'
#' @examples
#' # Generate a small random dataset
#' dummy_data <- data.table::data.table(
#'   V1 = rnorm(50),
#'   V2 = rnorm(50),
#'   V3 = rnorm(50)
#' )
#'
#' # Run without interactive rank-selection prompts
#' res <- Entropic.Scree(dummy_data, interactive_mode = FALSE, num_cores = 1)
#' res$measurement_loadings
#' # All measurements on one selected axis, with no loading cutoff:
#' res$measurement_loadings[, 1, drop = FALSE]
#' # Pole summaries use the same informational-axis names:
#' res$bipolar_modules$Axis_1
#' res$NSMR_percent
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
                           , extract_eigenvectors = TRUE
                           , extract_bipolar_modules = TRUE
                           , bipolar_top_n = NULL
                           , return_processed_data = FALSE
                           , categorical_missing = c("error", "category")
                           , numeric_missing = c("error", "median_indicator")) {

  # Enforce data.table requirement
  if (!data.table::is.data.table(data)) {
    stop("Error: 'data' must be a data.table. Please convert your dataset using data.table::as.data.table() before running Entropic.Scree.")
  }

  # Enforce eigenvector extraction if bipolar modules are requested
  if (extract_bipolar_modules) {
    extract_eigenvectors <- TRUE
  }

  .resolve_bipolar_top_n(bipolar_top_n, ncol(data))

  categorical_missing <- match.arg(categorical_missing)
  numeric_missing <- match.arg(numeric_missing)
  .scalar_number(triple_tap_window, "triple_tap_window", 1, integer = TRUE)
  .scalar_number(fwer_alpha, "fwer_alpha", .Machine$double.eps, 1 - .Machine$double.eps)
  .scalar_number(num_cores, "num_cores", -Inf, integer = TRUE)
  preparation <- .prepare_scree_data(data, num_bins, bin_multiplier, low_entropy_thresh,
    purge_constants, check_collinearity, categorical_missing, numeric_missing,
    return_analysis_data = return_processed_data)
  dt <- preparation$data
  H_vec <- preparation$entropy
  num_bins <- preparation$num_bins
  preparation_info <- preparation$information
  data_export <- .analysis_data_export(data, preparation$analysis_data, return_processed_data)
  missingness_indicators_out <- if (return_processed_data) preparation$raw_indicators else NULL
  rm(preparation)
  valid_vars <- names(dt)
  p <- ncol(dt)

  message(sprintf("[6/10] Computing %d x %d Mutual Information Matrix (C++ OpenMP)...", p, p))

  processed_data_out <- data_export$data

  mat_data <- as.matrix(dt)
  bin_sample_sizes <- apply(mat_data, 2, tabulate)
  rm(dt)
  gc(verbose = FALSE)

  safe_target_cores <- max(1, num_cores)
  MI_mat <- fast_parallel_MI(mat_data, num_bins = num_bins, cores = safe_target_cores)
  rownames(MI_mat) <- colnames(MI_mat) <- valid_vars

  message("[7/10] Applying Joint-Entropy Normalization...")
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

  message("[9/10] Extracting Informational Axes (Eigen Decomposition)...")
  if (extract_eigenvectors) {
    eigen_res <- eigen(NMI_mat_c, symmetric = TRUE)
    raw_eig_vals <- eigen_res$values
    vectors_out <- eigen_res$vectors
    dimnames(vectors_out) <- list(valid_vars, paste0("Axis_", seq_len(ncol(vectors_out))))
  } else {
    eigen_res <- eigen(NMI_mat_c, symmetric = TRUE, only.values = TRUE)
    raw_eig_vals <- eigen_res$values
    vectors_out <- NULL
  }

  positive_info <- .positive_axis_diagnostics(raw_eig_vals)
  if (identical(positive_info$status, "degenerate"))
    stop("Degenerate centered spectrum: no positive eigenvalue exceeds the documented numerical tolerance. Flooring cannot establish resolved informational axes.")

  # Calculate the NSMR ratio and percentage before flooring
  m_plus <- sum(raw_eig_vals[raw_eig_vals > 0])
  m_minus <- sum(abs(raw_eig_vals[raw_eig_vals < 0]))
  NSMR <- if (m_plus > 0) m_minus / m_plus else NA_real_
  NSMR_percent <- NSMR * 100

  eig_vals <- pmax(raw_eig_vals, 1e-9)

  # Count true positive dimensions before the geometric zero to bound plots
  max_plot_idx <- sum(raw_eig_vals > 1e-12)

  # Total working spectral weight (sum of floored eigenvalues)
  S_epsilon <- sum(eig_vals)

  message("[10/10] Calculating Total Unique Probabilistic Volume (R_eff) and Estimating Spectral Boundaries...")
  sig_vals <- eig_vals[eig_vals > 0]
  if (length(sig_vals) > 0) {
    p_vals <- sig_vals / sum(sig_vals)
    H_spec <- -sum(p_vals * log(p_vals))
    R_eff <- exp(H_spec)
  } else {
    R_eff <- 1
  }

  # ==========================================================================
  # STEP 1: MACRO GAP BOUNDARY (AUXILIARY MACRO-GAP DETECTION)
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
  # This engine isolates the Primary Informational Rank (K_roots)
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

    # 2. Nominal multiplicity-adjusted local threshold; family-wise calibration is not established
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
      sub_title_text <- sprintf("Converged: Primary Informational Rank & Extended Signal Tail boundary = %d", K_log_gap)
    } else {
      lines_geom <- list(
        ggplot2::geom_vline(xintercept = K_log_gap, color = "#D55E00", linetype = "dashed", linewidth = 1.2),
        ggplot2::geom_vline(xintercept = K_triple_tap, color = "purple", linetype = "twodash", linewidth = 1.2)
      )
      sub_title_text <- sprintf("Primary Informational Rank (Orange) = %d | Extended Signal Tail boundary (Purple) = %d", K_log_gap, K_triple_tap)
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
      ggplot2::labs(title = "Zoomed View", x = "Eigenvalue Index", y = "Working Eigenvalue (log scale)") +
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
      ggplot2::labs(title = "Macro View", x = "Eigenvalue Index", y = "Working Eigenvalue (log scale)") +
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
    message(sprintf(" -> %-50s : %d", "Retained Measurements (m)", m_valid))
    message(sprintf(" -> %-50s : %.3f", "Centered Trace (Tr_Mc)", Tr_Mc))
    message(sprintf(" -> %-50s : %.3f%%", "Negative Spectral Mass Ratio (NSMR; %)", NSMR_percent))
    message("-----------------------------------------------------------------")
    message(sprintf(" -> %-50s : %.3f", "Total Unique Probabilistic Volume (R_eff)", R_eff))
    message(sprintf(" -> %-50s : %.3f%%", "%", pct_prob_volume))
    message("      (Entropy effective rank of the full working spectrum)")
    message("-----------------------------------------------------------------")
    message(sprintf(" -> %-50s : %.3f", "Redundant Signal Volume (m - R_eff)", redundant_signal_volume))
    message(sprintf(" -> %-50s : %.3f%%", "%", pct_redundant_signal))
    message("-----------------------------------------------------------------")
    message(sprintf(" -> %-50s : %d", "Working Eigenvalues > Tr(M_c)/m", n_eigen_gt_mean))
    message(sprintf(" -> %-50s : %d", "Working Eigenvalues <= Tr(M_c)/m", n_eigen_le_mean))
    message("===================================================================================\n")

    message("===================================================================================")
    message(" AUTOMATED ELBOW DETECTION (DUAL-DIAGNOSTIC ENSEMBLE)")
    message("===================================================================================")
    message(" [Diagnostic: Continuous Bulk Boundary (Macro Cliff)]")
    if (!is.na(top_of_bulk_idx)) {
      message(sprintf(" -> %-43s : Index %d", "Auxiliary Macro-Boundary Index", top_of_bulk_idx))
      message("    (First index after the detected macro gap)")
      message(sprintf(" -> %-43s : %.6f", "Macro Gap Baseline (Max Bulk Gap)", macro_max_bulk_gap))
      message(sprintf(" -> %-43s : %.2fx Baseline", "Actual Macro Gap Magnitude", macro_gap_ratio))
    } else {
      message(sprintf(" -> %-43s : %s", "Auxiliary Macro-Boundary Index", "No macro gap (using mean-trace count)"))
    }
    message("-----------------------------------------------------------------")
    message(" [Engine A: Log-Gap]")
    message(sprintf(" -> %-43s : %d%s", "Primary Informational Rank (K_roots)", K_log_gap, ifelse(log_gap_fallback, ifelse(is.na(top_of_bulk_idx), " [No eligible secondary gap; fallback]", " [Low-rank macro fallback]"), "")))
    if (!is.na(log_gap_ratio)) {
      message(sprintf(" -> %-43s : %.2fx (%.1f%% Drop)", "Primary-Boundary Eigenvalue Ratio", log_gap_ratio, log_gap_pct_drop))
    }
    message("-----------------------------------------------------------------")
    message(sprintf(" [Engine B: Triple-Tap (%d-Point Linear)]", triple_tap_window))
    message(sprintf(" -> %-43s : %d%s", "Extended Signal Tail boundary (K_extended)", K_triple_tap, ifelse(triple_tap_fallback, " [Failed to trigger; defaulted to Log-Gap]", "")))
    message(sprintf(" -> %-43s : %s", "Macro-Stitch Applied", ifelse(stitch_applied, "Yes (+0.5x Macro Gap to Bulk)", "No")))
    if (!is.na(triple_tap_expected_val)) {
      message(sprintf(" -> %-43s : %.5f", "Expected Local Eigenvalue", triple_tap_expected_val))
      message(sprintf(" -> %-43s : %.5f", "Actual Local Eigenvalue", triple_tap_actual_val))
      message(sprintf(" -> %-43s : %.6f", "Calculated Local Baseline Sigma (SD)", triple_tap_local_sd))
    }
    if (!is.na(triple_tap_multiplier)) {
      message(sprintf(" -> %-43s : %.2f (nominal quantile = %.6f)", "Required Local t-Multiplier", triple_tap_multiplier, prob_target))
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
  rank_selection <- "automatic"
  .validate_ranks(K_roots, K_extended, m_valid)

  if (interactive_mode) {
    if (!interactive()) {
      # CRAN SAFETY NET: Force bypass if run in a headless test environment
      message("[-] Non-interactive session detected. Bypassing manual rank confirmation.")
    } else {
      message("\n[WARNING]: The automated extractors rely on statistical heuristics and")
      message("may not perfectly align with the true structural elbow of your specific dataset.")
      message("Please visually examine the generated scree plot.\n")

      while (TRUE) {
        # 1. Prompt for Primary Informational Rank
        ans_r <- trimws(readline(prompt = sprintf("Enter the Primary Informational Rank (K_roots) [Press Enter to keep %d]: ", K_roots)))
        candidate_k <- if (ans_r == "") K_roots else suppressWarnings(as.numeric(ans_r))
        ans_ext <- trimws(readline(prompt = sprintf("Enter the Extended Signal Tail boundary [Press Enter to keep %d]: ", K_extended)))
        candidate_e <- if (ans_ext == "") K_extended else suppressWarnings(as.numeric(ans_ext))
        rank_error <- tryCatch({ .validate_ranks(candidate_k, candidate_e, m_valid); NULL },
                               error = function(e) conditionMessage(e))
        if (!is.null(rank_error)) {
          message("[-] ", rank_error, " Please choose valid ranks.")
          next
        }
        K_roots <- as.integer(candidate_k)
        K_extended <- as.integer(candidate_e)

        # ====================================================================
        # INTERACTIVE GRAPH PREVIEW & METRICS
        # ====================================================================
        if (requireNamespace("ggplot2", quietly = TRUE)) {

          # Draw the automatic detector recommendations in gray
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
            ggplot2::labs(title = "Zoomed View", x = "Eigenvalue Index", y = "Working Eigenvalue (log scale)") +
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
            ggplot2::labs(title = "Macro View", x = "Eigenvalue Index", y = "Working Eigenvalue (log scale)") +
            ggplot2::theme_minimal(base_size = 14) +
            ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 12))

          # Render Side-by-Side
          if (requireNamespace("patchwork", quietly = TRUE)) {
            combined_plot_upd <- (p_scree_macro_upd + p_scree_zoom_upd) +
              patchwork::plot_annotation(
                title = "Entropic Scree Results",
                subtitle = sprintf("Selection Preview: Primary Informational Rank (Solid) = %d | Extended Signal Tail boundary (Dashed) = %d", K_roots, K_extended),
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
        selected_spectral_weight_prev <- sum(eig_vals[1:K_extended])
        P_sig_prev <- selected_spectral_weight_prev / S_epsilon
        unique_signal_volume_prev <- R_eff * P_sig_prev
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
        message(sprintf(" -> Average Informational Gravity (AIG): %.3f Variable Equivalents per Primary Axis", AIG_prev))
        message(" -> Factor-Specific Informational Gravity (FSIG):")
        message(paste(capture.output(print(round(FSIG_prev, 3))), collapse = "\n"))
        message(" -> Structural Topology Profile (Relative to FSIG_1):")
        message(paste(capture.output(print(round(STP_prev, 3))), collapse = "\n"))
        message("===================================================================================\n")

        ans_confirm <- trimws(readline(prompt = "Do you want to finalize these ranks? (Y/N): "))
        if (tolower(ans_confirm) %in% c("y", "yes")) {
          rank_selection <- "user_confirmed"
          message("\n[+] Finalizing rank selections.")
          break
        }
      }
    }
  }

  # --- FINAL GRAVITY CALCULATIONS ---
  # 1. Capture Total Shared Signal Volume strictly up to the Extended Signal Tail
  selected_spectral_weight <- sum(eig_vals[1:K_extended])
  P_sig <- selected_spectral_weight / S_epsilon
  unique_signal_volume <- R_eff * P_sig
  total_signal_volume <- unique_signal_volume + redundant_signal_volume

  # 2. Rebundle AIG and FSIG strictly into the Primary Informational Rank (K_roots)
  AIG <- total_signal_volume / K_roots
  core_eigenvals <- eig_vals[1:K_roots]
  p_core <- core_eigenvals / sum(core_eigenvals)
  FSIG_final <- p_core * total_signal_volume

  elbow_label <- ifelse(rank_selection == "user_confirmed", "User-Confirmed", "Automatic Candidate")
  positive_info <- .check_selected_geometry(raw_eig_vals, K_roots)
  STP_final <- FSIG_final / FSIG_final[1]

  trailing_spectral_weight <- if (K_extended < m_valid) sum(eig_vals[(K_extended + 1):m_valid]) else 0
  idiosyncratic_weight <- trailing_spectral_weight / S_epsilon
  idiosyncratic_volume <- R_eff * idiosyncratic_weight

  if (interactive_mode) {
    message("\n===================================================================================")
    message(sprintf(" (FINAL) ENTROPIC SCREE METRICS (based on %s)", elbow_label))
    message("===================================================================================")
    message(" [Rank Estimates]")
    message(sprintf(" -> Primary Informational Rank (K_roots)     : %d", K_roots))
    message(sprintf(" -> Extended Signal Tail boundary (K_extended) : %d", K_extended))
    message("-----------------------------------------------------------------")
    message(" [Informational Gravity Metrics Estimates] (K_extended Volume Rebundled to K_roots)")
    message(sprintf(" -> Average Informational Gravity (AIG): %.3f Variable Equivalents per Primary Axis", AIG))
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
    message(sprintf(" -> %-50s : %d", "Retained Measurements (m)", m_valid))
    message(sprintf(" -> %-50s : %.2f", "Centered Trace (Tr_Mc)", Tr_Mc))
    message("-----------------------------------------------------------------")
    message(sprintf(" -> %-50s : %.3f", "Total Shared Signal Volume (Unique + Redundant)", total_signal_volume))
    message(sprintf(" -> %-40s : %.3f", "   (Unique Signal Volume)", unique_signal_volume))
    message(sprintf(" -> %-40s : %.3f", "   (Redundant Signal Volume)", redundant_signal_volume))
    message("-----------------------------")
    message(sprintf(" -> %-50s : %.3f", "Idiosyncratic Informational Volume", idiosyncratic_volume))
    message("      (Residual spectral allocation; not a measured error fraction)")
    message("-----------------------------------------------------------------")
    message(sprintf(" -> %-50s : %.3f%%", "% Total Shared Signal", (total_signal_volume / m_valid) * 100))
    message(sprintf(" -> %-40s : %.3f%%", "   (% Unique Signal)", (unique_signal_volume / m_valid) * 100))
    message(sprintf(" -> %-40s : %.3f%%", "   (% Redundant Signal)", (redundant_signal_volume / m_valid) * 100))
    message("-----------------------------")
    message(sprintf(" -> %-50s : %.3f%%", "% Idiosyncratic Informational Volume", (idiosyncratic_volume / m_valid) * 100))
    message("      (Residual spectral allocation; not a measured error fraction)")
    message("===================================================================================\n")

    # ============================================================================
    # WAVE 4: METHODOLOGICAL REFERENCE
    # ============================================================================
    message("===================================================================================")
    message(" ENTROPIC SCREE (v1.0.2) - METHODOLOGICAL REFERENCE & LICENSE")
    message("===================================================================================")
    message(" -> Framework developed by Terrence J. Lee-St. John (Enli)")
    message(" -> Released under the Apache License 2.0 (Open Source)")
    message(" -> For full methods and metric definitions, see:")
    message("    The Entropic Scree (2026) - https://doi.org/10.5281/zenodo.22028086")
    message("===================================================================================\n")
  }

  # ============================================================================
  # AUXILIARY Factor-Specific Informational Gravity (FSIG) CALCULATIONS
  # ============================================================================

  auxiliary <- .auxiliary_gravity(eig_vals, R_eff, K_roots, K_extended,
                                 mean_trace, top_of_bulk_idx)
  rank_diagnostics <- list(selection = rank_selection,
    primary_method = if (n_valid_search < 3) "short_spectrum_fallback" else
      if (log_gap_fallback) "boundary_fallback" else
      if (!is.na(top_of_bulk_idx) && top_of_bulk_idx < 4) "low_rank_secondary_gap" else "secondary_log_gap",
    automatic_primary_rank = K_log_gap,
    automatic_extended_boundary = max(K_log_gap, K_triple_tap),
    extended_fallback = triple_tap_fallback, search_endpoint = search_start_idx,
    local_window = triple_tap_window, nominal_alpha = fwer_alpha,
    boundary_eigenvalue_ratio = log_gap_ratio, boundary_drop_percent = log_gap_pct_drop,
    stitch_applied = stitch_applied)

  # ============================================================================
  # EXTRACT BIPOLAR MODULES (INFORMATIONAL-AXIS POLES)
  # ============================================================================
  bipolar_modules_out <- NULL
  if (extract_bipolar_modules && !is.null(vectors_out) && K_roots > 0) {
    bipolar_modules_out <- list()

    # Dynamically resolve top N (proportion vs absolute count)
    resolved_top_n <- .resolve_bipolar_top_n(bipolar_top_n, m_valid)

    for (k in 1:K_roots) {
      # Bind variable names to the specific eigenvector
      loadings <- setNames(vectors_out[, k], valid_vars)

      # Keep all positive loadings by default, ordered strongest first
      pos_pole <- sort(loadings[loadings > 0], decreasing = TRUE)
      pos_pole <- head(pos_pole, resolved_top_n)

      # Keep all negative loadings by default, ordered strongest first
      neg_pole <- sort(loadings[loadings < 0], decreasing = FALSE)
      neg_pole <- head(neg_pole, resolved_top_n)

      bipolar_modules_out[[paste0("Axis_", k)]] <- list(
        Positive_Anchor = pos_pole,
        Negative_Anchor = neg_pole
      )
    }
  }

  return(list(
    eigenvalues = eig_vals,
    raw_eigenvalues = raw_eig_vals,
    centered_trace = Tr_Mc,
    NSMR = NSMR,
    NSMR_percent = NSMR_percent,
    similarity_matrix = NMI_mat_c,
    retained_features = valid_vars,
    bin_distributions = bin_sample_sizes,
    R_eff = R_eff,
    K_log_gap = K_log_gap,
    K_triple_tap = K_triple_tap,
    triple_tap_multiplier = triple_tap_multiplier,
    K_roots = K_roots,
    K_extended = K_extended,
    total_signal_volume = total_signal_volume,
    unique_signal_volume = unique_signal_volume,
    redundant_signal_volume = redundant_signal_volume,
    idiosyncratic_volume = idiosyncratic_volume,
    AIG = AIG,
    FSIG_final = FSIG_final,
    structural_topology_profile = STP_final,
    measurement_loadings = .primary_measurement_loadings(vectors_out, valid_vars, K_roots),
    eigenvectors = vectors_out,
    bipolar_modules = bipolar_modules_out,
    processed_data = processed_data_out,
    processed_data_status = data_export$status,
    preparation_info = preparation_info,
    missingness_indicators = missingness_indicators_out,
    positive_axis_diagnostics = positive_info,
    rank_diagnostics = rank_diagnostics,
    auxiliary_diagnostics = auxiliary$diagnostics,
    FSIG_auxiliary_bulk = auxiliary$bulk,
    FSIG_auxiliary_mean_reference = auxiliary$mean_reference
  ))
}

# ==============================================================================
# Update.Entropic.Scree (Post-Hoc Rank Modification)
# ==============================================================================

#' Update Entropic Scree Spectral Boundaries
#'
#' @description Changes the selected Primary Informational Rank and/or Extended
#' Signal Tail boundary of an existing result. Recalculates spectral-volume and
#' informational-gravity allocations and the primary-axis loading table without
#' recomputing mutual information or the eigendecomposition.
#'
#' @param scree_obj A valid output list from \code{Entropic.Scree()}.
#' @param new_K_roots Integer. New Primary Informational Rank. NULL retains the current selection.
#' @param new_K_extended Integer. New inclusive Extended Signal Tail boundary. NULL retains the current selection. The resulting pair must satisfy 1 <= K_roots <= K_extended <= m; invalid pairs are rejected, not silently clamped.
#' @param bipolar_top_n NULL or numeric. NULL (default) keeps every nonzero loading per pole. A fraction strictly between zero and one limits each pole to that proportion of the retained measurement count; a positive integer sets a per-pole count. The complete measurement_loadings matrix is never truncated by this option.
#' @param interactive_mode Logical. If TRUE, prints updated metric dashboards to the console.
#'
#' @details The returned primary-axis loading matrix uses the same saved full
#' eigenvectors. All retained measurements are included on each selected axis.
#' If the input result has no saved eigenvectors, loading outputs remain NULL;
#' this function does not silently repeat the eigendecomposition. Pole lists
#' use Axis_1, Axis_2, etc., matching the primary loading-matrix columns.
#' See \code{\link{Entropic.Scree}} for definitions and interpretation.
#' Both explicitly named auxiliary allocations refresh using the
#' current ranks. Their explicit inclusive cutoffs and mean raw-eigenvalue
#' reference are recorded in auxiliary_diagnostics. Macro metadata are
#' calculated from the saved spectrum without an eigensolve. Geometry status
#' is unavailable when the input has no raw eigenvalues. Preparation receipts
#' and any saved original-unit analysis data remain unchanged.
#'
#' @section Output structure:
#' Axis-specific pole lists use \code{Axis_1}, \code{Axis_2}, and so on,
#' matching the columns of \code{measurement_loadings}. Each pole list contains
#' \code{Positive_Anchor} and \code{Negative_Anchor} coefficients.
#' \code{NSMR} is a ratio; \code{NSMR_percent} is 100 times that ratio.
#'
#' Within \code{auxiliary_diagnostics}, \code{macro_boundary_index} identifies
#' the first axis after a detected macro gap (NA when no gap is detected).
#' \code{macro_endpoint} is the inclusive endpoint, and \code{macro_endpoint_source}
#' identifies whether it comes from a detected gap or the mean-reference count.
#' \code{bulk_allocation_cutoff} is the inclusive auxiliary allocation endpoint.
#' \code{Update.Entropic.Scree()} returns a result without modifying its input
#' object in place.
#'
#' @return The updated result list, returned invisibly, including the selected
#' boundaries, primary measurement_loadings, complete pole lists by default,
#' spectral-volume allocations, informational gravity, and Structural Topology
#' Profile. The saved full eigenvectors and eigenvalues are unchanged.
#'
#' @export
Update.Entropic.Scree <- function(scree_obj, new_K_roots = NULL, new_K_extended = NULL, bipolar_top_n = NULL, interactive_mode = TRUE) {

  # 1. Input Validation
  if (!is.list(scree_obj) || is.null(scree_obj$eigenvalues)) {
    stop("Error: 'scree_obj' must be a valid output list from Entropic.Scree().")
  }

  eig_vals <- scree_obj$eigenvalues
  m_valid <- length(scree_obj$retained_features)

  .resolve_bipolar_top_n(bipolar_top_n, m_valid)

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

  # Validate the same ordered inclusive boundaries used by the main function.
  K_roots <- if (is.null(new_K_roots)) scree_obj$K_roots else new_K_roots
  K_extended <- if (is.null(new_K_extended)) scree_obj$K_extended else new_K_extended
  .validate_ranks(K_roots, K_extended, m_valid)
  K_roots <- as.integer(K_roots)
  K_extended <- as.integer(K_extended)
  positive_info <- .check_selected_geometry(scree_obj$raw_eigenvalues, K_roots)

  # 2. Extract Saved Spectral Quantities from Object
  R_eff <- scree_obj$R_eff
  S_epsilon <- sum(eig_vals)
  redundant_signal_volume <- scree_obj$redundant_signal_volume
  vectors_out <- scree_obj$eigenvectors
  valid_vars <- scree_obj$retained_features
  Tr_Mc <- sum(diag(scree_obj$similarity_matrix))

  # ============================================================================
  # 3. RECALCULATE GRAVITY & STRUCTURAL COMPOSITION
  # ============================================================================

  # Calculate Total Shared Signal Volume through the selected Extended Signal Tail endpoint
  selected_spectral_weight <- sum(eig_vals[1:K_extended])
  P_sig <- selected_spectral_weight / S_epsilon
  unique_signal_volume <- R_eff * P_sig
  total_signal_volume <- unique_signal_volume + redundant_signal_volume
  new_pct_signal <- (total_signal_volume / m_valid) * 100
  new_pct_unique <- (unique_signal_volume / m_valid) * 100
  new_pct_redundant <- (redundant_signal_volume / m_valid) * 100

  # Allocate AIG and FSIG across the selected primary informational axes
  AIG <- total_signal_volume / K_roots
  core_eigenvals <- eig_vals[1:K_roots]
  p_core <- core_eigenvals / sum(core_eigenvals)
  FSIG_final <- p_core * total_signal_volume
  STP_final <- FSIG_final / FSIG_final[1]

  # Idiosyncratic Informational Volume
  trailing_spectral_weight <- if (K_extended < m_valid) sum(eig_vals[(K_extended + 1):m_valid]) else 0
  idiosyncratic_weight <- trailing_spectral_weight / S_epsilon
  idiosyncratic_volume <- R_eff * idiosyncratic_weight
  new_pct_idio <- (idiosyncratic_volume / m_valid) * 100

  # ============================================================================
  # 4. RE-EXTRACT BIPOLAR MODULES
  # ============================================================================
  bipolar_modules_out <- NULL
  if (!is.null(vectors_out) && K_roots > 0) {
    bipolar_modules_out <- list()

    resolved_top_n <- .resolve_bipolar_top_n(bipolar_top_n, m_valid)

    for (k in 1:K_roots) {
      loadings <- setNames(vectors_out[, k], valid_vars)

      pos_pole <- sort(loadings[loadings > 0], decreasing = TRUE)
      pos_pole <- head(pos_pole, resolved_top_n)

      neg_pole <- sort(loadings[loadings < 0], decreasing = FALSE)
      neg_pole <- head(neg_pole, resolved_top_n)

      bipolar_modules_out[[paste0("Axis_", k)]] <- list(
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
    message(sprintf(" -> Primary Informational Rank (K_roots)     : %d", K_roots))
    message(sprintf(" -> Extended Signal Tail boundary (K_extended) : %d", K_extended))
    message("-----------------------------------------------------------------")
    message(" [Informational Gravity Metrics Estimates] (K_extended Volume Rebundled to K_roots)")
    message(sprintf(" -> Average Informational Gravity (AIG): %.3f Variable Equivalents per Primary Axis", AIG))
    message(" -> Factor-Specific Informational Gravity (FSIG):")
    message(paste(capture.output(print(round(FSIG_final, 3))), collapse = "\n"))
    message(" -> Structural Topology Profile (Relative to FSIG_1):")
    message(paste(capture.output(print(round(STP_final, 3))), collapse = "\n"))
    message("===================================================================================\n")

    message("===================================================================================")
    message(sprintf(" (UPDATED) TRIPARTITE STRUCTURAL COMPOSITION (Rebundling Based on K_extended = %d)", K_extended))
    message("===================================================================================")
    message(sprintf(" -> %-50s : %d", "Retained Measurements (m)", m_valid))
    message(sprintf(" -> %-50s : %.2f", "Centered Trace (Tr_Mc)", Tr_Mc))
    message("-----------------------------------------------------------------")
    message(sprintf(" -> %-50s : %.3f", "Total Shared Signal Volume (Unique + Redundant)", total_signal_volume))
    message(sprintf(" -> %-40s : %.3f", "   (Unique Signal Volume)", unique_signal_volume))
    message(sprintf(" -> %-40s : %.3f", "   (Redundant Signal Volume)", redundant_signal_volume))
    message("-----------------------------")
    message(sprintf(" -> %-50s : %.3f", "Idiosyncratic Informational Volume", idiosyncratic_volume))
    message("      (Residual spectral allocation; not a measured error fraction)")
    message("-----------------------------------------------------------------")
    message(sprintf(" -> %-50s : %.3f%%", "% Total Shared Signal", new_pct_signal))
    message(sprintf(" -> %-40s : %.3f%%", "   (% Unique Signal)", new_pct_unique))
    message(sprintf(" -> %-40s : %.3f%%", "   (% Redundant Signal)", new_pct_redundant))
    message("-----------------------------")
    message(sprintf(" -> %-50s : %.3f%%", "% Idiosyncratic Informational Volume", new_pct_idio))
    message("      (Residual spectral allocation; not a measured error fraction)")
    message("===================================================================================\n")

    message("===================================================================================")
    message(" (DELTA) METRIC SHIFT COMPARISON")
    message("===================================================================================")
    message(sprintf(" -> %-37s : %d -> %d", "Primary Informational Rank (K_roots)", orig_K_roots, K_roots))
    message(sprintf(" -> %-37s : %d -> %d", "Extended Signal Tail boundary (K_extended)", orig_K_extended, K_extended))
    message("-----------------------------------------------------------------")
    message(sprintf(" -> %-37s : %.3f -> %.3f", "Average Informational Gravity (AIG)", orig_AIG, AIG))
    message("-----------------------------------------------------------------")
    message(sprintf(" -> %-37s : %.3f -> %.3f", "Total Shared Signal Volume", orig_total_signal, total_signal_volume))
    message(sprintf(" -> %-37s : %.3f -> %.3f", "   (Unique Signal Volume)", orig_unique_signal, unique_signal_volume))
    message(sprintf(" -> %-37s : %.3f -> %.3f", "   (Redundant Signal Volume)", orig_redundant_signal, redundant_signal_volume))
    message("-----------------------------")
    message(sprintf(" -> %-37s : %.3f -> %.3f", "Idiosyncratic Informational Volume", orig_idio_vol, idiosyncratic_volume))
    message("-----------------------------------------------------------------")
    message(sprintf(" -> %-37s : %.3f%% -> %.3f%%", "% Total Shared Signal", orig_pct_signal, new_pct_signal))
    message(sprintf(" -> %-37s : %.3f%% -> %.3f%%", "   (% Unique Signal)", orig_pct_unique, new_pct_unique))
    message(sprintf(" -> %-37s : %.3f%% -> %.3f%%", "   (% Redundant Signal)", orig_pct_redundant, new_pct_redundant))
    message("-----------------------------")
    message(sprintf(" -> %-37s : %.3f%% -> %.3f%%", "% Idiosyncratic Informational Volume", orig_pct_idio, new_pct_idio))
    message("===================================================================================\n")

    message("===================================================================================")
    message(" ENTROPIC SCREE (v1.0.2) - METHODOLOGICAL REFERENCE & LICENSE")
    message("===================================================================================")
    message(" -> Framework developed by Terrence J. Lee-St. John (Enli)")
    message(" -> Released under the Apache License 2.0 (Open Source)")
    message(" -> For full methods and metric definitions, see:")
    message("    The Entropic Scree (2026) - https://doi.org/10.5281/zenodo.22028086")
    message("===================================================================================\n")
  }

  # Calculate auxiliary allocations at the selected ranks and derive macro
  # metadata from the saved spectrum without an eigensolve.
  mean_raw <- Tr_Mc / m_valid
  macro_info <- .macro_endpoint(eig_vals, mean_raw)
  auxiliary <- .auxiliary_gravity(eig_vals, R_eff, K_roots, K_extended,
    mean_raw, macro_info$macro_boundary_index)

  # 6. Mutate and Return Object
  scree_obj$K_roots <- K_roots
  scree_obj$K_extended <- K_extended
  scree_obj$total_signal_volume <- total_signal_volume
  scree_obj$unique_signal_volume <- unique_signal_volume
  scree_obj$idiosyncratic_volume <- idiosyncratic_volume
  scree_obj$AIG <- AIG
  scree_obj$FSIG_final <- FSIG_final
  scree_obj$structural_topology_profile <- STP_final
  scree_obj$FSIG_auxiliary_bulk <- auxiliary$bulk
  scree_obj$FSIG_auxiliary_mean_reference <- auxiliary$mean_reference
  scree_obj$auxiliary_diagnostics <- auxiliary$diagnostics
  scree_obj$positive_axis_diagnostics <- positive_info
  if (is.null(scree_obj$rank_diagnostics)) scree_obj$rank_diagnostics <- list()
  scree_obj$rank_diagnostics$selection <- "posthoc_user_selection"
  scree_obj["measurement_loadings"] <- list(.primary_measurement_loadings(vectors_out, valid_vars, K_roots))
  scree_obj$bipolar_modules <- bipolar_modules_out

  # Express negative spectral mass in ratio and percentage units, and keep
  # only the documented output fields.
  if (is.null(scree_obj[["NSMR"]]) && !is.null(scree_obj[["SCDR"]]))
    scree_obj$NSMR <- scree_obj[["SCDR"]] / 100
  scree_obj$NSMR_percent <- if (is.null(scree_obj[["NSMR"]])) NA_real_ else
    100 * scree_obj[["NSMR"]]
  scree_obj[c("SCDR", "FSIG_extended_bulk", "FSIG_extended_kaiser", "top_of_bulk")] <- NULL

  return(invisible(scree_obj))
}

# Internal helpers for loading matrices and sign-specific pole lists.
.resolve_bipolar_top_n <- function(top_n, measurement_count) {
  if (is.null(top_n)) return(measurement_count)
  if (!is.numeric(top_n) || length(top_n) != 1L || is.na(top_n) ||
      !is.finite(top_n) || top_n <= 0 || (top_n >= 1 && top_n != floor(top_n))) {
    stop("'bipolar_top_n' must be NULL, a proportion between 0 and 1, or a positive integer.")
  }
  if (top_n < 1) return(max(1, floor(top_n * measurement_count)))
  min(measurement_count, top_n)
}

.primary_measurement_loadings <- function(vectors, measurements, primary_rank) {
  if (is.null(vectors)) return(NULL)
  selected <- vectors[, seq_len(primary_rank), drop = FALSE]
  dimnames(selected) <- list(measurements, paste0("Axis_", seq_len(primary_rank)))
  selected
}
