#' Forests for generative modeling for high-dimensional omics data.
#'
#' Uses a pre-trained harf model to generate synthetic omics data.
#'
#' @param harf_obj A pre-trained harf model.
#' @param n_synth Number of synthetic samples to generate.
#' @param evidence Optional set of conditioning events. This will be further passed
#' to the \code{forde} function in each isolated regions. See \code{forde} for details.
#' @param omx_onset_data Optional data.frame of conditional onset omics features.
#' @param evidence_row_mode Interpretation of rows in multi-row evidence.
#' See \code{forde} for details.
#' @param round Round continuous variables to their respective maximum precision
#' in the real data set? See \code{forde} for details.
#' @param sample_NAs Sample NAs respecting the probability for missing values in
#'  the original data? See \code{forde} for details.
#' @param nomatch What to do if no leaf matches a condition in evidence?
#' Options are to force sampling from a random leaf ("force") or return NA ("na").
#'  The default is "force".
#' @param verbose What to do if no leaf matches a condition in \code{evidence}?
#' See \code{forde} for details.
#' @param stepsize How many rows of evidence should be handled at each step?
#' See \code{forde} for details.
#' @param parallel Compute in parallel? See \code{forde} for details.
#'
#' @returns A data.table containing the generated synthetic omics data.
#' @importFrom foreach %do% %dopar% foreach
#' @export
#'
# TODO: Add examples
h_forge <- function (
    harf_obj,
    n_synth,
    evidence = NULL,
    omx_onset_data = NULL,
    evidence_row_mode = c("separate", "or"),
    round = TRUE,
    sample_NAs = FALSE,
    nomatch = c("force", "na"),
    verbose = TRUE,
    stepsize = 0,
    parallel = TRUE
) {
  i <- NULL  # To avoid R CMD check note for foreach
  # Test that harf_obj is of class harf
  if (!inherits(harf_obj, "harf")) {
    stop("harf_obj must be an object of class 'isoARF'.")
  }
  # n_synth must be an integer greater than 0
  if (!is.numeric(n_synth) | n_synth <= 0 | n_synth != round(n_synth)) {
    stop("n_synth must be a positive integer.")
  }
  # Extract those evidence columns that are in clinical/lab meta data
  if (!is.null(evidence)) {
    if (any(colnames(evidence) %in% c(unique(harf_obj$meta_model$meta_forde$cnt$variable),
                                      unique(harf_obj$meta_model$meta_forde$cat$variable)))) {
      meta_evidence <- evidence[ , colnames(evidence) %in% c(unique(harf_obj$meta_model$meta_forde$cnt$variable),
                                                             unique(harf_obj$meta_model$meta_forde$cat$variable)), drop = FALSE]
    } else {
      meta_evidence <- NULL
    }
  } else {
    meta_evidence <- NULL
  }
  # Retrieve meta_model from isoARF
  meta_model <- harf_obj$meta_model
  # Retrieve the omx_onset_pcs from harf_obj if omx_onset_data is provided
  if (!is.null(omx_onset_data)) {
    # Stop if evidence does not have the same number of rows as omx_onset_data
    if (!is.null(evidence)) {
      if (nrow(evidence) != nrow(omx_onset_data)) {
        stop("Number of rows in evidence must match number of rows in omx_onset_data.")
      }
    }
    onset_loadings <- meta_model$onset_loadings
    if (!is.null(onset_loadings)) {
      # Projection of onset features onto pcs
      if (ncol(omx_onset_data) != nrow(onset_loadings)) {
        stop("Number of columns in omx_onset_data must match number of rows in onset_loadings.")
      }
      onset_features <- as.data.frame(
        crossprod(t(as.matrix(omx_onset_data)), as.matrix(onset_loadings))
      )
      # Update meta_evidence
      if (is.null(meta_evidence)) {
        meta_evidence <- onset_features
      } else {
        meta_evidence <- data.frame(meta_evidence,
                                    onset_features)
      }
    } else {
      stop("onset_loadings not found in harf_obj object.")
    }
  }
  # Forge meta data
  meta_data <- arf::forge(
    params = meta_model$meta_forde,
    n_synth = n_synth,
    evidence = meta_evidence,
    evidence_row_mode = evidence_row_mode,
    round = round,
    sample_NAs = sample_NAs,
    nomatch = nomatch,
    verbose = verbose,
    stepsize = stepsize,
    parallel = parallel
  )
  # For the remaining clusters, generate synthetic data conditioned on the current evidence
  synth_omx <- function (mdl_idx) {
    if (verbose) {
      message(paste0("Generating synthetic data for cluster ", mdl_idx, " out of ", length(harf_obj$models), "..."))
    }
    synth_data <- arf::forge(
      params = harf_obj$models[[mdl_idx]]$iso_forde,
      n_synth = 1, # Always equals to 1 since we include evidence; nrow(evidence) == n_synth.
      evidence = meta_data,
      evidence_row_mode = evidence_row_mode,
      round = round,
      sample_NAs = sample_NAs,
      nomatch = nomatch,
      verbose = verbose,
      stepsize = stepsize,
      parallel = parallel
    )
    # Dataset to be kept
    synth_data <- synth_data[ , harf_obj$models[[mdl_idx]]$ftr_in_cluster, drop = FALSE]
    return(synth_data)
  }
  # If parallel is TRUE, use dopar
  if (parallel) {
    synth_omx_data_list <- foreach(i = 1:length(harf_obj$models)) %dopar% synth_omx(i)
  } else {
    synth_omx_data_list <- foreach(i = 1:length(harf_obj$models)) %do% synth_omx(i)
  }
  # Combine all synthetic omics data
  synth_omx_data <- do.call(cbind, synth_omx_data_list)
  if (!is.null(harf_obj$cli_lab_feature)) {
    synth_data <- data.table::data.table(synth_omx_data,
                             meta_data[ , h_arf$cli_lab_feature,
                                        drop = FALSE])
  } else {
    synth_data <- data.table::data.table(synth_omx_data)
  }
  return(synth_data)
}
