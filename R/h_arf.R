#' Isolated Adversarial Random Forest for Omics Data
#'
#' @param omx_data A data.frame containing omics data where rows
#' represent samples (e.g. patients or cells) and columns represent features
#' (e.g., gene or protein expressions).
#' @param cli_lab_data A data.frame of clinical or laboratory data. For example,
#' this may be a data.frame where rows represent cell types (for single cell data)
#' or additional clinical patient information (e.g. disease status).
#' @param omx_onset_data Optional data.frame of conditional onset omics features.
#' @param num_trees Number of trees to grow in each Adversarial Random Forest (ARF) model.
#' @param min_node_size Minimum number of samples required to split an internal node in the ARF model.
#' @param encoder Method to use for encoding. Currently,
#' only \code{fastcor} and Random Forest Auto-Encoder (\code{RFAE}) supported. Default is "fastcor".
#' @param correlation_method Methods to compute correlation between features. Only useful for the \code{fastcor} encoder.
#' Options include "pearson", "spearman", and "kendall". Default is "spearman".
#' @param num_clusters Number of clusters to form. If NULL, the optimal number
#' of clusters is found by finding the elbow method.
#' @param num_btwn_pcs Number of principal components to use for between cluster variability.
#' Usually set to 1 or 2. Default is 2.
#' @param num_onset_pcs Number of principal components to use for onset features.
#' Default is 2.
#' @param chunck_size Size of feature chuncks to process at a time. Default is 10.
#' @param oob Logical indicating whether to use out-of-bag samples for density estimation.
#' Default is FALSE. Also see \code{forde}
#' @param family Distribution to use for density estimation of continuous features. See \code{forde} for options.
#' @param finite_bounds Impose finite bounds on all continuous variables?
#' If "local", infinite bounds are set to empirical extrema within leaves.
#' If "global", infinite bounds are set to global empirical extrema. if "no"
#' (the default), infinite bounds are left unchanged. See \code{forde} for more details.
#' @param alpha Optional pseudocount for Laplace smoothing in density estimation. See \code{forde} for more details.
#' @param epsilon Optional slack parameter on empirical bounds. See \code{forde} for more details.
#' @param parallel Logical indicating whether to use parallel processing. Default is TRUE.
#' @param verbose Logical indicating whether to print progress messages. Default is FALSE.
#'
#' @returns An isoARF object containing the fitted adversarial models, and clustering
#' information.
#' @importFrom cluster pam
#' @importFrom rsvd rpca
#' @importFrom foreach %do% %dopar% foreach
#' @importFrom fastPLS fastcor
#' @importFrom gmodels fast.prcomp
#' @importFrom stats cor dist prcomp
#' @export
# TODO: Add examples

h_arf <- function (
    omx_data,
    cli_lab_data = NULL,
    omx_onset_data = NULL,
    num_trees = 20,
    min_node_size = 5,
    encoder = "fastcor",
    correlation_method = "spearman",
    num_clusters = NULL,
    num_btwn_pcs = 2,
    num_onset_pcs = 2,
    chunck_size = 10,
    oob = FALSE,
    family = "truncnorm",
    finite_bounds = "no",
    alpha = 0.05,
    epsilon = 0.05,
    parallel = TRUE,
    verbose = FALSE
) {
  # Stop if omx_data is not a data.frame or matrix
  if (!is.data.frame(omx_data)) {
    stop("omx_data must be a data.frame.")
  }
  cln_omx_data <- colnames(omx_data)
  if (!is.null(omx_onset_data) & !is.data.frame(omx_onset_data)) {
    stop("omx_onset_data must be a data.frame.")
  }
  if (!is.null(omx_onset_data)) {
    if (nrow(omx_data) != nrow(omx_onset_data)) {
      stop("Number of rows in omx_data must match number of rows in omx_onset_data.")
    }
  }
  if (!is.null(omx_onset_data)) {
    if (!all(colnames(omx_onset_data) %in% colnames(omx_data))) {
      stop("All column names in omx_onset_data must be present in omx_data.")
    }
  }
  if (!is.null(cli_lab_data) & !is.data.frame(cli_lab_data)) {
    stop("cli_lab_data must be a data.frame.")
  }
  if (!is.null(cli_lab_data)) {
    cln_cli_lab_data <- colnames(cli_lab_data)
  } else {
    cln_cli_lab_data <- NULL
  }
  if (nrow(omx_data) != nrow(cli_lab_data)) {
    stop("Number of rows in omx_data must match number of rows in cli_lab_data.")
  }
  if (!encoder %in% c("fastcor", "RFAE")) {
    stop("encoder must be one of 'fastcor' or 'RFAE'.")
  }
  if (!correlation_method %in% c("pearson", "spearman", "kendall")) {
    stop("correlation_method must be one of 'pearson', 'spearman', or 'kendall'.")
  }
  if (!is.numeric(chunck_size) | chunck_size <= 0 | chunck_size != round(chunck_size)) {
    stop("chunck_size must be a positive integer.")
  }
  message("Clustering features...\n")
  # Encoding via fast correlation and PCA
  if (encoder == "fastcor") {
    if (correlation_method == "spearman") {
      cor_matrix <- fastcor(as.matrix(omx_data))
    } else {
      # TODO: Fix me by using a fast correlation function that supports pearson and kendall
      cor_matrix <- cor(as.matrix(omx_data), method = correlation_method)
    }

    # PCA and elbow point detection to determine optimal number of clusters
    pca_res <- fast.prcomp(cor_matrix, retx = TRUE, center = TRUE, scale. = TRUE)
    variances <- pca_res$sdev^2
    k_pcs <- find_elbow(variances)
    selected_pcs <- pca_res$rotation[, 1:k_pcs, drop = FALSE]
    # Projected omx_data onto pcs with crossprod
    projected_data <- as.data.frame(crossprod(as.matrix(omx_data), as.matrix(selected_pcs)))
  }
  # Encoding via RFAE
  if (encoder == "RFAE") {
    # TODO: Encoding via RFAE; review this with high-dimensional version
    # omx_rfae <- RFAE::encode(rf = adversarial_rf(x = omx_data,
    #                                              num_trees = num_trees,
    #                                              min_node_size = min_node_size,
    #                                              verbose = FALSE),
    #                          x = omx_data,
    #                          k = num_btwn_pcs
    # )
    # projected_data <- as.data.frame(omx_rfae$Z)
    stop("RFAE encoder not yet implemented.")
  }
  # kmeans clustering around medoids
  dist_features <- dist(projected_data)
  if (!is.null(num_clusters)) {
    k_pcs <- num_clusters
  }
  pam_fit <- pam(dist_features, k = k_pcs, diss = TRUE)
  feature_clusters <- pam_fit$clustering
  # Re-split clusters that are larger than chunck_size
  if (!is.null(chunck_size)) {
    max_cluster_size <- chunck_size
    new_cluster_id <- max(feature_clusters) + 1
    for (cluster in unique(feature_clusters)) {
      ftr_in_cluster <- which(feature_clusters == cluster)
      if (length(ftr_in_cluster) > max_cluster_size) {
        num_subclusters <- ceiling(length(ftr_in_cluster) / max_cluster_size)
        sub_pam_fit <- pam(
          projected_data[ftr_in_cluster, ],
          k = num_subclusters,
          diss = FALSE
        )
        sub_clusters <- sub_pam_fit$clustering
        for (sub_cluster in unique(sub_clusters)) {
          ftr_in_subcluster <- ftr_in_cluster[sub_clusters == sub_cluster]
          feature_clusters[ftr_in_subcluster] <- new_cluster_id
          new_cluster_id <- new_cluster_id + 1
        }
      }
    }
  }
  # Compute meta cluster features for the joint structure between cluster.
  meta_features <- sapply(
    unique(feature_clusters),
    function(cluster) {
      cluster_data <- omx_data[ , which(feature_clusters == cluster), drop = FALSE]
      if (ncol(cluster_data) < 2) {
        return(as.matrix(cluster_data))
      }

      # fastcor encoder
      if (encoder == "fastcor") {
        # PCA of cluster features
        pca_meta <- rpca(t(as.matrix(cluster_data)),
                         k = num_btwn_pcs,
                         center = TRUE,
                         scale = TRUE)
        # Projection of cluster features onto pcs
        return(crossprod(t(as.matrix(cluster_data)), as.matrix(pca_meta$x)))
      }

      # Random forest encoder
      if (encoder == "RFAE") {
        cluster_rfae <- do.call(
          what = "RFAE::encode",
          args = list(
            rf = arf::adversarial_rf(
              x = cluster_data,
              num_trees = num_trees,
              min_node_size = min_node_size,
              verbose = FALSE
            ),
            x = cluster_data,
            k = num_btwn_pcs
          )
        )
        return(as.matrix(cluster_rfae$Z))
      }
    }
  )
  meta_features <- as.data.frame(meta_features)
  colnames(meta_features) <- sprintf("cluster_%s",
                                          ncol(meta_features))
  # Dimension reduction of meta features if num_btwn_pcs < ncol(meta_features)
  if (num_btwn_pcs < ncol(meta_features)) {
    pca_meta <- rpca(as.matrix(t(meta_features)),
                     k = num_btwn_pcs,
                     center = TRUE,
                     scale = TRUE)
    # Projection of meta features onto pcs
    meta_features <- as.data.frame(
      crossprod(t(as.matrix(meta_features)), as.matrix(pca_meta$x))
    )
    colnames(meta_features) <- sprintf("meta_pc_%s",
                                            1:num_btwn_pcs)
  }
  # Dimension reduction of omx_onset_data if provided
  onset_loadings <- NULL
  if (!is.null(omx_onset_data)) {
    if (ncol(omx_onset_data) < 2) {
      stop("For omx_onset_data, at least 2 features are required.")
    }
    pca_onset <- rpca(as.matrix(t(omx_onset_data)),
                      k = num_onset_pcs,
                      center = TRUE,
                      scale = TRUE)
    onset_loadings <- pca_onset$x
    # Projection of onset features onto pcs
    onset_features <- as.data.frame(
      crossprod(t(as.matrix(omx_onset_data)), as.matrix(onset_loadings))
    )
    colnames(onset_features) <- sprintf("onset_pc_%s",
                                        1:ncol(onset_features))
    # Add onset features to meta features
    meta_features <- data.frame(meta_features,
                                     onset_features)
  }
  # Add clinical/laboratory data if provided
  if (!is.null(cli_lab_data)) {
    meta_features <- data.frame(meta_features, cli_lab_data)
  }
  # Run the meta adversarial game
  if (isTRUE(verbose)) {
    message("Fitting meta-cluster ARF model...\n")
  }
  meta_arf <- arf::adversarial_rf(
    x = meta_features,
    num_trees = num_trees,
    min_node_size = min_node_size,
    verbose = verbose
  )
  # Density estimation for meta model
  meta_forde <- arf::forde(
    arf = meta_arf,
    x = meta_features,
    oob = oob,
    family = family,
    finite_bounds = finite_bounds,
    alpha = alpha,
    epsilon = epsilon,
    parallel = TRUE # TODO: Set parallel to parallel argument later.
  )

  meta_model <- list(meta_arf = meta_arf,
                     meta_forde = meta_forde,
                     onset_loadings = onset_loadings)
  # For each cluster fit ARF model
  arf_clusters <- function(cluster) {
    if (isTRUE(verbose)) {
      message(paste0("Fitting ARF model for cluster ",
                     cluster, " of ", length(unique(feature_clusters)),
                     " clusters...\n"))
    }
    ftr_in_cluster <- which(feature_clusters == cluster)
    ftr_data_subset <- omx_data[, ftr_in_cluster, drop = FALSE]
    other_clusters <- unique(feature_clusters[feature_clusters != cluster])
    ftr_data_subset <- data.frame(ftr_data_subset,
                                  meta_features)
    # Adversarial game in isolated regions
    iso_arf <- arf::adversarial_rf(
      x = ftr_data_subset,
      num_trees = num_trees,
      min_node_size = min_node_size,
      verbose = verbose
    )

    # Density estimation
    iso_forde <- arf::forde(
      arf = iso_arf,
      x = ftr_data_subset,
      oob = oob,
      family = family,
      finite_bounds = finite_bounds,
      alpha = alpha,
      epsilon = epsilon,
      parallel = TRUE # TODO: Set parallel to parallel argument later.
    )
    # Export the iso_arf object
    return(list(cluster_id = cluster,
                ftr_in_cluster = names(ftr_in_cluster),
                iso_arf = iso_arf,
                iso_forde = iso_forde))
  }
  if (isFALSE(parallel)) {
    # Use for loop if parallel is FALSE
    arf_models <- foreach(cluster = sort(unique(feature_clusters))) %do% arf_clusters(cluster)
  } else {
    # Use dopar if parallel is TRUE
    arf_models <- foreach(cluster = sort(unique(feature_clusters))) %dopar% arf_clusters(cluster)
  }
  # Export feature_cluster_df and meta to arf_models as well.
  hd_arf <- list(meta_model = meta_model,
                  models = arf_models,
                  cluster = data.frame(
                    feature = colnames(omx_data),
                    cluster = feature_clusters
                  ),
                  meta_features = meta_features,
                  omx_features = cln_omx_data,
                  cli_lab_features = cln_cli_lab_data
  )
  class(hd_arf) <- "harf"
  return(hd_arf)
}
