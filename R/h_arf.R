#' High-dimensional Adversarial Random Forest
#'
#' The algorithm partitions high-dimensional data into isolated regions and fits
#' Adversarial Random Forest (ARF) models within each region to capture local
#' dependencies.
#'
#' @param omx_data A data.frame containing omics data where rows
#' represent samples (e.g. patients or cells) and columns represent features
#' (e.g., gene or protein expressions).
#' @param cli_lab_data A data.frame of clinical or laboratory data. For example,
#' this may be a data.frame where rows represent cell types (for single cell data)
#' or additional clinical patient information (e.g. disease status).
#' @param target Optional name of the target variable in \code{cli_lab_data} for supervised
#' dimension reduction. If NULL, unsupervised dimension reduction is performed.
#' @param omx_onset_data Optional data.frame of conditional onset omics features.
#' @param num_trees Number of trees to grow in each Adversarial Random Forest (ARF) model.
#' @param min_node_size Minimum number of samples required to split an internal node in the ARF model.
#' @param correlation_method Methods to compute correlation between features.
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
#' @param always_split_meta Logical indicating whether to always split on meta features. Default is FALSE.
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
#' @importFrom ClusterR KMeans_rcpp
#' @importFrom fastPLS fastcor
#' @importFrom rsvd rpca
#' @importFrom stats cor dist prcomp
#' @importFrom ranger ranger
#' @importFrom RGCCA rgcca
#' @importFrom mixOmics block.spls block.splsda
#' @export
#' @seealso [h_forge], [arf::adversarial_rf], [arf::forde]
#' @examples
#' \dontrun{
#' data(single_cell)
#' harf_model <- h_arf(
#'   omx_data = single_cell[ , - which(colnames(single_cell)  == "cell_type")],
#'   cli_lab_data = data.frame(cell_type = single_cell$cell_type),
#'   parallel = TRUE,
#'   verbose = FALSE
#' )
#' # Unconditional sampling from harf_model
#' set.seed(123)
#' synth_single_cell <- h_forge(
#'  harf_obj = harf_model,
#'  n_synth = nrow(single_cell),
#'  evidence = NULL,
#'  parallel = TRUE,
#'  verbose = FALSE
#'  )
#'  # Conditional resampling from harf_model
#'  set.seed(142)
#'  lung_single_cell <- h_forge(
#'      harf_obj = harf_model,
#'      n_synth = sum(single_cell$cell_type == "lung"),
#'      evidence = data.frame(cell_type = "lung"),
#'      verbose = FALSE,
#'      parallel = TRUE
#'     )
#' }

h_arf <- function (
    omx_data,
    cli_lab_data = NULL,
    target = NULL,
    omx_onset_data = NULL,
    num_trees = 10,
    min_node_size = 5,
    # encoder = "fastcor",
    correlation_method = "spearman",
    num_clusters = NULL,
    num_btwn_pcs = 2,
    num_onset_pcs = 2,
    chunck_size = 10,
    oob = FALSE,
    always_split_meta = FALSE,
    family = "truncnorm",
    finite_bounds = "no",
    alpha = 0,
    epsilon = 0.0,
    parallel = TRUE,
    verbose = FALSE
) {
  # Stop if omx_data is not a data.frame or data.table
  if (!inherits(omx_data, "data.frame")) {
    stop("omx_data must be a data.frame.")
  } else {
    omx_data <- as.data.frame(omx_data)
  }
  cln_omx_data <- colnames(omx_data)
  if (!is.null(omx_onset_data)) {
    if (!is.data.frame(omx_onset_data)) {
      stop("omx_onset_data must be a data.frame.")
    }
    omx_onset_data <- as.data.frame(omx_onset_data)
    if (nrow(omx_data) != nrow(omx_onset_data)) {
      stop("Number of rows in omx_data must match number of rows in omx_onset_data.")
    }
  }
  if (!is.null(omx_onset_data)) {
    if (!all(colnames(omx_onset_data) %in% colnames(omx_data))) {
      stop("All column names in omx_onset_data must be present in omx_data.")
    }
  }
  if (!is.null(cli_lab_data)) {
    if(!inherits(cli_lab_data, "data.frame")) {
      stop("cli_lab_data must be a data.frame.")
    }
  } else {
    cli_lab_data <- as.data.frame(cli_lab_data)
  }
  if (!is.null(cli_lab_data)) {
    cln_cli_lab_data <- colnames(cli_lab_data)
  } else {
    cln_cli_lab_data <- NULL
  }
  if (nrow(omx_data) != nrow(cli_lab_data)) {
    stop("Number of rows in omx_data must match number of rows in cli_lab_data.")
  }
  # if (!encoder %in% c("fastcor", "RFAE")) {
  #   stop("encoder must be one of 'fastcor' or 'RFAE'.")
  # }
  if (!correlation_method %in% c("pearson", "spearman", "kendall")) {
    stop("correlation_method must be one of 'pearson', 'spearman', or 'kendall'.")
  }
  if (!is.numeric(chunck_size) | chunck_size <= 0 | chunck_size != round(chunck_size)) {
    stop("chunck_size must be a positive integer.")
  }
  if (isTRUE(verbose)) {
    message("Clustering features...\n")
  }
  # Encoding via fast correlation and PCA
  if (correlation_method == "spearman") {
    cor_matrix <- fastcor(t(as.matrix(omx_data)))
  } else {
    # TODO: Optmize me by using a fast correlation function that supports pearson and kendall
    cor_matrix <- cor(as.matrix(omx_data), method = correlation_method)
  }
  # PCA with encoded data
  projected_data <- fast.pca(
    X = cor_matrix,
    K = min(100, ncol(cor_matrix) - 1)
  )
  # kmeans clustering around medoids
  kmeanpp_fit <- KMeans_rcpp(projected_data,
                             clusters = floor(ncol(omx_data) / chunck_size),
                             num_init = 5,
                             max_iters = 100,
                             initializer = 'kmeans++')
  # feature_clusters <- clara_fit$clustering
  feature_clusters <- kmeanpp_fit$clusters
  # Re-split clusters that are larger than chunck_size
  # if (!is.null(chunck_size)) {
  cluster_splitting <- TRUE
  while (cluster_splitting) {
    max_cluster_size <- chunck_size
    new_cluster_id <- max(feature_clusters) + 1
    cluster_splitting <- FALSE
    for (cluster in unique(feature_clusters)) {
      ftr_in_cluster <- which(feature_clusters == cluster)
      if (length(ftr_in_cluster) > max_cluster_size + 1) {
        cluster_splitting <- TRUE
        num_subclusters <- ceiling(length(ftr_in_cluster) / max_cluster_size)

        kmeanpp_fit <- KMeans_rcpp(projected_data[ftr_in_cluster, ],
                                   clusters = num_subclusters,
                                   num_init = 5,
                                   max_iters = 100,
                                   initializer = 'kmeans++')

        sub_clusters <- kmeanpp_fit$clusters + max(feature_clusters)
        feature_clusters[ftr_in_cluster] <- sub_clusters
      }
    }
  }
  # Merge clusters with less than num_btwn_pcs to nearest cluster
  for (cluster in unique(feature_clusters)) {
    ftr_in_cluster <- which(feature_clusters == cluster)
    if (length(ftr_in_cluster) <= num_btwn_pcs) {
      other_clusters <- unique(feature_clusters[feature_clusters != cluster])
      cluster_center <- colMeans(projected_data[ftr_in_cluster, , drop = FALSE])
      dists <- sapply(
        other_clusters,
        function(other_cluster) {
          ftr_in_other_cluster <- which(feature_clusters == other_cluster)
          other_cluster_center <- colMeans(projected_data[ftr_in_other_cluster, , drop = FALSE])
          return(sqrt(sum((cluster_center - other_cluster_center)^2)))
        }
      )
      nearest_cluster <- other_clusters[which.min(dists)]
      feature_clusters[ftr_in_cluster] <- nearest_cluster

    }
  }
  # Prepare omix regions as list of matrices for CCA.
  meta_features <- lapply(
    unique(feature_clusters),
    function(cluster) {
      cluster_data <- omx_data[ , which(feature_clusters == cluster),
                                drop = FALSE]
      return(as.matrix(cluster_data))
    }
  )
  names(meta_features) <- sprintf("cluster_%s",
                                  unique(feature_clusters))
  # Prepare clinical region RFAE of clinical/laboratory data if provided
  if (!is.null(cli_lab_data) & (ncol(cli_lab_data) > num_btwn_pcs)) {
    if (isTRUE(verbose)) {
      message("Encoding clinical/laboratory data via RFAE...\n")
    }
    if (is.null(target)) {
      # Unsupervised RFAE
      cli_lab_rfae <- do.call(
        what = "RFAE::encode",
        args = list(
          rf = arf::adversarial_rf(
            x = cli_lab_data,
            num_trees = num_trees,
            min_node_size = min_node_size,
            verbose = FALSE
          ),
          x = cli_lab_data,
          k = num_btwn_pcs
        )
      )
    } else {
      # Supervised RFAE
      cli_lab_rf <- ranger::ranger(
        data = cli_lab_data,
        dependent.variable.name = target,
        num_trees = num_trees,
        min_node_size = min_node_size,
        verbose = verbose
      )
      cli_lab_rfae <- do.call(
        what = "RFAE::encode",
        args = list(
          rf = cli_lab_rf,
          x = cli_lab_data,
          k = num_btwn_pcs
        )
      )
    }
    cli_lab_features <- as.matrix(cli_lab_rfae$Z)
    meta_features <- c(meta_features,
                       list(cli_lab = cli_lab_features))
  }
  # (Un)supervised CCA on all blocks to get meta features.
  if (is.null(target)) {
    # Unsupervised RGCCA to capture between-cluster variability
    Z <- RGCCA::rgcca(
      blocks = meta_features,
      ncomp = num_btwn_pcs,
      matrix(1, length(meta_features),
             length(meta_features)) - diag(length(meta_features)),
      scheme = "centroid"
    )
    meta_features <- as.data.frame(do.call(cbind, Z$Y))
  } else {
    if (is.null(cli_lab_data) | !(target %in% colnames(cli_lab_data))) {
      stop("target must be a column name in cli_lab_data.")
    }
    if (is.numeric(cli_lab_data[[target]])) {
      target_type <- "numeric"
    } else {
      target_type <- "categorical"
    }
    # Supervised CCA depending on target type
    if (target_type == "categorical") {
      splsda_fit <- mixOmics::block.splsda(
        X = meta_features,
        Y = cli_lab_data[[target]],
        ncomp = num_btwn_pcs,
        keepX = lapply(meta_features,
                       function(x) {
                         rep(ncol(x), num_btwn_pcs)
                       })
      )
      meta_features <- as.data.frame(do.call(cbind, splsda_fit$variates))
    } else {
      spls_fit <- mixOmics::block.spls(
        X = meta_features,
        Y = cli_lab_data[[target]],
        ncomp = num_btwn_pcs,
        keepX = lapply(meta_features,
                       function(x) {
                         rep(ncol(x), num_btwn_pcs)
                       }),
        keepY = rep(ncol(cli_lab_data[[target]]), num_btwn_pcs)
      )
      meta_features <- as.data.frame(do.call(cbind, spls_fit$variates))
    }
  }
  # Rename meta features as cluster_i.1 and cluster_i.2 etc.
  colnames(meta_features) <- sprintf("cluster_%s.%s",
                                     rep(sort(unique(feature_clusters)),
                                         each = num_btwn_pcs),
                                     rep(1:num_btwn_pcs,
                                         times = length(unique(feature_clusters))))
  row.names(meta_features) <- row.names(omx_data)
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
  # Run the meta adversarial game on meta features
  acc <- vector(mode = "numeric", length = length(unique(feature_clusters)) + 1L)
  names(acc) <- c("meta_model",
                  paste0("cluster_", sort(unique(feature_clusters))))
  if (isTRUE(verbose)) {
    message("Fitting the meta ARF model...\n")
  }
  meta_arf <- arf::adversarial_rf(
    x = meta_features,
    num_trees = num_trees,
    min_node_size = min_node_size,
    verbose = verbose
  )
  acc["meta_model"] <- meta_arf$acc[length(meta_arf$acc)]
  # Density estimation for meta model
  meta_forde <- arf::forde(
    arf = meta_arf,
    x = meta_features,
    oob = oob,
    family = family,
    finite_bounds = finite_bounds,
    alpha = alpha,
    epsilon = epsilon,
    parallel = parallel
  )

  meta_model <- list(meta_arf = meta_arf,
                     meta_forde = meta_forde,
                     onset_loadings = onset_loadings)
  # For each cluster fit ARF model
  arf_clusters <- function(cluster) {
    if (isTRUE(verbose)) {
      message(paste0("Fitting ARF model for cluster ",
                     which(sort(unique(feature_clusters)) == cluster),
                     " of ", length(unique(feature_clusters)),
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
      always.split.variables = if (isTRUE(always_split_meta)) {
        colnames(meta_features)
      } else {
        NULL
      },
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
      parallel = parallel
    )
    # Export the iso_arf object
    return(list(cluster_id = cluster,
                ftr_in_cluster = colnames(omx_data)[ftr_in_cluster],
                iso_arf = iso_arf,
                iso_forde = iso_forde))
  }

  arf_models <- lapply(
    sort(unique(feature_clusters)),
    arf_clusters
  )
  # Save accuracy for each cluster model
  for (i in seq_along(arf_models)) {
    acc[paste0("cluster_", arf_models[[i]]$cluster_id)] <-
      arf_models[[i]]$iso_arf$acc[length(arf_models[[i]]$iso_arf$acc)]
  }
  # Export feature_cluster_df and meta to arf_models as well.
  hd_arf <- list(meta_model = meta_model,
                 models = arf_models,
                 cluster = data.frame(
                   feature = colnames(omx_data),
                   cluster = feature_clusters
                 ),
                 meta_features = meta_features,
                 clin_lab_rfae = if (exists("cli_lab_rfae")) cli_lab_rfae else NULL,
                 omx_features = cln_omx_data,
                 cli_lab_features = cln_cli_lab_data,
                 accuracy = acc
  )
  class(hd_arf) <- "harf"
  return(hd_arf)
}
