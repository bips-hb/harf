## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  dev = "png",
  type = "cairo"
)

## ----install_libraries, warning = FALSE, message = FALSE, eval = FALSE--------
# install.packages("data.table")
# install.packages("rsvd")
# install.packages("Rtsne")
# install.packages("cowplot")
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("SingleCellExperiment")
# BiocManager::install("scater")
# install.packages("pROC")
# install.packages("caret")
# install.packages("ggplot2")
# install.packages("corrplot")
# install.packages("ranger")
# install.packages("doParallel")

## ----set_par, echo=FALSE, include=FALSE---------------------------------------
# Save graphics parameters properly
old_par <- par(no.readonly = TRUE)
on.exit(par(old_par), add = TRUE)

# Save environment variables
old_env <- Sys.getenv(c("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"))

on.exit({
  Sys.setenv(
    OMP_NUM_THREADS = old_env[1],
    MKL_NUM_THREADS = old_env[2],
    OPENBLAS_NUM_THREADS = old_env[3]
  )
}, add = TRUE)

# Set safe CRAN-recommended thread limits
chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
num_cores <- if (nzchar(chk)) 1 else 2

Sys.setenv(
  OMP_NUM_THREADS = num_cores,
  MKL_NUM_THREADS = num_cores,
  OPENBLAS_NUM_THREADS = num_cores
)

# Ensure default layout (optional but safe)
par(mfrow = c(1, 1))

## ----load_libraries, warning = FALSE, message = FALSE-------------------------
library(harf)
library(data.table)
library(rsvd)
library(Rtsne)
library(cowplot)
library(SingleCellExperiment)
library(ggplot2)
library(corrplot)
library(scater)
library(pROC)
library(caret)
library(ranger)
library(doParallel)

## ----numcore_lib, echo=FALSE, include=FALSE-----------------------------------
data.table::setDTthreads(num_cores)

## ----single_cell_example, include=TRUE, eval=TRUE, message=FALSE, warning=FALSE----
data("single_cell")
chunk_size <- 5

## ----harf_training, include = TRUE, eval = TRUE, message=FALSE, warning = FALSE----
my_omx_data <- single_cell[ , - which(colnames(single_cell)  == "cell_type")]
my_cli_lab_data <- data.frame(cell_type = single_cell$cell_type)
harf_model <- h_arf(
 omx_data = my_omx_data,
 cli_lab_data = my_cli_lab_data,
 feature_ordering = colnames(single_cell),
 parallel = FALSE,
 chunk_size = chunk_size,
 verbose = TRUE
)
str(harf_model,max.level = 1)

## ----harf_accuracy, include = TRUE, eval = TRUE, message=FALSE, fig.width = 7, fig.height = 3----
acc_df <- data.frame(
  Region = names(harf_model$accuracy),
  Accuracy = harf_model$accuracy
)
acc_plot <- ggplot2::ggplot(acc_df, ggplot2::aes(x = Region, y = Accuracy)) +
  ggplot2::geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +
  ggplot2::geom_bar(stat = "identity", fill = "steelblue") +
  ggplot2::ylim(0, 1) +
  ggplot2::labs(title = "h-ARF Convergence Accuracy",
                x = "Regions",
                y = "Accuracy") +
  ggplot2::theme_minimal() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
acc_plot

## ----harf_synthetic_data, include = TRUE, eval = TRUE, message=FALSE----------
set.seed(32)
synth_single_cell <- h_forge(
  harf_obj = harf_model,
  n_synth = nrow(single_cell), 
  evidence = NULL,
  parallel = FALSE,
  verbose = FALSE
  )

## ----harf_correlation_matrices, include = TRUE, eval = TRUE, message = FALSE, fig.width = 3.5, fig.height = 3----
# Re-arrange data by grouping gene by clusters
cluster_feature <- copy(harf_model$cluster)
setorder(cluster_feature, cluster)
orig_clustered <- single_cell[ , c("cell_type", cluster_feature$feature)]
synth_clustered <- as.data.frame(synth_single_cell)[ , c("cell_type", cluster_feature$feature)]
plot_corr <- function(dt, title) {
  corr_matrix <- cor(dt[ , 2:21], method = "spearman")
  corrplot(corr_matrix,
           method = "circle",
           tl.col = "black",
           tl.pos = "n",
           title = title,
           mar = c(0, 0, 1, 0))
}

## ----harf_tsne, include = TRUE, eval = TRUE, message = FALSE, fig.width = 7.1, fig.height = 3.5----
tsne_it <- function (sc_data, perp = 30, title = "") {
  # Create SingleCellExperiment object
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = t(as.matrix(sc_data[ , - which(colnames(sc_data)  == "cell_type")])))
  )
  SingleCellExperiment::logcounts(sce) <- SingleCellExperiment::counts(sce) # Log-normalization
  sce$cell_type <- sc_data$cell_type
  pc_sce <- rpca(t(SingleCellExperiment::counts(sce)))
  # tSNE with rotated pcs
  ts_sce <- Rtsne::Rtsne(
    pc_sce$x %*% pc_sce$rotation,
    perplexity = perp,
    verb = FALSE,
    pca = FALSE,
    check_duplicates = FALSE
  )
  SingleCellExperiment::reducedDim(sce, "tsne") = ts_sce$Y
  sce_plot <- scater::plotReducedDim(sce, "tsne", colour_by = "cell_type") + 
  ggplot2::ggtitle(title) +
  ggplot2::theme(legend.position = "bottom")
  return(sce_plot)
}  
orig_plot <- tsne_it(single_cell,
                     perp = 30,
                     title = "Original")
synth_plot <- tsne_it(as.data.frame(synth_single_cell),
                      perp = 30,
                      title = "Synthetic")
legend <- cowplot::get_legend(
  orig_plot + theme(legend.position = "bottom")
)
orig_plot <- orig_plot + theme(legend.position = "none")
synth_plot <- synth_plot + theme(legend.position = "none")
all_plots <- cowplot::plot_grid(orig_plot,
                                synth_plot,
                                ncol = 2)
par(mfrow = c(1,2))
plot_corr(orig_clustered, "Original")
plot_corr(synth_clustered, "Synthetic")
par(mfrow = c(1, 1))
print(all_plots)

## ----harf_tsne_plot, include = FALSE, eval = FALSE, message = FALSE, fig.width = 7, fig.height = 6.5----
# par(mfrow = c(1, 1))
# plot_grid(
#   all_plots,
#   legend,
#   ncol = 1,
#   rel_heights = c(1, 0.2)
# )

## ----harf_conditional_expectation, include = TRUE, eval = TRUE, message=FALSE, fig.width = 7, fig.height = 3.5----
sub_cell_type <- c("lung", "liver", "esophagus_muc")
single_cell_list <- lapply(sub_cell_type, function (ct) {
  ct_synth <- h_forge(
        harf_obj = harf_model,
        n_synth = sum(single_cell$cell_type == ct),
        evidence = data.frame(cell_type = ct),
        verbose = FALSE,
        parallel = FALSE
      )
  return(ct_synth)
})
cond_synth_single_cell <- do.call(rbind, single_cell_list)
cond_synth_clustered <- as.data.frame(cond_synth_single_cell)[ , c("cell_type", cluster_feature$feature)]
cond_synth_plot <- tsne_it(as.data.frame(cond_synth_single_cell),
                           perp = 30,
                           title = "Conditional resampling")
sub_legend <- cowplot::get_legend(
  cond_synth_plot + theme(legend.position = "none")
)
cond_synth_plot <- cond_synth_plot + theme(legend.position = "none")
sub_single_cell <- orig_clustered[orig_clustered$cell_type %in% sub_cell_type , ]
sub_orig_plot <- tsne_it(sub_single_cell,
                     perp = 30,
                     title = "Original")
legend <- cowplot::get_legend(
  orig_plot + theme(legend.position = "bottom")
)
sub_all_plots <- cowplot::plot_grid(sub_orig_plot,
                                cond_synth_plot,
                                ncol = 2)
par(mfrow = c(1,2))
plot_corr(sub_single_cell, "Original")
plot_corr(cond_synth_clustered, "Synthetic")
par(mfrow = c(1, 1))
plot_grid(sub_all_plots, sub_legend, ncol = 1, rel_heights = c(1, 0.2))
par(mfrow = c(1, 1))

## ----kich_example, include=TRUE, eval=TRUE, message=FALSE, warning=FALSE------
data("kich")
seed <- 123
set.seed(seed)
train_idx <- caret::createDataPartition(
  kich$tumor_stage,
  p = 0.7,
  list = FALSE
)
train_idx <- train_idx[ , "Resample1"]

## ----kich_rf, include = TRUE, eval = TRUE, message=FALSE, warning = FALSE-----
set.seed(seed)
rf_model <- ranger(tumor_stage ~ .,
                   data = kich[train_idx, ],
                   num.trees = 500,
                   probability = FALSE)
# Estimate AUC on the test set
test_pred <- predict(rf_model, data = kich[-train_idx, ])$predictions
test_labels <- kich$tumor_stage[-train_idx]
auc_original <- roc(test_labels, as.numeric(test_pred))$auc
print(paste("AUC:", auc_original))

## ----kich_harf_training, include = TRUE, eval = TRUE, message=FALSE, warning = FALSE----
set.seed(seed)
kich_harf <- h_arf(
  omx_data = kich[train_idx , !(colnames(kich) %in% c("tumor_stage",
                                                      "age", "gender"))],
                   cli_lab_data = kich[train_idx, c("tumor_stage", "age", "gender")],
                   chunk_size = 10,
                   target = "tumor_stage",
                   verbose = TRUE
  )

## ----kich_synthetic_data, include = TRUE, eval = TRUE, message=FALSE----------
# Synthetic data
set.seed(seed)
synth_kich <- h_forge(
  harf_obj = kich_harf,
  n_synth = length(train_idx),
  evidence = NULL,
  parallel = FALSE
)
synth_kich <- as.data.frame(synth_kich)

## ----kich_rf_synth, include = TRUE, eval = TRUE, message=FALSE, warning = FALSE----
set.seed(seed)
rf_model_synth <- ranger(tumor_stage ~ .,
                          data = synth_kich,
                          num.trees = 500,
                          probability = FALSE)
# Estimate AUC on the test set
test_pred_synth <- predict(rf_model_synth, data = kich[-train_idx, ])$predictions
auc_synth <- roc(test_labels, as.numeric(test_pred_synth))$auc
auc_comparison <- data.frame(
  Model = c("Original", "Synthetic"),
  AUC = c(auc_original, auc_synth)
)
print(auc_comparison)

## ----reset_par, echo=FALSE, include=FALSE-------------------------------------
par(old_par)

