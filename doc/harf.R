## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
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
# install.packages("ggplot2")
# install.packages("corrplot")
# install.packages("doParallel")

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
# Register cores - Unix
library(doParallel)
registerDoParallel(cores = 2)
# Set seed
set.seed(123, "L'Ecuyer-CMRG")

## ----data_example, include=TRUE, eval=TRUE, message=FALSE, warning=FALSE------
data("single_cell")
# Do not parellelize on CRAN
parallel <- ifelse(Sys.getenv("NOT_CRAN") == "true", FALSE, FALSE)

## ----harf_training, include = TRUE, eval = TRUE, message=FALSE, warning = FALSE----
my_omx_data <- single_cell[ , - which(colnames(single_cell)  == "cell_type")]
my_cli_lab_data <- data.frame(cell_type = single_cell$cell_type)
harf_model <- h_arf(
 omx_data = my_omx_data,
 cli_lab_data = my_cli_lab_data,
 target = "cell_type",
 parallel = parallel,
 chunck_size = 4,
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
  ggplot2::labs(title = "HARF Convergence Accuracy",
                x = "Regions",
                y = "Accuracy") +
  ggplot2::theme_minimal() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
acc_plot

## ----harf_synthetic_data, include = TRUE, eval = TRUE, message=FALSE----------
synth_single_cell <- h_forge(
  harf_obj = harf_model,
  n_synth = nrow(single_cell), 
  evidence = NULL,
  parallel = parallel,
  verbose = FALSE
  )

## ----harf_correlation_matrices, include = TRUE, eval = TRUE, message = FALSE, fig.width = 3.5, fig.height = 3----
# Re-arrange data by grouping gene by clusters
cluster_feature <- copy(harf_model$cluster)
setorder(cluster_feature, cluster)
orig_clustered <- single_cell[ , c("cell_type", cluster_feature$feature)]
synth_clustered <- as.data.frame(synth_single_cell)[ , c("cell_type", cluster_feature$feature)]
plot_corr <- function(dt, title) {
  corr_matrix <- cor(dt[ , 2:51], method = "spearman")
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
        parallel = parallel
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

