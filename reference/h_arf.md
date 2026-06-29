# High-dimensional adversarial random forest (h-ARF).

This function extends the adversarial random forest (ARF) algorithm to
high-dimensional settings. It partitions high-dimensional data into
isolated regions and fits ARF models within each region and on a latent
space representing region joint distribution to capture within and
between region feature dependencies.

## Usage

``` r
h_arf(
  omx_data,
  cli_lab_data = NULL,
  target = NULL,
  feature_ordering = NULL,
  omx_onset_data = NULL,
  num_trees = 10,
  min_node_size = 5,
  correlation_method = "spearman",
  correlation_mat = NULL,
  num_btwn_pcs = 2,
  num_onset_pcs = 2,
  chunk_size = 10,
  oob = FALSE,
  family = "truncnorm",
  finite_bounds = "no",
  alpha = 0,
  epsilon = 0,
  parallel = FALSE,
  export_cor_mat = FALSE,
  verbose = FALSE
)
```

## Arguments

- omx_data:

  A data.frame containing omics data where rows represent samples (e.g.
  patients or cells) and columns represent features (e.g., gene or
  protein expressions).

- cli_lab_data:

  A data.frame of clinical or laboratory data. For example, this may be
  a data.frame where rows represent cell types (for single cell data or
  additional clinical patient information, e.g. disease status, age,
  sex, BMI, etc).

- target:

  Optional name of the target variable in `cli_lab_data` for supervised
  downstream analysis. If NULL, training is conducted for unsupervised
  downstream analysis.

- feature_ordering:

  Optional vector of feature names specifying the order of features in
  the synthesized data. If NULL, features are ordered according to their
  original order in `omx_data`, followed by order in `clin_lab_data`.

- omx_onset_data:

  Optional data.frame of conditional onset omics features. This can be
  used to condition the synthesis on specific omics features (e.g.
  expression of a gene of interest). If provided, dimension reduction
  will be performed on these features and the resulting components will
  be included as additional meta features for training the meta
  adversarial model.

- num_trees:

  Number of trees to grow in each Adversarial Random Forest (ARF) model.

- min_node_size:

  Minimum number of samples required to split an internal node in the
  ARF model.

- correlation_method:

  Methods to compute correlation between features. Options include
  "pearson", "spearman", and "kendall". Default is "spearman".

- correlation_mat:

  Optional pre-computed correlation matrix between features. If
  provided, this will be used instead of computing correlations from
  `omx_data`. This can be usefull for tuning hyperparameters of the
  clustering step without having to recompute the correlation matrix
  each time.

- num_btwn_pcs:

  Number of principal components to use for between cluster variability.
  Default is 2.

- num_onset_pcs:

  Number of principal components to use for onset features. Default is
  2.

- chunk_size:

  Size of feature chunks to process at a time. Default is 10.

- oob:

  Logical indicating whether to use out-of-bag samples for density
  estimation. Default is FALSE. Also see `forde`

- family:

  Distribution to use for density estimation of continuous features. See
  `forde` for options.

- finite_bounds:

  Impose finite bounds on all continuous variables? If "local", infinite
  bounds are set to empirical extrema within leaves. If "global",
  infinite bounds are set to global empirical extrema. if "no" (the
  default), infinite bounds are left unchanged. See `forde` for more
  details.

- alpha:

  Optional pseudocount for Laplace smoothing in density estimation. See
  `forde` for more details.

- epsilon:

  Optional slack parameter on empirical bounds. See `forde` for more
  details.

- parallel:

  Logical indicating whether to use parallel processing. Default is
  TRUE.

- export_cor_mat:

  Logical indicating whether to export the correlation matrix used for
  clustering. Default is FALSE.

- verbose:

  Logical indicating whether to print progress messages. Default is
  FALSE.

## Value

An isoARF object containing the fitted adversarial models, and
clustering information.

## References

- Watson et al. (2023). Adversarial Random Forests. Proceedings of the
  International Conference on Machine Learning (PMLR 206).
  <https://proceedings.mlr.press/v206/watson23a.html>

- Fouodo et al. (2026). High-Dimensional Adversarial Random Forests.
  Submitted / under review.

## See also

[h_forge](https://bips-hb.github.io/harf/reference/h_forge.md),
[arf::adversarial_rf](https://bips-hb.github.io/arf/reference/adversarial_rf.html),
[arf::forde](https://bips-hb.github.io/arf/reference/forde.html)

## Examples

``` r
# \donttest{
data(single_cell)
harf_model <- h_arf(
  omx_data = single_cell[ , - which(colnames(single_cell)  == "cell_type")],
  cli_lab_data = data.frame(cell_type = single_cell$cell_type)
)
# Unconditional sampling from harf_model
set.seed(123)
synth_single_cell <- h_forge(
 harf_obj = harf_model,
 n_synth = nrow(single_cell)
 )
#> Generating synthetic data for cluster 1 out of 12...
#> Generating synthetic data for cluster 2 out of 12...
#> Generating synthetic data for cluster 3 out of 12...
#> Generating synthetic data for cluster 4 out of 12...
#> Generating synthetic data for cluster 5 out of 12...
#> Generating synthetic data for cluster 6 out of 12...
#> Generating synthetic data for cluster 7 out of 12...
#> Generating synthetic data for cluster 8 out of 12...
#> Generating synthetic data for cluster 9 out of 12...
#> Generating synthetic data for cluster 10 out of 12...
#> Generating synthetic data for cluster 11 out of 12...
#> Generating synthetic data for cluster 12 out of 12...
 # Conditional resampling from harf_model
 set.seed(142)
 lung_single_cell <- h_forge(
     harf_obj = harf_model,
     n_synth = sum(single_cell$cell_type == "lung"),
     evidence = data.frame(cell_type = "lung")
    )
#> Generating synthetic data for cluster 1 out of 12...
#> Generating synthetic data for cluster 2 out of 12...
#> Generating synthetic data for cluster 3 out of 12...
#> Generating synthetic data for cluster 4 out of 12...
#> Generating synthetic data for cluster 5 out of 12...
#> Generating synthetic data for cluster 6 out of 12...
#> Generating synthetic data for cluster 7 out of 12...
#> Generating synthetic data for cluster 8 out of 12...
#> Generating synthetic data for cluster 9 out of 12...
#> Generating synthetic data for cluster 10 out of 12...
#> Generating synthetic data for cluster 11 out of 12...
#> Generating synthetic data for cluster 12 out of 12...
# }
```
