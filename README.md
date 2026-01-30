<!-- badges: start -->

[![R-CMD-check](https://github.com/bips-hb/harf/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/bips-hb/harf/actions/workflows/R-CMD-check.yaml)
[![Lifecycle:
Unstable](https://img.shields.io/badge/lifecycle-Unstable-red.svg)](https://lifecycle.r-lib.org/articles/stages.html#Stable)
[![CRAN
Status](https://img.shields.io/badge/CRAN-harf-orange)](https://cran.r-project.org/package=harf)
<!-- [![codecov](https://codecov.io/github/bips-hb/harf/graph/badge.svg?token=SZU7NGK8G8)](https://app.codecov.io/github/bips-hb/harf/) -->

<!-- badges: end -->

## harf: High-Dimensional Adversarial Random Forests for Omics Data

## Introduction

Adversarial random forests (ARFs) have recently been introduced as data
synthetizer. Based on unsupervised RFs, ARFs rely on a recursive
adversarial procedure in which trees progressively learn the structural
properties of the data through alternating rounds of data generation and
discrimination. The unsupervised classificatin task is achieved by
introducing a synthetic response variable $`y \in \{0,1\}`$, where
$`y = 0`$ denotes synthetic data and $`y = 1`$ denotes original data.
During the discrimination phase, the objective is to distinguish
original observations from synthetic ones obtained via marginal
resampling of the original features in terminal nodes. The adversarial
process stops when the prediction accuracy for $`y`$ falls below a
predefined threshold.

While adversarial random forests (ARFs) have demonstrated strong
performance in various low-dimensional settings, their behavior in
high-dimensional contexts, such as omics data, remains less well
understood. A key assumption of ARFs — that feature distributions are
independent within terminal nodes — may be violated in high-dimensional
settings, where a small subset of features can be highly predictive of
the synthetic response variable $`y`$, while features not used for
splitting may remain correlated within terminal nodes.

The high-dimensional adversarial random forests (HARFs) package
addresses this limitation by identifying isolated regions of the feature
space in which the independence assumption is more likely to hold. In
addition to these isolated regions, the HARF procedure constructs a
low-dimensional meta-space that captures the relationships among
regions. Within each isolated region, separate ARFs are trained to
better capture local data structures. To preserve dependencies between
regions, meta-features representing region membership are used in the
meta-space.

## Package installation

``` r
devtools::install_github("bips-hp/harf")
```

## Training a HARF model

Using the built-in `single_cell` dataset, we illustrate how to train a
HARF model.

``` r
data(single_cell)
harf_model <- h_arf(
 omx_data = single_cell[ , - which(colnames(single_cell)  == "cell_type")],
 cli_lab_data = data.frame(cell_type = single_cell$cell_type)
)
```

## Generating synthetic data

``` r
synth_single_cell <- h_forge(
  harf_obj = harf_model,
  n_synth = nrow(single_cell)
  )
```

## Conditional resampling

``` r
  lung_single_cell <- h_forge(
      harf_obj = harf_model,
      n_synth = sum(single_cell$cell_type == "lung"),
      evidence = data.frame(cell_type = "lung")
     )
```

We refer to the package vignette for detailed examples and explanations.

``` r
vignette("harf")
```

## References

- Fouodo, C. J. K., et al. (2026). High-dimensional adversarial random
  forests. Submission. Link [don’t
  click](https://arxiv.org/abs/2405.12345).

- Watson, D. S., Blesch, K., Kapar, J. & Wright, M. N. (2023).
  Adversarial random forests for density estimation and generative
  modeling. In Proceedings of the 26th International Conference on
  Artificial Intelligence and Statistics. Link
  [here](https://proceedings.mlr.press/v206/watson23a.html).
