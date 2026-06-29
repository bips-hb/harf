## harf: High-dimensional Adversarial Random Forests for Omics Data

## Introduction

Adversarial random forests (ARFs) have recently been introduced as a
well-performing generative method for low-dimensional datasets. Based on
unsupervised RFs, ARFs rely on a recursive adversarial procedure in
which trees progressively learn the structural properties of the data
through alternating rounds of data generation and discrimination. The
unsupervised classification task is achieved by introducing a synthetic
response variable $`y \in \{0,1\}`$, where $`y = 0`$ denotes synthetic
data and $`y = 1`$ denotes original data. During the discrimination
phase, the objective is to distinguish original observations from
synthetic ones. Synthetic data are created during the the adversarial
game by marginal resampling of the original features in terminal nodes.
The adversarial process stops when the prediction accuracy for $`y`$
falls below a predefined threshold, $`0.5`$, for example.

While ARFs have demonstrated strong performance in various
low-dimensional settings, their behavior in high-dimensional contexts,
such as omics data, remains non investigated. The key assumption of ARFs
— that feature distributions are independent within terminal nodes — may
be violated in high-dimensional settings. For example, if a small subset
of features is highly predictive of the synthetic response variable
$`y`$, the remaining features not used to reach the terminal node may
still be correlated, not allowing for marginal resampling. A typical
resulting behaviour of ARF in such a situation is that the algorithm
terminates without converging.

The high-dimensional adversarial random forests (*h*-ARFs) package
addresses the convergence issue of ARF by identifying isolated regions
of the feature space in which the independence assumption is more likely
to hold in terminal nodes. In addition, the *h*-ARF algorithm constructs
a low-dimensional meta-space that captures the relationships among
regions. Within each isolated region, separate ARFs are trained to
better capture local data structures, followed by training additional
ARF in a meta-region that capture the joint structure between isolated
regions. To synthesize observations, the *h*-ARF algorithm first samples
a region from the meta-space and then condition the resampling in the
other isolated region by synthesized meta-observations. The *h*-ARF
algorithm offers a flexible framework for unconditional and conditional
data generation for downstream clustering and prediction analyses. The
package includes a built-in single-cell RNA-seq datasets and TCGA-KICH
dataset to illustrate its usage. We refer to the package vignette for
detailed examples and explanations.

## Package installation

``` r

devtools::install_github("bips-hb/harf", build_vignettes  = TRUE)
```

Visit our vignette
[here](https://bips-hb.github.io/harf/articles/harf.html) for detailed
examples and explanations.

## References

- Fouodo, C. J. K., Kapar, J. Huels A., Qin, S. Z. & Wright, M. N.
  (2026). High-dimensional adversarial random forests. submitted for
  peer review. Link [don’t click](https://arxiv.org/abs/2405.12345).

- Watson, D. S., Blesch, K., Kapar, J. & Wright, M. N. (2023).
  Adversarial random forests for density estimation and generative
  modeling. In Proceedings of the 26th International Conference on
  Artificial Intelligence and Statistics. Link
  [here](https://proceedings.mlr.press/v206/watson23a.html).
