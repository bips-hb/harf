# Partial Least Squares (PLS) regression for supervised dimension reduction.

This function performs supervised dimension reduction using PLS
regression to find components that capture the covariance between the
omics data and the target variable. The resulting components can be used
as meta features for the meta adversarial model in `h_arf`.

## Usage

``` r
supervised_pls(omx_data, y, num_btwn_pcs = 2)
```

## Arguments

- omx_data:

  A data.frames or matrices representing different regions.

- y:

  A numeric or factor vector representing the target variable.

- num_btwn_pcs:

  Number of principal components to use for between region variability.

## Value

A matrix of latent supervised scores.
