# Forest generator for high-dimensional data.

The algorithm uses the high-dimensional arf model to generate synthetic
data

## Usage

``` r
h_forge(
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
)
```

## Arguments

- harf_obj:

  A pre-trained harf model.

- n_synth:

  Number of synthetic samples to generate.

- evidence:

  Optional set of conditioning events. This will be further passed to
  the `forde` function in each isolated regions. See `forde` for
  details.

- omx_onset_data:

  Optional data.frame of conditional onset omics features.

- evidence_row_mode:

  Interpretation of rows in multi-row evidence. See `forde` for details.

- round:

  Round continuous variables to their respective maximum precision in
  the real data set? See `forde` for details.

- sample_NAs:

  Sample NAs respecting the probability for missing values in the
  original data? See `forde` for details.

- nomatch:

  What to do if no leaf matches a condition in evidence? Options are to
  force sampling from a random leaf ("force") or return NA ("na"). The
  default is "force".

- verbose:

  What to do if no leaf matches a condition in `evidence`? See `forde`
  for details.

- stepsize:

  How many rows of evidence should be handled at each step? See `forde`
  for details.

- parallel:

  Compute in parallel? See `forde` for details.

## Value

A data.table containing the generated synthetic omics data.

## See also

[`forge`](https://bips-hb.github.io/arf/reference/forge.html) for
details on the forging process.

## Author

Césaire Fouodo

## Examples

``` r
if (FALSE) { # \dontrun{
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
 # Conditional resampling from harf_model
 set.seed(142)
 lung_single_cell <- h_forge(
     harf_obj = harf_model,
     n_synth = sum(single_cell$cell_type == "lung"),
     evidence = data.frame(cell_type = "lung")
    )
} # }
```
