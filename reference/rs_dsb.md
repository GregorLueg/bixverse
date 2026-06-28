# Run DSB normalisation on raw ADT counts

**\[experimental\]** This function applies the DSB algorithm to
normalise antibody-derived tag (ADT) counts from CITE-seq experiments.
Two variants are supported. When `background_counts` is provided,
per-protein ambient background is estimated from empty droplets ("Step
I" of the original paper). When `background_counts` is `NULL`,
per-protein background is estimated by a two-component k-means on the
log-transformed cell counts, with the lower centroid taken as the
background level. An optional second step removes cell-to-cell technical
noise by regressing out PC1 of a noise matrix built from isotype
controls (if available) and the per-cell background mean.

## Usage

``` r
rs_dsb(
  raw_counts,
  background_counts,
  isotype_indices,
  dsb_params,
  scale_factor,
  seed,
  verbose
)
```

## Arguments

- raw_counts:

  Numeric matrix. Cells x proteins matrix of raw ADT counts.

- background_counts:

  Optional numeric matrix. Cells x proteins matrix of empty-droplet ADT
  counts used to estimate per-protein ambient background. Must have the
  same number of columns as `raw_counts`. If `NULL`, the function falls
  back to k-means-based background estimation.

- isotype_indices:

  Integer vector. 0-based column indices into `raw_counts` identifying
  isotype control proteins. Used in Step II if
  `dsb_params$use_isotype_controls = TRUE`. Pass an empty integer vector
  if there are no isotype controls.

- dsb_params:

  List. DSB parameters.

- scale_factor:

  String. One of `"standardise"` or `"mean_subtract"`. Only used when
  `background_counts` is provided. `"standardise"` subtracts the
  per-protein background mean and divides by the per-protein background
  SD. `"mean_subtract"` subtracts the mean only.

- seed:

  Integer. Random seed for k-means initialisation.

- verbose:

  Integer. `0L` - quiet; `1L` - normal verbosity; `2L` - detailed
  verbosity.

## Value

A list with the following elements:

- norm_counts - Numeric matrix. The DSB-normalised cells x proteins
  matrix.

- protein_background_mean - Numeric vector of length `n_proteins`.
  Per-protein background mean used in Step I.

- protein_background_sd - Numeric vector of length `n_proteins`, or
  empty vector if `background_counts` was `NULL`. Per-protein background
  SD used in Step I.

- technical_component - Numeric vector of length `n_cells`, or empty
  vector if `dsb_params$denoise_counts = FALSE`. Per-cell technical
  component regressed out in Step II.

- cellwise_background_mean - Numeric vector of length `n_cells`, or
  empty vector if `dsb_params$denoise_counts = FALSE`. Per-cell
  background mean from the 2-component k-means clustering.
