# Identify the empty droplets from the per-barcode library sizes

**\[experimental\]** Kept on the Rust side so the knee detector exists
once: it has to match `scipy.ndimage.gaussian_filter1d` and
`numpy.gradient` closely enough to land on the same rank as the
CellSweep reference, and a second copy in R would drift.

## Usage

``` r
rs_sc_infer_empty_droplets(lib_size, empty_params)
```

## Arguments

- lib_size:

  Integer vector. Library size per barcode, in store order.

- empty_params:

  List. Parameter list, see
  [`params_sc_empty_droplets()`](https://gregorlueg.github.io/bixverse/reference/params_sc_empty_droplets.md).
  `"supplied"` is rejected, since there is nothing to infer.

## Value

A logical vector that is `TRUE` where the barcode is an empty droplet.
