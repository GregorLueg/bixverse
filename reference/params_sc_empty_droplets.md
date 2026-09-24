# Parameters for identifying empty droplets

CellSweep trains its ambient profile on the empty droplets, so it needs
them called before it runs. Four ways to do that, in descending order of
how much you should trust them. `"supplied"` takes an existing logical
column from the obs table. This is the common case: if you have
CellRanger's filtered barcode list you already know which barcodes are
empty, and nothing here will beat that. `"umi_cutoff"` calls everything
below an absolute library size empty. `"expected_cells"` turns a cell
count into that cutoff via the sorted library sizes. `"knee"` finds the
cutoff from the curvature of the rank / log-count curve; it is
experimental in the reference too, and smoothing puts the curvature
minimum a few ranks ahead of the actual cliff, so it sweeps up the real
barcodes nearest the transition.

## Usage

``` r
params_sc_empty_droplets(
  method = c("supplied", "umi_cutoff", "expected_cells", "knee"),
  is_empty_column = NULL,
  umi_cutoff = NULL,
  expected_cells = NULL
)
```

## Arguments

- method:

  String. How the empty droplets are called, see the description. One of
  `c("supplied", "umi_cutoff", "expected_cells", "knee")`. Defaults to
  `"supplied"`.

- is_empty_column:

  String or `NULL`. Name of the logical obs column holding the mask.
  Required for `method = "supplied"`, ignored otherwise. Defaults to
  `NULL`.

- umi_cutoff:

  Integer or `NULL`. Barcodes with a library size strictly below this
  are empty. Required for `method = "umi_cutoff"`, ignored otherwise.
  Defaults to `NULL`.

- expected_cells:

  Integer or `NULL`. Number of real cells expected in the run. Required
  for `method = "expected_cells"`, ignored otherwise. Defaults to
  `NULL`.

## Value

A named list with the following elements:

- method - String. How the empty droplets are called, see the
  description. One of
  `c("supplied", "umi_cutoff", "expected_cells", "knee")`. Defaults to
  `"supplied"`.

- is_empty_column - String or `NULL`. Name of the logical obs column
  holding the mask. Required for `method = "supplied"`, ignored
  otherwise. Defaults to `NULL`.

- umi_cutoff - Integer or `NULL`. Barcodes with a library size strictly
  below this are empty. Required for `method = "umi_cutoff"`, ignored
  otherwise. Defaults to `NULL`.

- expected_cells - Integer or `NULL`. Number of real cells expected in
  the run. Required for `method = "expected_cells"`, ignored otherwise.
  Defaults to `NULL`.
