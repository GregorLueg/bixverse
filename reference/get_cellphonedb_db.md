# Download the CellPhoneDB ligand-receptor database

Fetches the CellPhoneDB v5.0.0 database from the
[cellphonedb-data](https://github.com/ventolab/cellphonedb-data)
repository and flattens it to one row per interaction. Each partner is
resolved to the gene symbols of its subunits: one gene for a single
protein, several for a complex. The download is cached in `dir`.

## Usage

``` r
get_cellphonedb_db(dir = tempdir(), .verbose = TRUE)
```

## Arguments

- dir:

  String. Directory to store and look for `cellphonedb.zip`. Defaults to
  [`tempdir()`](https://rdrr.io/r/base/tempfile.html).

- .verbose:

  Boolean. Controls verbosity of the download.

## Value

A data.table with one row per interaction and the columns:

- interaction_id - The CellPhoneDB interaction identifier.

- partner_a, partner_b - The partner names: the gene symbol for a single
  protein, the complex name otherwise.

- genes_a, genes_b - List columns with the gene symbols of the subunits.

- is_complex_a, is_complex_b - Is the partner a complex.

- directionality - E.g. `"Ligand-Receptor"`.

- classification - The signalling classification.

## References

Troulé et al., Nat Protoc, 2025.
