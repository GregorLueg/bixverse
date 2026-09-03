# Download the Baran-Gale, et al. ageing thymus data

Downloads the mouse ageing thymus droplet experiment as a
`SingleCellExperiment`. Roughly 69k thymic epithelial cells across three
ages, with real change in cell type proportions between them.

That compositional change is what makes it the example for
[`get_miloR_abundances_sc()`](https://gregorlueg.github.io/bixverse/reference/get_miloR_abundances_sc.md)
and
[`meld_sc()`](https://gregorlueg.github.io/bixverse/reference/meld_sc.md).
A stimulation experiment moves cells in embedding space without moving
the proportions, so almost every neighbourhood comes out significant and
the result says nothing. Load it with
[`load_sce()`](https://gregorlueg.github.io/bixverse/reference/load_sce.md).

## Usage

``` r
download_thymus_ageing(quiet = FALSE)
```

## Arguments

- quiet:

  Boolean. If the download shall be quiet.

## Value

String. The path to the qs2 file holding the `SingleCellExperiment`.

## Details

Sourced from the `MouseThymusAgeing` Bioconductor package. See
`data-raw/zenodo_sce_datasets.R` for the preparation. The object carries
no gene names, so
[`load_sce()`](https://gregorlueg.github.io/bixverse/reference/load_sce.md)
falls back to the Ensembl identifiers in the `rowData`, and 336 genes
with no annotation at all get a generated identifier.

Served from the `bixverse-data` GitHub release. There is no Zenodo
fallback for this one yet. At 430 MB it is the largest of the example
datasets, so expect the first call to take a while.

## References

Baran-Gale, et al., Development, 2020
