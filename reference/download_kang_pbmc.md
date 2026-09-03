# Download the Kang, et al. IFN-beta stimulated PBMC data

Downloads the Kang, et al. demultiplexed PBMC experiment as a
`SingleCellExperiment`. 8 lupus donors, PBMCs, control against IFN-beta
stimulated in vitro, roughly 29k cells before quality control.

The design is paired within donor, which is what makes it the example
for
[`nebula_sc()`](https://gregorlueg.github.io/bixverse/reference/nebula_sc.md):
the donor is a genuine random effect and every donor contributes to both
arms. Load it with
[`load_sce()`](https://gregorlueg.github.io/bixverse/reference/load_sce.md).

## Usage

``` r
download_kang_pbmc(quiet = FALSE)
```

## Arguments

- quiet:

  Boolean. If the download shall be quiet.

## Value

String. The path to the qs2 file holding the `SingleCellExperiment`.

## Details

Sourced from the `muscData` Bioconductor package, itself from GEO
accession `GSE96583`. See `data-raw/zenodo_sce_datasets.R` for the
preparation. The `multiplets` column marks doublets and ambiguous
droplets, which are worth dropping before any modelling.

Served from the `bixverse-data` GitHub release. There is no Zenodo
fallback for this one yet.

## References

Kang, et al., Nat. Biotechnol., 2018
