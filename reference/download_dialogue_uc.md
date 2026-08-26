# Download the ulcerative colitis example data for DIALOGUE

This function downloads the ulcerative colitis subset used by the
DIALOGUE paper into the temporary directory. 5,374 cells across five
cell subtypes and 30 donors, with the inflammation status per biopsy.
Enough samples per cell type for
[`dialogue_sc()`](https://gregorlueg.github.io/bixverse/reference/dialogue_sc.md),
which most single sample example data sets are not.

## Usage

``` r
download_dialogue_uc(quiet = FALSE)
```

## Arguments

- quiet:

  Boolean. If the download shall be quiet.

## Value

String. The path to the ulcerative colitis h5ad file.

## Details

The published matrix holds `log2(TPM / 10 + 1)`, which is what Smillie,
et al. released; the raw counts were never published. The file served
here has been linearised back to the TPM/10 scale and rounded, so it
loads through
[`load_h5ad()`](https://gregorlueg.github.io/bixverse/reference/load_h5ad.md)
like any count matrix and gets bixverse's own log CPM on the way in. See
`data-raw/dialogue_uc.R` for the preparation.

## References

Smillie, et al., Cell, 2019; Jerby-Arnon and Regev, Nat. Biotechnol.,
2022
