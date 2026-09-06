# Derive sparse module membership from a loading matrix

Turns a `gene x k` loading matrix from a matrix factorisation (ICA, NMF,
DGRDL) into a membership table by keeping the tails of each component's
loading distribution.

Two thresholding rules, selected via `membership_params`:

- `"zscore"` - standardisation per component, keeping `abs(z) > cutoff`.
  No distributional assumption beyond rough symmetry.

- `"fdr"` - two-sided p-values against a Normal null fitted per
  component, Benjamini-Hochberg adjusted, keeping `padj < fdr`.

The standardisation itself is controlled by `membership_params$scaling`:
`"robust"` centres and scales by the median and MAD, `"standard"` by the
mean and standard deviation. The latter is stricter and keeps fewer
genes on skewed loadings, which is common for NMF.

## Usage

``` r
modules_from_loadings(loadings, membership_params = params_module_membership())
```

## Arguments

- loadings:

  Numeric matrix. `gene x k`, with row and column names.

- membership_params:

  List. See
  [`params_module_membership()`](https://gregorlueg.github.io/bixverse/reference/params_module_membership.md).

## Value

A data.table with columns `gene`, `module_id`, `loading`, `sign` and the
per-component score (`z` or `padj` depending on the method). One row per
surviving (gene, component) pair, ordered by component then by
descending absolute loading.

## Examples

``` r
# turn NMF gene loadings into module membership
syn <- generate_gene_module_data(n_samples = 24L, n_genes = 60L)
# NMF needs a non-negative matrix
mat <- syn$data - min(syn$data)
obj <- BulkCoExp(mat, syn$meta_data)
obj <- preprocess_bulk_coexp(
  obj, hvg = NULL, scaling = FALSE, .verbose = FALSE
)
obj <- nmf_bulk(obj, k = 4L, .verbose = FALSE)
head(modules_from_loadings(get_nmf_gene_loadings(obj)))
#>          gene module_id  loading   sign        z
#>        <char>    <char>    <num> <char>    <num>
#> 1: feature_44   comp_01 3.675330    pos 38.28581
#> 2: feature_40   comp_01 3.643421    pos 37.94433
#> 3: feature_42   comp_01 3.633641    pos 37.83967
#> 4: feature_31   comp_01 3.630990    pos 37.81130
#> 5: feature_34   comp_01 3.602253    pos 37.50377
#> 6: feature_45   comp_01 3.590235    pos 37.37517
```
