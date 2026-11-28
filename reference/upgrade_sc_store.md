# Upgrade the on-disk counts to the v5 file format

Stores written by older versions use the v3 or v4 binary format. v3
stores no longer open at all; v4 stores still load, with a warning. This
converts `counts_cells.bin` and `counts_genes.bin` in place. The result
is checked before it replaces the originals, so a failure leaves the
store untouched. Stores already in v5 are left alone. Do not run this
while anything else (another R session, a background job) uses the
files.

## Usage

``` r
upgrade_sc_store(
  object,
  codec = c("lz4", "zstd", "none"),
  keep_backup = TRUE,
  .verbose = TRUE
)
```

## Arguments

- object:

  `SingleCells` class.

- codec:

  String. Chunk compression, one of `c("lz4", "zstd", "none")`. `"lz4"`
  decodes fastest, `"zstd"` gives smaller files but decodes slower.

- keep_backup:

  Boolean. Keep the original files as `<file>.v<n>.bak`, `n` being their
  version. They are as large as the original store, delete them yourself
  once you are happy with the result.

- .verbose:

  Boolean. Controls verbosity of the function.

## Value

The object, invisibly.

## Examples

``` r
# demo stores are written in v5 already, so this is a no-op
sc <- demo_single_cells()
sc <- upgrade_sc_store(sc, .verbose = FALSE)

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
