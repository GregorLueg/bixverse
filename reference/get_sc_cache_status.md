# Status of everything held in a single cell object's caches

One row per cached artefact, saying whether it still agrees with the
current cell set and with the artefacts it was derived from. Artefacts
computed before provenance stamping existed report `stamped = FALSE` and
are never flagged stale, because nothing is known about them either way.

## Usage

``` r
get_sc_cache_status(object)
```

## Arguments

- object:

  `SingleCells`, `SingleCellsSubset`, `MetaCells` or
  `SingleCellsMultiModal` class.

## Value

A `data.table` with the columns

- modality - The modality the artefact lives in.

- artefact - One of `pca`, `embedding`, `knn`, `snn`, `magic`.

- name - The embedding name, `NA` for the others.

- stamped - Whether the artefact carries a provenance stamp.

- stale - Whether it disagrees with the current state.

- reason - Why it is stale, `NA` otherwise.

- id - The artefact's stamp id.

- from - List column of the parent stamp ids.

## Examples

``` r
# what the object holds and whether it still agrees with the cells
sc <- demo_single_cells()
get_sc_cache_status(sc)
#>    modality artefact   name stamped  stale reason               id
#>      <char>   <char> <char>  <lgcl> <lgcl> <char>           <char>
#> 1:      rna      pca   <NA>    TRUE  FALSE   <NA> cf29f84484d82c85
#> 2:      rna      knn   <NA>    TRUE  FALSE   <NA> 578d3d7b717b8abe
#> 3:      rna      snn   <NA>    TRUE  FALSE   <NA> 82e64ac3f6879dc5
#>                from
#>              <list>
#> 1:                 
#> 2: cf29f84484d82c85
#> 3: 578d3d7b717b8abe

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
