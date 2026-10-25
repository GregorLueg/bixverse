# Load the human (hg38) transcription factor list used by SCENIC

Returns the human transcription factor list from the Aerts lab
(`allTFs_hg38.txt`), the candidate regulator list pySCENIC uses for GRN
inference. A copy ships with the package because the original server is
not always reachable. All credit to the Aerts lab; the original lives at
`https://resources.aertslab.org/cistarget/tf_lists/`.

## Usage

``` r
load_scenic_tfs_hg38()
```

## Value

Character vector of transcription factor gene symbols.

## References

Aibar, et al., Nat Methods, 2017

## Examples

``` r
# the candidate regulators for SCENIC
tfs <- load_scenic_tfs_hg38()
head(tfs)
#> [1] "ZNF354C" "KLF12"   "ZNF143"  "ZIC2"    "ZNF274"  "SP2"    
```
