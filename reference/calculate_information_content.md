# Calculate the information content for each ontology term

this function will calculate the information content of each provided
term based on a list of ancestors, which is a named list of terms with
ancestor identifiers as their values. Can be calculated using
[`get_ontology_ancestry()`](https://gregorlueg.github.io/bixverse/reference/get_ontology_ancestry.md).
The information content is calculated as
`-log2(number descendant/total terms in the ontology)`. More information
can be found
[here](https://yulab-smu.top/biomedical-knowledge-mining-book/semantic-similarity-overview.html).

## Usage

``` r
calculate_information_content(ancestor_list)
```

## Arguments

- ancestor_list:

  List. Named list of terms with ancestor identifiers as their values

## Value

A named list of each term and their information content as values.

## Examples

``` r
# information content from the descendant counts of a toy ontology
onto <- data.table::data.table(
  parent = c("a", "b", "b", "b", "c"),
  child = c("b", "c", "d", "e", "f")
)
ancestry <- get_ontology_ancestry(onto)
unlist(calculate_information_content(ancestry$descendants))
#>         a         b         c         d         e         f 
#> 0.0000000 0.1823216 1.0986123 1.7917595 1.7917595 1.7917595 
```
