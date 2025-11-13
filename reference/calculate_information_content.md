# Calculate the information content for each ontology term

this function will calculate the information content of each provided
term based on a list of ancestors, which is a named list of terms with
ancestor identifiers as their values. Can be calculated using
[`get_ontology_ancestry()`](get_ontology_ancestry.md). The information
content is calculated as
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
