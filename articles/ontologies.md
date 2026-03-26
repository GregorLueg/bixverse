# Ontologies

## Ontologies

`bixverse` provides methods for computing semantic similarities over
ontologies — disease ontologies, gene ontologies, phenotype ontologies,
and so on. The heavy lifting is done yet again in Rust with heavy
parallelism being exploited. This makes this (usually) faster than the
alternative packages.

This vignette walks through the core workflow: building the ancestry
from a parent–child table, deriving information content, and computing
different flavours of semantic similarity.

``` r
library(bixverse)
library(data.table)
library(magrittr)
library(zeallot)
```

## Data

We will use the Human Phenotype Ontology (HPO) bundled with the
`ontologyIndex` package. The only thing `bixverse` needs is a
`data.table` with three columns: `parent`, `child`, and `type`. These
are describing the directed edges of the ontology graph.

``` r
library("ontologyIndex")

data("hpo")

hpo_data <- as.data.table(stack(hpo$children)) %>%
  setnames(old = c("values", "ind"), new = c("child", "parent")) %>%
  .[, c("parent", "child"), with = FALSE] %>%
  .[, type := "parent_of"]
```

## Semantic similarities

### Ancestry and information content

Before computing any similarity we need two ingredients: the full set of
ancestors for every term, and the information content (IC) derived from
the number of descendants each term has. Both are obtained in one step:

``` r
c(ancestors, descendants) %<-% get_ontology_ancestry(hpo_data)

information_content <- calculate_information_content(descendants)
```

### Full similarity matrices

Armed with ancestors and IC values we can build pairwise similarity
matrices over the entire ontology. `bixverse` supports three IC-based
measures:

- **Resnik** similarity — the IC of the most informative common
  ancestor.
- **Lin** similarity — Resnik similarity normalised by the IC of the two
  terms.
- **Combined** — the average of the (normalised) Resnik and Lin scores.

All matrices are rescaled to the $\lbrack 0,1\rbrack$ interval.

> **Note**
>
> The full matrix generations are not run during Vignette generation.
> While fast, the GitHub runners are VERY memory limited making this
> impossible.

``` r
resnik_similarity <- calculate_semantic_sim_mat(
  similarity_type = "resnik",
  ancestor_list = ancestors,
  ic_list = information_content
)

resnik_similarity[1:5, 1:5]
```

Switching to the Lin similarity is just a parameter change:

``` r
lin_similarity <- calculate_semantic_sim_mat(
  similarity_type = "lin",
  ancestor_list = ancestors,
  ic_list = information_content
)

lin_similarity[1:5, 1:5]
```

### Similarities for a subset of terms

Computing the full matrix is not always necessary or desirable. If you
only care about a handful of terms you can calculate pairwise
similarities directly, which returns a long-format `data.table` instead
of a dense matrix. Here we draw a random sample of ten terms and use the
combined similarity:

``` r
set.seed(10101L)

random_term_sample <- sort(sample(names(information_content), 10L))

combined_similarity <- calculate_semantic_sim(
  terms = random_term_sample,
  similarity_type = "combined",
  ancestor_list = ancestors,
  ic_list = information_content
)

head(combined_similarity)
#>         term1      term2        sims
#>        <char>     <char>       <num>
#> 1: HP:0000452 HP:0003986 0.003174386
#> 2: HP:0000452 HP:0004921 0.003211837
#> 3: HP:0000452 HP:0009655 0.003211837
#> 4: HP:0000452 HP:0031206 0.003114082
#> 5: HP:0000452 HP:0032411 0.003114082
#> 6: HP:0000452 HP:0034891 0.003114082
```

## Wang similarities

In addition to the IC-based measures above, `bixverse` implements the
[Wang
similarity](https://academic.oup.com/bioinformatics/article/23/10/1274/197095).
This approach operates directly on the DAG structure of the ontology and
lets you assign an edge weight to each relationship type. For ontologies
with multiple relationship types — the Gene Ontology, for instance —
this gives you fine-grained control over how different edges contribute
to similarity.

``` r
# The HPO only has one relationship type; for the GO you would specify
# several, e.g. c(is_a = 0.8, part_of = 0.6)

wang_weights <- setNames(0.8, "parent_of")

wang_similarity <- calculate_wang_sim_mat(
  parent_child_dt = hpo_data,
  weight = wang_weights
)

wang_similarity[1:5, 1:5]
```

As with the IC-based measures, you can restrict the computation to a
subset of terms:

``` r
# needs to be provided as previous chunk is not evaluated
wang_weights <- setNames(0.8, "parent_of")

wang_similarity_subset <- calculate_wang_sim(
  terms = random_term_sample,
  parent_child_dt = hpo_data,
  weight = wang_weights
)

head(wang_similarity_subset)
#>         term1      term2       sims
#>        <char>     <char>      <num>
#> 1: HP:0000452 HP:0003986 0.06945263
#> 2: HP:0000452 HP:0004921 0.11339664
#> 3: HP:0000452 HP:0009655 0.04711510
#> 4: HP:0000452 HP:0031206 0.06484329
#> 5: HP:0000452 HP:0032411 0.07044337
#> 6: HP:0000452 HP:0034891 0.12414544
```
