# Graphs and ontologies

Three independent graph classes and two routes into ontology semantic
similarity. None of them talk to each other, pick by the question.

## Network diffusion

Propagate a seed set over a network and see what lights up. Personalised
PageRank under the hood.

```r
diff_obj <- NetworkDiffusions(edge_data, weighted = FALSE, directed = FALSE)

diff_obj <- diffuse_seed_nodes(
  diff_obj,
  diffusion_vector = seed_scores,   # named numeric
  summarisation = "max"             # how to collapse duplicate names
)
```

Add a null so the scores mean something:

```r
diff_obj <- permute_seed_nodes(diff_obj, perm_iters = 1000L)
get_diffusion_vector(diff_obj)
get_diffusion_perms(diff_obj)
```

Benchmark against a gold standard:

```r
auc <- calculate_diffusion_auc(diff_obj, hit_nodes = gold_standard_genes,
                               permutation_test = TRUE)   # gives Z-scores
```

Pull out communities in the high-heat region:

```r
diff_obj <- community_detection(
  diff_obj,
  community_params = params_community_detection(
    min_nodes = 5L,
    threshold_type = "prop_based",   # or "pval_based" with pval_threshold
    min_seed_nodes = 1L
  )
)
res <- get_results(diff_obj)
```

`threshold_type = "pval_based"` needs `permute_seed_nodes()` to have run first,
otherwise there are no p-values to threshold on.

**Tied diffusion**, for two seed sets where you want what both reach:

```r
diff_obj <- tied_diffusion(
  diff_obj,
  diffusion_vector_1 = vec_1,
  diffusion_vector_2 = vec_2,
  summarisation = "max",
  score_aggregation = "min"   # how the two diffusion results combine
)
```

Standalone helpers: `generate_personalisation_vec()`,
`constrained_page_rank()`, `constrained_page_rank_ls()` (list version) and
`knn_graph_label_propagation()`.

## Reciprocal best hit graphs

For finding conserved module correspondence across datasets or species.

```r
obj <- RbhGraph(...)
obj <- generate_rbh_graph(obj, ...)
obj <- find_rbh_communities(obj, ...)
res <- get_rbh_res(obj)
```

## Similarity network fusion

Combining several data modalities into one patient (or sample) similarity
network.

```r
obj <- SimilarityNetworkFusion(data = mat_1, data_name = "rna",
                               snf_params = params_snf())
obj <- add_snf_data_modality(obj, data = mat_2, data_name = "methylation")
obj <- add_snf_data_modality(obj, data = mat_3, data_name = "protein")

obj <- run_snf(obj)

get_snf_final_mat(obj)
get_snf_adjcacency_mat(obj)
get_snf_params(obj)
```

Add every modality before calling `run_snf()`. Note the spelling of
`get_snf_adjcacency_mat()`, it's a typo that's now part of the API.

## Ontology semantic similarity

Two routes. The flat functions take a parent-child data.table and are fine for
one-off work. The `OntologySim` class transfers the structure into Rust once and
is the right choice when you're doing many queries or filtering the result.

### Flat functions

```r
c(ancestors, descendants) %<-% get_ontology_ancestry(parent_child_dt)
ic <- calculate_information_content(descendants)
```

Three IC-based measures: **Resnik** (IC of the most informative common
ancestor), **Lin** (Resnik normalised by the two terms' IC) and **combined**.

```r
# full matrix over the whole ontology
sim_mat <- calculate_semantic_sim_mat(
  similarity_type = "resnik",
  ancestor_list   = ancestors,
  ic_list         = ic
)

# just a subset of terms
sim <- calculate_semantic_sim(
  terms           = term_subset,
  similarity_type = "combined",
  ancestor_list   = ancestors,
  ic_list         = ic
)
```

**Wang** similarity is structure-based rather than IC-based, so it takes the
parent-child table directly plus per-relationship weights:

```r
wang <- calculate_wang_sim(terms = term_subset, parent_child_dt = hpo_data,
                          weight = c(is_a = 0.8, part_of = 0.6))
wang_mat <- calculate_wang_sim_mat(parent_child_dt = hpo_data,
                                   weight = wang_weights)
```

The `%<-%` destructuring comes from `zeallot`, which bixverse imports.

### The class route

```r
obj <- OntologySim(parent_child_dt)
obj <- pre_process_sim_onto(obj)

obj <- calculate_semantic_sim_onto(obj, similarity_type = "resnik")
#  or calculate_wang_sim_onto(obj, ...)

crit <- calculate_critical_value(...)
obj  <- filter_similarities(obj, ...)

get_sim_matrix(obj)
```

`calculate_critical_value()` gives you a data-driven threshold rather than a
guessed one, then `filter_similarities()` sparsifies on it. Worth doing, a full
ontology similarity matrix is large and mostly noise.

Human GO ships with the package (`get_go_data_human()`). For HPO or anything
else, bring your own parent-child table with `parent`, `child` and
`relationship` columns.

Semantic similarity also backs `simplify_hypergeom_res()` in `enrichment.md`,
which is the most common reason to want it.
