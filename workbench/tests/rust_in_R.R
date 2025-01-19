# Load in and extend the Rust documentation...
rextendr::document()
devtools::document()
devtools::load_all()

# devtools::install()

library(magrittr)

protein_coding_genes = data.table::fread("~/Desktop/protein_coding_genes.csv")

seed = 123
set.seed(seed)

universe = protein_coding_genes$id
gene_sets_no = 5000
target_gene_sets_no = 2500

gene_sets = purrr::map(1:gene_sets_no, ~{
  set.seed(seed + .x + 1)
  size = sample(20:100, 1)
  sample(universe, size, replace = FALSE)
})

names(gene_sets) <- purrr::map_chr(1:gene_sets_no, ~{
  set.seed(seed + .x + 1)
  paste(sample(LETTERS, 3), collapse = "")
})

target_gene_sets = purrr::map(1:target_gene_sets_no, ~{
  set.seed(.x * seed)
  size = sample(50:100, 1)
  sample(universe, size, replace = FALSE)
})

names(target_gene_sets) <- purrr::map_chr(1:target_gene_sets_no, ~{
  set.seed(seed + .x + 1)
  paste(sample(letters, 3), collapse = "")
})

target_genes = target_gene_sets[[1]]

tictoc::tic()
test = rs_hypergeom_test(
  target_genes = target_genes,
  gene_sets = gene_sets,
  gene_universe = universe
)
tictoc::toc()

GSE_hypergeometric(
  target_genes = target_genes,
  gene_set_list = gene_sets,
  gene_universe = universe,
  minimum_overlap = 0L,
  threshold = 1
)

tictoc::tic()
test_2 = rs_hypergeom_test_list(
  target_genes = target_gene_sets,
  gene_sets = gene_sets,
  gene_universe = universe
)
tictoc::toc()

GSE_hypergeometric_list(
  target_gene_sets,
  gene_set_list = gene_sets
)

min_genes = 3

go_data_dt = biomind_to_go_data("~/Desktop/biomind_downloads/processed_data/")

go_data_s7 = gene_ontology_enrich_data(go_data_dt, 3L)

go_to_genes = S7::prop(go_data_s7, "go_to_genes")
ancestry = S7::prop(go_data_s7, "ancestry")
levels = S7::prop(go_data_s7, "levels")

gene_universe_length = length(unique(unlist(go_to_genes)))

results_go = rs_gse_geom_elim(
  target_genes = target_genes,
  go_to_genes = go_to_genes,
  ancestors = ancestry,
  levels = levels,
  gene_universe_length = gene_universe_length,
  min_genes = 3,
  elim_threshold = 0,
  debug = TRUE
)

test_1 = data.table(do.call(cbind, results_go[-1])) %>%
  .[, go_id := results_go$go_ids]

test_2 = data.table(do.call(cbind, results_go[-1])) %>%
  .[, go_id := results_go$go_ids]

combined = merge(
  test_1,
  test_2,
  by = 'go_id'
) %>%
  setorder(pvals.x)

head(combined)

plot(
  combined$gene_set_lengths.x,
  combined$gene_set_lengths.y
)


plot(
  combined$pvals.x,
  combined$pvals.y
)
