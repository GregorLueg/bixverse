# Load in and extend the Rust documentation...
rextendr::document()
devtools::document()
devtools::load_all()
# devtools::check()

# devtools::install()

library(magrittr)

protein_coding_genes = data.table::fread("~/Desktop/protein_coding_genes.csv")

seed = 123
set.seed(seed)

universe = protein_coding_genes$id
gene_sets_no = 5000
target_gene_sets_no = 50

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

go_data_s7 = gene_ontology_data(go_data_dt, min_genes = 3L)

?GSE_GO_elim_method

x = GSE_GO_elim_method(go_data_s7, target_genes, min_genes = 3L, fdr_threshold = 1)

?GSE_GO_elim_method

head(x)

go_info = S7::prop(go_data_s7, "go_info")
go_to_genes = S7::prop(go_data_s7, "go_to_genes")
ancestry = S7::prop(go_data_s7, "ancestry")
levels = S7::prop(go_data_s7, "levels")

gene_universe_length = length(unique(unlist(go_to_genes)))

tictoc::tic()
results_go = rs_gse_geom_elim(
  target_genes = target_genes,
  go_to_genes = go_to_genes,
  ancestors = ancestry,
  levels = levels,
  gene_universe_length = gene_universe_length,
  min_genes = 3,
  elim_threshold = 0.1,
  debug = FALSE
)
tictoc::toc()

results_go_dt = data.table(do.call(cbind, results_go[-1])) %>%
  .[, `:=`(go_id = results_go$go_ids,
           FDR = p.adjust(pvals, method = 'BH'))] %>%
  merge(.,
        go_info,
        by = 'go_id') %>%
  data.table::setcolorder(
    .,
    c(
      'go_id',
      'go_name',
      'odds_ratios',
      'pvals',
      'FDR',
      'hits',
      'gene_set_lengths'
    )
  )


tictoc::tic()
results_go_list = rs_gse_geom_elim_list(
  target_genes = target_gene_sets,
  go_to_genes = go_to_genes,
  ancestors = ancestry,
  levels = levels,
  gene_universe_length = gene_universe_length,
  min_genes = 3,
  elim_threshold = 0.1,
  debug = FALSE
)
tictoc::toc()

devtools::document()

GSE_GO_elim_method_list(
  S7_obj = go_data_s7,
  target_gene_list = target_gene_sets
)





test =

plot(
  combined$gene_set_lengths.x,
  combined$gene_set_lengths.y
)


plot(
  combined$pvals.x,
  combined$pvals.y
)
