# Load in and extend the Rust documentation...
rextendr::document()
devtools::document()
devtools::load_all()
# devtools::check()

# devtools::install()

print('x')

library(magrittr)

protein_coding_genes <- data.table::fread("~/Desktop/protein_coding_genes.csv")

seed <- 123
set.seed(seed)

universe <- protein_coding_genes$id
gene_sets_no <- 5000
target_gene_sets_no <- 100

gene_sets <- purrr::map(1:gene_sets_no, ~ {
  set.seed(seed + .x + 1)
  size <- sample(20:100, 1)
  sample(universe, size, replace = FALSE)
})

names(gene_sets) <- purrr::map_chr(1:gene_sets_no, ~ {
  set.seed(seed + .x + 1)
  paste(sample(LETTERS, 3), collapse = "")
})

target_gene_sets <- purrr::map(1:target_gene_sets_no, ~ {
  set.seed(.x * seed)
  size <- sample(50:100, 1)
  sample(universe, size, replace = FALSE)
})

names(target_gene_sets) <- purrr::map_chr(1:target_gene_sets_no, ~ {
  set.seed(seed + .x + 1)
  paste(sample(letters, 3), collapse = "")
})

target_genes <- target_gene_sets[[1]]

rs_set_sim_list(target_genes, target_gene_sets, similarity_index = T)

tictoc::tic()
t1 <- GSE_hypergeometric(
  target_genes = target_genes,
  gene_set_list = gene_sets,
  gene_universe = universe,
  minimum_overlap = 0L,
  threshold = 1
)
tictoc::toc()

tictoc::tic()
t2 <- GSE_hypergeometric_list(
  target_gene_sets,
  gene_set_list = gene_sets,
  minimum_overlap = 0L,
  threshold = 1
)
tictoc::toc()

devtools::document()

go_data_dt <- biomind_to_go_data("~/Desktop/biomind_downloads/processed_data/")

go_data_s7 <- gene_ontology_data(go_data_dt, min_genes = 3L)

print(go_data_s7)


S7_obj <- go_data_s7

print("test")

print(go_data_s7)

go_data_s7

number_levels <- length(S7::prop(S7_obj, "levels"))
number_gene_sets <- length(S7::prop(S7_obj, "go_to_genes"))
min_genes <- S7::prop(S7_obj, "min_genes")

cat(paste(
  "Gene ontology class:",
  sprintf(" Contains %i gene ontology terms.", number_gene_sets),
  sprintf(" Total of %i levels represented in the ontology.", number_levels),
  sprintf(" Minimum genes per term set to %i.", min_genes),
  sep = "\n"
))


?network_diffusions

arrow::write_parquet(go_data_dt, sink = "~/Desktop/go_data.parquet")

tictoc::tic()
t3 <- GSE_GO_elim_method(go_data_s7,
                         target_genes,
                         minimum_overlap = 0L,
                         fdr_threshold = 1)
tictoc::toc()

tictoc::tic()
t4 <- GSE_GO_elim_method_list(
  go_data_s7,
  target_gene_sets,
  minimum_overlap = 0L,
  fdr_threshold = 1
)
tictoc::toc()
