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
