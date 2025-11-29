# cell cycle scoring -----------------------------------------------------------

## speed improvement -----------------------------------------------------------

dir_data <- path.expand("~/Downloads/single_cell_data/")

sc_object <- single_cell_exp(
  dir_data = dir_data
)

sc_object <- load_existing(sc_object)

data(cell_cycle_genes)

cell_cycle_genes_ls <- split(
  cell_cycle_genes$ensembl_gene_id,
  cell_cycle_genes$set
)

devtools::check(vignettes = FALSE)

tictoc::tic()
module_scores <- module_scores_sc(
  object = sc_object,
  gs_list = cell_cycle_genes_ls,
  streaming = TRUE,
  .verbose = TRUE
)
tictoc::toc()
