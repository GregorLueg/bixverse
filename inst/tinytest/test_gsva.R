# gsva tests -------------------------------------------------------------------

library(magrittr)

set.seed(42L)

no_genes <- 12000L
no_samples <- 500L
no_gs <- 200L

data <- matrix(
  data = rnorm(no_genes * no_samples, sd = 2),
  nrow = no_genes,
  ncol = no_samples
)

rownames(data) <- sprintf("gene_%i", 1:no_genes)
colnames(data) <- sprintf("sample_%i", 1:no_samples)

random_gene_sets <- purrr::map(
  1:no_gs,
  ~ {
    sample(rownames(data), 10)
  }
)
names(random_gene_sets) <- sprintf("gs_%i", 1:no_gs)

rextendr::document()

tictoc::tic()
rust_version <- rs_gsva(
  exp = data,
  gs_list = random_gene_sets,
  tau = 1.0,
  gaussian = TRUE,
  max_diff = TRUE,
  abs_rank = FALSE
)
tictoc::toc()

library(GSVA)

tictoc::tic()
gsvaPar <- gsvaParam(data, random_gene_sets)
gsva.es <- gsva(
  gsvaPar,
  verbose = FALSE,
  BPPARAM = BiocParallel::MulticoreParam()
)
tictoc::toc()

?gsva

rust_version[1:5, 1:5]

gsva.es[1:5, 1:5]

plot(rust_version[1, ], gsva.es[1, ])
