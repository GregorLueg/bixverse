library(devtools)
library(ggplot2)
library(magrittr)

rextendr::document()
devtools::document()
devtools::load_all()

# Test on real data ------------------------------------------------------------

gtex_brain <- recount3::create_rse_manual(
  project = "BRAIN",
  project_home = "data_sources/gtex",
  organism = "human",
  annotation = "gencode_v29",
  type = "gene"
)

coldata <- SummarizedExperiment::colData(gtex_brain) |> as.data.frame()
rowdata <- SummarizedExperiment::rowData(gtex_brain) |> as.data.frame()

d <- edgeR::DGEList(SummarizedExperiment::assay(gtex_brain))
d <- edgeR::calcNormFactors(d, method = "upperquartile")
to_keep <- suppressWarnings(edgeR::filterByExpr(d))
d <- d[to_keep, ]
d <- edgeR::cpm(d, log = TRUE)

d <- as.matrix(d)

new_meta_data <- data.table::data.table(
  sample_id = rownames(coldata),
  case_control = "case",
  gtex_subgrp = coldata$gtex.smtsd
)

samples_to_keep <- new_meta_data[
  gtex_subgrp == "Brain - Hippocampus",
  sample_id
]
samples_to_keep_2 <- new_meta_data[gtex_subgrp == "Brain - Amygdala", sample_id]
data_1 <- t(d)[samples_to_keep, ]
data_2 <- t(d)[samples_to_keep_2, ]
meta_data_1 <- new_meta_data[gtex_subgrp == "Brain - Hippocampus"]
meta_data_2 <- new_meta_data[gtex_subgrp == "Brain - Amygdala"]

cor_test <- bulk_coexp(raw_data = data_1, meta_data = meta_data_1) %>%
  preprocess_bulk_coexp(., mad_threshold = 1) %>%
  cor_module_processing(., correlation_method = "spearman")

devtools::document()

cor_test <- cor_module_check_epsilon(cor_test, rbf_func = "gaussian")

plot_epsilon_res(cor_test)

options(future.globals.maxSize = 2000 * 1024^2)
cor_test <- cor_module_check_res(
  cor_test,
  graph_params = params_cor_graph(epsilon = 1.5)
)

plot_resolution_res(cor_test)

devtools::load_all()

cor_test <- cor_module_final_modules(cor_test)

cor_test <- bulk_coexp(raw_data = data_2, meta_data = meta_data_2) %>%
  preprocess_bulk_coexp(., mad_threshold = 1) %>%
  cor_module_processing(., correlation_method = "spearman")

# Differential correlation -----------------------------------------------------

devtools::load_all()
devtools::document()

dim(data_1)

tictoc::tic()
cor_test_2 <- bulk_coexp(raw_data = data_1, meta_data = meta_data_1) %>%
  preprocess_bulk_coexp(., mad_threshold = 1) %>%
  diffcor_module_processing(
    .,
    background_mat = data_2,
    correlation_method = "spearman"
  ) %>%
  cor_module_check_res(.) %>%
  cor_module_final_modules(.)

cor_test_2 <- cor_module_final_modules(cor_test_2)
tictoc::toc()

# write TOM into Rust ----------------------------------------------------------

# Claude explanation
set.seed(123)
n <- 100

adj_matrix <- abs(cor(X))[1:300, 1:300]

dim(adj_matrix)

# Set diagonal to 0 (no self-connections)
diag(adj_matrix) <- 0

adj <- adj_matrix
TOMDenom <- "min"

tictoc::tic()
n <- nrow(adj)
TOM <- matrix(0, nrow = n, ncol = n)

connectivity <- rowSums(adj)

# Calculate TOM for each pair of nodes
for (i in 1:n) {
  for (j in 1:n) {
    if (i != j) {
      # Calculate the numerator: aij + sum(aik * akj) for k != i,j
      numerator <- adj[i, j]

      # Sum over all k except i and j
      shared_neighbors <- 0
      for (k in 1:n) {
        if (k != i && k != j) {
          shared_neighbors <- shared_neighbors + adj[i, k] * adj[k, j]
        }
      }
      numerator <- numerator + shared_neighbors

      # Calculate the denominator: f(ki, kj) + 1 - aij
      if (TOMDenom == "min") {
        f_ki_kj <- min(connectivity[i], connectivity[j])
      } else if (TOMDenom == "mean") {
        f_ki_kj <- (connectivity[i] + connectivity[j]) / 2
      }

      denominator <- f_ki_kj + 1 - adj[i, j]

      # Calculate TOM
      TOM[i, j] <- numerator / denominator
    }
  }
}

tictoc::toc()


tictoc::tic()
rs_res = rs_tom(adj)
tictoc::toc()

rextendr::document()

summary(c(rs_res))

rs_res[1:10, 1:10]

TOM[1:10, 1:10]
