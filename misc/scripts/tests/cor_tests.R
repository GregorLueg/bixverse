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

table(new_meta_data$gtex_subgrp)

samples_to_keep <- new_meta_data[
  gtex_subgrp == "Brain - Caudate (basal ganglia)",
  sample_id
]
samples_to_keep_2 <- new_meta_data[
  gtex_subgrp == "Brain - Frontal Cortex (BA9)",
  sample_id
]
data_1 <- t(d)[samples_to_keep, ]
data_2 <- t(d)[samples_to_keep_2, ]
meta_data_1 <- new_meta_data[gtex_subgrp == "Brain - Caudate (basal ganglia)"]
meta_data_2 <- new_meta_data[gtex_subgrp == "Brain - Frontal Cortex (BA9)"]

cor_test <- bulk_coexp(raw_data = data_1, meta_data = meta_data_1) %>%
  preprocess_bulk_coexp(., mad_threshold = 1) %>%
  cor_module_processing(., cor_method = "spearman")

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
  cor_module_processing(., cor_method = "spearman")

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
    cor_method = "spearman"
  ) %>%
  cor_module_check_res(.) %>%
  cor_module_final_modules(.)

cor_test_2 <- cor_module_final_modules(cor_test_2)
tictoc::toc()

# CoReMo -----------------------------------------------------------------------

devtools::load_all()


tictoc::tic()
cor_test <- bulk_coexp(raw_data = data_1, meta_data = meta_data_1) %>%
  preprocess_bulk_coexp(., mad_threshold = 1) %>%
  cor_module_processing(., cor_method = "spearman")

plot_hvgs(cor_test)

cor_test <- cor_module_check_epsilon(cor_test, rbf_func = "gaussian")

plot_epsilon_res(cor_test)

cor_test <- cor_module_coremo_clustering(
  object = cor_test,
  coremo_params = params_coremo(epsilon = 2)
)

cor_test <- cor_module_coremo_stability(object = cor_test)
tictoc::toc()

devtools::document()

object = cor_test
data_mat <- S7::prop(object, "processed_data")[["processed_data"]]
coremo_params <- get_params(object)[["coremo"]]
chunk_size <- 30L

epsilon <- coremo_params$epsilon

rextendr::document()

devtools::document()

dim(data_mat)

devtools::load_all()

tictoc::tic()
chunk_size <- 15L

total_samples <- seq_len(nrow(data_mat))
no_total_samples <- length(total_samples)
groups <- ceiling(seq_along(total_samples) / chunk_size)
chunks <- split(total_samples, groups)

all_results <- vector(mode = "list", length = length(chunks))

for (i in seq_along(chunks)) {
  indices <- chunks[[i]]
  chunk_res <- rs_coremo_stability(
    data = data_mat,
    indices = indices,
    epsilon = coremo_params$epsilon,
    rbf_type = coremo_params$rbf_func,
    spearman = TRUE
  )

  leave_one_out_clustering <- purrr::map(chunk_res, \(chunk) {
    dist_obj <- create_dist_obj(
      x = chunk,
      size = ncol(data_mat)
    )

    new_tree <- fastcluster::hclust(dist_obj, method = "ward.D")

    clusters <- cutree(new_tree, k = coremo_params$inflection_idx)

    clusters
  })

  all_results[[i]] <- leave_one_out_clustering

  message_txt <- sprintf(
    "Chunk %i out of %i: Processed %i samples out of a total of %i samples.",
    i,
    length(chunks),
    ifelse(chunk_size * i < no_total_samples, chunk_size * i, no_total_samples),
    no_total_samples
  )

  message(message_txt)
}

tictoc::toc()

all_results <- purrr::flatten(all_results)

chunk_res <- rs_coremo_stability(
  data = data_mat,
  indices = 1:25L,
  epsilon = coremo_params$epsilon,
  rbf_type = coremo_params$rbf_func,
  spearman = TRUE
)

leave_one_out_clustering <- purrr::map(chunk_res, \(chunk) {
  dist_obj <- create_dist_obj(
    x = chunk,
    size = ncol(data_mat)
  )

  new_tree <- fastcluster::hclust(dist_obj, method = "ward.D")

  clusters <- cutree(new_tree, k = coremo_params$inflection_idx)

  clusters
})
tictoc::toc()
