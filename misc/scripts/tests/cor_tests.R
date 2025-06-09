library(devtools)
library(ggplot2)
library(magrittr)

rextendr::document()
# devtools::document()
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

params_coremo(epsilon = 2)

cluster_data <- cor_module_coremo_stability(object = cor_test)

object <- cor_test

if (purrr::is_empty(S7::prop(object, "processed_data")[["processed_data"]])) {
  warning("No pre-processed data found. Defaulting to the raw data.")
  data_mat <- S7::prop(object, "raw_data")
} else {
  data_mat <- S7::prop(object, "processed_data")[["processed_data"]]
}


final_modules <- S7::prop(object = object, name = "outputs")[["final_modules"]]

cluster_mat <- do.call(cbind, cluster_data) %>% `rownames<-`(colnames(data_mat))

cluster_mat_red <- cluster_mat[final_modules$gene, ]

stability <- rs_cluster_stability(cluster_mat_red)

final_modules[, c("stability", "std_stability") := stability]

hist(stability$mean_jaccard)

rextendr::document()

stability_test <- apply(cluster_mat_red, 1, FUN = function(x) {
  sort(table(x), decreasing = TRUE)[[1]]
})

cluster_mat[1:10, 1:10]

tictoc::toc()

devtools::document()

object <- cor_test


plot_df <- data.table::copy(S7::prop(object = object, name = "outputs")[[
  "optimal_cuts"
]]) %>%
  data.table::setorder(gradient_change)

if (is.null(plot_df)) {
  warning(paste(
    "No optimal_cuts data.table found.",
    "Did you run cor_module_coremo_clustering()?",
    "Returning NULL."
  ))
}

optimal_cuts <- S7::prop(object = object, name = "params")[["coremo"]][[
  "inflection_idx"
]]

ggplot2::ggplot(
  data = plot_df,
  mapping = ggplot2::aes(x = k, y = R2_weighted_median)
) +
  ggplot2::geom_point(
    mapping = ggplot2::aes(fill = gradient_change),
    shape = 21,
    alpha = 0.7,
    size = 3
  ) +
  ggplot2::scale_fill_viridis_c() +
  ggplot2::theme_bw() +
  ggplot2::xlab("k cuts") +
  ggplot2::ylab("Median of median weighted R2") +
  ggplot2::labs(fill = "Gradient change") +
  ggplot2::geom_vline(
    xintercept = optimal_cuts,
    color = "darkgrey",
    linetype = "dashed"
  ) +
  ggplot2::ggtitle("k cuts vs. change in R2")
