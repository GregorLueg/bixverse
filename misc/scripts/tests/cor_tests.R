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

cor_test <- bulk_coexp(raw_data = data_1, meta_data = meta_data_1) %>%
  preprocess_bulk_coexp(., mad_threshold = 1) %>%
  cor_module_processing(., cor_method = "spearman")

plot_hvgs(cor_test)

cor_test <- cor_module_check_epsilon(cor_test, rbf_func = "gaussian")

plot_epsilon_res(cor_test)

devtools::load_all()

cor_test <- cor_module_coremo_clustering(object = cor_test, epsilon = 2)

checkCoReMoParams(coremo_params)

x <- coremo_params

test_choice_rules <- list(
  rbf_func = c("gaussian", "inverse_quadratic", "bump"),
  cor_method = c("spearman", "pearson")
)
test_choice_res <- purrr::imap_lgl(x, \(x, name) {
  if (name %in% names(test_choice_rules)) {
    checkmate::testChoice(x, test_choice_rules[[name]])
  } else {
    TRUE
  }
})

qtest_rules <- list(
  k_min = "I1",
  k_max = "I1",
  junk_module_threshold = "N1",
  min_size = c("I1", "0")
)
q_test_res <- purrr::imap_lgl(x, \(x, name) {
  if (name %in% names(qtest_rules)) {
    checkmate::qtest(x, qtest_rules[[name]])
  } else {
    TRUE
  }
})

object = cor_test
epsilon = 2
coremo_params = params_coremo()
.verbose = FALSE
.seed = 10101L

devtools::load_all()

cor_res <- S7::prop(object, "processed_data")$correlation_res
cor_mat <- cor_res$get_cor_matrix()

rextendr::document()

tom_mat <- rs_tom(
  cor_mat,
  "v2",
  TRUE
)

aff_mat <- with(
  params_coremo,
  rs_rbf_function_mat(
    x = 1 - abs(cor_mat),
    epsilon = epsilon,
    rbf_type = "gaussian"
  )
)

dist_mat <- 1 - aff_mat

dist_mat[1:10, 1:10]

tictoc::tic()
tree <- stats::hclust(as.dist(dist_mat), method = "ward.D")
tictoc::toc()

tictoc::tic()
ftree <- fastcluster::hclust(as.dist(dist_mat), method = "ward.D")
tictoc::toc()

class(ftree)

class(tree)

plot(tree, labels = FALSE)

# devtools::load_all()
rextendr::document()

tictoc::tic()
optimal_cuts <- with(
  params_coremo,
  tree_cut_iter(
    tree = ftree,
    cor_mat = cor_mat,
    dist_mat = dist_mat,
    k_min = k_min,
    k_max = k_max,
    min_size = min_size,
    cor_method = cor_method
  )
)
tictoc::toc()


c(inflection_idx, gradient_change) %<-%
  get_inflection_point(
    optimal_cuts$k,
    optimal_cuts$R2_weighted_median,
    span = 0.5
  )

optimal_cuts[, gradient_change := c(0, gradient_change)]

final_clusters <- with(
  params_coremo,
  coremo_tree_cut(
    tree = tree,
    k = as.integer(inflection_idx),
    dist_mat = dist_mat,
    cor_method = cor_method
  )
) %>%
  `names<-`(rownames(cor_mat))

cluster_list <- split(
  names(final_clusters),
  final_clusters
)

final_quality <- rs_coremo_quality(
  cluster_genes = cluster_list,
  cor_mat = cor_mat,
  row_names = rownames(cor_mat),
  seed = .seed
) %>%
  data.table::setDT() %>%
  .[, cluster_id := names(cluster_list)]

junk_modules <- final_quality[r2med <= 0.05, cluster_id]

module_dt <- data.table::as.data.table(
  stack(cluster_list)
) %>%
  data.table::setnames(
    old = c("values", "ind"),
    new = c("gene", "cluster_id")
  ) %>%
  .[!cluster_id %in% junk_modules] %>%
  merge(., final_quality, by = "cluster_id")


ggplot(data = optimal_cuts, mapping = aes(x = k, y = R2_weighted_median)) +
  geom_point(
    mapping = aes(fill = gradient_change),
    size = 4,
    alpha = .5,
    shape = 21
  ) +
  geom_smooth(method = "loess", se = TRUE, span = 0.25) +
  theme_bw() +
  xlab("k cuts") +
  ylab("Weighted R2 median") +
  scale_fill_viridis_c()

k = 75L
tictoc::tic()
modules <-
  coremo_tree_cut(
    tree = tree,
    k = k,
    min_size = min_size,
    dist_mat = dist_mat,
    cor_method = cor_method
  ) %>%
  `names<-`(rownames(cor_mat))
tictoc::toc()

cluster_list <- split(names(modules), modules)

devtools::load_all()

# rextendr::document()

tictoc::tic()
qc <- coremo_cluster_quality(modules = modules, cor_mat = cor_mat)
tictoc::toc()

tictoc::tic()
qc_2 <- rs_coremo_quality(
  cluster_genes = cluster_list,
  cor_mat = cor_mat,
  row_names = rownames(cor_mat),
  seed = 10101L
)
tictoc::toc()

coremo_cluster_quality_v2 <- function(modules, cor_mat, random_seed = 10101L) {
  # Checks
  checkmate::qassert(modules, c("S+", "I+"))
  checkmate::assertNamed(modules)
  checkmate::assertMatrix(cor_mat, mode = "numeric")
  checkmate::qassert(random_seed, "I1")
  # Function body
  cluster_list <- split(names(modules), modules)

  n_clusters <- length(cluster_list)
  result_list <- vector("list", n_clusters)

  for (i in seq_len(n_clusters)) {
    genes <- cluster_list[[i]]
    n <- length(genes)
    if (n < 2) {
      result_list[[i]] <- data.table::data.table(
        size = n,
        r2med = 1,
        r2mad = 0
      )
    } else {
      # Sample genes if needed
      scluster <- if (n > 1000) {
        sample(genes, 1000, replace = FALSE)
      } else {
        genes
      }

      # Extract submatrix more efficiently
      sub_cor_sq <- cor_mat[scluster, scluster, drop = FALSE]^2

      r2_values <- rs_coremo_quality(sub_cor_sq)

      result_list[[i]] <- data.table::data.table(
        size = n,
        r2med = r2_values$median,
        r2mad = r2_values$mad
      )
    }
  }

  data.table::rbindlist(result_list)
}
