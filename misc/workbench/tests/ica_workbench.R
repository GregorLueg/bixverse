# Check fastICA implementation

library(zeallot)

devtools::document()
rextendr::clean()
rextendr::document()

devtools::check()

syn_data = synthetic_signal_matrix()



meta_data = data.table::data.table(
  sample_id = names(syn_data$group),
  case_control = 'case',
  grp = syn_data$group
)


ica_test = bulk_coexp(X, meta_data)

ica_test = ica_processing(ica_test)

c(X_norm, K) %<-% rs_prepare_whitening(X)

dim(X_norm)

dim(K)

ica_res <- fast_ica_rust(
  X_norm,
  K,
  n_icas = 10L,
  ica_fun = "logcosh",
  seed = 246L,
  .verbose = FALSE
)

# Implement a Rust version that does the approach from the MSTD paper

rextendr::document()
rextendr::clean()

syn_data = synthetic_signal_matrix()

X = t(syn_data$mat) # Artificial count matrix; transposed for rows = samples, cols = features

c(X_norm, K) %<-% rs_prepare_whitening(X)

n_ica_vec <- c(2, 3, 4, 5, 10, 15, 20, 25, 30, 35, 45, 50, 55, 65, 70)

all_scores <- c()

pb = txtProgressBar(initial = 0, max = length(n_ica_vec), style = 3)

stepi = 1

for (no_comp in n_ica_vec) {
  ica_res <- rs_ica_iters(
    x_whiten = X_norm,
    k = K,
    no_comp = no_comp,
    no_iters = no_randomisations,
    maxit = 200L,
    alpha = 1,
    tol = 1e-04,
    ica_type = "logcosh",
    random_seed = 123L,
    verbose = FALSE
  )

  pearson_r <- abs(rs_cor(ica_res$s_combined, spearman = FALSE))

  dist <- as.dist(1 - pearson_r)

  clusters <- hclust(dist)

  clusterCut <- cutree(clusters, no_comp)

  scores <- vector(mode = 'double', length = no_comp)

  for (cluster in seq_len(no_comp)) {
    cluster_indx <- which(clusterCut == cluster)
    not_cluster_indx <- which(clusterCut != cluster)
    within_cluster <- sum(pearson_r[cluster_indx, cluster_indx]) / length(cluster_indx) ^ 2
    outside_cluster <- sum(pearson_r[cluster_indx, not_cluster_indx]) / (length(cluster_indx) * length(not_cluster_indx))
    scores[cluster] <- within_cluster - outside_cluster
  }

  all_scores <- append(all_scores, sort(scores, decreasing = TRUE))

  setTxtProgressBar(pb, stepi)

  stepi <- stepi + 1
}

close(pb)

ica_comps_rep <- unlist(purrr::map(n_ica_vec, \(x) {
  rep(x, x)
}))
ica_comps_no <- unlist(purrr::map(n_ica_vec, \(x) {
  seq_len(x)
}))

ica_stability_res <- list(
  component_rank = ica_comps_no,
  no_components = ica_comps_rep,
  stability = all_scores
) %>% data.table::setDT()

ica_stability_res[, .(mean_stab = mean(stability)), no_components]

library(ggplot2)

ggplot(data = ica_stability_res,
       mapping = aes(x = component_rank,
                     y = stability)) +
  geom_line(mapping = aes(color = factor(no_components)), linewidth = 1) +
  scale_color_viridis_d(option = "C") +
  theme_minimal() +
  labs(color = "No ICAs")

plot(x = ica_stability_res$component_ranl,
     y = ica_stability_res$stability)


# Test real data ----

library(devtools)
library(ggplot2)
library(magrittr)

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
d <- edgeR::calcNormFactors(d, method = 'upperquartile')
to_keep <- suppressWarnings(edgeR::filterByExpr(d))
d <- d[to_keep, ]
d <- edgeR::cpm(d, log = TRUE)

d <- as.matrix(d)

devtools::document()
devtools::load_all()

new_meta_data <- data.table::data.table(
  sample_id = rownames(coldata),
  case_control = 'case',
  gtex_subgrp = coldata$gtex.smtsd
)

samples_to_keep <- new_meta_data[gtex_subgrp == "Brain - Putamen (basal ganglia)", sample_id]
data_1 = t(d)[samples_to_keep, ]
meta_data = new_meta_data[gtex_subgrp == "Brain - Putamen (basal ganglia)"]

ica_test = bulk_coexp(raw_data = data_1, meta_data = meta_data)
ica_test = preprocess_bulk_coexp(ica_test, mad_threshold = 1)
ica_test = ica_processing(ica_test)

# Rust version
tictoc::tic()
c(X_norm, K) %<-% rs_prepare_whitening(X, TRUE, 123L, rank = 100, NULL, NULL)
tictoc::toc()

rextendr::document()

dim(test_2$s_combined)

length(test_2$converged)

max_components <- 100

n_ica_vec <- c(2, 3, 4, seq(from = 5, to = max_components, by = 5))

length(n_ica_vec)

all_scores <- c()
all_convergence <- c()
# all_loadings <- vector(mode = 'list', length = length(n_ica_vec))

tictoc::tic()
pb = txtProgressBar(initial = 0, max = length(n_ica_vec), style = 3)

stepi = 1

?rs_ica_iters()

no_comp = 2L


get_stability_scores <- function(s, ncomp) {

}


c(s_combined, converged) %<-% rs_ica_iters_cv(
  x_raw = X,
  no_comp = no_comp,
  no_folds = 10L,
  no_random_init = 5L,
  ica_type = "exp",
  random_seed = 123L,
  ica_params = list(
    maxit = 200L,
    alpha = 1.0,
    max_tol = 0.0001,
    verbose = FALSE
  )
)

n_comp = 2L
s <- s_combined
return_centrotypes = TRUE

abs_cor <- abs(rs_cor(s, spearman = FALSE))
dist <- as.dist(1 - abs_cor)

clusters <- hclust(dist)

clusterCut <- cutree(tree = clusters, k = no_comp)

scores <- vector(mode = 'double', length = no_comp)

for (cluster in seq_len(no_comp)) {
  cluster_indx <- which(clusterCut == cluster)
  not_cluster_indx <- which(clusterCut != cluster)
  within_cluster <- sum(abs_cor[cluster_indx, cluster_indx]) / length(cluster_indx) ^ 2
  outside_cluster <- sum(abs_cor[cluster_indx, not_cluster_indx]) / (length(cluster_indx) * length(not_cluster_indx))
  scores[cluster] <- within_cluster - outside_cluster
}

cluster <- 2L

cluster_indx <- which(clusterCut == cluster)
not_cluster_indx <- which(clusterCut != cluster)

max_similarity <- which(rowSums(abs_cor[cluster_indx, cluster_indx]) == max(rowSums(abs_cor[cluster_indx, cluster_indx])))

abs_cor[cluster_indx, cluster_indx][1:5, 1:5]



s_combined[, cluster_indx[max_similarity], drop = FALSE]

for (no_comp in n_ica_vec) {
  c(s_combined, converged) %<-% rs_ica_iters_cv(
    x_raw = X,
    no_comp = no_comp,
    no_folds = 10L,
    no_random_init = 5L,
    ica_type = "exp",
    random_seed = 123L,
    ica_params = list(
      maxit = 200L,
      alpha = 1.0,
      max_tol = 0.0001,
      verbose = FALSE
    )
  )

  pearson_r <- abs(rs_cor(ica_res$s_combined, spearman = FALSE))

  dist <- as.dist(1 - pearson_r)

  clusters <- hclust(dist)

  clusterCut <- cutree(clusters, no_comp)

  scores <- vector(mode = 'double', length = no_comp)

  for (cluster in seq_len(no_comp)) {
    cluster_indx <- which(clusterCut == cluster)
    not_cluster_indx <- which(clusterCut != cluster)
    within_cluster <- sum(pearson_r[cluster_indx, cluster_indx]) / length(cluster_indx) ^ 2
    outside_cluster <- sum(pearson_r[cluster_indx, not_cluster_indx]) / (length(cluster_indx) * length(not_cluster_indx))
    scores[cluster] <- within_cluster - outside_cluster
  }

  all_scores <- append(all_scores, sort(scores, decreasing = TRUE))
  all_convergence <- append(all_convergence, ica_res$converged)
  # all_loadings[[stepi]] <- ica_res$s_combined

  setTxtProgressBar(pb, stepi)

  stepi <- stepi + 1
}
tictoc::toc()

sum(n_ica_vec)

length(all_scores)

sum(all_convergence) / length(all_convergence)

ica_comps_rep <- do.call(c, purrr::map(n_ica_vec, \(x) {
  rep(x, x)
}))

ica_comps_no <- do.call(c, purrr::map(n_ica_vec, \(x) {
  seq_len(x)
}))

ica_stability_res <- list(
  component_rank = ica_comps_no,
  no_components = ica_comps_rep,
  stability = all_scores
) %>% data.table::setDT()


ica_res_sum = ica_stability_res[, .(mean_stability = mean(stability), sd_stability = sd(stability)), no_components] %>%
  .[, effect := mean_stability / sd_stability]

ggplot(data = ica_stability_res,
       mapping = aes(x = component_rank,
                     y = stability)) +
  geom_line(mapping = aes(color = factor(no_components)), linewidth = 1) +
  scale_color_viridis_d(option = "C") +
  theme_minimal() +
  labs(color = "No ICAs") +
  ylim(0, 1) +
  xlab("Component rank") +
  ylab("Stability index")


ggplot(ica_stability_res,
       aes(x = component_rank,
           y = stability)) +
  geom_point(mapping = aes(colour = factor(no_components))) +
  scale_color_viridis_d(option = "C") +
  theme_minimal() +
  labs(color = "No ICAs") +
  ylim(0, 1) +
  xlab("Component rank") +
  ylab("Stability index")


all_loadings_mat <- purrr::reduce(all_loadings, cbind)

## hc method -------

pearson_r <- abs(rs_cor(all_loadings_mat, spearman = FALSE))

dist <- as.dist(1 - pearson_r)

clusters <- hclust(dist)

clusterCut <- cutree(clusters, max(n_ica_vec))

scores <- vector(mode = 'double', length = max(n_ica_vec))

for (cluster in seq_len(max(n_ica_vec))) {
  cluster_indx <- which(clusterCut == cluster)
  not_cluster_indx <- which(clusterCut != cluster)
  within_cluster <- sum(pearson_r[cluster_indx, cluster_indx]) / length(cluster_indx) ^ 2
  outside_cluster <- sum(pearson_r[cluster_indx, not_cluster_indx]) / (length(cluster_indx) * length(not_cluster_indx))
  scores[cluster] <- within_cluster - outside_cluster
}

all_scores <- append(all_scores, sort(scores, decreasing = TRUE))

## graph based method ----------

dim(all_loadings_mat)

pearson_r_upper_triangle = rs_cor_upper_triangle(all_loadings_mat, shift = 1L, spearman = FALSE)

feature_names <- purrr::map2(n_ica_vec, no_randomisations, \(comp, iters) {
  sprintf(
    "IC_%i_ncomp_%i_iter_%i",
    seq_len(comp),
    comp,
    rep(seq_len(no_randomisations), each = comp)
  )
}) %>%
  do.call(c, .)

ica_cor_results = upper_triangular_cor_mat$new(
  cor_coef = pearson_r_upper_triangle,
  features = feature_names,
  shift = 1L
)

.verbose = TRUE
kernel_bandwidth = .15
min_affinity = 0.01

graph_df <- ica_cor_results$get_data_table(.verbose = .verbose) %>%
  .[, cor_abs := abs(cor)] %>%
  .[, dist := 1 - cor_abs] %>%
  .[, dist := data.table::fifelse(dist < 0, 0, dist)] %>%
  .[, affinity := rs_gaussian_affinity_kernel(x = dist, bandwidth = kernel_bandwidth)] %>%
  .[affinity >= min_affinity] %>%
  .[, c("feature_a", "feature_b", "affinity")] %>%
  data.table::setnames(
    .,
    old = c("feature_a", "feature_b", "affinity"),
    new = c("from", "to", "weight")
  )

correlation_examples <- 1 - seq(from = 0, to = 1, length.out = 101)

kernel <- rs_gaussian_affinity_kernel(correlation_examples, bandwidth = .15)

1 - max(correlation_examples[which(kernel > 0.01)])

par(mfcol = c(1, 1))

plot(correlation_examples, kernel)

ica_graph <- igraph::graph_from_data_frame(graph_df, directed = FALSE)
ica_graph <- igraph::simplify(ica_graph)

resolution_params = list(min_res = 0.1,
                         max_res = 10,
                         number_res = 15L)

resolutions <- with(resolution_params, exp(seq(log(min_res), log(max_res), length.out = number_res)))

.temp_workers <- 5L

future::plan(future::multisession(workers = .temp_workers))

modularity_scores <- furrr::future_map_dbl(
  resolutions,
  \(res) {
    set.seed(123L)
    ica_communities <- igraph::cluster_leiden(
      ica_graph,
      objective_function = 'modularity',
      resolution = res,
      n_iterations = 5L
    )

    modularity <- igraph::modularity(x = ica_graph, membership = ica_communities$membership)

    return(modularity)
  },
  .progress = TRUE,
  .options = furrr::furrr_options(seed = TRUE)
)

modularity_scores

final_res = resolutions[which(modularity_scores == max(modularity_scores))]

ica_communities <- igraph::cluster_leiden(
  ica_graph,
  objective_function = 'modularity',
  resolution = final_res,
  n_iterations = 5L
)

table(ica_communities$membership)

sum(n_ica_vec) * 50L

final_data <- list(
  node_name = ica_communities$names,
  membership = ica_communities$membership
) %>% setDT()

final_data[membership == 20]

modularity <- igraph::modularity(x = ica_graph, membership = ica_communities$membership)
modularity

# Is my ICA implementation correct ... ? ---------------------------------------

## FastICA ---------------------------------------------------------------------

S <- cbind(sin((1:1000) / 20), rep((((
  1:200
) - 100) / 100), 5))
A <- matrix(c(0.291, 0.6557, -0.5439, 0.5572), 2, 2)
X <- S %*% A

a <- fastICA::fastICA(
  X,
  2,
  alg.typ = "parallel",
  fun = "logcosh",
  alpha = 1,
  method = "R",
  row.norm = FALSE,
  maxit = 200,
  tol = 0.0001,
  verbose = TRUE
)

par(mfcol = c(2, 3))
plot(
  1:1000,
  S[, 1],
  type = "l",
  main = "Original Signals",
  xlab = "",
  ylab = ""
)
plot(1:1000,
     S[, 2],
     type = "l",
     xlab = "",
     ylab = "")
plot(
  1:1000,
  X[, 1],
  type = "l",
  main = "Mixed Signals",
  xlab = "",
  ylab = ""
)
plot(1:1000,
     X[, 2],
     type = "l",
     xlab = "",
     ylab = "")
plot(
  1:1000,
  a$S[, 1],
  type = "l",
  main = "ICA source estimates",
  xlab = "",
  ylab = ""
)
plot(1:1000,
     a$S[, 2],
     type = "l",
     xlab = "",
     ylab = "")

## Rust ------------------------------------------------------------------------

S <- cbind(sin((1:1000) / 20), rep((((
  1:200
) - 100) / 100), 5))
A <- matrix(c(0.291, 0.6557, -0.5439, 0.5572), 2, 2)
X <- S %*% A

c(X_norm, K) %<-% rs_prepare_whitening(X, TRUE, 123L, NULL, NULL, NULL)

rextendr::document()
rextendr::clean()
devtools::load_all()

?fast_ica_rust

ica_res_rs <- fast_ica_rust(
  X_norm,
  K,
  n_icas = 2L,
  ica_fun = "logcosh",
  seed = 42L,
  ica_params = list("x" = 515)
)

par(mfcol = c(2, 3))
plot(
  1:1000,
  S[, 1],
  type = "l",
  main = "Original Signals",
  xlab = "",
  ylab = ""
)
plot(1:1000,
     S[, 2],
     type = "l",
     xlab = "",
     ylab = "")
plot(
  1:1000,
  X[, 1],
  type = "l",
  main = "Mixed Signals",
  xlab = "",
  ylab = ""
)
plot(1:1000,
     X[, 2],
     type = "l",
     xlab = "",
     ylab = "")
plot(
  1:1000,
  ica_res_rs$S[1, ],
  type = "l",
  main = "ICA source estimates",
  xlab = "",
  ylab = ""
)
plot(1:1000,
     ica_res_rs$S[2, ],
     type = "l",
     xlab = "",
     ylab = "")

