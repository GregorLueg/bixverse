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

rextendr::clean()
rextendr::document()
devtools::document()
devtools::check()
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
ica_test = ica_evaluate_comp(
  ica_test,
  ica_type = 'logcosh',
  ica_params = list(
    maxit = 200L,
    alpha = 1.0,
    max_tol = 0.0001,
    verbose = FALSE
  ),
  iter_params = list(
    cross_validate = FALSE,
    random_init = 50L,
    folds = 10L
  )
)

?ica_evaluate_comp

get_params(ica_test, TRUE, TRUE)

object = ica_test
ica_type = 'logcosh'
ncomp_params = list(
  max_comp = 100L,
  steps = 5L,
  custom_seq = NULL
)
iter_params = list(
  bootstrap = FALSE,
  random_init = 50L,
  folds = 10L
)
ica_params = list(
  maxit = 200L,
  alpha = 1.0,
  max_tol = 0.0001,
  verbose = FALSE
)
random_seed = 42L
.verbose = TRUE

detection_method <- S7::prop(object, "params")[["detection_method"]]

X_raw <- S7::prop(object, "processed_data")[['processed_data']]
X1 <- S7::prop(object, "processed_data")[["X1"]]
K <- S7::prop(object, "processed_data")[["K"]]

n_comp_vector <- if (is.null(ncomp_params$custom_seq)) {
  with(ncomp_params, c(2, 3, 4, seq(
    from = 5, to = max_comp, by = steps
  )))
} else {
  ncomp_params$custom_seq
}

if (.verbose)
  message(sprintf("Using a total of %i different n_comp parameters", length(n_comp_vector)))

# Set up the loop
if (iter_params$bootstrap) {
  total_randomisations <- iter_params$random_init * iter_params$folds
  no_ica_runs <- total_randomisations * length(n_comp_vector)
  if (.verbose)
    message(
      sprintf(
        "Using bootstrapping with %i folds and %i random initialisations for a total of %i ICA runs",
        iter_params$folds,
        iter_params$random_init,
        no_ica_runs
      )
    )
} else {
  total_randomisations <- iter_params$random_init
  no_ica_runs <- iter_params$random_init * length(n_comp_vector)
  if (.verbose)
    message(
      sprintf(
        "Using %i random initialisations for a total of %i ICA runs",
        iter_params$random_init,
        no_ica_runs
      )
    )
}

all_scores <- c()
all_convergence <- c()

pb = txtProgressBar(initial = 0, max = length(n_comp_vector), style = 3)

for (i in seq_along(n_comp_vector)) {
  no_comp <- n_comp_vector[[i]]
  # Get the combined S matrix and convergence information
  c(s_combined, converged) %<-% with(iter_params, switch(
    as.integer(iter_params$bootstrap) + 1,
    rs_ica_iters(
      x_processed = X1,
      k = K,
      no_comp = no_comp,
      no_random_init = random_init,
      ica_type = ica_type,
      random_seed = random_seed,
      ica_params = ica_params
    ),
    rs_ica_iters_cv(
      x_raw = X_raw,
      no_comp = no_comp,
      no_folds = folds,
      no_random_init = random_init,
      ica_type = ica_type,
      random_seed = random_seed,
      ica_params = ica_params
    )
  ))

  c(stability_scores, centrotype) %<-% .community_stability(
    no_comp = as.integer(no_comp),
    s = s_combined,
    return_centrotype = FALSE
  )

  all_scores <- append(all_scores, sort(stability_scores, decreasing = TRUE))
  all_convergence <- append(all_convergence, converged)

  setTxtProgressBar(pb, i)
}

close(pb)

ica_comps_rep <- unlist(purrr::map(n_comp_vector, \(x) {
  rep(x, x)
}))
ica_comps_no <- unlist(purrr::map(n_comp_vector, \(x) {
  seq_len(x)
}))

convergence_split <- split(all_convergence, ceiling(seq_along(all_convergence)/total_randomisations))
names(convergence_split) <- n_comp_vector

ica_stability_res <- list(
  component_rank = ica_comps_no,
  no_components = ica_comps_rep,
  stability = all_scores
) %>% data.table::setDT()

prop_converged <- purrr::imap_dfr(convergence_split, \(bool, x) {
  total_converged <- sum(bool) / length(bool)
  data.table(no_components = as.integer(x), converged = total_converged)
})

ica_res_sum = ica_stability_res[, .(mean_stability = mean(stability)), no_components] %>%
  merge(., prop_converged, by = 'no_components')



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

