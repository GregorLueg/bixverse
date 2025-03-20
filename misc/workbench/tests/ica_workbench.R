# Check ICA stuff ---

## Test real data ----

library(devtools)
library(ggplot2)
library(magrittr)
library(zeallot)
devtools::document()

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

tictoc::tic()
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
tictoc::toc()

plot_ica_stability(ica_test)

# Write a plotting function

plot_df <- ica_test@outputs$ica_stability_res

p1 <- ggplot(data = plot_df,
             mapping = aes(x = component_rank, y = stability)) +
  geom_line(mapping = aes(color = factor(no_components)), linewidth = 1) +
  scale_color_viridis_d(option = "C") +
  theme_minimal() +
  labs(color = "No ICAs") +
  ylim(0, 1) +
  xlab("Component rank") +
  ylab("Stability index")


p2 <- ggplot(data = plot_df, aes(x = component_rank, y = stability)) +
  geom_point(mapping = aes(colour = factor(no_components))) +
  scale_color_viridis_d(option = "C") +
  theme_minimal() +
  labs(color = "No ICAs") +
  ylim(0, 1) +
  xlab("Component rank") +
  ylab("Stability index")

p1 + p2 + patchwork::plot_annotation(
  title = "Stability of independent components",
  subtitle = "Over different ncomps and randomisations"
)


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

