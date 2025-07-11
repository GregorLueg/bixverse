# ICA RBH ----------------------------------------------------------------------

library(devtools)
library(ggplot2)
library(magrittr)
library(zeallot)

## data ------------------------------------------------------------------------

devtools::load_all()

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

### data set 1 -----------------------------------------------------------------

samples_to_keep_1 <- new_meta_data[
  gtex_subgrp == "Brain - Putamen (basal ganglia)",
  sample_id
]
data_1 <- t(d)[samples_to_keep_1, ]
meta_data_1 <- new_meta_data[gtex_subgrp == "Brain - Putamen (basal ganglia)"]

ica_test_1 <- bulk_coexp(raw_data = data_1, meta_data = meta_data_1)
ica_test_1 <- preprocess_bulk_coexp(ica_test_1, mad_threshold = 1)

plot_hvgs(ica_test_1)

ica_test_1 <- ica_processing(ica_test_1)

ica_test_1 <- ica_evaluate_comp(ica_test_1, ica_type = "logcosh")

ica_test_1 <- ica_optimal_ncomp(ica_test_1, span = 0.5)

plot_ica_ncomp_params(ica_test_1)

ica_test_1 <- ica_stabilised_results(ica_test_1, ica_type = "logcosh")

ica_results_1 <- get_results(ica_test_1)

s_1 <- t(ica_results_1$S)

s_1[1:5, 1:5]

### data set 2 -----------------------------------------------------------------

samples_to_keep_2 <- new_meta_data[
  gtex_subgrp == "Brain - Spinal cord (cervical c-1)",
  sample_id
]
data_2 <- t(d)[samples_to_keep_2, ]
meta_data_2 <- new_meta_data[
  gtex_subgrp == "Brain - Spinal cord (cervical c-1)"
]

ica_test_2 <- bulk_coexp(raw_data = data_2, meta_data = meta_data_2)
ica_test_2 <- preprocess_bulk_coexp(ica_test_2, mad_threshold = 1)

plot_hvgs(ica_test_2)

ica_test_2 <- ica_processing(ica_test_2)

ica_test_2 <- ica_evaluate_comp(ica_test_2, ica_type = "logcosh")

ica_test_2 <- ica_optimal_ncomp(ica_test_2, span = 0.5)

plot_ica_ncomp_params(ica_test_2)

ica_test_2 <- ica_stabilised_results(ica_test_2, ica_type = "logcosh")

ica_results_2 <- get_results(ica_test_2)

s_2 <- t(ica_results_2$S)

s_2[1:5, 1:5]

rbh_data <- list(
  module_1 = s_1,
  module_2 = s_2
)

rextendr::document()

rs_rbh_cor(rbh_data, FALSE, 0.0)
