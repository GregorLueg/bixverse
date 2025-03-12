library(devtools)
library(ggplot2)
library(magrittr)

rextendr::document()
devtools::document()
devtools::load_all()
devtools::check()

syn_data = synthetic_signal_matrix()

X = t(syn_data$mat)

meta_data = data.table::data.table(
  sample_id = names(syn_data$group),
  case_control = 'case',
  grp = syn_data$group
)

cor_test = bulk_coexp(X, meta_data)

cor_test = preprocess_bulk_coexp(cor_test)

cor_test = cor_module_processing(cor_test, correlation_method = 'spearman')

rm(list = c("assertCorGraphParams", "checkCorGraphParams"))

devtools::document()

cor_test@processed_data$correlation_res$get_data_table()

cor_test = cor_module_identification(cor_test)

cor_test@outputs$cluser_quality

plot_df <- cor_test@outputs$cluser_quality

head(plot_df)

ggplot(data = plot_df,
       mapping = aes(x = res, y = r_median_of_adjust)) +
  geom_point(mapping = aes(size = log10(median_size), fill = r_weighted_median),
             shape = 21) +
  theme_minimal() +
  xlab("Leiden resolution") +
  ylab("Median adjusted r")

plot_df

# Test on real data ----

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

new_meta_data <- data.table::data.table(
  sample_id = rownames(coldata),
  case_control = 'case',
  gtex_subgrp = coldata$gtex.smtsd
)

samples_to_keep <- new_meta_data[gtex_subgrp == "Brain - Hippocampus", sample_id]
samples_to_keep_2 <- new_meta_data[gtex_subgrp == "Brain - Amygdala", sample_id]
data_1 = t(d)[samples_to_keep, ]
data_2 = t(d)[samples_to_keep_2, ]
meta_data = new_meta_data[gtex_subgrp == "Brain - Hippocampus"]

cor_test = bulk_coexp(raw_data = data_1, meta_data = meta_data)

devtools::load_all()
devtools::document()

tictoc::tic()
cor_test = bulk_coexp(raw_data = data_1, meta_data = meta_data) %>%
  preprocess_bulk_coexp(., mad_threshold = 1) %>%
  cor_module_processing(., correlation_method = 'spearman') %>%
  cor_module_check_res(.) %>%
  cor_module_final_modules(.)
tictoc::toc()

?cor_module_check_res

# Differential correlation -----

devtools::load_all()
devtools::document()

dim(data_1)

tictoc::tic()
cor_test_2 = bulk_coexp(raw_data = data_1, meta_data = meta_data) %>%
  preprocess_bulk_coexp(., mad_threshold = 1) %>%
  diffcor_module_processing(., background_mat = data_2, correlation_method = 'spearman') %>%
  cor_module_check_res(.) %>%
  cor_module_final_modules(.)
tictoc::toc()

cor_test_2 <- cor_module_final_modules(cor_test_2)
