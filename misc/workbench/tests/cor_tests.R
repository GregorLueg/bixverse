library(devtools)
library(ggplot2)
library(magrittr)

devtools::document()
devtools::load_all()
rextendr::document()
devtools::check()

?S7::new_generic

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

tictoc::tic()
cor_test = bulk_coexp(t(d)[samples_to_keep, ], new_meta_data[gtex_subgrp == "Brain - Hippocampus"]) %>%
  preprocess_bulk_coexp(., mad_threshold = 1) %>%
  cor_module_processing(., correlation_method = 'spearman') %>%
  cor_module_check_res(.)

get_resolution_res(cor_test)
tictoc::toc()

devtools::load_all()

cor_test <- cor_module_check_res(cor_test)

get_params(cor_test)

S7::prop(x, "params")[['correlation_graph']][['no_nodes']]

x <- cor_test

x

