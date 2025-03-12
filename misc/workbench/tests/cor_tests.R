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

x <- get_results(cor_test_2)

table(x$cluster_id)

cor_test <- diffcor_module_processing(bulk_coexp = cor_test_2, background_mat = data_2, correlation_method = 'spearman')

cor_test@params$detection_method


x <- cor_test_2
background_data <- data_2
correlation_method = 'spearman'
.verbose = TRUE

if (purrr::is_empty(S7::prop(x, "processed_data")[['processed_data']])) {
  warning("No pre-processed data found. Defaulting to the raw data")
  target_mat <- S7::prop(x, "raw_data")
} else {
  target_mat <- S7::prop(x, "processed_data")[['processed_data']]
}

spearman <- if (correlation_method == 'pearson') {
  if (.verbose)
    message("Using Pearson correlations.")
  FALSE
} else {
  if (.verbose)
    message("Using Spearman correlations.")
  TRUE
}

features <- colnames(target_mat)

background_mat <- background_data[, features]

combined_mad_df <- list(
  feature = features,
  MAD_target = matrixStats::colMads(target_mat),
  MAD_background = matrixStats::colMads(background_mat)
) %>% data.table::setDT()


final_features <- features

target_mat <- target_mat[, final_features]
background_mat <- background_mat[, final_features]

diff_cor_res <- rs_differential_cor(target_mat, background_mat, spearman = TRUE)

test <- upper_triangle_diffcor_mat$new(diff_cor_res = diff_cor_res, features = final_features)

x <- test$get_data_table()

head(x)

pvals <- x$p_val

length(pvals)

tictoc::tic()
fdr_r <- p.adjust(pvals, method = 'fdr')
tictoc::toc()

tictoc::tic()
fdr_rust <- rs_fdr_adjustment(pvals)
tictoc::toc()

all(fdr_r == fdr_rust)

fdr_rust[1:25]

rextendr::document()

min_cor <- 0.2
fdr_threshold <- 0.05



x_red <- data.table::copy(x)[, delta_cor := cor_a - cor_b] %>%
  .[, `:=`(
    cor_a = abs(cor_a),
    cor_b = abs(cor_b),
    fdr = p.adjust(p_val, method = 'fdr')
  )]

x_red <- x_red[fdr <= fdr_threshold & (cor_a >= min_cor | cor_b >= min_cor)] %>%
  .[, weight := rs_range_norm(abs(delta_cor), max_val = 1, min_val = 0.05)] %>%
  .[, c("feature_a", "feature_b", 'weight')]

vec <- abs(x$z_score)

rextendr::document()

r_norm <- function(x, max_val, min_val) {
  x_max = max(x)
  x_min = min(x)
  x_norm = (x - x_min) / (x_max - x_min) * (max_val - min_val) + min_val
}

tictoc::tic()
vec_norm <- r_norm(vec, 1, 0.05)
tictoc::toc()

vec_norm[1:5]

tictoc::tic()
vec_norm_rs <- rs_min_max(vec, 1, 0.05)
tictoc::toc()

vec_norm_rs[1:5]

head(x)

rextendr::document()

mat <- test$get_cor_matrix(to_ret = "pvals")

mat[1:5, 1:5]

hist(x$z_score)

summary(x$z_score)

head(x)

rust_cor_a <- rs_cor(target_mat, spearman = TRUE)
rust_cor_a.2 <- rs_cor_upper_triangle(target_mat, spearman = TRUE, shift = 0L)
rownames(rust_cor_a) <- colnames(rust_cor_a) <- colnames(target_mat)


a = 0.999980603
b = 1

a.1 = atanh(a)
a.2 = atanh(b)



max(rust_cor_a)

max(rust_cor_a.2)

a <- target_mat[, "ENSG00000130720.12"]
b <- target_mat[, "ENSG00000263690.2"]

plot(a, b)

cor(a, b, method = 'spearman')

plot(rank(a), rank(b))

head(x)

length(unique(purrr::map_dbl(diff_cor_res, length))) == 1
names(diff_cor_res)

plot(MAD_target, MAD_background)
