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

dim(d)

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

plot_resolution_res(cor_test)

?plot_resolution_res

plot_df <- get_resolution_res(cor_test) %>%
  data.table::setorder(., -modularity)



ggplot(data = plot_df,
       mapping =  aes(x = resolution, y = modularity)) +
  geom_point(mapping = aes(size = log10(good_clusters), fill = log10(avg_size)),
             shape = 21,
             alpha = .7) +
  xlab("Leiden cluster resolution") +
  ylab("Modularity") +
  theme_minimal() +
  scale_fill_viridis_c() +
  scale_size_continuous(range = c(2, 8)) +
  labs(size = "Number of good clusters (log10)", fill = "Average cluster size (log10)") +
  ggtitle("Resolution vs. modularity", subtitle = 'With cluster number and size')


vector_to_corr_matrix <- function(vec, n = NULL) {
  # If n (matrix dimension) is not provided, calculate it from vector length
  if (is.null(n)) {
    # For a correlation matrix of size n×n, the upper triangle (excluding diagonal)
    # has n*(n-1)/2 elements
    # Solve quadratic equation: n^2 - n - 2*length(vec) = 0
    n <- (1 + sqrt(1 + 8 * length(vec))) / 2

    # Check if n is an integer
    if (n != floor(n)) {
      stop("Vector length is not compatible with a triangular matrix structure")
    }
    n <- floor(n)  # Ensure n is an integer
  }

  # Initialize correlation matrix with 1s on diagonal
  corr_matrix <- matrix(0, nrow = n, ncol = n)
  diag(corr_matrix) <- 1

  # Fill upper triangle (excluding diagonal)
  idx <- 1
  for (i in 1:(n-1)) {  # Only go up to n-1 for row index
    for (j in (i+1):n) {  # Start from i+1 to skip diagonal
      corr_matrix[i, j] <- vec[idx]
      # For correlation matrices, we fill the lower triangle as well
      corr_matrix[j, i] <- vec[idx]
      idx <- idx + 1
    }
  }

  return(corr_matrix)
}



devtools::document()

rextendr::document()

rextendr::clean()

?upper_triangular_cor_mat

vec <-  rs_cor_upper_triangle(cor_test@processed_data$processed_data, spearman = TRUE, shift = 0L)

test = upper_triangular_cor_mat$new(
  cor_coef = vec,
  features = colnames(cor_test@processed_data$processed_data),
  shift = 0L
)

matrix <- test$get_cor_matrix()

dt_version <- test$get_data_table()

original_size <- ncol(cor_test@processed_data$processed_data)

tictoc::tic()
corr_matrix <- matrix(0, nrow = original_size, ncol = original_size)

idx <- 1

for (i in 1:original_size) {
  for (j in i:original_size) {
    corr_matrix[i, j] <- vec[idx]
    corr_matrix[j, i] <- vec[idx]
    idx <- idx + 1
  }
}
tictoc::toc()

tictoc::tic()
corr_matrix_v2 = rs_upper_triangle_to_dense(vec, shift = 0L, n = original_size)
tictoc::toc()

corr_matrix_v2[1:5, 1:5]

tictoc::tic()
corr_matrix <- vector_to_corr_matrix(example_vec, n = 9414)
tictoc::toc()

dim(corr_matrix)

correlations <- rs_cor(cor_test@processed_data$processed_data, spearman = TRUE)

dim(correlations)



