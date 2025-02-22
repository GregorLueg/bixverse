library(magrittr)
library(zeallot)

## Synthetic data ----

devtools::document()
devtools::load_all()
rextendr::document()

cPCA_data = synthetic_cPCA_data()

cPCA_data$target[1:5, 1:5]

## CoVar tests ----

ncol = 10000
nrow = 100

set.seed(123)
x = matrix(rnorm(ncol * nrow), nrow, ncol)

tictoc::tic()
y_1 = cov(x)
tictoc::toc()

tictoc::tic()
y_2 = rs_covariance(x)
tictoc::toc()

## Class tests ----

raw_data = t(cPCA_data$target)
background_mat = t(cPCA_data$background)

sample_meta = data.table(
  sample_id = rownames(raw_data),
  grp = cPCA_data$target_labels,
  case_control = "case"
)

bulk_coexp_class = bulk_coexp(raw_data = raw_data, meta_data = sample_meta)

bulk_coexp_class = preprocess_bulk_coexp(bulk_coexp_class)

# If this is run without pre-processing it will throw a warning
bulk_coexp_class = contrastive_pca_processing(bulk_coexp_class, background_mat = background_mat)

c_pca_plot_alphas(bulk_coexp_class, label_column = 'grp', n_alphas = 10L, max_alpha = 100)

bulk_coexp_class <- apply_contrastive_pca(bulk_coexp_class, alpha = 2.5, no_pcs	= 10L)

# Compare against some old code

scale <- FALSE

target_matrix <- scale(raw_data, scale = scale)
background_matrix <- scale(background_mat, scale = scale)

target_matrix_covar <- coop::covar(target_matrix)
background_matrix_covar <- coop::covar(background_matrix)

.calculate_cPCA <- function(target_covar,
                            background_covar,
                            target_matrix,
                            alpha,
                            nPCs,
                            return_loadings = T) {
  # Checks
  # checkmate::assert(.check_cPCAdim_covar(target_covar, background_covar))
  checkmate::assertMatrix(target_matrix, mode = "numeric")
  checkmate::qassert(alpha, "N1[0,)")
  checkmate::qassert(nPCs, "I1")
  checkmate::qassert(return_loadings, "B1")
  # Final co-variance matrix and calculate cPCAs
  final_covar <- target_covar - alpha * background_covar
  # Calculate the Eigen values
  cPCA_results <- irlba::partial_eigen(final_covar, nPCs)
  print(cPCA_results$vectors)
  # Get loadings and factors
  cPCA_loadings <- cPCA_results$vectors %>%
    `colnames<-`(paste0("cPC_loading", 1:nPCs)) %>%
    `rownames<-`(rownames(target_covar))
  cPCA_factors <- target_matrix %*% cPCA_loadings %>%
    `colnames<-`(paste0("cPC", 1:nPCs))

  if (return_loadings) {
    res <- list(factors = cPCA_factors, loadings = cPCA_loadings)
  } else {
    res <- cPCA_factors
  }

  res
}


test <- .calculate_cPCA(
  target_covar = target_matrix_covar,
  background_covar = background_matrix_covar,
  target_matrix = target_matrix,
  alpha = 2.5,
  nPCs = 2L
)


rs_test = rs_contrastive_pca(
  target_covar = target_matrix_covar,
  background_covar = background_matrix_covar,
  target_mat = target_matrix,
  alpha = 2.5,
  n_pcs = 2L,
  return_loadings = FALSE
)

library(data.table)
library(ggplot2)

test_factors <- data.table::data.table(test$factors)
test_factors[, grp := bulk_coexp_class@meta_data$grp]

head(test_factors)

ggplot(data = test_factors,
       mapping = aes(x = cPC1,
                     y = cPC2)) +
  geom_point(aes(color = grp))

head(test_factors)

