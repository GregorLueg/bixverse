library(magrittr)
library(zeallot)

## Synthetic data ----

devtools::document()
devtools::load_all()

cPCA_data = create_synthetic_cPCA_data()

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

get_params(bulk_coexp_class, TRUE, TRUE)

devtools::document()

test = bulk_coexp(raw_data = raw_data, meta_data = sample_meta)

mat <- S7::prop(bulk_coexp_class, "raw_data")

feature_meta <- data.table::data.table(
  feature_name = colnames(mat),
  mean_exp = colMeans(mat),
  MAD = matrixStats::colMads(mat),
  var_exp = matrixStats::colVars(mat)
) %>%
  data.table::setorder(-MAD)

hvg = NULL
mad_threshold = NULL
scaling = TRUE
scaling_type = 'robust'

if(is.null(hvg) & is.null(mad_threshold)) {
  hvg = 1
}


if (!is.null(hvg)) {
  no_genes_to_take <-
    ifelse(is.integer(hvg), hvg, ceiling(hvg * ncol(mat)))
  hvg_genes <- feature_meta[1:no_genes_to_take, feature_name]
} else {
  hvg_genes <- feature_meta[MAD >= mad_threshold, feature_name]
}

feature_meta[, hvg := feature_name %in% hvg_genes]

# Process the matrix
matrix_processed <- mat[, hvg_genes]

if (scaling) {
  fun <-
    ifelse(scaling_type == "normal",
           "scale",
           "bixverse::robust_scale"
    )
  matrix_processed <- rlang::eval_tidy(rlang::quo(apply(
    matrix_processed, 1, !!!rlang::parse_exprs(fun)
  )))
}


S7::method(format, bulk_coexp) <- function(x, ...) {
  pre_processed <- !purrr::is_empty(x@processed_data)
  co_exp_method <- if (purrr::is_empty(x@params[["detection_method"]])) {
    'not defined yet'
  } else {
    S7::prop(x, "params")[["detection_method"]]
  }
  outputs_available = !purrr::is_empty(S7::prop(bulk_coexp, "params")[["outputs"]])

  sprintf(
    "Bulk co-expression class object:\n  Pre-processsed: %b\n  Method: %s\n  Outputs available: %b",
    pre_processed,
    co_exp_method,
    outputs_available
  )
}

bulk_coexp_class

print(bulk_coexp_class)

purrr::is_empty(bulk_coexp_class@params[["outputs"]])



bulk_coexp_class = contrastive_pca_processing(bulk_coexp_class, background_mat = background_mat)

?c_pca_plot_alphas

c_pca_plot_alphas(bulk_coexp_class, label_column = 'grp', n_alphas = 15L, max_alpha = 1000)

bulk_coexp_class <- apply_contrastive_pca(bulk_coexp_class, alpha = 1, no_pcs	= 10L)

get_params(bulk_coexp_class, TRUE, TRUE)

get_outputs(bulk_coexp_class)
