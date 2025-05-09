library(magrittr)
library(zeallot)

## Synthetic data ----

devtools::document()
devtools::load_all()
rextendr::document()

cPCA_data <- synthetic_cPCA_data()

cPCA_data$target[1:5, 1:5]

## Class tests ----

raw_data <- t(cPCA_data$target)
background_mat <- t(cPCA_data$background)

sample_meta <- data.table(
  sample_id = rownames(raw_data),
  grp = cPCA_data$target_labels,
  case_control = "case"
)

bulk_coexp_class <- bulk_coexp(raw_data = raw_data, meta_data = sample_meta)

bulk_coexp_class <- preprocess_bulk_coexp(bulk_coexp_class)


# If this is run without pre-processing it will throw a warning
bulk_coexp_class <- contrastive_pca_processing(
  bulk_coexp_class,
  background_mat = background_mat
)


c_pca_plot_alphas(
  bulk_coexp_class,
  label_column = "grp",
  n_alphas = 10L,
  max_alpha = 100
)

bulk_coexp_class <- apply_contrastive_pca(
  bulk_coexp_class,
  alpha = 2.5,
  no_pcs = 10L
)
