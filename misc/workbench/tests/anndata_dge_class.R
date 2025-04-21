library("anndata")
library(data.table)
library(magrittr)

library(reticulate)

reticulate::py_list_packages()

reticulate::py_install('anndata')

ad <- read_h5ad("~/Desktop/geo_data/GSE154773_anndata.h5ad")

all(rownames(ad$var) == h5_obj$get_var_info()$var_id)

counts <- t(ad$X)

devtools::load_all()
devtools::document()

obs_original = ad$obs

h5_file <- "~/Desktop/geo_data/GSE154773_anndata.h5ad"

h5_path <- h5_file

h5_obj <- anndata_parser$new(h5_file)


if (.verbose) message("Loading data from the h5ad object")
c(meta_data, var_info, counts)
h5_contentc(meta_data, var_info, counts) %<-% h5_obj$get_key_data()


h5_obj$get_var_info()

bulk_dge(raw_counts = counts, meta_data = meta_data, variable_info = var_info)

devtools::load_all()
devtools::document()


object <- bulk_dge_from_h5ad(h5_file) %>%
  change_gene_identifier(alternative_gene_id = "Symbol") %>%
  calculate_pca_bulk_dge()

plot_pca_res(object)

# DGE code ----

genes_to_filter <- object@variable_info[GeneType == 'protein-coding', Symbol]
meta_data <- object@meta_data
dge_list <- object@outputs$dge_list
main_contrast <- 'contrast_info'
co_variates <- NULL
source_col <- 'sample_source'
.verbose <- TRUE

unique(meta_data$sample_source)


meta_data_red <- meta_data[sample_source == "whole blood"]
dge_list_red <- dge_list[genes_to_filter, meta_data_red$sample_id]
variables <- c(main_contrast, co_variates)

class(dge_list)

edgeR::DGEList()

run_limma_voom <- function(
  meta_info,
  main_contrast,
  dge_list,
  co_variates = NULL,
  ...,
  .verbose = TRUE
) {
  variabales <- c(main_contrast, co_variates)
  # Checks
  checkmate::assertDataFrame(meta_info)
  checkmate::qassert(main_contrast, "S1")
  checkmate::assertClass(dge_list, "DGEList")
  checkmate::qassert(co_variates, c("S+", "0"))
  checkmate::assertNames(
    names(meta_info),
    must.include = variabales
  )
  checkmate::qassert(.verbose, "B1")

  # Fix any names
  meta_info[,
    (variables) := lapply(.SD, fix_contrast_names),
    .SDcols = variables
  ]
  if (.verbose)
    message(paste(
      "Fixing any naming issues for the selected main contrast",
      "and any co-variates."
    ))

  model_matrix <- model.matrix(as.formula(model_formula), data = meta_data_red)
  colnames(model_matrix) <- gsub(main_contrast, "", colnames(model_matrix))

  # Filter lowly expressed genes
  to_keep <- edgeR::filterByExpr(
    y = dge_list,
    design = model_matrix,
    ...
  )
  if (.verbose) message(sprintf("A total of %i genes are kept.", sum(to_keep)))

  voom_obj <- limma::voom(
    counts = dge_list[to_keep, ],
    design = model_matrix,
    normalize.method = "quantile",
    plot = TRUE
  )
  limma_fit <- limma::lmFit(voom_obj, model_matrix)

  contrasts <- all_contrasts(limma_fit)

  final_fit <- limma::contrasts.fit(limma_fit, contrasts)
  final_fit <- limma::eBayes(final_fit)

  tested_contrasts <- attributes(contrasts)$dimnames$Contrasts

  all_dge_res <- purrr::map(tested_contrasts, \(coef) {
    top.table <- as.data.table(
      limma::topTable(fit = fit2, coef = coef, sort.by = "P", n = Inf),
      keep.rownames = TRUE
    ) %>%
      setnames(old = 'rn', new = 'gene_id') %>%
      .[, contrast := gsub("-", "_vs_", coef)]
  }) %>%
    rbindlist()

  return(all_dge_res)
}


run_limma_voom(
  meta_info = meta_data_red,
  main_contrast = main_contrast,
  dge_list = dge_list_red
)

meta_data_red[,
  (variables) := lapply(.SD, fix_contrast_names),
  .SDcols = variables
]
if (.verbose)
  message(paste(
    "Fixing any naming issues for the selected main contrast",
    "and any co-variates."
  ))

model_formula <- sprintf("~ 0 + %s", paste(variables, collapse = " + "))

model_matrix <- model.matrix(as.formula(model_formula), data = meta_data_red)
colnames(model_matrix) <- gsub(main_contrast, "", colnames(model_matrix))

to_keep <- edgeR::filterByExpr(
  dge_list_red,
  design = model_matrix,
  ...
)


voom_obj <- limma::voom(
  counts = dge_list_red[to_keep, ],
  design = model_matrix,
  normalize.method = "quantile",
  plot = TRUE
)
limma_fit <- limma::lmFit(voom_obj, model_matrix)

contrasts <- all_contrasts(limma_fit)

fit2 <- limma::contrasts.fit(limma_fit, contrasts)
fit2 <- limma::eBayes(fit2)

tested_contrasts <- attributes(contrasts)$dimnames$Contrasts


?limma::topTable
