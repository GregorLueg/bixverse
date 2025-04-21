library("anndata")
library(data.table)
library(magrittr)

ad <- read_h5ad("~/Desktop/geo_data/GSE65832/final/GSE154773_anndata.h5ad")

ad$var

counts <- t(ad$X)

devtools::load_all()
devtools::document()

obs_original = ad$obs

h5_file <- "~/Desktop/geo_data/GSE154773_anndata.h5ad"

h5_path <- h5_file

h5_obj <- anndata_parser$new(h5_path)

h5_obj$.__enclos_env__$private$h5_content


if (.verbose) message("Loading data from the h5ad object")
c(meta_data, var_info, counts)
h5_contentc(meta_data, var_info, counts) %<-% h5_obj$get_key_data()


h5_obj$get_var_info()

bulk_dge(raw_counts = counts, meta_data = meta_data, variable_info = var_info)

object <- bulk_dge_from_h5ad(h5_file)

devtools::load_all()
devtools::document()

object <- calculate_pca_bulk_dge(object)

plot_pca_res(object)

get_metadata(object)

# DGE code ----

meta_data <- get_metadata(test)

colnames(meta_data)

dge_list <- edgeR::DGEList(counts)
dge_list <- edgeR::calcNormFactors(dge_list, method = "TMM")
model_matrix <- model.matrix(~Sample_characteristics_ch1_4, data = meta_data)
voom_obj <- limma::voom(counts = dge_list, design = model_matrix, plot = TRUE)
limma_fit <- limma::lmFit(voom_obj, model_matrix)


all_contrasts <- function(group, delim = "vs") {
  group <- unique(as.character(group))
  cb <- combn(group, 2, FUN = function(x) {
    paste0(sprintf("`%s`", x[1]), "-", sprintf("`%s`", x[2]))
  })
  contrasts <- limma::makeContrasts(contrasts = cb, levels = group)
  colnames(contrasts) <- gsub("-", delim, colnames(contrasts))
  return(contrasts)
}

group <- meta_data$Sample_characteristics_ch1_4

group <- sprintf("`%s`", unique(as.character(group)))

cb <- combn(group, 2, FUN = function(x) {
  paste0(sprintf("`%s`", x[1]), "-", sprintf("`%s`", x[2]))
})

contrasts <- limma::makeContrasts(contrasts = cb, levels = group)

all_contrasts(c("a", "b", "c"))


limma_fit <- limma::eBayes(limma_fit)

targ_coefs = setdiff(colnames(model_matrix), "(Intercept)")

dge_res <- limma::topTable(
  limma_fit,
  coef = targ_coefs,
  number = Inf,
  sort.by = "t",
  genelist = rownames(counts)
) %>%
  setDT() %>%
  .[, ensembl_id := entrez_ensembl[ID]] %>%
  setorder(logFC)


dge_res[ensembl_id == "ENSG00000153563"]


head(coef(limma_fit))
