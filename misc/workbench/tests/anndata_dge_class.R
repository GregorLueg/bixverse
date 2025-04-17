library("anndata")
library(data.table)
library(magrittr)

ad <- read_h5ad("~/Desktop/geo_data/GSE65832/final/GSE65832_anndata.h5ad")

ad$var

counts <- t(ad$X)

obs_original = ad$obs

h5_file <- "~/Desktop/geo_data/GSE157194/final/GSE157194_anndata.h5ad"

devtools::document()

devtools::load_all()

test <- bulk_dge_from_h5ad(h5_file)

counts <- test@filtered_counts

dim(counts)

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

group <- unique(as.character(group))

cb <- combn(group, 2, FUN = function(x) {
  paste0(x[1], "-", x[2])
})

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
