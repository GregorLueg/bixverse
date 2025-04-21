library("anndata")
library(data.table)
library(magrittr)

ad <- read_h5ad("~/Desktop/geo_data/GSE65832/final/GSE65832_anndata.h5ad")

ad$var

counts <- t(ad$X)

devtools::load_all()
devtools::document()

obs_original = ad$obs

h5_file <- "~/Desktop/geo_data/GSE154773_anndata.h5ad"

h5_path <- h5_file

h5_obj <- anndata_parser$new(h5_path)
if (.verbose) message("Loading data from the h5ad object")
c(meta_data, var_info, counts) %<-% h5_obj$get_key_data()

bulk_dge(raw_counts = counts, meta_data = meta_data, variable_info = var_info)

object <- bulk_dge_from_h5ad(h5_file)

# PCA code ----

object <- object
top_hgv <- 2500L

mat <- S7::prop(object, "filtered_counts")

dge_list <- edgeR::DGEList(mat)
dge_list <- edgeR::calcNormFactors(dge_list, method = "TMM")

cpm <- edgeR::cpm(dge_list, log = TRUE)

voom_obj <- limma::voom(dge_list, design = NULL, )

hvg_data <- data.table::setDT(
  list(
    gene_id = rownames(cpm),
    mad = matrixStats::rowMads(cpm)
  )
) %>%
  data.table::setorder(-mad)

hvg_genes <- hvg_data[1:top_hgv, gene_id]

input_genes <- t(cpm[hvg_genes, ])

rextendr::document()

tictoc::tic()
pca_results <- prcomp(input_genes, scale. = FALSE)
tictoc::toc()

tictoc::tic()
rs_pca_results <- rs_svd(input_genes, scale = FALSE)
tictoc::toc()


rs_pca_results$scores[1:5, 1:5]


pca_results$x[1:5, 1:5]

x <- input_genes

x <- as.matrix(x)
x <- scale(x, center = TRUE, scale = FALSE)

s <- svd(x, nu = 0)

scores <- x %*% s$v

scores[1:10, 1:10]

pca_results$x[1:10, 1:10]


scores_2 <- x %*% rs_pca_results$v


scores_2[1:10, 1:10]

rs_pca_results$s

test_x <- input_genes %*% rs_pca_results$v

test_x <- rs_pca_results$v %*% diag(rs_pca_results$s)

plot(test_x[1, ], pca_results$x[1, ])

pca_results$x[1:5, 1:5]

test_x[1:5, 1:5]

plot(pca_results$sdev, rs_pca_results$s)

rs_pca_results$v[1:10, 1:10]

pca_results$rotation[1:10, 1:10]

rs_pca_results$s

pca_results$sdev


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
