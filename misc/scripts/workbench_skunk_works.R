# cell cycle scoring -----------------------------------------------------------

library(Seurat)

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes


BiocManager::install("biomaRt")


library(biomaRt)

# Connect to Ensembl
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

ensembl_set1 <- getBM(
  attributes = c("hgnc_symbol", "ensembl_gene_id"),
  filters = "hgnc_symbol",
  values = s.genes,
  mart = mart
) %>%
  as.data.table()

ensembl_set1[, set := "S phase"]

ensembl_set2 <- getBM(
  attributes = c("hgnc_symbol", "ensembl_gene_id"),
  filters = "hgnc_symbol",
  values = g2m.genes,
  mart = mart
) %>%
  as.data.table()

ensembl_set2[, set := "G2/M phase"]

cell_cycle_genes <- rbindlist(list(ensembl_set1, ensembl_set2))

usethis::use_data(my_pkg_data)

usethis::use_data_raw("seurat_cell_cycle")
