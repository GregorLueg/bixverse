library("anndata")
library(data.table)
library(magrittr)

ad <- read_h5ad("~/Desktop/geo_data/GSE65832/final/GSE65832_anndata.h5ad")

counts <- t(ad$X)

h5_file <- "~/Desktop/geo_data/GSE65832/final/GSE65832_anndata.h5ad"

h5_content <- rhdf5::h5ls(
  h5_file
) %>%
  setDT()


obs_grps <- h5_content[group %like% "/obs/"] %>%
  .[, full_path := paste(group, name, sep = "/")]

categories_paths <- obs_grps[name == "categories"]
codes_paths <- obs_grps[name == "codes"]

rhdf5::h5read(file = h5_file, name = file.path(group, "feature_meta"))

test_obj <- bulk_dge(raw_counts = counts, obs)
