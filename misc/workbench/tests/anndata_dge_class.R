library("anndata")
library(data.table)
library(magrittr)

ad <- read_h5ad("~/Desktop/geo_data/GSE65832/final/GSE65832_anndata.h5ad")

counts <- t(ad$X)

obs_original = ad$obs

h5_file <- "~/Desktop/geo_data/GSE65832/final/GSE65832_anndata.h5ad"

h5_class = anndata_parser$new(h5_file)

h5_class$get_obs_data()

h5_content <- rhdf5::h5ls(
  h5_file
) %>%
  setDT()

devtools::document()

devtools::load_all()


test <- bulk_dge_from_h5ad(h5_file)


# Obs

obs_index = h5_content[group == "/obs" & otype == "H5I_DATASET"] %>%
  .[, path := paste(group, name, sep = "/")] %>%
  .[, path]

sample_ids <- c(rhdf5::h5read(
  file = h5_file,
  name = obs_index
))

obs_grps <- h5_content[group %like% "/obs/"] %>%
  .[, full_path := paste(group, name, sep = "/")]

categories_paths <- obs_grps[name == "categories", full_path]
codes_paths <- obs_grps[name == "codes", full_path]

categories <- purrr::map(categories_paths, \(path) {
  as.character(rhdf5::h5read(file = h5_file, name = path))
})
codes <- purrr::map(codes_paths, \(path) {
  as.integer(rhdf5::h5read(file = h5_file, name = path)) + 1
})

colnames <- unique(gsub("/obs/", "", obs_grps[["group"]]))

obs = purrr::map2(categories, codes, \(cat, code) {
  factor(cat[code], levels = cat)
}) %>%
  setDT() %>%
  `colnames<-`(colnames) %>%
  .[, sample_ids := sample_ids] %>%
  .[, c("sample_ids", colnames), with = FALSE]

# Vars

var_index = h5_content[group == "/var" & otype == "H5I_DATASET"] %>%
  .[, path := paste(group, name, sep = "/")] %>%
  .[, path]

var_ids <- c(rhdf5::h5read(
  file = h5_file,
  name = var_index
))

# X

X <- rhdf5::h5read(
  file = h5_file,
  name = "/X"
)

dim(t(X))
