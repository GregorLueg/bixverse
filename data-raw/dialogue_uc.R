# DIALOGUE ulcerative colitis example data -------------------------------------

# Prepares the Smillie, et al. ulcerative colitis subset used by the original
# DIALOGUE paper (and by pertpy) for ingestion into `SingleCells`.
#
# Source: https://figshare.com/ndownloader/files/43462662 (pertpy's
# `dialogue_example.h5ad`), itself derived from
# https://github.com/livnatje/DIALOGUE/wiki/Example.
#
# The source `X` holds `log2(TPM / 10 + 1)`, which is what Smillie, et al.
# published. `load_h5ad()` wants counts and applies its own log CPM on the way
# in, so feeding it the published matrix would log twice. We linearise back to
# the TPM/10 scale and round, which gives an integer pseudo-count matrix on a
# ~92k-per-cell scale. The raw counts are not recoverable from what was
# published, and DIALOGUE only ever reads the normalised layer, so this loses
# nothing that matters.
#
# The `CD8+ IL17+` cells are kept here. They are absent from twelve donors and
# have to be dropped before DIALOGUE runs, but that is the vignette's job to
# show rather than something to bake into the artefact.

library(Matrix)
library(data.table)

src <- commandArgs(trailingOnly = TRUE)[1]
out <- commandArgs(trailingOnly = TRUE)[2]

stopifnot(file.exists(src))

# counts -----------------------------------------------------------------------

dims <- rhdf5::h5readAttributes(src, "X")$shape

counts <- new(
  "dgCMatrix",
  x = as.numeric(rhdf5::h5read(src, "X/data")),
  i = as.integer(rhdf5::h5read(src, "X/indices")),
  p = as.integer(rhdf5::h5read(src, "X/indptr")),
  Dim = as.integer(dims)
)

counts@x <- round(2^counts@x - 1)
stopifnot(all(counts@x >= 1))

cell_ids <- as.character(rhdf5::h5read(src, "obs/index"))
gene_ids <- as.character(rhdf5::h5read(src, "var/index"))

dimnames(counts) <- list(cell_ids, gene_ids)

# obs --------------------------------------------------------------------------

# categoricals come back as a `categories` / 0-based `codes` pair
read_cat <- function(name) {
  x <- rhdf5::h5read(src, paste0("obs/", name))
  as.character(x$categories[x$codes + 1L])
}

obs <- data.table(
  cell_id = cell_ids,
  cell_subtype = read_cat("cell.subtypes"),
  sample = read_cat("sample"),
  clinical_status = read_cat("clinical.status"),
  pathology = as.logical(rhdf5::h5read(src, "obs/pathology")),
  cell_q = as.numeric(rhdf5::h5read(src, "obs/cellQ")),
  location = read_cat("location"),
  gender = read_cat("gender"),
  subset = read_cat("subset")
)

var <- data.table(
  gene_id = gene_ids,
  gene_symbol = as.character(rhdf5::h5read(src, "var/name"))
)

# write ------------------------------------------------------------------------

bixverse::write_h5ad_sc(
  f_path = out,
  counts = counts,
  obs = obs,
  var = var,
  overwrite = TRUE
)

message(sprintf("Wrote %s (%s cells x %s genes)", out, nrow(obs), nrow(var)))
