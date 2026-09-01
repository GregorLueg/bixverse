# Zenodo example datasets ------------------------------------------------------

# Pulls the two SingleCellExperiments the differential expression and
# differential abundance vignettes run on, and writes them out for upload.
#
# Both come from Bioconductor ExperimentHub packages, so provenance and
# versioning are handled upstream. What goes on Zenodo is a convenience cache:
# one file per dataset that the vignettes fetch instead of pulling the whole
# Bioconductor stack at knit time.
#
# Saved as qs2 rather than MTX on purpose. MTX-style output is the counts plus
# barcodes and features, which drops `colData` entirely, and `colData` is where
# the donor and condition columns live. Those are the whole point of both
# datasets, so the format has to preserve them.

## parameters ------------------------------------------------------------------

out_dir <- "~/Desktop/zenodo_data"

# qs2 compression. 3 is the package default and already close to the knee of
# the size/time curve for these.
compress_level <- 3L

# Overwrite files that are already there.
overwrite <- FALSE

# Which of the two thymus assays to take. The droplet one is the larger and the
# one miloR's own vignette works from.
thymus_assay <- c("droplet", "smartseq")[1]

datasets <- c("kang", "thymus")

## libraries -------------------------------------------------------------------

library(SingleCellExperiment)
library(SummarizedExperiment)
library(data.table)

stopifnot(
  "qs2 is needed to write the output" = requireNamespace("qs2", quietly = TRUE),
  "muscData is needed for the Kang data" = requireNamespace(
    "muscData",
    quietly = TRUE
  ),
  "MouseThymusAgeing is needed for the thymus data" = requireNamespace(
    "MouseThymusAgeing",
    quietly = TRUE
  )
)

out_dir <- path.expand(out_dir)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## helpers ---------------------------------------------------------------------

#' Report the shape and the grouping columns of an SCE
#'
#' @param sce SingleCellExperiment.
#' @param cols Character vector. The `colData` columns worth tabulating, i.e.
#' the donor and condition ones the vignettes model over.
describe_sce <- function(sce, cols = character()) {
  cat(sprintf(
    "  %i genes x %i cells | assays: %s\n",
    nrow(sce),
    ncol(sce),
    paste(assayNames(sce), collapse = ", ")
  ))
  for (col in intersect(cols, colnames(colData(sce)))) {
    counts_tab <- table(colData(sce)[[col]], useNA = "ifany")
    cat(sprintf(
      "  %-14s %i level(s): %s\n",
      col,
      length(counts_tab),
      paste(
        sprintf("%s (%i)", names(counts_tab), as.integer(counts_tab)),
        collapse = ", "
      )
    ))
  }
  invisible(NULL)
}

#' Write an SCE out, unless it is already there
#'
#' @param sce SingleCellExperiment.
#' @param name String. File stem.
write_sce <- function(sce, name) {
  path <- file.path(out_dir, sprintf("%s.qs2", name))
  if (file.exists(path) && !overwrite) {
    cat(sprintf(
      "  %s exists, skipping. Set overwrite <- TRUE to replace.\n",
      path
    ))
    return(invisible(path))
  }
  qs2::qs_save(sce, path, compress_level = compress_level)
  cat(sprintf(
    "  wrote %s (%.1f MB)\n",
    path,
    file.size(path) / 1024^2
  ))
  invisible(path)
}

## Kang, et al. 2018 -----------------------------------------------------------

# 8 lupus donors, PBMCs, control vs IFN-beta stimulated in vitro. Paired within
# donor, which is what makes it the NEBULA example: the donor is a genuine
# random effect and every donor contributes to both arms.
#
# Kang, et al., Nat Biotechnol, 2018. GEO: GSE96583.

if ("kang" %in% datasets) {
  cat("Kang18_8vs8\n")
  kang <- muscData::Kang18_8vs8()
  describe_sce(kang, cols = c("ind", "stim", "cell", "multiplets"))
  write_sce(kang, "kang18_8vs8")
}

## Baran-Gale, et al. 2020 -----------------------------------------------------

# Mouse thymic epithelial cells across five ages. Real compositional change with
# age, which is what Milo needs and what the Kang data does not have: in vitro
# stimulation moves cells in embedding space without moving cell type
# proportions, so nearly every neighbourhood would come out significant.
#
# Baran-Gale, et al., Development, 2020.

if ("thymus" %in% datasets) {
  cat("\nMouseThymusAgeing (", thymus_assay, ")\n", sep = "")
  thymus <- switch(
    thymus_assay,
    "droplet" = MouseThymusAgeing::MouseDropletData(),
    "smartseq" = MouseThymusAgeing::MouseSMARTseqData(),
    stop("`thymus_assay` must be 'droplet' or 'smartseq'.")
  )
  describe_sce(thymus, cols = c("Age", "SampID", "SortType", "ClusterAnnot"))
  write_sce(thymus, sprintf("thymus_ageing_%s", thymus_assay))
}

cat("\nDone. Files in:", out_dir, "\n")
