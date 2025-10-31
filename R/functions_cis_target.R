# cistarget --------------------------------------------------------------------

## motif to tf annotations -----------------------------------------------------

#' Read in the motif annotation file
#'
#' @description
#' This function loads in the motif2tf information that you can get from
#' \code{https://resources.aertslab.org/cistarget/motif2tf/}.
#' The function will generate a data.table that can be subsequently used.
#'
#' @param annot_file String. Path to the motif2tf file that you downloaded.
#'
#' @returns data.table with the motif to transcription factor information.
#'
#' @export
read_motif_annotation_file <- function(annot_file) {
  # checks
  checkmate::assertFileExists(annot_file)

  # function body
  motif_annotations <- data.table::fread(
    annotFile,
    colClasses = setNames("character", "source_version")
  )

  data.table::setnames(
    motif_annotations,
    old = c("#motif_id", "gene_name"),
    new = c("motif", "TF"),
    skip_absent = TRUE
  )

  motif_annotations[, `:=`(
    direct_annotation = description == "gene is directly annotated",
    inferred_orthology = orthologous_gene_name != "None",
    inferred_motif_sim = similar_motif_id != "None"
  )]

  # fcase is cool!
  motif_annotations[,
    annotationSource := data.table::fcase(
      inferred_orthology & inferred_motif_sim , "inferredBy_MotifSimilarity_n_Orthology" ,
      inferred_motif_sim                      , "inferredBy_MotifSimilarity"             ,
      inferred_orthology                      , "inferredBy_Orthology"                   ,
      direct_annotation                       , "direct_Annotation"                      ,
      default = ""
    )
  ]
  motif_annotations[, annotationSource := factor(annotationSource)]

  selectedColumns <- c(
    "motif",
    "TF",
    "direct_annotation",
    "inferred_orthology",
    "inferred_motif_sim",
    "annotationSource",
    "description"
  )
  motif_annotations <- motif_annotations[, ..selectedColumns]

  data.table::setkeyv(motif_annotations, c("motif", "TF"))

  return(motif_annotations)
}

## prepare motif rankings ------------------------------------------------------

#' Read in the motif rankings and transform them into a matrix
#'
#' @description
#' This function loads in the .feather files with the motif to target gene
#' rankings. These can be found here:
#' \code{https://resources.aertslab.org/cistarget/databases/}
#'
#' @param ranking_file String. The file path to the .feather file
#'
#' @returns An integer matrix that has been transposed for easier use in the
#' underlying Rust code.
#'
#' @export
#'
#' @importFrom magrittr `%>%`
read_motif_ranking <- function(ranking_file) {
  # checks
  checkmate::assertFileExists(ranking_file)

  # transform into matrix for easier use subsequently
  rankings <- setDT(arrow::read_feather(
    ranking_file
  ))

  motifs <- unlist(rankings[, ncol(rankings), with = FALSE])

  transposed <- as.matrix(rankings[,
    -ncol(rankings),
    with = FALSE
  ]) %>%
    t() %>%
    `colnames<-`(motifs)

  return(transposed)
}
