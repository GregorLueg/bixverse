# cistarget --------------------------------------------------------------------

## helpers ---------------------------------------------------------------------

#' Helper to process CisTarget results
#'
#' @param cs_ls List. The result list from the Rust wrapper.
#' @param gs_name String. Name of the tested gene set.
#' @param represented_motifs Character vector. The represented motifs in the
#' rankings.
#' @param represented_genes Character vector. The represented genes in the
#' rankings.
#'
#' @return A data.table with the results if there were any significant motifs.
process_cistarget_res <- function(
  cs_ls,
  gs_name,
  represented_motifs,
  represented_genes
) {
  # checks
  checkmate::qassert(gs_name, "S1")
  checkmate::assertList(cs_ls)
  checkmate::assertNames(
    names(cs_ls),
    must.include = c(
      "motif_idx",
      "nes",
      "auc",
      "rank_at_max",
      "n_enriched",
      "leading_edge"
    )
  )
  checkmate::qassert(represented_motifs, "S+")
  checkmate::qassert(represented_genes, "S+")

  # early return
  if (length(cs_ls$nes) == 0) {
    return(NULL)
  }

  gs_res <- data.table::data.table(
    gs_name = gs_name,
    motif = represented_motifs[cs_ls$motif_idx],
    nes = cs_ls$nes,
    auc = cs_ls$auc,
    rank_at_max = cs_ls$rank_at_max,
    n_enriched = cs_ls$n_enriched,
    leading_edge_genes = purrr::map_chr(
      cs_ls$leading_edge,
      \(leading_edge_idx) {
        paste(represented_genes[leading_edge_idx], collapse = ";")
      }
    )
  ) %>%
    data.table::setorder(., -nes)

  gs_res
}

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
    annot_file,
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
      direct_annotation                       , "directAnnotation"                       ,
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

## main ------------------------------------------------------------------------

#' Main function to run CisTarget
#'
#' @description
#' The `bixverse` implementation of the RCisTarget workflow, one of the
#' algorithms used in SCENIC, see Aibar, et al. You will need motif to target
#' gene rankings, see [bixverse::read_motif_ranking()] and the motif to TF
#' annotations, see [bixverse::read_motif_annotation_file()].
#'
#' @param gs_list Named list of character vectors. Each element is a gene set
#' containing gene identifiers that must match row names in `rankings`.
#' @param rankings Integer matrix. Motif rankings for genes. Row names are gene
#' identifiers, column names are motif identifiers. Lower values indicate
#' higher regulatory potential.
#' @param annot_data data.table. Motif annotation database mapping motifs to
#' transcription factors. Must contain columns: `motif`, `TF`, and
#' `annotationSource`.
#' @param cis_target_params List. Output of [bixverse::params_cistarget()]:
#' \itemize{
#'   \item{auc_threshold - Numeric. Proportion of genes to use for AUC
#'   threshold calculation. Default 0.05 means top 5 percent of genes.}
#'   \item{nes_threshold - Numeric. Normalised Enrichment Score threshold for
#'   determining significant motifs. Default is 3.0.}
#'   \item{rcc_method - Character. Recovery curve calculation method: "approx"
#'   (approximate, faster) or "icistarget" (exact, slower).}
#'   \item{high_conf_cats - Character vector. Annotation categories considered
#'   high confidence (e.g., "directAnnotation", "inferredBy_Orthology").}
#'   \item{low_conf_cats - Character vector. Annotation categories considered
#'   lower confidence (e.g., "inferredBy_MotifSimilarity").}
#' }
#'
#' @return data.table with enriched motifs and corresponding statistics and
#' high & low confidence TFs for each gene set.
#'
#' @references Aibar, et al., Nat Methods, 2017
#'
#' @export
run_cistarget <- function(
  gs_list,
  rankings,
  annot_data,
  cis_target_params = params_cistarget()
) {
  # checks
  checkmate::assertList(gs_list, names = "named", types = "character")
  checkmate::assertMatrix(
    rankings,
    mode = "integer",
    row.names = "named",
    col.names = "named"
  )
  checkmate::assertDataTable(annot_data)
  checkmate::assertNames(
    names(annot_data),
    must.include = c("motif", "TF", "annotationSource")
  )
  assertCistargetParams(cis_target_params)

  # function body
  annot_red <- annot_data[motif %in% colnames(rankings)]
  gs_indices <- purrr::map(gs_list, \(gene) {
    which(rownames(rankings) %in% gene)
  })

  no_represented_genes <- purrr::map_dbl(gs_indices, length)

  if (
    !all(
      c(cis_target_params$low_conf_cats, cis_target_params$high_conf_cats) %in%
        annot_red$annotationSource
    )
  ) {
    warning("Not all of the high and low confidence categories were found")
  }

  if (any(no_represented_genes == 0)) {
    warning("Some of the gene sets have zero overlap with the ranking.")
  }

  rs_res <- with(
    cis_target_params,
    rs_cistarget(
      rankings = rankings,
      gs_list = gs_indices,
      auc_threshold = as.integer(auc_threshold * nrow(rankings)),
      nes_threshold = nes_threshold,
      max_rank = nrow(rankings),
      method = "approx",
      n_mean = 10L
    )
  )

  cis_res <- purrr::map2(
    .x = rs_res,
    .y = names(gs_indices),
    process_cistarget_res,
    represented_motifs = colnames(rankings),
    represented_genes = rownames(rankings)
  ) %>%
    purrr::keep(
      .,
      ~ {
        !is.null(.x)
      }
    ) %>%
    data.table::rbindlist()

  tf_high <- with(
    cis_target_params,
    annot_red[
      annotationSource %in% high_conf_cats,
      .(TF_highConf = paste(sort(unique(TF)), collapse = ";")),
      by = motif
    ]
  )

  tf_low <- with(
    cis_target_params,
    annot_red[
      annotationSource %in% low_conf_cats,
      .(TF_lowConf = paste(sort(unique(TF)), collapse = ";")),
      by = motif
    ]
  )

  cis_res_final <- cis_res %>%
    merge(tf_high, by = "motif", all.x = TRUE) %>%
    merge(tf_low, by = "motif", all.x = TRUE)

  setorder(cis_res_final, gs_name, -nes)

  return(cis_res_final)
}
