# cistarget --------------------------------------------------------------------

## data ------------------------------------------------------------------------

#' Download CisTarget reference files for human (hg38)
#'
#' @param cache_dir String. Directory to cache the files. Defaults to a
#' package-specific user cache directory.
#' @param overwrite Logical. Re-download even if files already exist.
#'
#' @returns Named list with paths: `rankings` and `motif_annotations`.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # fetch the hg38 rankings and motif annotations (large, network)
#' paths <- download_cistarget_hg38()
#' paths$rankings
#' }
download_cistarget_hg38 <- function(
  cache_dir = tools::R_user_dir("bixverse", which = "cache"),
  overwrite = FALSE
) {
  checkmate::qassert(cache_dir, "S1")
  checkmate::qassert(overwrite, "B1")

  urls <- list(
    rankings = paste0(
      "https://resources.aertslab.org/cistarget/databases/",
      "homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/",
      "hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"
    ),
    motif_annotations = paste0(
      "https://resources.aertslab.org/cistarget/motif2tf/",
      "motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"
    )
  )

  if (!dir.exists(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE)
  }

  paths <- purrr::imap(urls, \(url, name) {
    dest <- file.path(cache_dir, basename(url))
    if (!file.exists(dest) || overwrite) {
      message(sprintf("Downloading %s ...", name))
      utils::download.file(url, destfile = dest, mode = "wb")
    } else {
      message(sprintf("Using cached %s: %s", name, dest))
    }
    dest
  })

  paths
}

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
#' @returns A data.table with the results if there were any significant motifs.
#'
#' @keywords internal
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
#' `https://resources.aertslab.org/cistarget/motif2tf/`.
#' The function will generate a data.table that can be subsequently used.
#'
#' @param annot_file String. Path to the motif2tf file that you downloaded.
#'
#' @returns data.table with the motif to transcription factor information.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # motif to transcription factor table from the downloaded reference
#' paths <- download_cistarget_hg38()
#' annot <- read_motif_annotation_file(paths$motif_annotations)
#' head(annot)
#' }
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
      inferred_orthology & inferred_motif_sim  ,
      "inferredBy_MotifSimilarity_n_Orthology" ,
      inferred_motif_sim                       ,
      "inferredBy_MotifSimilarity"             ,
      inferred_orthology                       ,
      "inferredBy_Orthology"                   ,
      direct_annotation                        ,
      "directAnnotation"                       ,
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
#' `https://resources.aertslab.org/cistarget/databases/`
#'
#' @param ranking_file String. The file path to the .feather file
#'
#' @returns An integer matrix that has been transposed for easier use in the
#' underlying Rust code.
#'
#' @export
#'
#' @importFrom magrittr %>%
#'
#' @examples
#' \dontrun{
#' # transposed motif rankings from the downloaded feather file
#' paths <- download_cistarget_hg38()
#' rankings <- read_motif_ranking(paths$rankings)
#' dim(rankings)
#' }
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
#' @param .verbose Boolean. Controls verbosity of the function.
#'
#' @returns data.table with enriched motifs and corresponding statistics and
#' high & low confidence TFs for each gene set.
#'
#' @references Aibar, et al., Nat Methods, 2017
#'
#' @export
#'
#' @examples
#' # motif enrichment against a tiny synthetic ranking database
#' rankings <- matrix(
#'   c(1L, 5L, 4L, 2L, 2L, 2L, 1L, 5L, 4L, 3L, 2L, 1L, 3L, 1L, 5L, 4L,
#'     5L, 4L, 3L, 3L),
#'   nrow = 5,
#'   byrow = TRUE,
#'   dimnames = list(sprintf("gene_%i", 1:5), sprintf("motif_%i", 1:4))
#' )
#' annot <- data.table::data.table(
#'   motif = sprintf("motif_%i", 1:4),
#'   TF = sprintf("TF%i", 1:4),
#'   annotationSource = factor(c(
#'     "directAnnotation",
#'     "inferredBy_Orthology",
#'     "inferredBy_MotifSimilarity",
#'     "inferredBy_MotifSimilarity_n_Orthology"
#'   ))
#' )
#' res <- run_cistarget(
#'   gs_list = list(set_a = c("gene_1", "gene_2", "gene_3")),
#'   rankings = rankings,
#'   annot_data = annot,
#'   cis_target_params = params_cistarget(
#'     auc_threshold = 1,
#'     nes_threshold = 0.2
#'   ),
#'   .verbose = FALSE
#' )
#' res[, c("gs_name", "motif", "nes")]
run_cistarget <- function(
  gs_list,
  rankings,
  annot_data,
  cis_target_params = params_cistarget(),
  .verbose = TRUE
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
      max_rank = min(max_rank, nrow(rankings)),
      method = rcc_method,
      n_mean = n_mean,
      verbose = .verbose
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

# binarisation -----------------------------------------------------------------

#' Binarise regulon activity into on/off calls
#'
#' @description
#' The last SCENIC step. Each regulon gets its own threshold derived from the
#' shape of its AUC distribution across cells, and a cell counts as on when its
#' score sits strictly above that threshold.
#'
#' The thresholds come back alongside the calls, so you can inspect them,
#' override any that look wrong and re-apply the comparison yourself. SCENIC
#' does the same, it writes them to an editable file between scoring and
#' assignment.
#'
#' Regulons flagged as not bimodal fell back to `mean + 2 * sd`. If most of your
#' regulons land there, the AUC distributions are too flat to separate, which
#' usually points at the scoring statistic rather than the thresholding. Check
#' you used `auc_type = "recovery"`, see [params_sc_aucell()].
#'
#' @param auc_matrix Numeric matrix of cells x regulons, or the `ScMatrixRes`
#' returned by [aucell_sc()].
#' @param binarise_params List. Output of [params_scenic_binarise()].
#' @param .verbose Boolean. Controls verbosity of the function.
#'
#' @returns A list with:
#' \itemize{
#'   \item{binary - Logical matrix of cells x regulons, `TRUE` where the
#'   regulon is on.}
#'   \item{thresholds - data.table with one row per regulon, holding the
#'   `threshold`, whether it was `bimodal` and the number of cells called on.}
#' }
#'
#' @references Aibar, et al., Nat Methods, 2017
#'
#' @export
#'
#' @examples
#' # on/off calls for a bimodal and a unimodal regulon
#' set.seed(7L)
#' auc <- cbind(
#'   regulon_a = c(rnorm(50, 0.1, 0.02), rnorm(50, 0.4, 0.02)),
#'   regulon_b = rnorm(100, 0.2, 0.05)
#' )
#' rownames(auc) <- sprintf("cell_%i", 1:100)
#' binarise_regulon_activity(auc, .verbose = FALSE)$thresholds
binarise_regulon_activity <- function(
  auc_matrix,
  binarise_params = params_scenic_binarise(),
  .verbose = TRUE
) {
  # ScMatrixRes is a classed matrix, so this just strips the attributes
  if (inherits(auc_matrix, "ScMatrixRes")) {
    auc_matrix <- unclass(auc_matrix)
    attr(auc_matrix, "cell_indices") <- NULL
  }

  # checks
  checkmate::assertMatrix(
    auc_matrix,
    mode = "numeric",
    min.rows = 2,
    min.cols = 1,
    col.names = "named"
  )
  assertScenicBinariseParams(binarise_params)
  checkmate::qassert(.verbose, "B1")

  storage.mode(auc_matrix) <- "double"

  res <- rs_regulon_thresholds(
    auc_matrix = auc_matrix,
    binarise_params = binarise_params
  )

  binary <- sweep(auc_matrix, 2L, res$thresholds, FUN = ">")

  thresholds <- data.table::data.table(
    regulon = colnames(auc_matrix),
    threshold = res$thresholds,
    bimodal = res$bimodal,
    n_cells_on = colSums(binary)
  )

  if (.verbose) {
    message(sprintf(
      "Binarised %d regulons, %d of which were bimodal.",
      ncol(auc_matrix),
      sum(res$bimodal)
    ))
  }

  return(list(binary = binary, thresholds = thresholds))
}

# binary heatmaps --------------------------------------------------------------

#' Extract plot-ready data for a binary heatmap
#'
#' @description
#' Filters, orders and optionally bins a logical samples x features matrix,
#' e.g. the regulon on/off calls from [binarise_regulon_activity()], so it can
#' be drawn as a single raster with `bixverse.plots::plot_binary_heatmap()`.
#'
#' Features are clustered within their group on the Jaccard distance, so
#' shared absences do not count as similarity. Samples are clustered within
#' their group on the Hamming distance. Groups larger than `max_cluster_n` are
#' instead ordered by barycentre, the mean plot position of the features a
#' sample has on, which avoids the quadratic distance matrix. With more
#' samples than `max_cols`, consecutive samples within a group are collapsed
#' into roughly `max_cols` bins that hold the fraction of samples on. Every
#' group keeps at least one bin.
#'
#' @param binary_mat Logical matrix of samples x features with unique row and
#' column names.
#' @param sample_groups Optional named character vector or factor mapping
#' every row name of `binary_mat` to a group. Factor levels set the group
#' order, otherwise groups are sorted.
#' @param feature_groups Optional named character vector or factor mapping
#' every column name of `binary_mat` to a group. Factor levels set the group
#' order, otherwise groups are sorted.
#' @param heatmap_params List. Output of [params_binary_heatmap()].
#' @param .verbose Boolean. Controls verbosity of the function.
#'
#' @returns An object of class `BinaryHeatmapData`, a list with:
#' \itemize{
#'   \item{mat - Numeric matrix of features x columns in plot order. Values
#'   are 0/1, or the fraction of samples on if binned.}
#'   \item{col_annot - data.table with one row per column: `col_idx`, `group`
#'   (factor) and `n_samples` in that column.}
#'   \item{row_annot - data.table with one row per feature: `row_idx`,
#'   `feature` and `group` (factor).}
#'   \item{binned - Boolean. Whether the samples were binned.}
#' }
#' Without groups, `group` is a single level `"all"`.
#'
#' @references Aibar, et al., Nat Methods, 2017
#'
#' @export
#'
#' @examples
#' set.seed(7L)
#' binary_mat <- matrix(runif(200 * 30) > 0.7, nrow = 200)
#' dimnames(binary_mat) <- list(
#'   sprintf("cell_%i", 1:200),
#'   sprintf("regulon_%i", 1:30)
#' )
#' sample_groups <- setNames(rep(c("a", "b"), each = 100), rownames(binary_mat))
#' res <- extract_binary_heatmap_data(
#'   binary_mat,
#'   sample_groups = sample_groups,
#'   .verbose = FALSE
#' )
#' dim(res$mat)
extract_binary_heatmap_data <- function(
  binary_mat,
  sample_groups = NULL,
  feature_groups = NULL,
  heatmap_params = params_binary_heatmap(),
  .verbose = TRUE
) {
  # checks
  checkmate::assertMatrix(
    binary_mat,
    mode = "logical",
    any.missing = FALSE,
    min.rows = 1L,
    min.cols = 1L,
    row.names = "unique",
    col.names = "unique"
  )
  .assert_group_map(sample_groups, rownames(binary_mat))
  .assert_group_map(feature_groups, colnames(binary_mat))
  assertBinaryHeatmapParams(heatmap_params)
  checkmate::qassert(.verbose, "B1")

  # filter
  frac_on <- colMeans(binary_mat)
  keep <- frac_on >= heatmap_params$min_frac_on &
    frac_on <= heatmap_params$max_frac_on
  if (!any(keep)) {
    stop("No features pass the min_frac_on / max_frac_on filter.")
  }
  if (.verbose) {
    message(sprintf(
      "Kept %d of %d features after the on-fraction filter.",
      sum(keep),
      length(keep)
    ))
  }
  binary_mat <- binary_mat[, keep, drop = FALSE]

  sample_groups <- .as_group_factor(sample_groups, rownames(binary_mat))
  feature_groups <- .as_group_factor(feature_groups, colnames(binary_mat))

  # ordering
  feature_ord <- .order_within_groups(feature_groups, \(idx) {
    if (!heatmap_params$cluster_features || length(idx) < 3L) {
      return(idx)
    }
    d <- stats::dist(t(binary_mat[, idx, drop = FALSE]), method = "binary")
    idx[fastcluster::hclust(d, method = "ward.D2")$order]
  })

  # inverse permutation, i.e. plot position of each feature
  feature_pos <- order(feature_ord)

  sample_ord <- .order_within_groups(sample_groups, \(idx) {
    if (!heatmap_params$cluster_samples || length(idx) < 3L) {
      return(idx)
    }
    sub_mat <- binary_mat[idx, , drop = FALSE]
    if (length(idx) <= heatmap_params$max_cluster_n) {
      storage.mode(sub_mat) <- "integer"
      d <- stats::as.dist(rs_hamming_dist(t(sub_mat)))
      return(idx[fastcluster::hclust(d, method = "ward.D2")$order])
    }
    n_on <- rowSums(sub_mat)
    # samples with nothing on give NaN and go last
    barycentre <- as.vector(sub_mat %*% feature_pos) / n_on
    idx[order(barycentre, n_on, na.last = TRUE)]
  })

  ordered_mat <- binary_mat[sample_ord, feature_ord, drop = FALSE]
  storage.mode(ordered_mat) <- "double"
  ordered_groups <- sample_groups[sample_ord]

  # binning
  n_samples <- nrow(ordered_mat)
  binned <- n_samples > heatmap_params$max_cols

  if (binned) {
    grp_n <- as.vector(table(ordered_groups))
    # n_bins <= grp_n, so every bin id within a group gets hit
    n_bins <- pmax(1L, round(grp_n / n_samples * heatmap_params$max_cols))
    offsets <- cumsum(c(0, utils::head(n_bins, -1L)))
    bin_id <- unlist(
      purrr::map(seq_along(grp_n), \(i) {
        ceiling(seq_len(grp_n[i]) * n_bins[i] / grp_n[i]) + offsets[i]
      }),
      use.names = FALSE
    )
    bin_size <- as.vector(table(bin_id))
    ordered_mat <- rowsum(ordered_mat, bin_id, reorder = TRUE) / bin_size
    rownames(ordered_mat) <- NULL
    col_groups <- ordered_groups[!duplicated(bin_id)]
  } else {
    bin_size <- rep(1L, n_samples)
    col_groups <- ordered_groups
  }

  if (.verbose && binned) {
    message(sprintf(
      "Binned %d samples into %d columns.",
      n_samples,
      nrow(ordered_mat)
    ))
  }

  res <- list(
    mat = t(ordered_mat),
    col_annot = data.table::data.table(
      col_idx = seq_along(col_groups),
      group = unname(col_groups),
      n_samples = as.integer(bin_size)
    ),
    row_annot = data.table::data.table(
      row_idx = seq_along(feature_ord),
      feature = colnames(binary_mat)[feature_ord],
      group = unname(feature_groups[feature_ord])
    ),
    binned = binned
  )
  class(res) <- "BinaryHeatmapData"

  return(res)
}

## helpers ---------------------------------------------------------------------

#' Assert an optional group mapping
#'
#' @param groups Optional named character vector or factor.
#' @param ids Character vector. Ids that all need a group.
#'
#' @returns Invisibly `TRUE`, errors otherwise.
#'
#' @keywords internal
.assert_group_map <- function(groups, ids) {
  checkmate::qassert(ids, "S+")
  if (is.null(groups)) {
    return(invisible(TRUE))
  }
  checkmate::assert(
    checkmate::checkCharacter(groups, any.missing = FALSE),
    checkmate::checkFactor(groups, any.missing = FALSE)
  )
  checkmate::assertNames(names(groups), must.include = ids)
  invisible(TRUE)
}

#' Align a group mapping to ids as a factor
#'
#' @param groups Optional named character vector or factor. `NULL` puts
#' everything into a single group `"all"`.
#' @param ids Character vector. Ids to align to.
#'
#' @returns Factor of the same length as `ids`, without unused levels.
#'
#' @keywords internal
.as_group_factor <- function(groups, ids) {
  checkmate::qassert(ids, "S+")
  checkmate::assert(
    checkmate::checkNull(groups),
    checkmate::checkCharacter(groups),
    checkmate::checkFactor(groups)
  )
  if (is.null(groups)) {
    return(factor(rep("all", length(ids))))
  }
  droplevels(as.factor(groups[ids]))
}

#' Order indices group by group
#'
#' @param groups Factor. Groups in the original order.
#' @param order_fun Function. Takes the integer indices of one group and
#' returns them reordered.
#'
#' @returns Integer vector of indices, groups in level order.
#'
#' @keywords internal
.order_within_groups <- function(groups, order_fun) {
  checkmate::assertFactor(groups)
  checkmate::assertFunction(order_fun)
  unlist(
    purrr::map(levels(groups), \(lvl) order_fun(which(groups == lvl))),
    use.names = FALSE
  )
}
