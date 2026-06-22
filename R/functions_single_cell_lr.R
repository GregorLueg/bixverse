# ------------------------------------------------------------------------------
# Single cell ligand receptor analysis frameworks
# - NicheNet approach
# ------------------------------------------------------------------------------

# nichenet ---------------------------------------------------------------------

## class -----------------------------------------------------------------------

#' Constructor for ligand-target influence results
#'
#' @description
#' Stores the ligand x gene influence matrix from the NicheNet-style
#' computation, the ligand seed groupings, the gene universe, and the parameters
#' used.
#'
#' @param influence Numeric matrix. Rows are ligand seeds, columns are genes.
#' @param ligand_seeds List of character vectors. Each element is a single
#' ligand or a group treated as one complex.
#' @param ligand_names Character vector. Row names for the influence matrix.
#' @param gene_ids Character vector. Column names for the influence matrix
#' (full gene universe used in the run).
#' @param params List. Parameters used for the run.
#'
#' @returns An object of class `LigandTargetInfluence`.
#'
#' @keywords internal
new_ligand_target_influence <- function(
  influence,
  ligand_seeds,
  ligand_names,
  gene_ids,
  params
) {
  rownames(influence) <- ligand_names
  colnames(influence) <- gene_ids

  res <- list(
    influence = influence,
    ligand_seeds = ligand_seeds,
    ligand_names = ligand_names,
    gene_ids = gene_ids,
    params = params
  )

  class(res) <- "LigandTargetInfluence"
  res
}

### primitives -----------------------------------------------------------------

#' @export
#'
#' @keywords internal
print.LigandTargetInfluence <- function(x, ...) {
  cat("LigandTargetInfluence\n")
  cat(sprintf("  No ligand seeds:    %d\n", nrow(x$influence)))
  cat(sprintf("  No genes:           %d\n", ncol(x$influence)))
  cat(sprintf("  Damping factor:     %.3f\n", x$params$damping_factor))
  cat(sprintf("  Max iter:           %d\n", x$params$max_iter))
  cat(sprintf("  Secondary targets:  %s\n", x$params$secondary_targets))
  invisible(x)
}

#' @export
#'
#' @keywords internal
dim.LigandTargetInfluence <- function(x) {
  dim(x$influence)
}

### getters --------------------------------------------------------------------

#' Get the ligand-target influence matrix
#'
#' @param x An object holding ligand-target influence results.
#'
#' @export
get_influence <- function(x) {
  UseMethod("get_influence")
}

#' @rdname get_influence
#'
#' @export
get_influence.LigandTargetInfluence <- function(x) {
  checkmate::assertClass(x, "LigandTargetInfluence")
  x$influence
}

## ligand gene influence scoring -----------------------------------------------

#' Generate the ligand to target influence matrix
#'
#' @description
#' Computes the NicheNet-style ligand to target gene influence matrix. Builds
#' the gene universe as the union of all symbols across both networks, remaps to
#' 0-indexed integer node IDs for the Rust side, and wraps the resulting matrix
#' in a `LigandTargetInfluence` object.
#'
#' @param ligand_seeds List of character vectors. Each element is one
#' ligand symbol or a group treated as one complex (e.g.
#' `list(c("TGFB1", "TGFB2"))`). Optionally named; if unnamed, row names
#' default to the symbols joined with `+`.
#' @param ppi_network data.table with columns `from`, `to`, `weight`
#' (character, character, numeric). Protein-protein / signalling layer.
#' @param grn_network data.table with columns `from`, `to`, `weight`
#' (character, character, numeric). Gene regulatory layer.
#' @param params List. As returned by [params_ligand_target()].
#'
#' @returns A `LigandTargetInfluence` object.
#'
#' @export
generate_ligand_target_influence <- function(
  ligand_seeds,
  ppi_network,
  grn_network,
  params = params_ligand_target()
) {
  checkmate::assertList(ligand_seeds, types = "character", min.len = 1L)
  checkmate::assertDataTable(ppi_network)
  checkmate::assertNames(
    names(ppi_network),
    must.include = c("from", "to", "weight")
  )
  checkmate::assertDataTable(grn_network)
  checkmate::assertNames(
    names(grn_network),
    must.include = c("from", "to", "weight")
  )
  assertLigandTarget(params)

  # gene universe = nodes appearing in either network
  all_genes <- unique(c(
    ppi_network$from,
    ppi_network$to,
    grn_network$from,
    grn_network$to
  ))

  missing_ligands <- setdiff(
    unlist(ligand_seeds, use.names = FALSE),
    all_genes
  )
  if (length(missing_ligands) > 0L) {
    stop(sprintf(
      "Ligand seeds not present in either network: %s",
      paste(missing_ligands, collapse = ", ")
    ))
  }

  # 0-indexed for Rust
  gene_to_idx <- setNames(seq_along(all_genes) - 1L, all_genes)

  ppi_list <- list(
    from = unname(gene_to_idx[ppi_network$from]),
    to = unname(gene_to_idx[ppi_network$to]),
    weight = ppi_network$weight
  )
  grn_list <- list(
    from = unname(gene_to_idx[grn_network$from]),
    to = unname(gene_to_idx[grn_network$to]),
    weight = grn_network$weight
  )

  seeds_int <- lapply(ligand_seeds, \(x) unname(gene_to_idx[x]))

  influence_mat <- rs_generate_ligand_target_influence(
    ligand_seeds = seeds_int,
    ppi_network = ppi_list,
    grn_network = grn_list,
    n_nodes = length(all_genes),
    params = params
  )

  ligand_names <- names(ligand_seeds)
  if (is.null(ligand_names) || any(ligand_names == "")) {
    ligand_names <- vapply(
      ligand_seeds,
      \(x) paste(x, collapse = "+"),
      character(1)
    )
  }

  new_ligand_target_influence(
    influence = influence_mat,
    ligand_seeds = ligand_seeds,
    ligand_names = ligand_names,
    gene_ids = all_genes,
    params = params
  )
}

## scoring ---------------------------------------------------------------------

#' Compute ligand activity scores against gene sets
#'
#' @description
#' For each gene set, ranks ligands by how well their influence vector aligns
#' with set membership across a background. Returns AUROC, AUPR, AUPR corrected
#' against the random baseline, Pearson, and Spearman per ligand and per gene
#' set.
#'
#' The influence matrix is restricted to `background` columns before
#' scoring, so AUROC / AUPR reflect ranking within the background (typical
#' NicheNet practice: background = genes expressed in the receiver cells,
#' gene set = the DEGs).
#'
#' @param ligand_influence A `LigandTargetInfluence` object.
#' @param gene_sets Either a character vector (one gene set) or a list of
#' character vectors. Names of the list are propagated to the output.
#' @param background Character vector or `NULL`. Genes against which gene
#' sets are scored. Defaults to all genes in `ligand_influence`. Background
#' members not present in the influence matrix are silently dropped.
#'
#' @returns A `data.table` with one row per (gene set, ligand) pair and
#' columns `gene_set`, `ligand`, `auroc`, `aupr`, `aupr_corrected`,
#' `pearson`, `spearman`.
#'
#' @export
ligand_activity_scores <- function(
  ligand_influence,
  gene_sets,
  background = NULL
) {
  checkmate::assertClass(ligand_influence, "LigandTargetInfluence")
  if (is.character(gene_sets)) {
    gene_sets <- list(gene_sets)
  }
  checkmate::assertList(gene_sets, types = "character", min.len = 1L)

  inf_mat <- get_influence(ligand_influence)
  if (is.null(background)) {
    background <- colnames(inf_mat)
  }
  checkmate::assertCharacter(background, min.len = 1L, any.missing = FALSE)

  bg_present <- intersect(background, colnames(inf_mat))
  if (length(bg_present) == 0L) {
    stop("No background genes overlap with the influence matrix columns.")
  }
  inf_sub <- inf_mat[, bg_present, drop = FALSE]

  # logical membership per gene set, in bg_present order
  in_sets <- lapply(gene_sets, \(gs) bg_present %in% gs)

  res <- rs_ligand_activity_scores(inf_sub, in_sets)

  set_names <- names(gene_sets)
  if (is.null(set_names) || any(set_names == "")) {
    set_names <- as.character(seq_along(gene_sets))
  }
  ligand_names <- rownames(inf_sub)

  out <- lapply(seq_along(res), \(i) {
    data.table::data.table(
      gene_set = set_names[i],
      ligand = ligand_names,
      auroc = res[[i]]$auroc,
      aupr = res[[i]]$aupr,
      aupr_corrected = res[[i]]$aupr_corrected,
      pearson = res[[i]]$pearson,
      spearman = res[[i]]$spearman
    )
  })
  data.table::rbindlist(out)
}
