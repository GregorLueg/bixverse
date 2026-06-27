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

  res <- structure(
    list(
      influence = influence,
      ligand_seeds = ligand_seeds,
      ligand_names = ligand_names,
      gene_ids = gene_ids,
      params = params
    ),
    class = "LigandTargetInfluence"
  )

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

## prioritisation --------------------------------------------------------------

### helpers --------------------------------------------------------------------

#### scoring -------------------------------------------------------------------

#' Rank-based min-max scaling
#'
#' @description
#' Returns `rank(x) / max(rank(x))` with average ties and NAs sent to the
#' bottom. Used to map differential-expression statistics into `[0, 1]` before
#' aggregation.
#'
#' @param x Numeric vector.
#'
#' @returns Numeric vector in `[0, 1]`.
#'
#' @keywords internal
.rank_scale <- function(x) {
  r <- rank(x, ties.method = "average", na.last = FALSE)
  r / max(r)
}

#' Z-score scaling
#'
#' @description
#' Returns `(x - mean(x)) / sd(x)`. Returns zeros if `sd(x)` is zero or `NA`.
#'
#' @param x Numeric vector.
#'
#' @returns Numeric vector.
#'
#' @keywords internal
.scaling_zscore <- function(x) {
  s <- stats::sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) {
    return(rep(0, length(x)))
  }
  (x - mean(x, na.rm = TRUE)) / s
}

#' Quantile-clipped min-max scaling
#'
#' @description
#' Clips `x` to its `[outlier_cutoff, 1 - outlier_cutoff]` quantiles, then maps
#' linearly to `[0, 1]`. Returns 0.5 everywhere if the two quantiles coincide.
#'
#' @param x Numeric vector.
#' @param outlier_cutoff Quantile cutoff. Defaults to 0.05.
#'
#' @returns Numeric vector in `[0, 1]`.
#'
#' @keywords internal
.scale_quantile_adapted <- function(x, outlier_cutoff = 0.05) {
  q_low <- quantile(x, outlier_cutoff, na.rm = TRUE)
  q_high <- quantile(x, 1 - outlier_cutoff, na.rm = TRUE)
  if (q_high == q_low) {
    return(rep(0.5, length(x)))
  }
  clipped <- pmin(pmax(x, q_low), q_high)
  as.numeric((clipped - q_low) / (q_high - q_low))
}

#### assertions ----------------------------------------------------------------

#' Assert a DE table has the expected columns and types
#'
#' @param x The candidate `data.table`.
#' @param name Variable name for error messages.
#' @param needs_cluster If TRUE, also require a `cluster_id` column.
#'
#' @keywords internal
assertDegTable <- function(x, name, needs_cluster = TRUE) {
  checkmate::assertDataTable(x, .var.name = name)
  required <- c("gene", "lfc", "pval")
  if (needs_cluster) {
    required <- c("cluster_id", required)
  }
  checkmate::assertNames(names(x), must.include = required, .var.name = name)
  checkmate::assertNumeric(
    x$lfc,
    any.missing = FALSE,
    .var.name = paste0(name, "$lfc")
  )
  checkmate::assertNumeric(
    x$pval,
    lower = 0,
    upper = 1,
    any.missing = FALSE,
    .var.name = paste0(name, "$pval")
  )
  checkmate::assertCharacter(
    x$gene,
    any.missing = FALSE,
    .var.name = paste0(name, "$gene")
  )
  if (needs_cluster) {
    checkmate::assertCharacter(
      x$cluster_id,
      any.missing = FALSE,
      .var.name = paste0(name, "$cluster_id")
    )
  }
  invisible(x)
}

#' Assert an expression info table has the expected columns and types
#'
#' @param x The candidate `data.table`.
#' @param name Variable name for error messages.
#'
#' @keywords internal
assertExpTable <- function(x, name) {
  checkmate::assertDataTable(x, .var.name = name)
  checkmate::assertNames(
    names(x),
    must.include = c("cluster_id", "gene", "avg_expr"),
    .var.name = name
  )
  checkmate::assertNumeric(
    x$avg_expr,
    any.missing = FALSE,
    .var.name = paste0(name, "$avg_expr")
  )
  checkmate::assertCharacter(
    x$gene,
    any.missing = FALSE,
    .var.name = paste0(name, "$gene")
  )
  checkmate::assertCharacter(
    x$cluster_id,
    any.missing = FALSE,
    .var.name = paste0(name, "$cluster_id")
  )
  invisible(x)
}

#### others --------------------------------------------------------------------

#' Resolve prioritisation weights from the scenario or validate a user vector
#'
#' @param weights ...
#' @param scenario ...
#' @param has_condition_de ...
#'
#' @returns ...
#'
#' @keywords internal
resolve_weights <- function(weights, scenario, has_condition_de) {
  weight_names <- c(
    "de_ligand",
    "de_receptor",
    "activity_scaled",
    "exprs_ligand",
    "exprs_receptor",
    "ligand_condition_specificity",
    "receptor_condition_specificity"
  )
  if (is.null(weights)) {
    weights <- switch(
      scenario,
      "case_control" = setNames(rep(1, 7), weight_names),
      "one_condition" = setNames(c(1, 1, 1, 1, 1, 0, 0), weight_names)
    )
  } else {
    checkmate::assertNumeric(
      weights,
      lower = 0,
      finite = TRUE,
      names = "named"
    )
    missing_names <- setdiff(weight_names, names(weights))
    if (length(missing_names) > 0L) {
      stop(sprintf(
        "weights is missing: %s",
        paste(missing_names, collapse = ", ")
      ))
    }
    weights <- weights[weight_names]
  }
  cond_keys <- c(
    "ligand_condition_specificity",
    "receptor_condition_specificity"
  )
  if (!has_condition_de && any(weights[cond_keys] > 0)) {
    stop("Non-zero condition-specificity weights require `condition_de`.")
  }
  weights
}

### cluster l/r expression -----------------------------------------------------

#' Compute per-cluster mean expression and expressing fraction for a gene set
#'
#' @description
#' Thin R wrapper around the Rust `compute_cluster_expression_stats` routine.
#' Streams gene chunks from the on-disk store and aggregates expression across
#' user-supplied cell clusters. Cells outside any cluster are ignored.
#'
#' If `condition_colname` and `condition_oi` are supplied, only cells from
#' that condition contribute to the aggregation.
#'
#' @param object A `SingleCells` object.
#' @param celltype_colname Name of the cluster column in `obs`.
#' @param genes Character vector of gene IDs to aggregate over.
#' @param condition_colname Optional. Name of a condition column in `obs`.
#' @param condition_oi Optional. Value of `condition_colname` to subset to.
#'
#' @returns A long `data.table` with columns `cluster_id`, `gene`,
#' `avg_expr`, `frac_expr`.
#'
#' @export
compute_expression_info_sc <- function(
  object,
  celltype_colname,
  genes,
  condition_colname = NULL,
  condition_oi = NULL
) {
  checkmate::assertTRUE(S7::S7_inherits(object, SingleCells))
  checkmate::qassert(celltype_colname, "S1")
  checkmate::assertCharacter(genes, min.len = 1L, any.missing = FALSE)
  checkmate::qassert(condition_colname, c("0", "S1"))
  checkmate::qassert(condition_oi, c("0", "S1"))
  if (xor(is.null(condition_colname), is.null(condition_oi))) {
    stop("Provide both `condition_colname` and `condition_oi`, or neither.")
  }

  # gene id -> 0-indexed lookup
  var <- get_sc_var(object)
  checkmate::assertNames(names(var), must.include = "gene_id")
  gene_to_idx <- setNames(seq_along(var$gene_id) - 1L, var$gene_id)
  missing_genes <- setdiff(genes, var$gene_id)
  if (length(missing_genes) > 0L) {
    warning(sprintf(
      "Dropping %d gene(s) not in the var table.",
      length(missing_genes)
    ))
    genes <- intersect(genes, var$gene_id)
    if (length(genes) == 0L) stop("No requested genes are present in `object`.")
  }
  gene_indices <- unname(gene_to_idx[genes])

  obs <- get_sc_obs(object, filtered = TRUE)
  checkmate::assertNames(names(obs), must.include = celltype_colname)
  if (!is.null(condition_colname)) {
    checkmate::assertNames(names(obs_keep), must.include = condition_colname)
    obs_keep <- obs_keep[get(condition_colname) == condition_oi]
  }
  if (nrow(obs_keep) == 0L) {
    stop("No cells remain after filtering by `to_keep` and condition.")
  }

  cluster_labels <- obs_keep[[celltype_colname]]
  cluster_levels <- sort(unique(cluster_labels))
  clusters <- lapply(cluster_levels, function(cl) {
    obs_keep[cluster_labels == cl, .cell_idx]
  })

  res <- rs_compute_cluster_expr_stats(
    f_path_gene = bixverse:::get_rust_count_gene_f_path(object),
    gene_indices = gene_indices,
    clusters = clusters
  )

  rownames(res$mean) <- genes
  colnames(res$mean) <- cluster_levels
  rownames(res$frac) <- genes
  colnames(res$frac) <- cluster_levels

  data.table::data.table(
    cluster_id = rep(cluster_levels, each = length(genes)),
    gene = rep(genes, length(cluster_levels)),
    avg_expr = as.numeric(res$mean),
    frac_expr = as.numeric(res$frac)
  )
}

### main -----------------------------------------------------------------------

#' Prioritise sender-ligand-receiver-receptor interactions
#'
#' @description
#' For every (sender, ligand, receiver, receptor) tuple drawn from the LR
#' network, compute a weighted prioritisation score combining:
#'
#' \itemize{
#'   \item ligand DE in the sender (`de_ligand`)
#'   \item receptor DE in the receiver (`de_receptor`)
#'   \item ligand activity in the receiver (`activity_scaled`)
#'   \item ligand expression specificity across senders (`exprs_ligand`)
#'   \item receptor expression specificity across receivers (`exprs_receptor`)
#'   \item ligand condition specificity (`ligand_condition_specificity`)
#'   \item receptor condition specificity (`receptor_condition_specificity`)
#' }
#'
#' Components with weight zero are skipped (their joins are not performed).
#'
#' DE tables are supplied by the user. `find_markers_sc` and
#' `find_all_markers_sc` give a quick Wilcox; pseudo-bulk + DESeq2 / edgeR is
#' generally preferable for statistical inference.
#'
#' @param celltype_de A `data.table` with columns `cluster_id`, `gene`,
#' `lfc`, `pval`. One row per (cluster, gene). Must cover both senders and
#' receivers.
#' @param expression_info A `data.table` with columns `cluster_id`,
#' `gene`, `avg_expr`. As returned by `compute_expression_info_sc`.
#' @param ligand_activities A `data.table` with at least the columns
#' `ligand` and `aupr_corrected`. One row per ligand. Subset to a single
#' gene set before calling.
#' @param lr_network A `data.table` with columns `ligand`, `receptor`.
#' @param senders_oi Character vector of sender cluster IDs.
#' @param receivers_oi Character vector of receiver cluster IDs.
#' @param condition_de Optional `data.table` with columns `gene`, `lfc`,
#' `pval` from a condition contrast (e.g. case vs control). Required when
#' `scenario == "case_control"` and the corresponding weights are non-zero.
#' @param weights Optional named numeric vector. If `NULL`, defaults are
#' chosen by `scenario`. Names: `de_ligand`, `de_receptor`,
#' `activity_scaled`, `exprs_ligand`, `exprs_receptor`,
#' `ligand_condition_specificity`, `receptor_condition_specificity`.
#' @param scenario `"case_control"` (all weights 1) or `"one_condition"`
#' (condition-specificity weights 0). Ignored if `weights` is supplied.
#'
#' @returns A `data.table` with one row per surviving
#' (sender, ligand, receiver, receptor) tuple, sorted by descending
#' `prioritisation_score`, with a `prioritisation_rank` column.
#'
#' @export
prioritise_interactions <- function(
  celltype_de,
  expression_info,
  ligand_activities,
  lr_network,
  senders_oi,
  receivers_oi,
  condition_de = NULL,
  weights = NULL,
  scenario = c("case_control", "one_condition")
) {
  scenario <- match.arg(scenario)

  # checks
  assertDegTable(celltype_de, "celltype_de", needs_cluster = TRUE)
  assertExpTable(expression_info, "expression_info")
  checkmate::assertDataTable(ligand_activities)
  checkmate::assertNames(
    names(ligand_activities),
    must.include = c("ligand", "aupr_corrected")
  )
  checkmate::assertDataTable(lr_network)
  checkmate::assertNames(
    names(lr_network),
    must.include = c("ligand", "receptor")
  )
  checkmate::assertCharacter(senders_oi, min.len = 1L, any.missing = FALSE)
  checkmate::assertCharacter(receivers_oi, min.len = 1L, any.missing = FALSE)
  if (!is.null(condition_de)) {
    assertDegTable(condition_de, "condition_de", needs_cluster = FALSE)
  }
  if (anyDuplicated(ligand_activities$ligand)) {
    stop(
      paste(
        "`ligand_activities` must have one row per ligand;",
        "filter to a single gene set first."
      )
    )
  }

  w <- resolve_weights(weights, scenario, !is.null(condition_de))

  ligands <- unique(lr_network$ligand)
  receptors <- unique(lr_network$receptor)

  # 1. sender ligand DE
  sender_de <- celltype_de[
    cluster_id %in% senders_oi & gene %in% ligands,
    .(sender = cluster_id, ligand = gene, lfc_ligand = lfc, pval_ligand = pval)
  ]
  if (w["de_ligand"] > 0 && nrow(sender_de) > 0L) {
    sender_de[,
      p_val_adapted_ligand := -log10(pmax(pval_ligand, .Machine$double.xmin)) *
        sign(lfc_ligand)
    ]
    sender_de[, scaled_p_val_adapted_ligand := rank_scale(p_val_adapted_ligand)]
  }

  # 2. receiver receptor DE
  receiver_de <- celltype_de[
    cluster_id %in% receivers_oi & gene %in% receptors,
    .(
      receiver = cluster_id,
      receptor = gene,
      lfc_receptor = lfc,
      pval_receptor = pval
    )
  ]
  if (w["de_receptor"] > 0 && nrow(receiver_de) > 0L) {
    receiver_de[,
      p_val_adapted_receptor := -log10(pmax(
        pval_receptor,
        .Machine$double.xmin
      )) *
        sign(lfc_receptor)
    ]
    receiver_de[,
      scaled_p_val_adapted_receptor := rank_scale(p_val_adapted_receptor)
    ]
  }

  # 3. ligand expression per sender, scaled within each ligand across senders
  ligand_expr <- expression_info[
    cluster_id %in% senders_oi & gene %in% ligands,
    .(sender = cluster_id, ligand = gene, avg_ligand = avg_expr)
  ]
  if (w["exprs_ligand"] > 0 && nrow(ligand_expr) > 0L) {
    ligand_expr[,
      scaled_avg_exprs_ligand := scale_quantile_adapted(avg_ligand),
      by = ligand
    ]
  }

  # 4. receptor expression per receiver, scaled within each receptor
  receptor_expr <- expression_info[
    cluster_id %in% receivers_oi & gene %in% receptors,
    .(receiver = cluster_id, receptor = gene, avg_receptor = avg_expr)
  ]
  if (w["exprs_receptor"] > 0 && nrow(receptor_expr) > 0L) {
    receptor_expr[,
      scaled_avg_exprs_receptor := scale_quantile_adapted(avg_receptor),
      by = receptor
    ]
  }

  # 5. ligand activity
  la <- data.table::copy(ligand_activities)
  la_keep <- c("ligand", "aupr_corrected")
  if (w["activity_scaled"] > 0) {
    la[, activity_zscore := scaling_zscore(aupr_corrected)]
    la[,
      scaled_activity := scale_quantile_adapted(
        aupr_corrected,
        outlier_cutoff = 0.01
      )
    ]
    la_keep <- c(la_keep, "activity_zscore", "scaled_activity")
  }
  la <- la[, la_keep, with = FALSE]
  data.table::setnames(la, "aupr_corrected", "activity")

  # 6. condition DE
  ligand_cond <- NULL
  receptor_cond <- NULL
  if (!is.null(condition_de)) {
    if (w["ligand_condition_specificity"] > 0) {
      ligand_cond <- condition_de[
        gene %in% ligands,
        .(ligand = gene, lfc_ligand_group = lfc, pval_ligand_group = pval)
      ]
      ligand_cond[,
        p_val_adapted_ligand_group := -log10(pmax(
          pval_ligand_group,
          .Machine$double.xmin
        )) *
          sign(lfc_ligand_group)
      ]
      ligand_cond[,
        scaled_p_val_adapted_ligand_group := rank_scale(
          p_val_adapted_ligand_group
        )
      ]
    }
    if (w["receptor_condition_specificity"] > 0) {
      receptor_cond <- condition_de[
        gene %in% receptors,
        .(receptor = gene, lfc_receptor_group = lfc, pval_receptor_group = pval)
      ]
      receptor_cond[,
        p_val_adapted_receptor_group := -log10(pmax(
          pval_receptor_group,
          .Machine$double.xmin
        )) *
          sign(lfc_receptor_group)
      ]
      receptor_cond[,
        scaled_p_val_adapted_receptor_group := rank_scale(
          p_val_adapted_receptor_group
        )
      ]
    }
  }

  # 7. base table: sender x receiver x LR pair
  sr <- data.table::CJ(sender = senders_oi, receiver = receivers_oi)
  sr[, .merge_k := 1L]
  lr_dt <- data.table::copy(lr_network[, .(ligand, receptor)])
  lr_dt[, .merge_k := 1L]
  out <- merge(sr, lr_dt, by = ".merge_k", allow.cartesian = TRUE)
  out[, .merge_k := NULL]

  # 8. inner-join active components
  if (w["de_ligand"] > 0) {
    out <- sender_de[out, on = c("sender", "ligand"), nomatch = NULL]
  }
  if (w["de_receptor"] > 0) {
    out <- receiver_de[out, on = c("receiver", "receptor"), nomatch = NULL]
  }
  if (w["exprs_ligand"] > 0) {
    out <- ligand_expr[out, on = c("sender", "ligand"), nomatch = NULL]
  }
  if (w["exprs_receptor"] > 0) {
    out <- receptor_expr[out, on = c("receiver", "receptor"), nomatch = NULL]
  }
  if (w["activity_scaled"] > 0) {
    out <- la[out, on = "ligand", nomatch = NULL]
  }
  if (!is.null(ligand_cond)) {
    out <- ligand_cond[out, on = "ligand", nomatch = NULL]
  }
  if (!is.null(receptor_cond)) {
    out <- receptor_cond[out, on = "receptor", nomatch = NULL]
  }

  if (nrow(out) == 0L) {
    warning("No interactions remain after all joins. Check inputs.")
    return(out)
  }

  # 9. weighted score
  # Mirrors the nichenetr convention: DE and expression terms are halved
  # because each contributes both ligand and receptor sides.
  sum_w <- 0.5 *
    w["de_ligand"] +
    0.5 * w["de_receptor"] +
    w["activity_scaled"] +
    0.5 * w["exprs_ligand"] +
    0.5 * w["exprs_receptor"] +
    w["ligand_condition_specificity"] +
    w["receptor_condition_specificity"]

  score <- numeric(nrow(out))
  add_term <- function(col, wt) {
    if (wt > 0 && col %in% names(out)) {
      score <<- score + wt * out[[col]]
    }
  }
  add_term("scaled_p_val_adapted_ligand", 0.5 * w["de_ligand"])
  add_term("scaled_p_val_adapted_receptor", 0.5 * w["de_receptor"])
  add_term("scaled_activity", w["activity_scaled"])
  add_term("scaled_avg_exprs_ligand", 0.5 * w["exprs_ligand"])
  add_term("scaled_avg_exprs_receptor", 0.5 * w["exprs_receptor"])
  add_term(
    "scaled_p_val_adapted_ligand_group",
    w["ligand_condition_specificity"]
  )
  add_term(
    "scaled_p_val_adapted_receptor_group",
    w["receptor_condition_specificity"]
  )

  out[, prioritisation_score := score / as.numeric(sum_w)]
  out[,
    prioritisation_rank := data.table::frank(
      -prioritisation_score,
      ties.method = "min"
    )
  ]
  data.table::setorder(out, -prioritisation_score)

  id_cols <- c("sender", "receiver", "ligand", "receptor")
  data.table::setcolorder(out, c(id_cols, setdiff(names(out), id_cols)))

  out
}
