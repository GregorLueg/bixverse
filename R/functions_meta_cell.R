# meta cell functions ----------------------------------------------------------

## merging ---------------------------------------------------------------------

#' Merge meta cell objects into one
#'
#' @description
#' Row-binds several [bixverse::MetaCells()] objects into a single one. The
#' typical use case is sample-pure meta cells: generate meta cells per patient
#' (see [meta_cells_per_group()]), then merge them so that methods like SCENIC,
#' AUCell or NMF can run across the full set.
#'
#' The normalised counts are carried over as generated. They are normalised
#' per meta cell, so row-binding leaves them valid. Caches (PCA, kNN, sNN,
#' embeddings) are per source and are dropped; recompute them on the merged
#' object.
#'
#' @param inputs List of `MetaCells` objects.
#' @param source_ids Optional character vector of the same length as `inputs`
#' with the source (e.g. patient) identifiers. Defaults to `names(inputs)` and
#' falls back to `source_01`, `source_02`, ... Needs to be unique.
#' @param feature_space String. One of `c("intersect", "union")`. Controls how
#' differing gene spaces are resolved. With `"union"` genes missing from an
#' input become structural zeros for its meta cells. Irrelevant when all inputs
#' came from the same source object, as their gene spaces are then identical.
#' @param prefix_ids Boolean. Prefix the meta cell identifiers with the source
#' identifier. If `FALSE`, duplicated identifiers across inputs are an error.
#' @param .verbose Boolean. Controls verbosity of the function.
#'
#' @returns A single `MetaCells` object with all meta cells of the inputs. The
#' observation table gains a `source_id` column; `original_cell_idx` stays in
#' the index space of its own source.
#'
#' @export
merge_meta_cells <- function(
  inputs,
  source_ids = NULL,
  feature_space = c("intersect", "union"),
  prefix_ids = TRUE,
  .verbose = TRUE
) {
  feature_space <- match.arg(feature_space)

  # checks
  checkmate::assertList(inputs, min.len = 1L)
  purrr::walk(inputs, \(x) checkmate::assertClass(x, "bixverse::MetaCells"))
  checkmate::assertChoice(feature_space, c("intersect", "union"))
  checkmate::qassert(prefix_ids, "B1")
  checkmate::qassert(.verbose, "B1")

  if (is.null(source_ids)) {
    source_ids <- if (!is.null(names(inputs))) {
      names(inputs)
    } else {
      sprintf("source_%02d", seq_along(inputs))
    }
  }
  checkmate::assertCharacter(source_ids, len = length(inputs), unique = TRUE)

  methods <- unique(purrr::map_chr(
    inputs,
    \(x) S7::prop(x, "meta_cell_method")
  ))
  if (length(methods) > 1L) {
    stop(sprintf(
      "All inputs need the same meta cell method. Found: %s.",
      paste(methods, collapse = ", ")
    ))
  }

  # gene space
  gene_lists <- purrr::map(inputs, \(x) S7::prop(x, "var_table")$gene_id)
  target_genes <- .resolve_mc_gene_space(gene_lists, feature_space)
  gene_maps <- purrr::map(gene_lists, \(g) match(g, target_genes))

  if (.verbose) {
    message(sprintf(
      "Merging %i meta cell objects over %i genes (%s).",
      length(inputs),
      length(target_genes),
      feature_space
    ))
  }

  dropped <- purrr::map_int(gene_maps, \(m) sum(is.na(m)))
  if (any(dropped > 0L)) {
    warning(sprintf(
      paste(
        "Dropping up to %i genes from individual inputs. The normalised",
        "counts of the affected meta cells no longer sum to their original",
        "target library size."
      ),
      max(dropped)
    ))
  }

  # observations
  obs_merged <- .merge_mc_obs(inputs, source_ids, prefix_ids)

  # counts
  counts <- .merge_mc_counts(
    inputs = inputs,
    gene_maps = gene_maps,
    meta_cell_ids = obs_merged$meta_cell_id,
    gene_ids = target_genes
  )

  # provenance
  per_source <- purrr::map(inputs, \(x) S7::prop(x, "original_assignment"))
  names(per_source) <- source_ids

  # sources built from the same parent share one obs space, so summing would
  # count it several times over
  source_cells <- purrr::map_dbl(per_source, "n_cells")
  shared_parent <- length(unique(source_cells)) == 1L
  n_cells <- if (shared_parent) source_cells[[1]] else sum(source_cells)
  # meta cells can share originating cells, so count the union rather than the
  # per meta cell totals
  n_unassigned <- if (shared_parent) {
    n_cells - length(unique(unlist(obs_merged$original_cell_idx)))
  } else {
    sum(purrr::map_dbl(per_source, "n_unassigned"))
  }

  other_data <- purrr::map(inputs, \(x) S7::prop(x, "other_data"))
  dropped_other <- unique(unlist(purrr::map(other_data, names)))
  if (length(dropped_other) > 0L) {
    warning(sprintf(
      "Dropping per-source entries from `other_data`: %s.",
      paste(dropped_other, collapse = ", ")
    ))
  }

  # `S7::new_object()` can only be called from a constructor, so the merged
  # object goes through `MetaCells()` and the merge-specific props are set after
  merged <- MetaCells(
    meta_cell_data = list(
      aggregated = list(
        indptr = counts$raw@p,
        indices = counts$raw@j,
        raw_counts = counts$raw@x,
        norm_counts = counts$norm@x,
        nrow = nrow(counts$raw),
        ncol = ncol(counts$raw)
      ),
      assignments = list(
        # source-specific; there is no single index space to map back into
        assignments = NULL,
        metacells = obs_merged$original_cell_idx,
        n_cells = as.integer(n_cells),
        n_metacells = nrow(obs_merged),
        n_unassigned = as.integer(n_unassigned)
      )
    ),
    var_data = data.table::data.table(
      gene_idx = seq_along(target_genes),
      gene_id = target_genes
    ),
    meta_cell_method = methods
  )

  S7::prop(merged, "obs_table") <- obs_merged
  S7::prop(merged, "data") <- counts
  S7::prop(merged, "original_assignment") <- c(
    S7::prop(merged, "original_assignment"),
    list(per_source = per_source)
  )
  S7::prop(merged, "other_data") <- list(merged = TRUE, sources = source_ids)

  merged
}

## internal helpers ------------------------------------------------------------

#' Resolve the target gene space for a meta cell merge
#'
#' @param gene_lists List of character vectors with the gene identifiers.
#' @param feature_space String. One of `c("intersect", "union")`.
#'
#' @returns Character vector with the target gene identifiers, ordered by the
#' first input.
#'
#' @keywords internal
.resolve_mc_gene_space <- function(gene_lists, feature_space) {
  checkmate::assertList(gene_lists, min.len = 1L, types = "character")
  checkmate::assertChoice(feature_space, c("intersect", "union"))

  target <- if (feature_space == "intersect") {
    Reduce(intersect, gene_lists)
  } else {
    Reduce(\(acc, g) c(acc, setdiff(g, acc)), gene_lists)
  }

  if (length(target) == 0L) {
    stop("The gene intersection of the inputs is empty.")
  }

  target
}

#' Row-bind the observation tables of meta cell objects
#'
#' @param inputs List of `MetaCells` objects.
#' @param source_ids Character vector with the source identifiers.
#' @param prefix_ids Boolean. Prefix the meta cell identifiers.
#'
#' @returns A data.table with the merged observations, re-indexed.
#'
#' @keywords internal
.merge_mc_obs <- function(inputs, source_ids, prefix_ids) {
  obs_list <- purrr::map(
    inputs,
    \(x) data.table::copy(S7::prop(x, "obs_table"))
  )

  shared <- Reduce(intersect, purrr::map(obs_list, names))
  cols <- names(obs_list[[1]])[names(obs_list[[1]]) %in% shared]
  dropped <- setdiff(unique(unlist(purrr::map(obs_list, names))), cols)
  if (length(dropped) > 0L) {
    warning(sprintf(
      "Dropping observation columns not shared by all inputs: %s.",
      paste(dropped, collapse = ", ")
    ))
  }

  merged <- data.table::rbindlist(purrr::map2(
    obs_list,
    source_ids,
    \(dt, sid) {
      dt <- dt[, cols, with = FALSE]
      dt[, source_id := sid]
      dt
    }
  ))

  if (prefix_ids) {
    merged[, meta_cell_id := sprintf("%s__%s", source_id, meta_cell_id)]
  } else if (anyDuplicated(merged$meta_cell_id) > 0L) {
    stop(paste(
      "Duplicated meta cell identifiers across the inputs.",
      "Use `prefix_ids = TRUE` or rename them."
    ))
  }

  merged[, meta_cell_idx := seq_len(.N)]
  data.table::setcolorder(
    merged,
    c("meta_cell_idx", "meta_cell_id", "source_id")
  )

  merged
}

#' Row-bind the count matrices of meta cell objects
#'
#' @param inputs List of `MetaCells` objects.
#' @param gene_maps List of integer vectors mapping each input's gene positions
#' onto the target gene space. `NA` marks genes to drop.
#' @param meta_cell_ids Character vector with the merged meta cell identifiers.
#' @param gene_ids Character vector with the target gene identifiers.
#'
#' @returns A list with the merged `raw` and `norm` CSR matrices.
#'
#' @keywords internal
.merge_mc_counts <- function(inputs, gene_maps, meta_cell_ids, gene_ids) {
  n_rows <- purrr::map_int(inputs, \(x) nrow(S7::prop(x, "data")$raw))
  offsets <- cumsum(c(0L, n_rows))

  parts <- purrr::map(seq_along(inputs), \(k) {
    raw <- S7::prop(inputs[[k]], "data")$raw
    norm <- S7::prop(inputs[[k]], "data")$norm
    # raw and norm are built off one index structure; reuse it for both
    checkmate::assertTRUE(identical(raw@j, norm@j))

    j <- gene_maps[[k]][raw@j + 1L]
    keep <- !is.na(j)

    list(
      i = rep.int(seq_len(nrow(raw)), diff(raw@p))[keep] + offsets[k],
      j = j[keep],
      raw = raw@x[keep],
      norm = norm@x[keep]
    )
  })

  dims <- c(sum(n_rows), length(gene_ids))
  dim_names <- list(meta_cell_ids, gene_ids)
  i <- unlist(purrr::map(parts, "i"))
  j <- unlist(purrr::map(parts, "j"))

  purrr::map(
    c(raw = "raw", norm = "norm"),
    \(assay) {
      Matrix::sparseMatrix(
        i = i,
        j = j,
        x = unlist(purrr::map(parts, assay)),
        dims = dims,
        dimnames = dim_names,
        repr = "R"
      )
    }
  )
}
