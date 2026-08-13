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
#' @returns A single `MetaCells` object with all meta cells of the inputs and
#' `is_merged` set to `TRUE`. The observation table gains a `source_id` column;
#' `original_cell_idx` stays in the index space of its own source, which is why
#' methods that resolve it against the source single cell data
#' ([calc_diffusion_coordinates()], [calc_manifold_metrics()]) refuse to run on
#' the result. `other_data` holds the source identifiers under `sources`.
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
    n_genes = length(target_genes)
  )

  # provenance. Only the scalars: the full `assignments` element holds one
  # integer vector per cell of the *parent* count file, so keeping it per source
  # costs O(sources * cells) vectors that nothing ever reads.
  per_source <- purrr::map(
    inputs,
    \(x) {
      S7::prop(x, "original_assignment")[c(
        "n_cells",
        "n_metacells",
        "n_unassigned"
      )]
    }
  )
  names(per_source) <- source_ids

  # sources built from the same parent share one obs space, so summing would
  # count it several times over
  source_cells <- purrr::map_dbl(per_source, "n_cells")
  shared_parent <- length(unique(source_cells)) == 1L
  n_cells <- if (shared_parent) source_cells[[1]] else sum(source_cells)
  # meta cells can share originating cells, so count the union rather than the
  # per meta cell totals. A bitmap over the obs space beats unlist() + unique(),
  # whose intermediate is a multiple of `n_cells` when meta cells overlap.
  n_unassigned <- if (shared_parent) {
    seen <- logical(n_cells)
    for (idx in obs_merged$original_cell_idx) {
      seen[idx] <- TRUE
    }
    n_cells - sum(seen)
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
  # object goes through `MetaCells()` and the merge-specific props are set after.
  # `obs_ids` makes the constructor build the count matrices with their final row
  # names, so the merged pair is materialised exactly once.
  merged <- MetaCells(
    meta_cell_data = list(
      aggregated = counts,
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
    meta_cell_method = methods,
    obs_ids = obs_merged$meta_cell_id
  )

  S7::prop(merged, "obs_table") <- obs_merged
  S7::prop(merged, "original_assignment") <- c(
    S7::prop(merged, "original_assignment"),
    list(per_source = per_source)
  )
  S7::prop(merged, "other_data") <- list(sources = source_ids)
  S7::prop(merged, "is_merged") <- TRUE

  merged
}

## index space bridge ----------------------------------------------------------

#' Resolve meta cell memberships onto the row space of a source artefact
#'
#' @description
#' `original_cell_idx` holds 1-indexed positions in the *full, unfiltered* obs
#' table of the source object. Embeddings, kNN graphs and diffusion maps have
#' one row per QC-passing cell, in `cells_to_keep` order (the invariant the
#' cache stamps rely on). Indexing such an artefact with the memberships
#' directly is silently off by the filtering the moment QC dropped a cell, which
#' is why the source's `cells_to_keep` is recorded on the object at generation
#' time. This resolves the one space into the other.
#'
#' @param object `MetaCells` class.
#' @param n_rows Integer. Number of rows the artefact has, checked against the
#' recorded cell mapping.
#'
#' @returns A list of integer vectors, one per meta cell, holding 1-indexed row
#' positions in the artefact.
#'
#' @keywords internal
.mc_artefact_rows <- function(object, n_rows) {
  checkmate::assertTRUE(S7::S7_inherits(object, MetaCells))
  checkmate::qassert(n_rows, "X1[0,)")

  assignment <- S7::prop(object, "original_assignment")
  keep <- assignment$cells_to_keep
  n_cells <- assignment$n_cells
  memberships <- S7::prop(object, "obs_table")$original_cell_idx

  # meta cells generated before the mapping was recorded. Passing them through
  # is only safe when nothing was filtered, as the two spaces then coincide.
  if (is.null(keep)) {
    if (n_rows != n_cells) {
      stop(sprintf(
        paste(
          "This object does not record the source cell mapping, and the",
          "artefact covers %i of the %i original cells, so the meta cell",
          "memberships cannot be resolved against it. Regenerate the meta",
          "cells to record the mapping."
        ),
        n_rows,
        n_cells
      ))
    }

    return(memberships)
  }

  if (n_rows != length(keep)) {
    stop(sprintf(
      paste(
        "The artefact covers %i cells but the meta cells were built over the",
        "%i QC-passing cells of the source. Regenerate it so the two agree."
      ),
      n_rows,
      length(keep)
    ))
  }

  # reverse lookup, so this stays linear instead of a `match()` per meta cell
  lookup <- integer(n_cells)
  lookup[keep + 1L] <- seq_along(keep)

  rows <- purrr::map(memberships, \(idx) lookup[idx])

  # 0 marks a member that is not part of the artefact at all, which cannot
  # happen for memberships drawn from the kept cells
  if (any(purrr::map_lgl(rows, \(r) any(r == 0L)))) {
    stop(paste(
      "Some meta cell members are absent from the recorded cell mapping.",
      "The object is inconsistent; regenerate the meta cells."
    ))
  }

  rows
}

## purity ----------------------------------------------------------------------

#' Build the meta cell purity columns
#'
#' @description
#' Shared body of [bixverse::calc_meta_cell_purity()] and
#' [bixverse::get_meta_cell_purity()]. Resolves the label vector against the
#' right index space and tabulates the label distribution per meta cell.
#'
#' @param object `MetaCells` class.
#' @param original_cell_type Character vector. The original cell type
#' annotations, either in full obs order or in QC-passing order. Which one was
#' passed is inferred from the length.
#' @param add_additional_info String. One of
#' `c("none", "top_label", "top_two_labels")`. Already resolved by `match.arg()`
#' in the caller.
#' @param add_entropy Boolean. Add the normalised Shannon entropy.
#'
#' @returns A named list of atomic vectors, one element per meta cell in
#' `obs_table` row order:
#' \itemize{
#'   \item mc_purity - Fraction of the meta cell's cells that carry the most
#'   abundant label. Always present.
#'   \item mc_top_label - Name of the most abundant label. Present for
#'   `add_additional_info %in% c("top_label", "top_two_labels")`.
#'   \item mc_second_label - Name of the second most abundant label, `NA` for a
#'   pure meta cell. Present for `add_additional_info = "top_two_labels"`.
#'   \item mc_second_frac - Fraction of the meta cell's cells carrying that
#'   second label, `0` for a pure meta cell. Present for
#'   `add_additional_info = "top_two_labels"`.
#'   \item mc_entropy - Normalised Shannon entropy of the label distribution.
#'   Present for `add_entropy = TRUE`.
#' }
#'
#' @keywords internal
.mc_purity_cols <- function(
  object,
  original_cell_type,
  add_additional_info,
  add_entropy
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(object, MetaCells))
  checkmate::qassert(original_cell_type, "S+")
  checkmate::assertChoice(
    add_additional_info,
    c("none", "top_label", "top_two_labels")
  )
  checkmate::qassert(add_entropy, "B1")

  # memberships index into the full obs space; a short vector would silently
  # recycle or produce NA rather than error
  assignment <- S7::prop(object, "original_assignment")
  n_cells <- assignment$n_cells

  # a merged object only has one obs space if all its sources shared a parent
  if (isTRUE(S7::prop(object, "is_merged"))) {
    source_cells <- unique(purrr::map_dbl(assignment$per_source, "n_cells"))
    if (length(source_cells) > 1L) {
      stop(paste(
        "The meta cells were merged from sources with different obs spaces,",
        "so `original_cell_idx` cannot be resolved against a single column.",
        "Calculate the purity per source before merging."
      ))
    }
  }

  # a full-obs vector indexes straight, a QC-passing one needs the memberships
  # translated into the filtered row space
  n_kept <- length(assignment$cells_to_keep)

  mc_rows <- if (length(original_cell_type) == n_cells) {
    S7::prop(object, "obs_table")$original_cell_idx
  } else if (n_kept > 0L && length(original_cell_type) == n_kept) {
    .mc_artefact_rows(object, n_kept)
  } else {
    stop(sprintf(
      paste(
        "`original_cell_type` has %i entries but the meta cells were built over",
        "%i original cells%s. Pass the unfiltered obs column, or the",
        "QC-passing one."
      ),
      length(original_cell_type),
      n_cells,
      if (n_kept > 0L) sprintf(" of which %i pass QC", n_kept) else ""
    ))
  }

  # one factor pass up front, so each meta cell becomes a tabulate() over
  # integer codes instead of a table() rebuilding the level set every time
  labels_fct <- factor(original_cell_type)
  codes <- as.integer(labels_fct)
  lvls <- levels(labels_fct)
  n_labels <- length(lvls)

  counts <- purrr::map(
    mc_rows,
    \(idx) tabulate(codes[idx], nbins = n_labels)
  )
  sizes <- purrr::map_int(mc_rows, length)

  new_cols <- list(
    mc_purity = purrr::map2_dbl(counts, sizes, \(cnt, n) max(cnt) / n)
  )

  if (add_additional_info != "none") {
    new_cols$mc_top_label <- purrr::map_chr(counts, \(cnt) lvls[which.max(cnt)])
  }

  if (add_additional_info == "top_two_labels") {
    # zero the winner and take the argmax again. An all-zero remainder means the
    # meta cell only ever saw a single label
    second <- purrr::map(counts, \(cnt) {
      cnt[which.max(cnt)] <- 0L
      cnt
    })

    new_cols$mc_second_label <- purrr::map_chr(
      second,
      \(cnt) if (max(cnt) == 0L) NA_character_ else lvls[which.max(cnt)]
    )
    new_cols$mc_second_frac <- purrr::map2_dbl(
      second,
      sizes,
      \(cnt, n) max(cnt) / n
    )
  }

  if (add_entropy) {
    new_cols$mc_entropy <- if (n_labels < 2L) {
      # log(1) would turn the normalisation into a NaN
      rep(0, length(counts))
    } else {
      purrr::map2_dbl(counts, sizes, \(cnt, n) {
        p <- cnt[cnt > 0L] / n
        -sum(p * log(p)) / log(n_labels)
      })
    }
  }

  return(new_cols)
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
  # no `copy()` needed: `dt[, cols, with = FALSE]` below already copies every
  # column, so the `:=` never reaches the caller's table
  obs_list <- purrr::map(inputs, \(x) S7::prop(x, "obs_table"))

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
#' @description
#' The inputs are already CSR, so a row-bind is concatenation of their `@p`,
#' `@j` and `@x` slots. No COO round-trip and no `Matrix::sparseMatrix(i, j, x)`
#' call, which would sort and copy the whole index structure once per assay. The
#' output vectors are allocated once off a counting pass, so peak memory is the
#' result plus one input's scratch.
#'
#' @param inputs List of `MetaCells` objects.
#' @param gene_maps List of integer vectors mapping each input's gene positions
#' onto the target gene space, 1-indexed. `NA` marks genes to drop.
#' @param n_genes Integer. Number of genes in the target gene space.
#'
#' @returns A list in the shape `get_meta_cell_matrices()` consumes:
#' \itemize{
#'  \item indptr - Integer vector. 0-indexed CSR row pointer.
#'  \item indices - Integer vector. 0-indexed column indices.
#'  \item raw_counts - Numeric vector. Raw counts.
#'  \item norm_counts - Numeric vector. Normalised counts.
#'  \item nrow - Integer. Number of merged meta cells.
#'  \item ncol - Integer. Number of genes.
#' }
#'
#' @keywords internal
.merge_mc_counts <- function(inputs, gene_maps, n_genes) {
  checkmate::assertList(inputs, min.len = 1L)
  checkmate::assertList(gene_maps, len = length(inputs))
  checkmate::qassert(n_genes, "X1[1,)")

  n_rows <- purrr::map_int(inputs, \(x) nrow(S7::prop(x, "data")$raw))

  # sources off one parent share its gene space, so the map is the identity and
  # the column indices need no translation at all. That is the common case.
  identity_map <- purrr::map_lgl(
    gene_maps,
    \(m) length(m) == n_genes && !anyNA(m) && all(m == seq_len(n_genes))
  )

  # pass 1: what each input contributes, so the output is allocated rather than
  # grown. `map_dbl` keeps the sum out of integer range while we check it.
  nnz_kept <- purrr::map_dbl(seq_along(inputs), \(k) {
    raw <- S7::prop(inputs[[k]], "data")$raw
    if (identity_map[[k]]) {
      length(raw@x)
    } else {
      sum(!is.na(gene_maps[[k]][raw@j + 1L]))
    }
  })

  total_nnz <- sum(nnz_kept)
  if (total_nnz > .Machine$integer.max) {
    stop(sprintf(
      paste(
        "The merged matrices would hold %.0f stored values, which overflows",
        "the integer index of a sparse matrix. Merge fewer sources."
      ),
      total_nnz
    ))
  }

  # pass 2: fill each input's span
  indptr <- integer(sum(n_rows) + 1L)
  indices <- integer(total_nnz)
  raw_counts <- numeric(total_nnz)
  norm_counts <- numeric(total_nnz)

  row_offset <- 0L
  nnz_offset <- 0L

  for (k in seq_along(inputs)) {
    raw <- S7::prop(inputs[[k]], "data")$raw
    norm <- S7::prop(inputs[[k]], "data")$norm
    # raw and norm are built off one index structure; reuse it for both
    checkmate::assertTRUE(identical(raw@j, norm@j))

    if (identity_map[[k]]) {
      j_new <- raw@j
      p_new <- raw@p[-1L]
      x_raw <- raw@x
      x_norm <- norm@x
    } else {
      # `gene_maps` is 1-indexed, the CSR indices are 0-indexed
      j_new <- gene_maps[[k]][raw@j + 1L] - 1L
      keep <- !is.na(j_new)
      # kept values per row straight off the CSR pointer, no row index needed
      p_new <- c(0L, cumsum(keep))[raw@p[-1L] + 1L]
      j_new <- j_new[keep]
      x_raw <- raw@x[keep]
      x_norm <- norm@x[keep]

      # a non-monotone map breaks the within-row ordering of the column indices,
      # which only `feature_space = "union"` over differing gene orders can do
      map_kept <- gene_maps[[k]][!is.na(gene_maps[[k]])]
      if (is.unsorted(map_kept, strictly = TRUE)) {
        row_ids <- rep.int(seq_along(p_new), diff(c(0L, p_new)))
        ord <- order(row_ids, j_new)
        j_new <- j_new[ord]
        x_raw <- x_raw[ord]
        x_norm <- x_norm[ord]
      }
    }

    n_new <- length(j_new)
    if (n_new > 0L) {
      span <- nnz_offset + seq_len(n_new)
      indices[span] <- j_new
      raw_counts[span] <- x_raw
      norm_counts[span] <- x_norm
    }
    indptr[row_offset + seq_len(n_rows[[k]]) + 1L] <- p_new + nnz_offset

    row_offset <- row_offset + n_rows[[k]]
    nnz_offset <- nnz_offset + n_new
  }

  list(
    indptr = indptr,
    indices = indices,
    raw_counts = raw_counts,
    norm_counts = norm_counts,
    nrow = sum(n_rows),
    ncol = as.integer(n_genes)
  )
}

## hotspot ---------------------------------------------------------------------

#' Resolve the shared inputs of the meta cell HotSpot methods
#'
#' @description
#' Both HotSpot entry points need the same four things resolved off the object
#' before anything can go to Rust. Kept in one place so the two methods cannot
#' drift apart.
#'
#' @param object `MetaCells` class.
#' @param embd_to_use String. The embedding to use.
#' @param use_knn Boolean. Shall the cached kNN graph be used. Forced to `FALSE`
#' when `cells_to_take` is provided, since the cached graph covers all meta
#' cells.
#' @param no_embd_to_use Optional integer. Number of embedding dimensions to
#' use. If `NULL` all will be used.
#' @param cells_to_take Optional string vector. Meta cell identifiers. If `NULL`
#' all meta cells are used.
#' @param genes_to_take Optional string vector. If `NULL` all genes are used.
#'
#' @returns A list with the following elements:
#' \itemize{
#'   \item embd - The embedding matrix, subset to `cells_to_take`.
#'   \item cells_to_take - The resolved meta cell identifiers.
#'   \item genes_to_take - The resolved gene identifiers.
#'   \item knn_data - The cached kNN object or `NULL`.
#' }
#'
#' @keywords internal
.prep_mc_hotspot <- function(
  object,
  embd_to_use,
  use_knn,
  no_embd_to_use,
  cells_to_take,
  genes_to_take
) {
  checkmate::assertTRUE(S7::S7_inherits(object, MetaCells))
  checkmate::assertTRUE(embd_to_use %in% get_available_embeddings(object))

  embd <- get_embedding(x = object, embd_name = embd_to_use)

  if (!is.null(no_embd_to_use)) {
    to_take <- min(c(no_embd_to_use, ncol(embd)))
    embd <- embd[, 1:to_take]
  }

  if (is.null(cells_to_take)) {
    cells_to_take <- S7::prop(object, "obs_table")$meta_cell_id
  } else {
    # if the user overwrites this specifically, recreate the kNN data internally
    use_knn <- FALSE
  }

  if (is.null(genes_to_take)) {
    genes_to_take <- S7::prop(object, "var_table")$gene_id
  }

  if (length(intersect(rownames(embd), cells_to_take)) == 0) {
    stop(
      paste(
        "There is no intersection in the provided meta cell names and the",
        "embedding space. Please double check your parameters!"
      )
    )
  }

  list(
    embd = embd[cells_to_take, ],
    cells_to_take = cells_to_take,
    genes_to_take = genes_to_take,
    knn_data = if (use_knn) get_knn_obj(object) else NULL
  )
}
