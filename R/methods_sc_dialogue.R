# dialogue ---------------------------------------------------------------------

# generic found in base_generics_sc.R

## constants -------------------------------------------------------------------

# Mirrors MIN_SHARED_SAMPLES on the Rust side. Checked here so the error names
# the study design rather than surfacing from the decomposition.
DIALOGUE_MIN_SHARED_SAMPLES <- 5L

## helpers ---------------------------------------------------------------------

#' Resolve the per-cell-type inputs DIALOGUE needs
#'
#' @description
#' Turns the obs table plus the user's feature matrices into the flat, 0-indexed
#' vectors the Rust side takes. The cell types analysed are `names(features)`,
#' in that order: anything else in `cell_type_col` is ignored, and that order is
#' what the pair columns of `emp_p` follow.
#'
#' @details
#' `sample_ids` and `cell_quality` are indexed by *global* cell position, not by
#' position within a cell type, so both are scattered into vectors long enough
#' to cover the largest index in use. Cells belonging to no analysed cell type
#' leave holes, which Rust never reads.
#'
#' @param object `SingleCells`, `SingleCellsSubset` or `MetaCells` class.
#' @param obs data.table. The observation table to read the labels off.
#' @param cell_type_col String. Column holding the cell type labels.
#' @param sample_col String. Column holding the sample labels.
#' @param features Named list of numeric matrices, one per cell type.
#' @param quality_col Optional string. Column to use as the quality covariate.
#' @param cell_id_col String. Column holding the cell identifiers.
#' @param default_quality Numeric vector aligned to `obs` rows. Used when
#' `quality_col` is `NULL`, after a `log1p` and a z-score. Typically the library
#' size.
#'
#' @returns A list with the following items
#' \itemize{
#'   \item cell_type_indices - List of integer vectors. 0-indexed global cell
#'   positions per cell type.
#'   \item features - List of numeric matrices, rows reordered to match.
#'   \item sample_ids - Integer vector. 0-indexed sample code per global cell.
#'   \item cell_quality - Numeric vector. Quality covariate per global cell.
#'   \item cell_types - Character vector. The analysed cell types, in order.
#'   \item cell_ids - List of character vectors. Cell ids per cell type.
#'   \item sample_levels - Character vector. Sample labels, position `i` being
#'   sample code `i - 1`.
#' }
#'
#' @keywords internal
.prep_dialogue_inputs <- function(
  object,
  obs,
  cell_type_col,
  sample_col,
  features,
  quality_col,
  cell_id_col,
  default_quality
) {
  # checks
  checkmate::assertDataFrame(obs)
  checkmate::qassert(cell_type_col, "S1")
  checkmate::qassert(sample_col, "S1")
  checkmate::assertList(features, types = "matrix", names = "named")
  checkmate::qassert(quality_col, c("S1", "0"))
  checkmate::qassert(cell_id_col, "S1")
  checkmate::qassert(default_quality, "N+")

  checkmate::assertNames(
    colnames(obs),
    must.include = c(cell_id_col, cell_type_col, sample_col)
  )
  if (!is.null(quality_col)) {
    checkmate::assertNames(colnames(obs), must.include = quality_col)
  }
  checkmate::assertTRUE(length(default_quality) == nrow(obs))

  cell_types <- names(features)
  if (length(cell_types) < 2L) {
    stop(sprintf(
      paste(
        "DIALOGUE needs at least two cell types, `features` names %i.",
        "Nothing to correlate across."
      ),
      length(cell_types)
    ))
  }
  if (anyDuplicated(cell_types) > 0L) {
    stop("`features` has duplicated names.")
  }

  observed <- unique(as.character(obs[[cell_type_col]]))
  missing_types <- setdiff(cell_types, observed)
  if (length(missing_types) > 0L) {
    stop(sprintf(
      "`features` names cell types absent from `%s`: %s.",
      cell_type_col,
      paste(missing_types, collapse = ", ")
    ))
  }

  labels <- as.character(obs[[cell_type_col]])
  ids <- as.character(obs[[cell_id_col]])
  samples <- as.character(obs[[sample_col]])

  quality_raw <- if (is.null(quality_col)) {
    default_quality
  } else {
    obs[[quality_col]]
  }
  checkmate::qassert(quality_raw, "N+")

  cell_ids <- purrr::map(cell_types, \(ct) ids[labels == ct])
  names(cell_ids) <- cell_types

  cell_type_indices <- purrr::map(cell_types, \(ct) {
    as.integer(get_cell_indices(
      x = object,
      cell_ids = cell_ids[[ct]],
      rust_index = TRUE
    ))
  })
  names(cell_type_indices) <- cell_types

  features_ordered <- purrr::map(cell_types, \(ct) {
    .assert_dialogue_features(features[[ct]], cell_ids[[ct]], ct)
  })
  names(features_ordered) <- cell_types

  # samples have to be present in every cell type, otherwise there is nothing
  # to line the blocks up on
  per_type_samples <- purrr::map(
    cell_types,
    \(ct) unique(samples[labels == ct])
  )
  shared <- Reduce(intersect, per_type_samples)
  if (length(shared) < DIALOGUE_MIN_SHARED_SAMPLES) {
    stop(sprintf(
      paste(
        "Only %i sample(s) are present in every cell type, DIALOGUE needs at",
        "least %i. Check `%s` and whether the rarer cell types were captured",
        "in every sample."
      ),
      length(shared),
      DIALOGUE_MIN_SHARED_SAMPLES,
      sample_col
    ))
  }

  sample_levels <- sort(unique(samples))
  sample_codes <- match(samples, sample_levels) - 1L

  # global (0-indexed) space, so the vectors have to reach the largest index
  n_slots <- max(unlist(cell_type_indices)) + 1L
  sample_ids <- integer(n_slots)
  cell_quality <- numeric(n_slots)

  quality_final <- if (is.null(quality_col)) {
    as.numeric(scale(log1p(quality_raw)))
  } else {
    as.numeric(quality_raw)
  }

  for (ct in cell_types) {
    rows <- which(labels == ct)
    slots <- cell_type_indices[[ct]] + 1L
    sample_ids[slots] <- sample_codes[rows]
    cell_quality[slots] <- quality_final[rows]
  }

  list(
    cell_type_indices = unname(cell_type_indices),
    features = unname(features_ordered),
    sample_ids = sample_ids,
    cell_quality = cell_quality,
    cell_types = cell_types,
    cell_ids = cell_ids,
    sample_levels = sample_levels
  )
}

#' Validate one cell type's feature matrix and align its rows
#'
#' @param feature Numeric matrix. The cell type's features.
#' @param cell_ids Character. The cell type's cell identifiers, in the order the
#' indices were resolved in.
#' @param cell_type String. Name of the cell type, for the error message.
#'
#' @returns The matrix, rows reordered to `cell_ids`.
#'
#' @keywords internal
.assert_dialogue_features <- function(feature, cell_ids, cell_type) {
  checkmate::qassert(cell_ids, "S+")
  checkmate::qassert(cell_type, "S1")

  if (!is.matrix(feature) || !is.numeric(feature)) {
    stop(sprintf("`features[['%s']]` is not a numeric matrix.", cell_type))
  }
  if (ncol(feature) < 2L) {
    stop(sprintf(
      "`features[['%s']]` has %i column(s); DIALOGUE needs at least 2.",
      cell_type,
      ncol(feature)
    ))
  }
  if (is.null(rownames(feature))) {
    stop(sprintf(
      paste(
        "`features[['%s']]` has no row names. Rows are matched to cells by",
        "name, not by position."
      ),
      cell_type
    ))
  }

  missing_cells <- setdiff(cell_ids, rownames(feature))
  if (length(missing_cells) > 0L) {
    stop(sprintf(
      "`features[['%s']]` is missing %i of the %i cells of that cell type.",
      cell_type,
      length(missing_cells),
      length(cell_ids)
    ))
  }

  out <- feature[cell_ids, , drop = FALSE]
  if (anyNA(out)) {
    stop(sprintf("`features[['%s']]` holds NAs.", cell_type))
  }

  out
}

#' Resolve the genes DIALOGUE builds signatures from
#'
#' @param object `SingleCells`, `SingleCellsSubset` or `MetaCells` class.
#' @param gene_ids Optional character. Genes to use. If `NULL`, uses the HVGs.
#'
#' @returns Integer vector with the 0-indexed gene positions.
#'
#' @keywords internal
.resolve_dialogue_genes <- function(object, gene_ids) {
  checkmate::qassert(gene_ids, c("S+", "0"))

  indices <- if (is.null(gene_ids)) {
    get_hvg(object)
  } else {
    get_gene_indices(x = object, gene_ids = gene_ids, rust_index = TRUE)
  }

  if (length(indices) == 0L) {
    stop(paste(
      "No genes to run DIALOGUE on. Either run find_hvg_sc() first or pass",
      "`gene_ids` explicitly."
    ))
  }

  as.integer(indices)
}

## methods ---------------------------------------------------------------------

#' @method dialogue_sc ScOrScSubset
S7::method(dialogue_sc, ScOrScSubset) <- function(
  object,
  cell_type_col,
  sample_col,
  features,
  quality_col = NULL,
  gene_ids = NULL,
  pmd_params = params_dialogue_pmd(),
  hlm_params = params_dialogue_hlm(),
  refine_params = params_dialogue_refine(),
  .verbose = TRUE
) {
  # checks
  checkmate::assertTRUE(
    S7::S7_inherits(object, SingleCells) ||
      S7::S7_inherits(object, SingleCellsSubset)
  )
  checkmate::qassert(cell_type_col, "S1")
  checkmate::qassert(sample_col, "S1")
  checkmate::assertList(features, names = "named")
  checkmate::qassert(quality_col, c("S1", "0"))
  checkmate::qassert(gene_ids, c("S+", "0"))
  assertDialoguePmd(pmd_params)
  assertDialogueHlm(hlm_params)
  assertDialogueRefine(refine_params)
  checkmate::qassert(.verbose, c("B1", "I1[0,2]"))

  obs <- get_sc_obs(object, filtered = TRUE)

  # `lib_size` lands in obs at ingestion, but an object built by hand may not
  # carry it and the covariate is not optional downstream
  if (is.null(quality_col) && !("lib_size" %in% colnames(obs))) {
    stop(paste(
      "No `lib_size` column in the obs table to derive the cell quality",
      "covariate from. Pass `quality_col` instead."
    ))
  }

  prepped <- .prep_dialogue_inputs(
    object = object,
    obs = obs,
    cell_type_col = cell_type_col,
    sample_col = sample_col,
    features = features,
    quality_col = quality_col,
    cell_id_col = "cell_id",
    default_quality = as.numeric(obs$lib_size)
  )

  dialogue_res <- rs_dialogue_sc(
    f_path_gene = get_rust_count_gene_f_path(object),
    cell_type_indices = prepped$cell_type_indices,
    features = prepped$features,
    sample_ids = prepped$sample_ids,
    cell_quality = prepped$cell_quality,
    gene_indices = .resolve_dialogue_genes(object, gene_ids),
    dialogue_params = c(pmd_params, hlm_params, refine_params),
    verbose = parse_verbosity(.verbose)
  )

  new_dialogue_result(
    dialogue_res = dialogue_res,
    prepped = prepped,
    gene_names = get_gene_names(object),
    source_class = if (S7::S7_inherits(object, SingleCells)) {
      "SingleCells"
    } else {
      "SingleCellsSubset"
    },
    params = list(
      cell_type_col = cell_type_col,
      sample_col = sample_col,
      quality_col = quality_col,
      pmd_params = pmd_params,
      hlm_params = hlm_params,
      refine_params = refine_params
    )
  )
}
