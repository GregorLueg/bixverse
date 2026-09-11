# cellsweep --------------------------------------------------------------------

## helpers ---------------------------------------------------------------------

#' Largest smallest-library-size that still looks like raw ingest
#'
#' @description A raw, unfiltered barcode list always contains near-empty
#' droplets. If the minimum library size across the whole object is above this,
#' a QC cutoff was applied at load time and the empty droplets are gone.
#'
#' @keywords internal
CELLSWEEP_MAX_MIN_LIB_SIZE <- 100L

#' Resolve the empty droplet mask
#'
#' @description Either reads the mask out of obs or hands the library sizes to
#' Rust to infer it. Kept separate from [cellsweep_sc()] so the resolution is
#' testable on an obs table alone.
#'
#' @param obs data.table. The unfiltered obs table, with `lib_size`.
#' @param empty_params List. See [params_sc_empty_droplets()]. Required: there
#' is no safe default, since the recommended `method = "supplied"` needs the
#' name of the obs column holding the mask.
#' @param .verbose Logical. Controls verbosity.
#'
#' @returns Logical vector, `TRUE` where the barcode is an empty droplet.
#'
#' @keywords internal
.resolve_empty_droplets <- function(obs, empty_params, .verbose = TRUE) {
  if (identical(empty_params$method, "supplied")) {
    col <- empty_params$is_empty_column
    if (!col %in% names(obs)) {
      stop(sprintf("'%s' is not a column in the obs table.", col))
    }
    is_empty <- obs[[col]]
    if (!is.logical(is_empty)) {
      stop(sprintf("'%s' must be a logical column.", col))
    }
    if (anyNA(is_empty)) {
      stop(sprintf("'%s' contains NAs.", col))
    }
    return(is_empty)
  }

  is_empty <- rs_sc_infer_empty_droplets(
    lib_size = as.integer(obs$lib_size),
    empty_params = unclass(empty_params)
  )

  if (.verbose) {
    message(sprintf(
      "Empty droplets (method '%s'): %s of %s barcodes called empty.",
      empty_params$method,
      format(sum(is_empty), big.mark = ","),
      format(nrow(obs), big.mark = ",")
    ))
  }

  is_empty
}

#' Assemble the per-barcode CellSweep diagnostics
#'
#' @description Stitches the per-sample fits back into one table in output
#' order. Per-sample scalars are broadcast across that sample's barcodes, which
#' keeps them queryable in obs without needing a second table.
#'
#' @param res List. The Rust return value.
#' @param celltype_levels Character vector. Factor levels the `z_hat` codes
#' index into.
#'
#' @returns A data.table with one row per written barcode.
#'
#' @keywords internal
.cellsweep_obs_diagnostics <- function(res, celltype_levels) {
  parts <- lapply(seq_along(res$fits), function(i) {
    fit <- res$fits[[i]]
    data.table::data.table(
      cellsweep_alpha = as.numeric(fit$alpha),
      cellsweep_z = celltype_levels[as.integer(fit$z_hat)],
      cellsweep_beta = as.numeric(fit$beta),
      cellsweep_ll = as.numeric(fit$log_likelihood),
      cellsweep_converged = as.logical(fit$converged)
    )
  })

  out <- data.table::rbindlist(parts)
  stopifnot(nrow(out) == length(res$cell_order))
  out
}

#' Average the per-sample ambient profiles
#'
#' @description The ambient profile is per-emulsion, so a single var column can
#' only carry a summary. The unweighted mean across samples is that summary; the
#' per-sample profiles are not persisted.
#'
#' @param fits List. One entry per sample, each with an `ambient` vector.
#'
#' @returns Numeric vector, one entry per gene.
#'
#' @keywords internal
.cellsweep_mean_ambient <- function(fits) {
  profiles <- vapply(
    fits,
    function(fit) as.numeric(fit$ambient),
    numeric(length(fits[[1L]]$ambient))
  )
  if (length(fits) == 1L) {
    return(as.numeric(profiles))
  }
  rowMeans(profiles)
}

## generic ---------------------------------------------------------------------

#' Remove ambient and bulk contamination with CellSweep
#'
#' @description
#' Fits the CellSweep multinomial mixture model and writes a new `SingleCells`
#' object holding the denoised counts. Every observed count is split three ways
#' by EM: ambient contamination drawn from a per-emulsion profile, bulk
#' contamination drawn from a global profile, and true expression drawn from the
#' barcode's cell-type profile. Subtracting the first two gives the denoised
#' matrix.
#'
#' @details
#' Two things about where this sits in the workflow, because both are easy to
#' get wrong.
#'
#' **It runs after annotation, not before.** The model subtracts against
#' cell-type profiles, so it needs the labels on input. The chain is: ingest the
#' raw barcodes, mask and cluster and annotate as usual, then `cellsweep_sc()`,
#' then redo feature selection and reduction on the clean counts.
#'
#' **It needs the empty droplets.** The ambient profile is estimated from their
#' pooled counts, and they stay in the EM with their contamination fraction
#' pinned at 1. That means the object has to have been ingested permissively:
#'
#' ```r
#' obj <- load_mtx(
#'   obj,
#'   sc_mtx_io_param = get_cell_ranger_params(path),
#'   sc_qc_param = params_sc_min_quality(
#'     min_unique_genes = 0L, min_lib_size = 0L, min_cells = 0L
#'   )
#' )
#' ```
#'
#' The load-time cutoffs are irreversible, so the defaults
#' (`min_lib_size = 250L`) delete exactly the barcodes CellSweep trains on. This
#' function checks for that and errors rather than fitting a garbage profile.
#'
#' One EM fit per sample, since the ambient profile is a property of a single
#' emulsion. `sample_column` is required for that reason; pooling samples into
#' one ambient profile is wrong even when it runs.
#'
#' Barcodes partition three ways: empty droplets, annotated barcodes that passed
#' QC, and everything else. The third group is excluded from the fit and does
#' not appear in the output.
#'
#' The denoised counts land in a new directory. The raw layer takes the
#' stochastically rounded values so the negative binomial methods downstream
#' still see integers, and the normalised layer keeps the float magnitudes.
#' Entries whose denoised value rounds to zero are dropped from both layers:
#' the two share one index set, so keeping them would put explicit zeros in the
#' raw layer and `library_size` would stop being the sum of what is stored,
#' breaking every consumer that computes a fraction of the library.
#'
#' @param target `SingleCells`. A fresh object pointing at the directory the
#' denoised counts should be written to.
#' @param input `SingleCells`. The raw object, with the empty droplets still in
#' it.
#' @param celltype_column String. Obs column with the cell-type labels.
#' @param sample_column String. Obs column identifying the emulsion. One
#' independent fit per level.
#' @param empty_params List. See [params_sc_empty_droplets()]. Required: there
#' is no safe default, since the recommended `method = "supplied"` needs the
#' name of the obs column holding the mask.
#' @param cellsweep_params List. See [params_sc_cellsweep()].
#' @param sc_qc_param List. See [params_sc_min_quality()]. Only `target_size` is
#' used, to scale the new normalised layer.
#' @param streaming Integer. `0L` in-memory, `1L` light streaming (default) or
#' `2L` memory-bounded, for the CSR to CSC conversion.
#' @param batch_size Integer. Cells per batch for `streaming = 1L`.
#' @param max_genes_in_memory Integer. Genes held at once for `streaming = 2L`.
#' @param cell_batch_size Integer. Cells per batch for `streaming = 2L`.
#' @param .verbose Boolean or integer. Controls verbosity and returns run times.
#' `FALSE` -> quiet, `TRUE` or `1L` -> normal verbosity, `2L` -> detailed
#' verbosity.
#'
#' @returns The `target` object, populated with the denoised counts. The fit
#' lands in obs as `cellsweep_alpha`, `cellsweep_z`, `cellsweep_beta`,
#' `cellsweep_ll` and `cellsweep_converged`, and in var as `cellsweep_ambient`.
#'
#' @references Sullivan et al., CellSweep, 2025
#'
#' @export
cellsweep_sc <- S7::new_generic(
  name = "cellsweep_sc",
  dispatch_args = "target",
  fun = function(
    target,
    input,
    celltype_column,
    sample_column,
    empty_params,
    cellsweep_params = params_sc_cellsweep(),
    sc_qc_param = params_sc_min_quality(),
    streaming = 1L,
    batch_size = 1000L,
    max_genes_in_memory = 2000L,
    cell_batch_size = 100000L,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

## method ----------------------------------------------------------------------

#' @method cellsweep_sc SingleCells
#'
#' @export
S7::method(cellsweep_sc, SingleCells) <- function(
  target,
  input,
  celltype_column,
  sample_column,
  empty_params,
  cellsweep_params = params_sc_cellsweep(),
  sc_qc_param = params_sc_min_quality(),
  streaming = 1L,
  batch_size = 1000L,
  max_genes_in_memory = 2000L,
  cell_batch_size = 100000L,
  .verbose = TRUE
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(target, SingleCells))
  checkmate::assertTRUE(S7::S7_inherits(input, SingleCells))
  checkmate::qassert(celltype_column, "S1")
  checkmate::qassert(sample_column, "S1")
  assertScEmptyDroplets(empty_params)
  assertScCellsweep(cellsweep_params)
  assertScMinQC(sc_qc_param)
  checkmate::qassert(streaming, "I1")
  checkmate::assertTRUE(streaming %in% c(0L, 1L, 2L))
  checkmate::qassert(batch_size, "I1")
  checkmate::qassert(.verbose, c("B1", "I1[0,2]"))

  if (identical(S7::prop(target, "dir_data"), S7::prop(input, "dir_data"))) {
    stop("`target` and `input` must point at different directories.")
  }

  duckdb_in <- get_sc_duckdb(input)
  obs_cols <- duckdb_in$get_obs_cols()
  for (col in c(celltype_column, sample_column)) {
    if (!col %in% obs_cols) {
      stop(sprintf("'%s' is not a column in the obs table.", col))
    }
  }

  # unfiltered obs: the empty droplets are exactly what cells_to_keep excludes
  obs <- get_sc_obs(input, filtered = FALSE)
  data.table::setDT(obs)
  n_barcodes <- nrow(obs)

  # The load-time cutoffs are irreversible, so an object ingested with the
  # default QC has no empty droplets left and the ambient profile would be
  # fitted to noise.
  if (!"lib_size" %in% names(obs)) {
    stop("The obs table has no `lib_size` column; re-ingest the object.")
  }
  if (isTRUE(cellsweep_params$freeze_ambient_profile)) {
    min_lib <- min(obs$lib_size, na.rm = TRUE)
    if (min_lib > CELLSWEEP_MAX_MIN_LIB_SIZE) {
      stop(sprintf(
        paste(
          "The smallest library size in the object is %d, so the empty",
          "droplets were filtered out at ingest and the ambient profile",
          "cannot be estimated. Re-ingest with",
          "`params_sc_min_quality(min_unique_genes = 0L, min_lib_size = 0L,",
          "min_cells = 0L)`, or set `freeze_ambient_profile = FALSE`."
        ),
        min_lib
      ))
    }
  }

  is_empty <- .resolve_empty_droplets(obs, empty_params, .verbose)

  keep_mask <- logical(n_barcodes)
  keep_mask[as.integer(get_cells_to_keep(input)) + 1L] <- TRUE

  celltype <- as.factor(obs[[celltype_column]])
  celltype_levels <- levels(celltype)
  annotated <- !is.na(celltype)

  real_mask <- keep_mask & annotated & !is_empty
  dropped <- sum(!real_mask & !is_empty)

  if (.verbose && dropped > 0L) {
    message(sprintf(
      paste(
        "%d barcodes are neither empty nor annotated-and-passing-QC.",
        "They are excluded from the fit and from the output."
      ),
      dropped
    ))
  }
  if (!any(real_mask)) {
    stop("No barcode is both annotated and passing QC.")
  }

  sample_id <- as.character(obs[[sample_column]])
  sample_levels <- sort(unique(sample_id[real_mask]))

  samples <- lapply(sample_levels, function(sid) {
    in_sample <- sample_id == sid
    real_idx <- which(in_sample & real_mask)
    empty_idx <- which(in_sample & is_empty)
    list(
      sample_id = sid,
      real_cells = as.integer(real_idx - 1L),
      empty_cells = as.integer(empty_idx - 1L),
      celltype_idx = as.integer(celltype[real_idx]) - 1L,
      n_celltypes = length(celltype_levels)
    )
  })

  if (.verbose) {
    message(sprintf(
      "Running CellSweep over %d samples, %d barcodes, %d empty droplets.",
      length(samples),
      sum(real_mask),
      sum(is_empty)
    ))
  }

  rust_con <- get_sc_rust_ptr(target)
  res <- rust_con$cellsweep(
    f_path_source = get_rust_count_cell_f_path(input),
    samples = samples,
    cellsweep_params = unclass(cellsweep_params),
    target_size = sc_qc_param$target_size,
    verbose = parse_verbosity(.verbose)
  )

  if (.verbose) {
    message("Generating gene-based binary.")
  }

  if (streaming == 1L) {
    rust_con$generate_gene_based_data_streaming(
      batch_size = batch_size,
      verbose = as.logical(.verbose)
    )
  } else if (streaming == 2L) {
    rust_con$generate_gene_based_data_memory_bounded(
      max_genes_in_memory = max_genes_in_memory,
      cell_batch_size = cell_batch_size,
      verbose = as.logical(.verbose)
    )
  } else {
    rust_con$generate_gene_based_data(verbose = as.logical(.verbose))
  }

  if (.verbose) {
    message("Populating obs and var tables.")
  }

  duckdb_out <- get_sc_duckdb(target)
  source_db <- file.path(S7::prop(input, "dir_data"), "sc_duckdb.db")

  duckdb_out$populate_obs_from_duckdb_subset(
    source_db_path = source_db,
    cell_idx_to_keep = as.integer(res$cell_order) + 1L
  )

  duckdb_out$populate_vars_from_duckdb_reordered(
    source_db_path = source_db,
    final_gene_names = duckdb_in$get_vars_table()$gene_id
  )

  duckdb_out$add_data_obs(
    new_data = .cellsweep_obs_diagnostics(res, celltype_levels)
  )
  duckdb_out$add_data_var(
    new_data = data.table::data.table(
      no_cells_exp = rust_con$get_nnz_genes(gene_indices = NULL),
      cellsweep_ambient = .cellsweep_mean_ambient(res$fits)
    )
  )
  duckdb_out$set_to_keep_column()

  S7::prop(target, "dims") <- as.integer(rust_con$get_shape())
  target <- set_cell_mapping(
    x = target,
    cell_map = duckdb_out$get_obs_index_map()
  )
  target <- set_gene_mapping(
    x = target,
    gene_map = duckdb_out$get_var_index_map()
  )

  return(target)
}
