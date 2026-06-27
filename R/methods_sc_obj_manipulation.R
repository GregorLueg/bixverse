# methods for merging / splitting single cell objects --------------------------

## merging ---------------------------------------------------------------------

#' Merge multiple `SingleCells` experiments into one
#'
#' @description
#' Merges N existing `SingleCells` objects into a freshly constructed target
#' object. The feature space of the result is the **intersection** of the input
#' gene sets. Each input's `cells_to_keep` filter is honoured (i.e. only cells
#' with `to_keep = TRUE` in the input's obs are carried over).
#'
#' If `renormalise = FALSE`, the stored `data_norm` values are copied through
#' unchanged. This is valid only when all inputs were normalised against the
#' same `target_size`. If the gene intersection is much smaller than the
#' individual input gene sets, the inherited `data_norm` becomes a lossy
#' approximation (it was computed against the pre-intersection library size). In
#' that case set `renormalise = TRUE` to recompute `data_norm` against the
#' surviving raw counts using `sc_qc_param$target_size`.
#'
#' Obs columns are intersected across inputs. The result obs gains an `exp_id`
#' column. Inputs that already have an `exp_id` column are rejected. The
#' `sc_cache` and `sc_map` of the target are populated fresh; any PCA, kNN, sNN
#' or HVG state on the inputs is not carried over and must be re-run.
#'
#' @param target A freshly constructed `SingleCells` pointing at the output
#' directory.
#' @param inputs List of `SingleCells` objects to merge. Length >= 2.
#' @param exp_ids Character vector of experiment identifiers, one per input.
#' Must be unique.
#' @param renormalise Boolean. Whether to recompute `data_norm` against
#' `sc_qc_param$target_size`. Defaults to `FALSE`.
#' @param sc_qc_param List. Output of [bixverse::params_sc_min_quality()].
#' Only `target_size` is consulted here; no QC filtering is applied during
#' merge.
#' @param streaming Integer. `0` -> no streaming, `1` -> light streaming, `2` ->
#' heavy streaming with memory upper boundaries. This enables you to control
#' the memory pressure during ingestion.
#' @param batch_size Integer. Batch size when `streaming = 1L`.
#' @param max_genes_in_memory Integer. How many genes shall be held in memory
#' at a given point. Defaults to `2000L`. Only relevant if streaming is set to
#' `2`.
#' @param cell_batch_size Integer. How big are the batch sizes for the cells
#' in the transformation from the cell-based to gene-based format. Defaults to
#' `100000L`. Only relevant if streaming is set to `2`.
#' @param .verbose Boolean.
#'
#' @return The populated target `SingleCells`.
#'
#' @export
merge_sc_experiments <- S7::new_generic(
  name = "merge_sc_experiments",
  dispatch_args = "target",
  fun = function(
    target,
    inputs,
    exp_ids,
    renormalise = FALSE,
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

#' @method merge_sc_experiments SingleCells
#'
#' @export
S7::method(merge_sc_experiments, SingleCells) <- function(
  target,
  inputs,
  exp_ids,
  renormalise = FALSE,
  sc_qc_param = params_sc_min_quality(),
  streaming = 1L,
  batch_size = 1000L,
  max_genes_in_memory = 2000L,
  cell_batch_size = 100000L,
  .verbose = TRUE
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(target, SingleCells))
  checkmate::assertList(inputs, min.len = 2L)
  for (inp in inputs) {
    checkmate::assertTRUE(S7::S7_inherits(inp, SingleCells))
  }
  checkmate::assertCharacter(exp_ids, len = length(inputs), unique = TRUE)
  checkmate::qassert(renormalise, "B1")
  assertScMinQC(sc_qc_param)
  checkmate::qassert(streaming, "I1")
  checkmate::assertTRUE(streaming %in% c(0L, 1L, 2L))
  checkmate::qassert(batch_size, "I1")
  checkmate::qassert(.verbose, "B1")

  # error on inputs that were themselves merged
  for (i in seq_along(inputs)) {
    obs_cols <- get_sc_duckdb(inputs[[i]])$get_obs_cols()
    if ("exp_id" %in% obs_cols) {
      stop(sprintf(
        paste(
          "Input %d ('%s') already has an exp_id column.",
          "Merging already-merged objects is not supported."
        ),
        i,
        exp_ids[i]
      ))
    }
  }

  # gene intersection universe
  if (.verbose) {
    message("Building gene intersection universe.")
  }
  input_gene_sets <- lapply(inputs, function(inp) {
    get_sc_duckdb(inp)$get_vars_table()$gene_id
  })
  universe <- Reduce(intersect, input_gene_sets)
  if (length(universe) == 0L) {
    stop("Gene intersection across inputs is empty.")
  }
  if (.verbose) {
    message(sprintf(
      "Universe size: %s genes (intersection of %d inputs).",
      format(length(universe), big.mark = ","),
      length(inputs)
    ))
  }

  # build merge tasks
  merge_tasks <- vector("list", length(inputs))
  for (i in seq_along(inputs)) {
    inp <- inputs[[i]]
    local_genes <- input_gene_sets[[i]]
    match_idx <- match(local_genes, universe)
    gene_local_to_universe <- as.integer(match_idx - 1L)

    cells_to_keep_0idx <- sort(as.integer(get_cells_to_keep(inp)))

    merge_tasks[[i]] <- list(
      exp_id = exp_ids[i],
      bin_cells_path = get_rust_count_cell_f_path(inp),
      cells_to_keep = cells_to_keep_0idx,
      gene_local_to_universe = as.integer(gene_local_to_universe)
    )
  }

  # Rust merge
  if (.verbose) {
    message("Merging binary cell files.")
  }
  rust_con <- get_sc_rust_ptr(target)

  merge_res <- rust_con$merge_sc_files(
    merge_tasks = merge_tasks,
    universe_size = length(universe),
    renormalise = renormalise,
    target_size = sc_qc_param$target_size,
    verbose = .verbose
  )

  # gene-based regeneration
  if (.verbose) {
    message("Generating gene-based binary.")
  }
  if (streaming == 1L) {
    if (.verbose) {
      message(" Using light streaming for the CSR to CSC conversion.")
    }
    rust_con$generate_gene_based_data_streaming(
      batch_size = batch_size,
      verbose = .verbose
    )
  } else if (streaming == 2L) {
    if (.verbose) {
      message(paste(
        " Using heavy streaming with reduced memory",
        "pressure for the CSR to CSC conversion."
      ))
    }
    rust_con$generate_gene_based_data_memory_bounded(
      max_genes_in_memory = max_genes_in_memory,
      cell_batch_size = cell_batch_size,
      verbose = .verbose
    )
  } else {
    if (.verbose) {
      message(paste(
        " Loading data directly into memory for CSR to CSC conversion."
      ))
    }
    rust_con$generate_gene_based_data(verbose = .verbose)
  }

  gene_nnz <- rust_con$get_nnz_genes(gene_indices = NULL)
  gene_nnz_dt <- data.table::data.table(no_cells_exp = gene_nnz)

  # populate target DuckDB
  duckdb_con <- get_sc_duckdb(target)

  if (.verbose) {
    message("Merging obs tables.")
  }
  per_file_info <- lapply(seq_along(inputs), function(i) {
    list(
      db_path = file.path(S7::prop(inputs[[i]], "dir_data"), "sc_duckdb.db"),
      exp_id = exp_ids[i]
    )
  })
  duckdb_con$populate_obs_from_multi_duckdb(per_file_info = per_file_info)

  if (.verbose) {
    message("Populating var table.")
  }
  first_db_path <- file.path(
    S7::prop(inputs[[1L]], "dir_data"),
    "sc_duckdb.db"
  )
  duckdb_con$populate_vars_from_duckdb_reordered(
    source_db_path = first_db_path,
    final_gene_names = universe
  )

  # add QC columns from merge result
  per_file_qc <- lapply(merge_res$per_file, function(f) {
    data.table::data.table(nnz = f$nnz, lib_size = f$lib_size)
  })
  cell_res_dt <- data.table::rbindlist(per_file_qc)

  duckdb_con$add_data_obs(new_data = cell_res_dt)
  duckdb_con$add_data_var(new_data = gene_nnz_dt)
  duckdb_con$set_to_keep_column()

  cell_map <- duckdb_con$get_obs_index_map()
  gene_map <- duckdb_con$get_var_index_map()

  S7::prop(target, "dims") <- as.integer(rust_con$get_shape())
  target <- set_cell_mapping(x = target, cell_map = cell_map)
  target <- set_gene_mapping(x = target, gene_map = gene_map)

  return(target)
}
