# ------------------------------------------------------------------------------
# tree-based representations of single cells. unlike the trajectory methods
# these infer a hierarchy of cell states over every cell, not a process with a
# start and a direction, so they live here rather than next to Palantir and
# PAGA.
# ------------------------------------------------------------------------------

# tree methods -----------------------------------------------------------------

## bonsai ----------------------------------------------------------------------

#' Build a Bonsai tree over the cells
#'
#' @description
#' Bonsai reconstructs a tree over the cells in which every cell is a leaf and
#' the internal nodes are inferred ancestral states, with branch lengths that
#' carry the amount of change between them. Unlike a kNN graph it uses each
#' measurement's error bar, which is what Sanity provides: the raw counts go
#' through Sanity first, for posterior log fold changes with error bars, and
#' Bonsai builds the tree on those. The tree is then laid out in 2D. For
#' details, please refer to de Groot, et al. and Breda, et al.
#'
#' No HVG selection needed. By default every gene goes in, streamed through
#' Sanity in chunks, and only the genes with enough signal over their own noise
#' (`min_signal_to_noise` in [bixverse::params_sc_bonsai()]) are kept for the
#' tree. That is the gene selection of the Bonsai paper, and it keeps memory at
#' one chunk plus the survivors. Pass `hvg` to restrict the candidates.
#'
#' Runtime grows a little faster than linearly with the number of cells. In
#' bonsai-rs's own benchmarks the search took about 70 seconds at 10,000 cells
#' and 210 seconds at 25,000 on ten cores. Sanity comes on top, linear in the
#' number of genes it has to fit; on the CPU that is the larger share once all
#' genes go in.
#'
#' The object itself is not touched. Store the leaf coordinates with
#' [bixverse::set_bonsai_embedding()] if you want them next to the other
#' embeddings.
#'
#' @param object `SingleCells` class.
#' @param hvg Optional integer. Restrict the candidate genes to these, e.g. the
#' output of [bixverse::get_hvg()] plus one. Please provide 1-indexed genes
#' here! If `NULL`, every gene in the object is a candidate.
#' @param bonsai_params List. See [bixverse::params_sc_bonsai()].
#' @param .verbose Boolean or integer. Controls verbosity and returns run times.
#' `FALSE` -> quiet, `TRUE` or `1L` -> normal verbosity, `2L` -> detailed
#' verbosity.
#'
#' @returns A `BonsaiTree` S3 object with:
#' \itemize{
#'   \item nodes - data.table with `node`, `parent` (`NA` for the root),
#'   `branch`, `is_leaf`, `cell_id` (`NA` for the inferred ancestors), `x` and
#'   `y`. The leaves come first, in the order of the object's cells.
#'   \item loglik - The loglikelihood of the final tree.
#'   \item steps - data.table with the loglikelihood after each search step
#'   and the step's wall time in seconds.
#'   \item timings - data.table with the wall time in seconds of each stage:
#'   `sanity`, `ingest`, `bonsai` (the whole search), `layout`, and `total`, the
#'   whole call as R saw it.
#'   \item genes_used - The genes the tree was built on.
#'   \item genes_dropped - The candidate genes left out: no counts in the
#'   cells, ill-conditioned Sanity posteriors, or a signal-to-noise ratio below
#'   `min_signal_to_noise`.
#'   \item cell_idx - The cells the tree was built over (0-indexed).
#'   \item layout, hyperbolic - The current layout.
#'   \item params - The parameters of the run.
#' }
#'
#' @references de Groot, et al., Nat. Biotechnol., 2026; Breda, et al., Nat.
#' Biotechnol., 2021.
#'
#' @export
#'
#' @examples
#' # a tree over the demo cells, genes selected by their signal-to-noise
#' sc <- demo_single_cells(prepped = FALSE)
#' tree <- bonsai_sc(sc, .verbose = FALSE)
#' tree
#'
#' unlink(sc@dir_data, recursive = TRUE, force = TRUE)
bonsai_sc <- S7::new_generic(
  name = "bonsai_sc",
  dispatch_args = "object",
  fun = function(
    object,
    hvg = NULL,
    bonsai_params = params_sc_bonsai(),
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method bonsai_sc SingleCells
#'
#' @export
S7::method(bonsai_sc, SingleCells) <- function(
  object,
  hvg = NULL,
  bonsai_params = params_sc_bonsai(),
  .verbose = TRUE
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(object, SingleCells))

  .bonsai_sc_run(
    object = object,
    hvg = hvg,
    bonsai_params = bonsai_params,
    runner = rs_sc_bonsai,
    .verbose = .verbose
  )
}

#' Run Bonsai with a given Rust entry point
#'
#' @description
#' Everything [bixverse::bonsai_sc()] does around the Rust call: candidate
#' genes, cells, the `BonsaiTree` and the total timing. The Rust entry point is
#' an argument so `bixverse.gpu` can hand in its GPU Sanity one and get back the
#' identical class.
#'
#' @param object `SingleCells` class.
#' @param hvg Optional integer. 1-indexed candidate genes, `NULL` for all.
#' @param bonsai_params List. See [bixverse::params_sc_bonsai()].
#' @param runner Function with the signature of [bixverse::rs_sc_bonsai()].
#' @param .verbose Boolean or integer. Controls verbosity.
#'
#' @returns A `BonsaiTree`, see [bixverse::bonsai_sc()].
#'
#' @keywords internal
.bonsai_sc_run <- function(object, hvg, bonsai_params, runner, .verbose) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(object, SingleCells))
  checkmate::qassert(hvg, c("I+", "0"))
  assertScBonsaiParams(bonsai_params)
  checkmate::assertFunction(
    runner,
    args = c(
      "f_path_gene",
      "f_path_cell",
      "cell_indices",
      "gene_indices",
      "bonsai_params",
      "verbose"
    )
  )
  checkmate::qassert(.verbose, c("B1", "I1[0,2]"))

  genes_in <- if (!is.null(hvg)) {
    hvg - 1L
  } else {
    get_gene_indices(object, get_gene_names(object), rust_index = TRUE)
  }
  # ascending, so genes_used comes back in gene order
  genes_in <- sort(unique(as.integer(genes_in)))

  cell_idx <- get_cells_to_keep(object)

  if (.verbose) {
    message(sprintf(
      "Running Sanity and Bonsai over %i cells and %i candidate genes.",
      length(cell_idx),
      length(genes_in)
    ))
  }

  started <- Sys.time()
  rs_res <- runner(
    f_path_gene = get_rust_count_gene_f_path(object),
    f_path_cell = get_rust_count_cell_f_path(object),
    cell_indices = cell_idx,
    gene_indices = genes_in,
    bonsai_params = bonsai_params,
    verbose = parse_verbosity(.verbose)
  )

  tree <- new_bonsai_tree(
    rs_res = rs_res,
    cell_idx = cell_idx,
    cell_names = get_cell_names(object, filtered = TRUE),
    genes_in = genes_in,
    gene_ids = unname(get_gene_names_from_idx(object, genes_in)),
    bonsai_params = bonsai_params
  )

  # anything the Rust stages do not account for shows up as the difference
  tree$timings <- rbind(
    tree$timings,
    data.table::data.table(
      stage = "total",
      seconds = as.numeric(difftime(Sys.time(), started, units = "secs"))
    )
  )

  tree
}
