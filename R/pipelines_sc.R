# single cell pipeline ---------------------------------------------------------

## classes ---------------------------------------------------------------------

### step -----------------------------------------------------------------------

#' Constructor for a single pipeline step
#'
#' @description
#' Internal constructor used by the public `step_*()` functions. Wraps the
#' generic to dispatch on the `SingleCells`/`SingleCellsSubset` object together
#' with the captured arguments and a short human-readable name.
#'
#' @param name String. Short identifier shown when printing the pipeline.
#' @param fn Function. Typically an S7 generic such as [find_hvg_sc()] that
#' takes the object as the first argument.
#' @param args Named list. Arguments passed to `fn` at apply time (object is
#' prepended automatically).
#'
#' @return An `ScStep` object.
#'
#' @keywords internal
new_sc_step <- function(name, fn, args) {
  checkmate::qassert(name, "S1")
  checkmate::assertFunction(fn)
  checkmate::assertList(args, names = "unique")

  structure(
    list(name = name, fn = fn, args = args),
    class = "ScStep"
  )
}

### pipeline -------------------------------------------------------------------

#' Construct an empty single cell pipeline
#'
#' @description
#' Linear container of `ScStep`s. Append steps with [`%>>%`] and execute with
#' [apply_pipeline()]. Pipelines are inert until applied; steps can be
#' inspected via `pipeline$steps`.
#'
#' @return An empty `ScPipeline` object.
#'
#' @export
sc_pipeline <- function() {
  structure(list(steps = list()), class = "ScPipeline")
}

## operator --------------------------------------------------------------------

#' Append a step to a pipeline
#'
#' @description
#' `%>>%` chains pipeline steps. Either side can be a `ScStep` or a
#' `ScPipeline`; the result is always a `ScPipeline`.
#'
#' @param lhs `ScPipeline` or `ScStep`.
#' @param rhs `ScStep`.
#'
#' @return A `ScPipeline`.
#'
#' @export
`%>>%` <- function(lhs, rhs) UseMethod("%>>%")

#' @export
`%>>%.ScPipeline` <- function(lhs, rhs) {
  checkmate::assertClass(rhs, "ScStep")
  lhs$steps[[length(lhs$steps) + 1L]] <- rhs
  lhs
}

#' @export
`%>>%.ScStep` <- function(lhs, rhs) {
  checkmate::assertClass(rhs, "ScStep")
  structure(list(steps = list(lhs, rhs)), class = "ScPipeline")
}

## primitives ------------------------------------------------------------------

#' @export
length.ScPipeline <- function(x) length(x$steps)

#' @export
print.ScStep <- function(x, ...) {
  cat(sprintf("<ScStep> %s(%s)\n", x$name, format_step_args(x$args)))
  invisible(x)
}

#' @export
print.ScPipeline <- function(x, ...) {
  n <- length(x$steps)
  cat(sprintf("<ScPipeline> %d step%s\n", n, if (n == 1L) "" else "s"))
  if (n == 0L) {
    cat("  (empty)\n")
    return(invisible(x))
  }
  width <- max(nchar(vapply(x$steps, `[[`, character(1), "name")))
  for (i in seq_along(x$steps)) {
    s <- x$steps[[i]]
    cat(sprintf(
      "  %d. %-*s  %s\n",
      i,
      width,
      s$name,
      format_step_args(s$args)
    ))
  }
  invisible(x)
}

# Internal: compact one-line representation of a step's args. Atomic scalars
# rendered with deparse; everything else collapsed to <type>.
format_step_args <- function(args) {
  if (length(args) == 0L) {
    return("")
  }
  parts <- vapply(
    names(args),
    function(nm) {
      v <- args[[nm]]
      if (is.null(v)) {
        val <- "NULL"
      } else if (is.atomic(v) && length(v) == 1L) {
        val <- deparse(v, width.cutoff = 30L)[1]
      } else if (is.atomic(v) && length(v) <= 4L) {
        val <- deparse(v, width.cutoff = 40L)[1]
      } else {
        val <- sprintf("<%s>", class(v)[1])
      }
      paste0(nm, " = ", val)
    },
    character(1)
  )
  paste(parts, collapse = ", ")
}

## apply -----------------------------------------------------------------------

#' Apply a pipeline to a single cell object
#'
#' @description
#' Runs each step in order. The first step receives `object`; subsequent steps
#' receive the result of the previous step. Errors propagate; nothing is
#' caught.
#'
#' @param pipeline `ScPipeline`.
#' @param object `SingleCells` or `SingleCellsSubset`. Dispatch happens inside
#' each step's underlying generic, so the same pipeline works on either.
#'
#' @return The object after all steps have run.
#'
#' @export
apply_pipeline <- function(pipeline, object) {
  checkmate::assertClass(pipeline, "ScPipeline")
  if (length(pipeline$steps) == 0L) {
    return(object)
  }
  Reduce(
    function(obj, step) do.call(step$fn, c(list(obj), step$args)),
    pipeline$steps,
    init = object
  )
}

#' Apply a pipeline independently to each group of a `SingleCells` object
#'
#' @description
#' Splits `object` by `group_col`, applies `pipeline` to each subset, and
#' returns a named list of processed `SingleCellsSubset`s. Useful for
#' per-sample / per-cell-type re-analysis where the same chain (HVG, PCA,
#' neighbours, clusters, ...) is run on each group, e.g. sample-pure metacell
#' generation followed by an external merge.
#'
#' @param pipeline `ScPipeline`.
#' @param object `SingleCells`.
#' @param group_col String. Column in obs used to split.
#' @param groups Optional character vector. Restrict to these group values; if
#' `NULL`, all unique values of `group_col` are used.
#'
#' @return Named list of `SingleCellsSubset` objects, names being the group
#' values.
#'
#' @export
apply_pipeline_per_group <- function(
  pipeline,
  object,
  group_col,
  groups = NULL
) {
  checkmate::assertClass(pipeline, "ScPipeline")
  checkmate::assertClass(object, "bixverse::SingleCells")
  checkmate::qassert(group_col, "S1")
  checkmate::qassert(groups, c("S+", "0"))

  obs <- get_sc_obs(object, filtered = TRUE)
  checkmate::assertNames(colnames(obs), must.include = group_col)

  if (is.null(groups)) {
    groups <- unique(as.character(obs[[group_col]]))
  }

  out <- lapply(groups, function(g) {
    sub <- SingleCellsSubset(object, grouping_column = group_col, group = g)
    apply_pipeline(pipeline, sub)
  })
  names(out) <- groups
  out
}

## steps -----------------------------------------------------------------------

### hvg ------------------------------------------------------------------------

#' Pipeline step: identify highly variable genes
#'
#' @description
#' Wraps [find_hvg_sc()] as an `ScStep`.
#'
#' @inheritParams find_hvg_sc
#'
#' @return An `ScStep`.
#'
#' @export
step_hvg_sc <- function(
  hvg_no = 2000L,
  hvg_params = params_sc_hvg(),
  streaming = NULL,
  .verbose = TRUE
) {
  new_sc_step(
    "hvg",
    find_hvg_sc,
    list(
      hvg_no = hvg_no,
      hvg_params = hvg_params,
      streaming = streaming,
      .verbose = .verbose
    )
  )
}

### pca ------------------------------------------------------------------------

#' Pipeline step: PCA
#'
#' @description
#' Wraps [calculate_pca_sc()] as an `ScStep`.
#'
#' @inheritParams calculate_pca_sc
#'
#' @return An `ScStep`.
#'
#' @export
step_pca_sc <- function(
  no_pcs = 30L,
  pca_params = params_sc_pca(),
  sparse_svd = FALSE,
  hvg = NULL,
  seed = 42L,
  .verbose = TRUE
) {
  new_sc_step(
    "pca",
    calculate_pca_sc,
    list(
      no_pcs = no_pcs,
      pca_params = pca_params,
      sparse_svd = sparse_svd,
      hvg = hvg,
      seed = seed,
      .verbose = .verbose
    )
  )
}

### neighbours -----------------------------------------------------------------

#' Pipeline step: nearest neighbours
#'
#' @description
#' Wraps [find_neighbours_sc()] as an `ScStep`.
#'
#' @inheritParams find_neighbours_sc
#'
#' @return An `ScStep`.
#'
#' @export
step_neighbours_sc <- function(
  embd_to_use = "pca",
  no_embd_to_use = NULL,
  modality = c("rna", "adt"),
  neighbours_params = params_sc_neighbours(),
  seed = 42L,
  .verbose = TRUE
) {
  new_sc_step(
    "neighbours",
    find_neighbours_sc,
    list(
      embd_to_use = embd_to_use,
      no_embd_to_use = no_embd_to_use,
      modality = modality,
      neighbours_params = neighbours_params,
      seed = seed,
      .verbose = .verbose
    )
  )
}

### clusters -------------------------------------------------------------------

#' Pipeline step: graph-based clustering
#'
#' @description
#' Wraps [find_clusters_sc()] as an `ScStep`.
#'
#' @inheritParams find_clusters_sc
#'
#' @return An `ScStep`.
#'
#' @export
step_clusters_sc <- function(
  cluster_algorithm = c("leiden", "louvain"),
  res = 1.0,
  name = "leiden_clustering",
  modality = c("rna", "adt", "wnn"),
  seed = 42L
) {
  new_sc_step(
    "clusters",
    find_clusters_sc,
    list(
      cluster_algorithm = cluster_algorithm,
      res = res,
      name = name,
      modality = modality,
      seed = seed
    )
  )
}

### batch correction ----------------------------------------------------------

#' Pipeline step: Harmony batch correction
#'
#' @description
#' Wraps [harmony_sc()] as an `ScStep`.
#'
#' @inheritParams harmony_sc
#'
#' @return An `ScStep`.
#'
#' @export
step_harmony_sc <- function(
  batch_column,
  additional_batch_columns = NULL,
  modality = c("rna", "adt"),
  harmony_params = params_sc_harmony(),
  seed = 42L,
  .verbose = TRUE
) {
  new_sc_step(
    "harmony",
    harmony_sc,
    list(
      batch_column = batch_column,
      additional_batch_columns = additional_batch_columns,
      modality = modality,
      harmony_params = harmony_params,
      seed = seed,
      .verbose = .verbose
    )
  )
}

#' Pipeline step: Harmony v2 batch correction
#'
#' @description
#' Wraps [harmony_v2_sc()] as an `ScStep`.
#'
#' @inheritParams harmony_v2_sc
#'
#' @return An `ScStep`.
#'
#' @export
step_harmony_v2_sc <- function(
  batch_column,
  additional_batch_columns = NULL,
  modality = c("rna", "adt"),
  harmony_params = params_sc_harmony_v2(),
  seed = 42L,
  .verbose = TRUE
) {
  new_sc_step(
    "harmony_v2",
    harmony_v2_sc,
    list(
      batch_column = batch_column,
      additional_batch_columns = additional_batch_columns,
      modality = modality,
      harmony_params = harmony_params,
      seed = seed,
      .verbose = .verbose
    )
  )
}

#' Pipeline step: BBKNN batch correction
#'
#' @description
#' Wraps [bbknn_sc()] as an `ScStep`.
#'
#' @inheritParams bbknn_sc
#'
#' @return An `ScStep`.
#'
#' @export
step_bbknn_sc <- function(
  batch_column,
  no_neighbours_to_keep = 5L,
  embd_to_use = "pca",
  no_embd_to_use = NULL,
  bbknn_params = params_sc_bbknn(),
  seed = 42L,
  .verbose = TRUE
) {
  new_sc_step(
    "bbknn",
    bbknn_sc,
    list(
      batch_column = batch_column,
      no_neighbours_to_keep = no_neighbours_to_keep,
      embd_to_use = embd_to_use,
      no_embd_to_use = no_embd_to_use,
      bbknn_params = bbknn_params,
      seed = seed,
      .verbose = .verbose
    )
  )
}

#' Pipeline step: fastMNN batch correction
#'
#' @description
#' Wraps [fast_mnn_sc()] as an `ScStep`.
#'
#' @inheritParams fast_mnn_sc
#'
#' @return An `ScStep`.
#'
#' @export
step_fast_mnn_sc <- function(
  batch_column,
  batch_hvg_genes,
  fastmnn_params = params_sc_fastmnn(),
  use_precomputed_pca = FALSE,
  seed = 42L,
  .verbose = TRUE
) {
  new_sc_step(
    "fast_mnn",
    fast_mnn_sc,
    list(
      batch_column = batch_column,
      batch_hvg_genes = batch_hvg_genes,
      fastmnn_params = fastmnn_params,
      use_precomputed_pca = use_precomputed_pca,
      seed = seed,
      .verbose = .verbose
    )
  )
}
