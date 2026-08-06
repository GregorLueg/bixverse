# single cell pipeline ---------------------------------------------------------

# Classes a pipeline step can consume or produce.
SC_STEP_CLASSES <- c("SingleCells", "SingleCellsSubset", "MetaCells")

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
#' @param accepts Character vector. Class names the underlying generic has
#' methods for. Checked by [validate_pipeline()] before anything runs.
#' @param returns String. Class name the step returns, or `"input"` if it hands
#' back whatever it was given.
#'
#' @return An `ScStep` object.
#'
#' @keywords internal
new_sc_step <- function(
  name,
  fn,
  args,
  accepts = c("SingleCells", "SingleCellsSubset"),
  returns = "input"
) {
  checkmate::qassert(name, "S1")
  checkmate::assertFunction(fn)
  checkmate::assertList(args, names = "unique")
  checkmate::assertSubset(accepts, SC_STEP_CLASSES)
  checkmate::assertChoice(returns, c("input", SC_STEP_CLASSES))

  structure(
    list(
      name = name,
      fn = fn,
      args = args,
      accepts = accepts,
      returns = returns
    ),
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

# Internal: short class name of a single cell object, e.g. "bixverse::MetaCells"
# becomes "MetaCells".
.sc_object_class <- function(object) {
  cls <- sub("^bixverse::", "", class(object)[1])
  if (!cls %in% SC_STEP_CLASSES) {
    stop(sprintf(
      "Pipelines run on %s, not on a `%s`.",
      paste(SC_STEP_CLASSES, collapse = ", "),
      cls
    ))
  }
  cls
}

## apply -----------------------------------------------------------------------

#' Check that a pipeline can run on a given class
#'
#' @description
#' Walks the pipeline and tracks which class each step would receive. Steps
#' declare the classes their generic has methods for, and most steps hand back
#' what they were given; [step_metacells_sc()] does not, it turns a
#' `SingleCells`/`SingleCellsSubset` into a `MetaCells`. Without this you would
#' only find out about a mismatch when S7 fails to dispatch, i.e. after the
#' expensive steps already ran.
#'
#' @param pipeline `ScPipeline`.
#' @param class String. Class of the object the pipeline would start on. One of
#' `c("SingleCells", "SingleCellsSubset", "MetaCells")`.
#'
#' @return Invisibly, the class the pipeline would return.
#'
#' @export
validate_pipeline <- function(pipeline, class) {
  checkmate::assertClass(pipeline, "ScPipeline")
  checkmate::assertChoice(class, SC_STEP_CLASSES)

  current <- class
  for (i in seq_along(pipeline$steps)) {
    step <- pipeline$steps[[i]]
    if (!current %in% step$accepts) {
      stop(sprintf(
        paste(
          "Step %d ('%s') cannot run on a `%s`. It has methods for: %s.",
          "\nThe pipeline started on a `%s`."
        ),
        i,
        step$name,
        current,
        paste(step$accepts, collapse = ", "),
        class
      ))
    }
    if (step$returns != "input") {
      current <- step$returns
    }
  }

  invisible(current)
}

#' Apply a pipeline to a single cell object
#'
#' @description
#' Runs each step in order. The first step receives `object`; subsequent steps
#' receive the result of the previous step. The chain is validated against the
#' class of `object` before anything runs, see [validate_pipeline()].
#' Errors propagate; nothing is caught.
#'
#' @param pipeline `ScPipeline`.
#' @param object `SingleCells`, `SingleCellsSubset` or `MetaCells`. Dispatch
#' happens inside each step's underlying generic, so the same pipeline works on
#' any class its steps have methods for.
#'
#' @return The object after all steps have run.
#'
#' @export
apply_pipeline <- function(pipeline, object) {
  checkmate::assertClass(pipeline, "ScPipeline")
  if (length(pipeline$steps) == 0L) {
    return(object)
  }
  validate_pipeline(pipeline, .sc_object_class(object))
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
#' @param progress Boolean. Shall big progress messages be printed to the
#' console. Defaults to `FALSE`.
#'
#' @return Named list of processed objects, names being the group values.
#' Usually `SingleCellsSubset`, or `MetaCells` if the pipeline ends on
#' [step_metacells_sc()], in which case [merge_meta_cells()] puts them back
#' together.
#'
#' @export
apply_pipeline_per_group <- function(
  pipeline,
  object,
  group_col,
  groups = NULL,
  progress = FALSE
) {
  checkmate::assertClass(pipeline, "ScPipeline")
  checkmate::assertClass(object, "bixverse::SingleCells")
  checkmate::qassert(group_col, "S1")
  checkmate::qassert(groups, c("S+", "0"))
  checkmate::qassert(progress, "B1")

  # fail before any subset is built rather than on the first group
  validate_pipeline(pipeline, "SingleCellsSubset")

  group_vec <- unlist(get_sc_obs[[c("group_col")]], use.names = FALSE)
  checkmate::qassert(group_vec, "S+")

  if (is.null(groups)) {
    groups <- unique(group_vec)
  }

  out <- vector(mode = "list", length = length(groups))

  for (i in seq_along(groups)) {
    g <- groups[[i]]
    if (progress) {
      cat(sprintf("\n=== Applying pipeline to %s\n ===", g))
    }

    sub <- SingleCellsSubset(object, grouping_column = group_col, group = g)
    res <- apply_pipeline(pipeline, sub)

    out[[i]] <- res
  }

  names(out) <- groups

  out
}

## meta cells ------------------------------------------------------------------

#' Generate source-pure meta cells and merge them
#'
#' @description
#' Splits `object` by `group_col`, optionally runs `pipeline` on each subset,
#' generates meta cells per group and merges the results into a single
#' [bixverse::MetaCells()] object. This gives you meta cells that never mix
#' cells from two patients/samples, while still being one object you can run
#' SCENIC, AUCell or NMF over.
#'
#' The meta cell generators need an embedding, so `pipeline` will normally be
#' `step_hvg_sc() %>>% step_pca_sc() %>>% step_neighbours_sc()` unless every
#' subset already carries one.
#'
#' @param object `SingleCells`.
#' @param group_col String. Column in obs used to split.
#' @param method String. One of `c("bootstrapped", "seacells", "supercells")`.
#' Picks the meta cell generator.
#' @param mc_params Named list. Arguments passed on to the generator, e.g.
#' `list(sc_meta_cell_params = params_sc_bt_metacells(), target_size = 1e5)`.
#' @param pipeline Optional `ScPipeline` applied to each subset before the meta
#' cells are generated.
#' @param groups Optional character vector. Restrict to these group values; if
#' `NULL`, all unique values of `group_col` are used.
#' @param feature_space String. One of `c("intersect", "union")`. Passed to
#' [merge_meta_cells()]. Irrelevant here as all groups share the gene space of
#' the parent object.
#' @param .verbose Boolean. Controls verbosity of the function.
#'
#' @returns A merged [bixverse::MetaCells()] object with a `source_id` column
#' in its observation table.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' prep <- step_hvg_sc() %>>% step_pca_sc(no_pcs = 20L) %>>% step_neighbours_sc()
#' mc <- meta_cells_per_group(
#'   object = sc_obj,
#'   group_col = "patient_id",
#'   method = "bootstrapped",
#'   pipeline = prep
#' )
#' }
meta_cells_per_group <- function(
  object,
  group_col,
  method = c("bootstrapped", "seacells", "supercells"),
  mc_params = list(),
  pipeline = NULL,
  groups = NULL,
  feature_space = c("intersect", "union"),
  .verbose = TRUE
) {
  method <- match.arg(method)
  feature_space <- match.arg(feature_space)

  # checks
  checkmate::assertClass(object, "bixverse::SingleCells")
  checkmate::qassert(group_col, "S1")
  checkmate::assertChoice(method, c("bootstrapped", "seacells", "supercells"))
  checkmate::assertList(mc_params, names = "unique")
  if (!is.null(pipeline)) {
    checkmate::assertClass(pipeline, "ScPipeline")
  }
  checkmate::qassert(groups, c("S+", "0"))
  checkmate::qassert(.verbose, "B1")

  obs <- get_sc_obs(object, filtered = TRUE)
  checkmate::assertNames(colnames(obs), must.include = group_col)

  if (is.null(groups)) {
    groups <- unique(as.character(obs[[group_col]]))
  }

  if (!is.null(pipeline)) {
    validate_pipeline(pipeline, "SingleCellsSubset")
  }

  mc_fn <- .sc_metacell_generator(method)

  meta_cells <- lapply(groups, function(g) {
    if (.verbose) {
      message(sprintf("Generating meta cells for group '%s'.", g))
    }
    sub <- SingleCellsSubset(object, grouping_column = group_col, group = g)
    if (!is.null(pipeline)) {
      sub <- apply_pipeline(pipeline, sub)
    }
    do.call(mc_fn, c(list(sub), mc_params))
  })

  merge_meta_cells(
    inputs = meta_cells,
    source_ids = groups,
    feature_space = feature_space,
    .verbose = .verbose
  )
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
    ),
    accepts = SC_STEP_CLASSES
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
    ),
    accepts = SC_STEP_CLASSES
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
    ),
    accepts = SC_STEP_CLASSES
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
    ),
    accepts = SC_STEP_CLASSES
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

### meta cells -----------------------------------------------------------------

#' Pipeline step: generate meta cells
#'
#' @description
#' Wraps the meta cell generators as an `ScStep`. Unlike the other steps this
#' one changes the class of the object: it takes a `SingleCells` or
#' `SingleCellsSubset` and returns a [bixverse::MetaCells()]. Steps that follow
#' it need `MetaCells` methods, which [validate_pipeline()] checks up front.
#'
#' Combined with [apply_pipeline_per_group()] this gives you per-group
#' pre-processing (HVG, PCA, batch correction within a patient, kNN) followed by
#' source-pure meta cells, which you then hand to [merge_meta_cells()].
#'
#' @param method String. One of `c("bootstrapped", "seacells", "supercells")`.
#' @param ... Arguments passed on to the generator, e.g.
#' `sc_meta_cell_params`, `target_size` or `.verbose`.
#'
#' @return An `ScStep`.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' pipeline <- step_hvg_sc() %>>%
#'   step_pca_sc(no_pcs = 20L) %>>%
#'   step_harmony_sc(batch_column = "plate") %>>%
#'   step_neighbours_sc(embd_to_use = "harmony") %>>%
#'   step_metacells_sc("bootstrapped")
#'
#' per_patient <- apply_pipeline_per_group(pipeline, sc_obj, "patient_id")
#' mc <- merge_meta_cells(per_patient)
#' }
step_metacells_sc <- function(
  method = c("bootstrapped", "seacells", "supercells"),
  ...
) {
  method <- match.arg(method)
  checkmate::assertChoice(method, c("bootstrapped", "seacells", "supercells"))

  new_sc_step(
    sprintf("metacells (%s)", method),
    .sc_metacell_generator(method),
    list(...),
    accepts = c("SingleCells", "SingleCellsSubset"),
    returns = "MetaCells"
  )
}

# Internal: the generic behind a meta cell method.
.sc_metacell_generator <- function(method) {
  switch(
    method,
    "bootstrapped" = generate_bt_meta_cells_sc,
    "seacells" = generate_seacells_sc,
    "supercells" = generate_supercells_sc
  )
}
