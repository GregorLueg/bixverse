# single cell cache state ------------------------------------------------------

# every artefact that lands in an `ScCache` (PCA factors, embeddings, the kNN
# object, the sNN graph) carries a provenance stamp as an attribute on the
# payload itself. the stamp records the cell set the artefact was computed on
# and the ids of the artefacts it was derived from, which is what lets us tell
# a stale PCA from a fresh one after the cell filter moved.
#
# the stamp rides on the payload rather than in a parallel registry so that
# deleting a payload deletes its stamp. there is no way for the two to drift.

## globals ---------------------------------------------------------------------

# attribute name under which every stamp is stored
SC_STAMP_ATTR <- "bixverse_stamp"

# artefact kinds a cache can hold
SC_ARTEFACTS <- c("pca", "embedding", "knn", "snn")

# embedding names that would collide with an artefact label
SC_RESERVED_EMBEDDINGS <- c("pca", "knn", "snn")

# monotonic counter so two stamps minted in the same clock tick still differ
sc_stamp_env <- new.env(parent = emptyenv())
sc_stamp_env$counter <- 0L

## stamp primitives ------------------------------------------------------------

#' Generate a unique artefact id
#'
#' @description
#' Builds a 16 character id from the wall clock, the process id and a package
#' local counter. Deliberately avoids `runif()`/`sample()`: consuming the RNG
#' stream would shift the results of every seeded test in the package.
#'
#' @returns String. A 16 character identifier.
#'
#' @keywords internal
.new_sc_id <- function() {
  sc_stamp_env$counter <- sc_stamp_env$counter + 1L

  substr(
    rlang::hash(list(Sys.time(), Sys.getpid(), sc_stamp_env$counter)),
    1,
    16
  )
}

#' Attach a provenance stamp to an artefact
#'
#' @param obj The artefact payload (matrix, kNN object, igraph).
#' @param stamp List. The stamp as built by [.new_sc_stamp()].
#'
#' @returns The payload with the stamp attached as an attribute.
#'
#' @keywords internal
.set_stamp <- function(obj, stamp) {
  # checks
  checkmate::assertList(stamp)

  attr(obj, SC_STAMP_ATTR) <- stamp

  return(obj)
}

#' Read the provenance stamp of an artefact
#'
#' @param obj The artefact payload.
#'
#' @returns The stamp list, or `NULL` when the artefact carries none (which is
#' the case for anything computed before stamping existed).
#'
#' @keywords internal
.read_stamp <- function(obj) {
  attr(obj, SC_STAMP_ATTR, exact = TRUE)
}

#' Strip the provenance stamp from an artefact
#'
#' @description
#' Used by the getters on the way out. The stamp is an internal consistency
#' device and has no business leaking into user code, where it would show up in
#' `str()` and make `identical()` on two numerically equal matrices `FALSE`.
#'
#' @param obj The artefact payload.
#'
#' @returns The payload without the stamp attribute.
#'
#' @keywords internal
.drop_stamp <- function(obj) {
  if (!is.null(obj)) {
    attr(obj, SC_STAMP_ATTR) <- NULL
  }

  return(obj)
}

#' Build a provenance stamp
#'
#' @param cells String. Hash of the cell state at write time.
#' @param from Character vector. Ids of the artefacts this one derives from.
#' Empty for a root artefact such as the PCA.
#'
#' @returns A list with `cells`, `id` and `from`.
#'
#' @keywords internal
.new_sc_stamp <- function(cells, from = character()) {
  # checks
  checkmate::qassert(cells, "S1")
  checkmate::qassert(from, c("S*", "0"))

  list(
    cells = cells,
    id = .new_sc_id(),
    from = as.character(from)
  )
}

## cell state hash -------------------------------------------------------------

#' Hash of the cell state an artefact must agree with
#'
#' @description
#' Order sensitive on purpose. [set_cells_to_keep()] does not sort, and the PCA
#' rows come back from Rust in exactly `cell_indices` order, so row order is
#' what aligns a cached matrix to a cell list.
#'
#' @param x The object owning the cache.
#'
#' @returns String. The hash of the current cell state.
#'
#' @keywords internal
.sc_state_hash <- S7::new_generic(
  name = ".sc_state_hash",
  dispatch_args = "x",
  fun = function(x) {
    S7::S7_dispatch()
  }
)

#' @noRd
S7::method(.sc_state_hash, SingleCells) <- function(x) {
  rlang::hash(get_cells_to_keep(x))
}

#' @noRd
S7::method(.sc_state_hash, SingleCellsSubset) <- function(x) {
  rlang::hash(get_cells_to_keep(x))
}

#' @noRd
S7::method(.sc_state_hash, MetaCells) <- function(x) {
  # meta cells have no `cells_to_keep`; their cell set is fixed at construction
  # and identified by the meta cell ids, which `merge_meta_cells()` rewrites
  # with source prefixes. hash the plain character vector, never the
  # data.table, whose internal self reference pointer is not stable.
  rlang::hash(as.character(S7::prop(x, "obs_table")[["meta_cell_id"]]))
}

## cache routing ---------------------------------------------------------------

#' The modalities an object can hold a cache for
#'
#' @param x The object owning the cache(s).
#'
#' @returns Character vector of modality names.
#'
#' @keywords internal
.sc_modalities <- function(x) {
  # multi modal must be tested first, it has `parent = SingleCells`
  if (S7::S7_inherits(x, SingleCellsMultiModal)) {
    c("rna", "adt", "atac", "wnn")
  } else {
    "rna"
  }
}

#' Read the cache for a given modality
#'
#' @param x The object owning the cache.
#' @param modality String. One of `c("rna", "adt", "atac", "wnn")`.
#'
#' @returns The cache list, or `NULL` when the object has no such cache.
#'
#' @keywords internal
.sc_cache_get <- function(x, modality = "rna") {
  if (S7::S7_inherits(x, SingleCellsMultiModal)) {
    if (identical(modality, "wnn")) {
      # tolerant on purpose; `.integration_slot()` errors when wnn is absent
      return(S7::prop(x, "other_data")[["wnn"]])
    }
    return(S7::prop(x, .cache_slot_from_modality(modality)))
  }

  if (!identical(modality, "rna")) {
    return(NULL)
  }

  S7::prop(x, "sc_cache")
}

#' Write the cache for a given modality back onto the object
#'
#' @param x The object owning the cache.
#' @param cache The cache list to write.
#' @param modality String. One of `c("rna", "adt", "atac", "wnn")`.
#'
#' @returns The object with the cache written back.
#'
#' @keywords internal
.sc_cache_put <- function(x, cache, modality = "rna") {
  if (S7::S7_inherits(x, SingleCellsMultiModal)) {
    if (identical(modality, "wnn")) {
      other_data <- S7::prop(x, "other_data")
      other_data[["wnn"]] <- cache
      S7::prop(x, "other_data") <- other_data
      return(x)
    }
    S7::prop(x, .cache_slot_from_modality(modality)) <- cache
    return(x)
  }

  S7::prop(x, "sc_cache") <- cache

  return(x)
}

#' Read an artefact payload straight out of a cache
#'
#' @description
#' Raw slot access with no freshness checking and no stamp stripping. This is
#' what the provenance machinery reads; user facing code goes through the
#' getters.
#'
#' The `wnn` pseudo modality stores its graph under `snn` and its embeddings
#' under `embeddings`, where a real `ScCache` uses `snn_graph` and
#' `other_embeddings`, hence the flag.
#'
#' @param cache The cache list.
#' @param artefact String. One of `c("pca", "embedding", "knn", "snn")`.
#' @param name String. Embedding name, only used when `artefact = "embedding"`.
#' @param wnn Boolean. Whether `cache` is the `wnn` slot rather than an
#' `ScCache`.
#'
#' @returns The payload, or `NULL`.
#'
#' @keywords internal
.sc_payload <- function(cache, artefact, name = NULL, wnn = FALSE) {
  if (is.null(cache)) {
    return(NULL)
  }

  switch(
    artefact,
    pca = if (wnn) NULL else cache[["pca_factors"]],
    knn = cache[["knn"]],
    snn = if (wnn) cache[["snn"]] else cache[["snn_graph"]],
    embedding = cache[[
      if (wnn) "embeddings" else "other_embeddings"
    ]][[name]],
    NULL
  )
}

#' Write an artefact payload back into a cache
#'
#' @inheritParams .sc_payload
#' @param value The payload to write.
#'
#' @returns The cache with the payload written back.
#'
#' @keywords internal
.sc_payload_set <- function(cache, value, artefact, name = NULL, wnn = FALSE) {
  slot <- switch(
    artefact,
    pca = "pca_factors",
    knn = "knn",
    snn = if (wnn) "snn" else "snn_graph",
    embedding = if (wnn) "embeddings" else "other_embeddings",
    stop(sprintf("Unknown artefact '%s'.", artefact))
  )

  if (identical(artefact, "embedding")) {
    cache[[slot]][[name]] <- value
  } else {
    cache[[slot]] <- value
  }

  return(cache)
}

#' Is an artefact present in a cache?
#'
#' @description
#' Pure presence probe with no freshness checking. For the call sites that ask
#' whether something is there before overwriting it, where erroring on a stale
#' artefact that is about to be replaced would be actively wrong.
#'
#' @param x The object owning the cache.
#' @inheritParams .sc_payload
#' @param modality String. The modality to look in.
#'
#' @returns Boolean.
#'
#' @keywords internal
.sc_has_artefact <- function(x, artefact, name = NULL, modality = "rna") {
  # checks
  checkmate::assertChoice(artefact, SC_ARTEFACTS)

  !is.null(.sc_payload(
    cache = .sc_cache_get(x, modality),
    artefact = artefact,
    name = name,
    wnn = identical(modality, "wnn")
  ))
}

## artefact enumeration --------------------------------------------------------

#' Label an artefact for reporting and provenance lookups
#'
#' @inheritParams .sc_payload
#' @param modality String. The modality the artefact lives in.
#'
#' @returns String of the form `"modality:artefact"`, where embeddings use
#' their name (`"rna:umap"`).
#'
#' @keywords internal
.sc_label <- function(modality, artefact, name = NULL) {
  sprintf(
    "%s:%s",
    modality,
    if (identical(artefact, "embedding")) name else artefact
  )
}

#' Enumerate every artefact present on an object
#'
#' @param x The object owning the cache(s).
#'
#' @returns A list of lists, each with `modality`, `artefact`, `name`, `label`
#' and `stamp` (which is `NULL` for an unstamped, i.e. legacy, artefact).
#'
#' @keywords internal
.sc_artefacts <- function(x) {
  out <- list()

  for (modality in .sc_modalities(x)) {
    cache <- .sc_cache_get(x, modality)
    if (is.null(cache)) {
      next
    }
    wnn <- identical(modality, "wnn")

    embd_names <- names(cache[[if (wnn) "embeddings" else "other_embeddings"]])

    entries <- c(
      list(list(artefact = "pca", name = NA_character_)),
      purrr::map(embd_names, \(nm) list(artefact = "embedding", name = nm)),
      list(list(artefact = "knn", name = NA_character_)),
      list(list(artefact = "snn", name = NA_character_))
    )

    for (entry in entries) {
      payload <- .sc_payload(cache, entry$artefact, entry$name, wnn)
      if (is.null(payload)) {
        next
      }
      label <- .sc_label(modality, entry$artefact, entry$name)
      out[[label]] <- list(
        modality = modality,
        artefact = entry$artefact,
        name = entry$name,
        label = label,
        stamp = .read_stamp(payload)
      )
    }
  }

  return(out)
}

#' Resolve parent artefact names to their current stamp ids
#'
#' @description
#' Producers pass parent *names* (`"pca"`, `"umap"`) rather than ids, so they
#' never have to reach into the stamps themselves. A name may be modality
#' qualified (`"adt:pca"`); otherwise `modality` is assumed. Names that do not
#' resolve to a stamped artefact are dropped, which is what keeps a chain built
#' on top of a legacy artefact from claiming a provenance it does not have.
#'
#' @param x The object owning the cache(s).
#' @param from Character vector of parent artefact names.
#' @param modality String. Modality to resolve unqualified names against.
#'
#' @returns Character vector of stamp ids.
#'
#' @keywords internal
.stamp_ids_from <- function(x, from, modality = "rna") {
  if (length(from) == 0) {
    return(character())
  }

  artefacts <- .sc_artefacts(x)

  ids <- purrr::map_chr(from, \(spec) {
    label <- if (grepl(":", spec, fixed = TRUE)) {
      spec
    } else {
      sprintf("%s:%s", modality, spec)
    }
    stamp <- artefacts[[label]][["stamp"]]
    if (is.null(stamp)) NA_character_ else stamp[["id"]]
  })

  return(ids[!is.na(ids)])
}

## stamping --------------------------------------------------------------------

#' Parents of a manifold embedding (UMAP, t-SNE, PHATE)
#'
#' @description
#' These three read their source embedding from `cache_modality` but write the
#' result under `modality`, which differ for `"wnn"`. Both parents are therefore
#' spelled out modality qualified rather than left to the write modality.
#'
#' @param embd_to_use String. Name of the source embedding.
#' @param cache_modality String. Modality the source embedding was read from.
#' @param modality String. Modality the result is written to.
#' @param has_knn Boolean. Whether a cached kNN fed the manifold.
#'
#' @returns Character vector of parent artefact names.
#'
#' @keywords internal
.manifold_from <- function(embd_to_use, cache_modality, modality, has_knn) {
  c(
    sprintf("%s:%s", cache_modality, embd_to_use),
    if (has_knn) sprintf("%s:knn", modality) else NULL
  )
}

#' Pull the parent artefact names out of a setter's dots
#'
#' @description
#' The setter generics all carry `...`, so producers opt into provenance by
#' passing `from = `. Producers that do not get a root stamp and keep working,
#' which is what makes this backwards compatible with downstream packages
#' calling the setters positionally.
#'
#' @param ... The setter's dots.
#'
#' @returns Character vector of parent artefact names, possibly empty.
#'
#' @keywords internal
.stamp_from <- function(...) {
  as.character(list(...)[["from"]] %||% character())
}

#' Stamp an artefact that has just been written into a cache
#'
#' @description
#' Called by the S7 setter forwarders, which are the only layer with access to
#' both the cache and the cell state. A no-op when the artefact is absent.
#'
#' @param x The object owning the cache.
#' @param artefact String. One of `c("pca", "embedding", "knn", "snn")`.
#' @param name String. Embedding name, only used when `artefact = "embedding"`.
#' @param modality String. The modality the artefact was written to.
#' @param from Character vector of parent artefact names.
#'
#' @returns The object with the artefact stamped.
#'
#' @keywords internal
.stamp_artefact <- function(
  x,
  artefact,
  name = NULL,
  modality = "rna",
  from = character()
) {
  # checks
  checkmate::assertChoice(artefact, SC_ARTEFACTS)
  checkmate::qassert(from, c("S*", "0"))

  cache <- .sc_cache_get(x, modality)
  if (is.null(cache)) {
    return(x)
  }

  wnn <- identical(modality, "wnn")
  payload <- .sc_payload(cache, artefact, name, wnn)
  if (is.null(payload)) {
    return(x)
  }

  stamp <- .new_sc_stamp(
    cells = .sc_state_hash(x),
    from = .stamp_ids_from(x, from, modality)
  )

  cache <- .sc_payload_set(
    cache = cache,
    value = .set_stamp(payload, stamp),
    artefact = artefact,
    name = name,
    wnn = wnn
  )

  .sc_cache_put(x, cache, modality)
}

## state evaluation ------------------------------------------------------------

#' Why an artefact is stale, if it is
#'
#' @description
#' Walks the provenance chain. An artefact is stale when the cell set moved
#' under it, when a parent it was built on no longer exists (i.e. was
#' recomputed and minted a new id), or when a parent is itself stale. The last
#' case is why this recurses: without it, re-running the PCA marks the kNN
#' stale but the sNN still resolves its untouched parent and passes.
#'
#' @param label String. The artefact label to evaluate.
#' @param artefacts List as returned by [.sc_artefacts()].
#' @param cells_hash String. The object's current cell state hash.
#' @param seen Character vector. Labels already on the recursion stack.
#'
#' @returns `NA_character_` when fresh (or unstamped), otherwise the reason.
#'
#' @keywords internal
.sc_stale_reason <- function(label, artefacts, cells_hash, seen = character()) {
  stamp <- artefacts[[label]][["stamp"]]

  # unstamped artefacts predate this machinery; unknown provenance passes
  if (is.null(stamp) || label %in% seen) {
    return(NA_character_)
  }

  if (!identical(stamp[["cells"]], cells_hash)) {
    return("the cell set changed since it was computed")
  }

  if (length(stamp[["from"]]) == 0) {
    return(NA_character_)
  }

  known <- purrr::map_chr(artefacts, \(a) {
    a[["stamp"]][["id"]] %||% NA_character_
  })

  for (parent_id in stamp[["from"]]) {
    parent <- names(known)[which(known == parent_id)]

    if (length(parent) == 0) {
      return("the artefact it was derived from was re-computed or removed")
    }

    upstream <- .sc_stale_reason(
      label = parent[1],
      artefacts = artefacts,
      cells_hash = cells_hash,
      seen = c(seen, label)
    )

    if (!is.na(upstream)) {
      return(sprintf("its upstream `%s` is stale", parent[1]))
    }
  }

  return(NA_character_)
}

#' Full state report for an object's caches
#'
#' @param x The object owning the cache(s).
#'
#' @returns A `data.table` with one row per present artefact and the columns
#' \itemize{
#'   \item modality - The modality the artefact lives in.
#'   \item artefact - One of `pca`, `embedding`, `knn`, `snn`.
#'   \item name - The embedding name, `NA` for the others.
#'   \item stamped - Whether the artefact carries a provenance stamp.
#'   \item stale - Whether it disagrees with the current state.
#'   \item reason - Why it is stale, `NA` otherwise.
#'   \item id - The artefact's stamp id.
#'   \item from - List column of the parent ids.
#' }
#'
#' @keywords internal
.sc_state_report <- function(x) {
  artefacts <- .sc_artefacts(x)

  if (length(artefacts) == 0) {
    return(data.table::data.table(
      modality = character(),
      artefact = character(),
      name = character(),
      stamped = logical(),
      stale = logical(),
      reason = character(),
      id = character(),
      from = list()
    ))
  }

  # hashed once for the whole table, never per row
  cells_hash <- .sc_state_hash(x)

  rows <- purrr::map(artefacts, \(entry) {
    reason <- .sc_stale_reason(entry$label, artefacts, cells_hash)
    data.table::data.table(
      modality = entry$modality,
      artefact = entry$artefact,
      name = as.character(entry$name),
      stamped = !is.null(entry$stamp),
      stale = !is.na(reason),
      reason = reason,
      id = entry$stamp[["id"]] %||% NA_character_,
      from = list(entry$stamp[["from"]] %||% character())
    )
  })

  data.table::rbindlist(rows)
}

#' Labels of every stale artefact on an object
#'
#' @param x The object owning the cache(s).
#'
#' @returns Character vector of `"modality:artefact"` labels, possibly empty.
#' Never signals: this backs a warning and a `print` method.
#'
#' @keywords internal
.sc_stale_artefacts <- function(x) {
  report <- tryCatch(.sc_state_report(x), error = function(e) NULL)

  if (is.null(report) || nrow(report) == 0) {
    return(character())
  }

  # base subsetting rather than data.table NSE, which would need a
  # `globalVariables()` declaration the package does not otherwise use
  stale <- report[which(report$stale), ]

  if (nrow(stale) == 0) {
    return(character())
  }

  purrr::pmap_chr(
    list(stale$modality, stale$artefact, stale$name),
    \(modality, artefact, name) .sc_label(modality, artefact, name)
  )
}

#' Stale artefacts rendered for a `print` method
#'
#' @param x The object owning the cache(s).
#'
#' @returns String. Comma separated labels, or `"none"`. Never signals: a
#' `print` method that can fail is worse than one that says "unknown".
#'
#' @keywords internal
.print_stale_str <- function(x) {
  stale <- tryCatch(.sc_stale_artefacts(x), error = function(e) NULL)

  if (is.null(stale)) {
    return("unknown")
  }

  if (length(stale) == 0) {
    return("none")
  }

  paste(stale, collapse = ", ")
}

## severity --------------------------------------------------------------------

#' Resolve the configured strictness
#'
#' @description
#' `bixverse.cache_check` accepts `"error"`, `"warn"` and `"none"`. Unset means
#' the two tier default: the getters warn, [assert_sc_state()] errors.
#' `"error"` promotes the getter warnings to errors, `"warn"` is the explicit
#' spelling of the default getter behaviour, and `"none"` disables both tiers.
#'
#' @returns String. One of `c("default", "error", "warn", "none")`.
#'
#' @keywords internal
.sc_check_mode <- function() {
  mode <- getOption("bixverse.cache_check", "default")

  if (!checkmate::testChoice(mode, c("default", "error", "warn", "none"))) {
    warning(sprintf(
      "Unknown `bixverse.cache_check` value '%s'. Falling back to the default.",
      as.character(mode)[1]
    ))
    return("default")
  }

  return(mode)
}

#' Warn (or error) when a getter hands back a stale artefact
#'
#' @description
#' The soft tier. Every consumer reads through the getters, so this catches
#' everything, but it stays a warning because several call sites use the
#' getters as presence probes on artefacts they are about to overwrite.
#'
#' @param x The object owning the cache.
#' @param artefact String. One of `c("pca", "embedding", "knn", "snn")`.
#' @param name String. Embedding name, only used when `artefact = "embedding"`.
#' @param modality String. The modality read from.
#'
#' @returns Invisibly `TRUE`. Called for the side effect.
#'
#' @keywords internal
.warn_sc_state <- function(x, artefact, name = NULL, modality = "rna") {
  mode <- .sc_check_mode()

  if (identical(mode, "none")) {
    return(invisible(TRUE))
  }

  res <- check_sc_state(
    x = x,
    artefacts = if (identical(artefact, "embedding")) name else artefact,
    modality = modality
  )

  if (isTRUE(res)) {
    return(invisible(TRUE))
  }

  if (identical(mode, "error")) {
    stop(res, call. = FALSE)
  }

  warning(res, call. = FALSE)

  invisible(TRUE)
}

## resetting -------------------------------------------------------------------

#' Ask before wiping a populated cache
#'
#' @description
#' Returns `TRUE` when the reset should go ahead. Skips the prompt entirely
#' when there is nothing to destroy, and refuses to prompt in a non-interactive
#' session, where `readline()` returns `""` immediately and would silently pick
#' a branch nobody chose.
#'
#' @param x The object about to be reset.
#' @param force Boolean. Whether the caller already opted in.
#'
#' @returns Boolean. Whether to proceed.
#'
#' @keywords internal
.confirm_sc_reset <- function(x, force) {
  # checks
  checkmate::qassert(force, "B1")

  sc_map <- S7::prop(x, "sc_map")

  destroyed <- c(
    names(.sc_artefacts(x)),
    if (!is.null(sc_map[["hvg_gene_indices"]])) "the HVG selection"
  )

  # nothing cached: the reset destroys nothing and needs no confirmation
  if (length(destroyed) == 0) {
    return(TRUE)
  }

  if (force) {
    return(TRUE)
  }

  if (!interactive()) {
    stop(sprintf(
      paste(
        "Resetting would drop %s. There is no interactive session to confirm",
        "in. Pass `force = TRUE` if that is what you want."
      ),
      paste(destroyed, collapse = ", ")
    ))
  }

  cat(sprintf(
    paste0(
      "Resetting restores every cell in the binary and drops:\n  %s\n",
      "A PCA or kNN computed on a filtered subset does not describe the full\n",
      "cell set, so it cannot be kept. This cannot be undone.\n"
    ),
    paste(destroyed, collapse = ", ")
  ))

  tolower(trimws(readline("Proceed? (y/N): "))) %in% c("y", "yes")
}

#' Restore every cell and clear the gene selection
#'
#' @description
#' Sets `cells_to_keep_idx` back to `NULL` rather than materialising the full
#' index vector: `get_cells_to_keep()` already falls back to all cells on an
#' empty slot, so this is the genuine pristine state and costs nothing at a
#' million cells. The HVG selection goes too, since [find_hvg_sc()] reads the
#' cell set and a selection made on a subset does not describe the whole.
#'
#' @param x `SingleCells` or `SingleCellsMultiModal` class.
#'
#' @returns The object with its map and DuckDB reset.
#'
#' @keywords internal
.reset_sc_cells <- function(x) {
  n_cells <- as.integer(get_sc_rust_ptr(x)$get_shape()[1])

  if (!identical(n_cells, S7::prop(x, "dims")[1])) {
    stop(sprintf(
      paste(
        "The object reports %i cells but the binary holds %i. Refusing to",
        "reset an object that disagrees with its own count file."
      ),
      S7::prop(x, "dims")[1],
      n_cells
    ))
  }

  sc_map <- S7::prop(x, "sc_map")
  # `[<-` with list(NULL) nulls the value; `[[<-` would drop the element
  sc_map["cells_to_keep_idx"] <- list(NULL)
  sc_map["hvg_gene_indices"] <- list(NULL)
  S7::prop(x, "sc_map") <- sc_map

  get_sc_duckdb(x)$reset_cells_to_keep()

  return(x)
}

## user facing -----------------------------------------------------------------

#' Check that cached artefacts still match the object's state
#'
#' @description
#' Returns `TRUE` when every requested artefact is fresh, otherwise a message
#' describing the first problem. An artefact that is absent, or that predates
#' provenance stamping, passes: unknown provenance is not the same as known
#' bad.
#'
#' Named in snake case rather than the `checkXxx` camel case of the other
#' checkmate extensions in this package, because it validates object state
#' rather than a `params_*()` list and sits alongside [get_sc_cache_status()].
#'
#' @param x `SingleCells`, `SingleCellsSubset`, `MetaCells` or
#' `SingleCellsMultiModal` class.
#' @param artefacts Character vector. Artefacts to check, either an artefact
#' kind (`"pca"`, `"knn"`, `"snn"`) or an embedding name (`"umap"`). May be
#' modality qualified (`"adt:pca"`).
#' @param modality String. Modality to resolve unqualified names against. One
#' of `c("rna", "adt", "atac", "wnn")`.
#'
#' @returns `TRUE` if the check was successful, otherwise a string describing
#' the failure.
#'
#' @export
check_sc_state <- function(x, artefacts, modality = "rna") {
  # checks
  checkmate::qassert(artefacts, "S+")
  checkmate::assertChoice(modality, c("rna", "adt", "atac", "wnn"))

  if (identical(.sc_check_mode(), "none")) {
    return(TRUE)
  }

  present <- .sc_artefacts(x)

  if (length(present) == 0) {
    return(TRUE)
  }

  labels <- purrr::map_chr(artefacts, \(spec) {
    if (grepl(":", spec, fixed = TRUE)) {
      spec
    } else {
      sprintf("%s:%s", modality, spec)
    }
  })

  # absent artefacts are handled by the caller's own presence guards, and
  # unstamped ones predate this machinery. bail before hashing the cell state
  # when there is nothing to compare against; this is what keeps the getters
  # cheap on legacy objects and at a million cells.
  relevant <- present[intersect(labels, names(present))]

  if (length(relevant) == 0) {
    return(TRUE)
  }

  if (all(purrr::map_lgl(relevant, \(a) is.null(a[["stamp"]])))) {
    return(TRUE)
  }

  cells_hash <- .sc_state_hash(x)

  for (label in labels) {
    if (is.null(present[[label]])) {
      next
    }

    reason <- .sc_stale_reason(label, present, cells_hash)

    if (!is.na(reason)) {
      # no trailing full stop: `assert_sc_state()` appends its own
      return(sprintf(
        paste(
          "Cached `%s` is stale: %s. Re-run the step that produces it, or",
          "drop it, before using it. See `get_sc_cache_status()`"
        ),
        label,
        reason
      ))
    }
  }

  return(TRUE)
}

#' Assert that cached artefacts still match the object's state
#'
#' @description
#' The hard tier. Called at the entry of anything that hands cached indices to
#' Rust or writes a derived result back into the object, where a stale artefact
#' means either a crash in the binding layer or silently mis-aligned biology.
#'
#' @inheritParams check_sc_state
#'
#' @returns Invisibly `x` when the assertion holds, otherwise an error.
#'
#' @export
assert_sc_state <- checkmate::makeAssertionFunction(check_sc_state)

#' Status of everything held in a single cell object's caches
#'
#' @description
#' One row per cached artefact, saying whether it still agrees with the current
#' cell set and with the artefacts it was derived from. Artefacts computed
#' before provenance stamping existed report `stamped = FALSE` and are never
#' flagged stale, because nothing is known about them either way.
#'
#' @param object `SingleCells`, `SingleCellsSubset`, `MetaCells` or
#' `SingleCellsMultiModal` class.
#'
#' @returns A `data.table` with the columns
#' \itemize{
#'   \item modality - The modality the artefact lives in.
#'   \item artefact - One of `pca`, `embedding`, `knn`, `snn`.
#'   \item name - The embedding name, `NA` for the others.
#'   \item stamped - Whether the artefact carries a provenance stamp.
#'   \item stale - Whether it disagrees with the current state.
#'   \item reason - Why it is stale, `NA` otherwise.
#'   \item id - The artefact's stamp id.
#'   \item from - List column of the parent stamp ids.
#' }
#'
#' @export
get_sc_cache_status <- function(object) {
  # checks
  checkmate::assertTRUE(
    S7::S7_inherits(object, SingleCells) ||
      S7::S7_inherits(object, SingleCellsSubset) ||
      S7::S7_inherits(object, MetaCells)
  )

  .sc_state_report(object)
}
