# methods for meta cells and their processing/metrics --------------------------

## metrics ---------------------------------------------------------------------

### meta cell purity -----------------------------------------------------------

#' Calculate meta cell purity
#'
#' @description
#' A potential metric to see how well the meta cells are aggregated is their
#' cell type purity. This helper function helps to plot the meta-cell purity
#' based on annotated cell types. These can be also just unsupervised
#' memberships to graph-based clustering, etc.
#'
#' @details
#' Ties for the top (or second) label resolve to whichever label sorts first,
#' as the labels are factorised internally. The entropy is the Shannon entropy
#' of the label distribution within a meta cell, divided by
#' `log(<number of distinct labels in original_cell_type>)`, so it sits in
#' `[0, 1]` and stays comparable between meta cells. It is `0` if there is only
#' a single label in the data.
#'
#' @param object `MetaCells` class.
#' @param original_cell_type Character vector. The original cell type
#' annotations of the object the meta cells came from. Either in the row order
#' of its full (unfiltered) obs table, i.e. `get_sc_obs(x)$<column>`, or of the
#' QC-passing cells only, i.e. `get_sc_obs(x, filtered = TRUE)$<column>`. Which
#' one you passed is inferred from the length, so a vector matching neither is
#' an error rather than a silently wrong purity.
#' @param add_additional_info String. Which label information to add on top of
#' the purity. One of `c("none", "top_label", "top_two_labels")`. Defaults to
#' `"none"`, i.e. the purity only.
#' @param add_entropy Boolean. Shall the normalised Shannon entropy of the
#' label distribution be added as a diversity measure. Defaults to `FALSE`.
#'
#' @returns The `MetaCells` with added columns to the observation table:
#' \itemize{
#'   \item mc_purity - Fraction of the meta cell's cells that carry the most
#'   abundant label. Always added.
#'   \item mc_top_label - Name of the most abundant label. Added for
#'   `add_additional_info %in% c("top_label", "top_two_labels")`.
#'   \item mc_second_label - Name of the second most abundant label, `NA` for a
#'   pure meta cell. Added for `add_additional_info = "top_two_labels"`.
#'   \item mc_second_frac - Fraction of the meta cell's cells carrying that
#'   second label, `0` for a pure meta cell. Added for
#'   `add_additional_info = "top_two_labels"`.
#'   \item mc_entropy - Normalised Shannon entropy of the label distribution,
#'   see details. Added for `add_entropy = TRUE`.
#' }
#'
#' @export
calc_meta_cell_purity <- S7::new_generic(
  name = "calc_meta_cell_purity",
  dispatch_args = "object",
  fun = function(
    object,
    original_cell_type,
    add_additional_info = c("none", "top_label", "top_two_labels"),
    add_entropy = FALSE
  ) {
    S7::S7_dispatch()
  }
)

#' @method calc_meta_cell_purity MetaCells
S7::method(calc_meta_cell_purity, MetaCells) <- function(
  object,
  original_cell_type,
  add_additional_info = c("none", "top_label", "top_two_labels"),
  add_entropy = FALSE
) {
  add_additional_info <- match.arg(add_additional_info)

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

  # the single-column branch of `[[<-` asserts an atomic value
  object[[names(new_cols)]] <- if (length(new_cols) == 1L) {
    new_cols[[1L]]
  } else {
    new_cols
  }

  return(object)
}

### meta cell diffusion coords -------------------------------------------------

#' Calculate diffusion coordinates
#'
#' @description
#' To leverage the quality metrics from Persad, et al., we need the diffusion
#' coordinates to then calculate if a cell is a dense or sparse region of the
#' manifold, its compactness and separation to other meta cells. To do so,
#' generate a diffusion map on the original data based on the approach of
#' SEACells and add the data to the object, see Persad, et al.
#'
#' @param object `MetaCells` class.
#' @param knn_data `SingleCellNearestNeighbour` class. Contains the kNN graph
#' from the original cells.
#' @param n_dcs Integer. Number of diffusion coordinates to use. Defaults to
#' `10L`.
#' @param k_density Integer. The k-th neighbour to use for the density region
#' estimation. Defaults to `150L`.
#' @param seed Integer. Seed for reproducibility
#' @param .verbose Boolean or integer. Controls verbosity and returns run times.
#' `FALSE` -> quiet, `TRUE` or `1L` -> normal verbosity, `2L` -> detailed
#' verbosity.
#'
#' @returns The class with the diffusion map coordinates, density distance and
#' region attached.
#'
#' @export
#'
#' @references
#' Persad, et al. Nat Biotechnol, 2023
calc_diffusion_coordinates <- S7::new_generic(
  name = "calc_diffusion_coordinates",
  dispatch_args = "object",
  fun = function(
    object,
    knn_data,
    n_dcs = 10L,
    k_density = 150L,
    seed = 42L,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method calc_diffusion_coordinates MetaCells
S7::method(calc_diffusion_coordinates, MetaCells) <- function(
  object,
  knn_data,
  n_dcs = 10L,
  k_density = 150L,
  seed = 42L,
  .verbose = TRUE
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(object, MetaCells))
  checkmate::assertClass(knn_data, "SingleCellNearestNeighbour")
  checkmate::qassert(n_dcs, "I1")
  checkmate::qassert(k_density, "I1")
  checkmate::qassert(seed, "I1")
  checkmate::qassert(.verbose, c("B1", "I1[0,2]"))

  # every source keeps `original_cell_idx` in its own index space, so a single
  # kNN over one source's cells cannot be resolved against the merged memberships
  if (isTRUE(S7::prop(object, "is_merged"))) {
    warning(
      paste(
        "The meta cells were merged, so `original_cell_idx` cannot be resolved",
        "against a single kNN graph.",
        "Run calc_diffusion_coordinates() per source before merging.",
        "Returning object as is."
      )
    )

    return(object)
  }

  # the kNN rows cover the QC-passing cells, the memberships index the full obs
  # space. Resolve before the density computation so a mismatch fails fast.
  mc_rows <- .mc_artefact_rows(object, length(knn_data$used_cells))

  # deal with knn params here
  knn_params <- params_knn_defaults()
  knn_params$k <- 1L

  res <- rs_metacell_density(
    knn_data = knn_data,
    n_dcs = n_dcs,
    k_density = k_density,
    knn_params = knn_params,
    verbose = parse_verbosity(.verbose),
    seed = seed
  )

  density_region <- purrr::map_chr(
    mc_rows,
    function(rows) {
      region <- res$regions[rows]
      names(which.max(table(region)))
    }
  )

  object[["density_region"]] <- density_region

  S7::prop(object = object, name = "other_data")[["dcs"]] <- res$dcs
  S7::prop(object = object, name = "other_data")[[
    "density_dist"
  ]] <- res$density_distances
  S7::prop(object = object, name = "other_data")[["regions"]] <- res$regions

  return(object)
}

### manifold metrics -----------------------------------------------------------

#' Calculate manifold metrics
#'
#' @description
#' This function will calculate the compactness and separation of your metacells
#' on the manifold (defined by the diffusion map). You must have run
#' [calc_diffusion_coordinates()] before calling this function. The idea is
#' that compactness indicates how tight the metacell spans the manifold, whereas
#' separation indicates how well the different metacells span the manifold.
#'
#' @param object `MetaCells` class for which to calculate the different metrics.
#'
#' @returns The class with the compactness and separation scores added.
#'
#' @export
#'
#' @references
#' Persad, et al. Nat Biotechnol, 2023
calc_manifold_metrics <- S7::new_generic(
  name = "calc_manifold_metrics",
  dispatch_args = "object",
  fun = function(
    object
  ) {
    S7::S7_dispatch()
  }
)

#' @method calc_manifold_metrics MetaCells
S7::method(calc_manifold_metrics, MetaCells) <- function(
  object
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(object, MetaCells))

  # the diffusion map lives in the source's cell index space, which a merged
  # object no longer has. Say so, rather than pointing at a function that will
  # refuse as well.
  if (isTRUE(S7::prop(object, "is_merged"))) {
    warning(
      paste(
        "The meta cells were merged, so there is no single diffusion map to",
        "resolve `original_cell_idx` against.",
        "Run calc_manifold_metrics() per source before merging.",
        "Returning object as is."
      )
    )

    return(object)
  }

  dcs <- S7::prop(object, "other_data")[["dcs"]]

  if (is.null(dcs)) {
    warning(
      paste(
        "No diffusion coordinates found.",
        "Please run calc_diffusion_coordinates().",
        "Returning object as is."
      )
    )

    return(object)
  }

  # the diffusion map has one row per QC-passing cell, the memberships index the
  # full obs space. Rust shifts these to 0-based itself, so hand it 1-based rows.
  mc_rows <- .mc_artefact_rows(object, nrow(dcs))

  # calculate compactness
  compactness <- rs_metacell_compactness(
    dc = dcs,
    mc_rows
  )

  # separation
  separation <- rs_metacell_separation(
    dc = dcs,
    mc_rows
  )

  object[["compactness"]] <- compactness
  object[["separation"]] <- separation

  return(object)
}

## processing ------------------------------------------------------------------

### hvg ------------------------------------------------------------------------

#### with object state ---------------------------------------------------------

# generic found in R/base_generics_sc.R

#' @method find_hvg_sc MetaCells
#'
#' @export
#'
#' @importFrom zeallot %<-%
#' @importFrom magrittr %>%
S7::method(find_hvg_sc, MetaCells) <- function(
  object,
  hvg_no = 2000L,
  hvg_params = params_sc_hvg(),
  streaming = NULL,
  .verbose = TRUE
) {
  checkmate::assertTRUE(S7::S7_inherits(object, MetaCells))
  checkmate::qassert(hvg_no, "I1")
  assertScHvg(hvg_params)
  checkmate::qassert(streaming, c("B1", "0"))
  checkmate::qassert(.verbose, c("B1", "I1[0,2]"))

  assay <- if (hvg_params$method == "vst") {
    "raw"
  } else {
    "norm"
  }

  count_list <- mc_counts_to_list(object = object, assay = assay)

  res <- with(
    hvg_params,
    rs_mc_hvg(
      sparse_data = count_list,
      hvg_method = method,
      loess_span = loess_span,
      binning = bin_method,
      n_bins = num_bin,
      clip_max = NULL,
      verbose = parse_verbosity(.verbose)
    )
  )

  # `:=` on the property alone silently no-ops once the data.table's
  # over-allocation is gone, e.g. after a saveRDS/readRDS round trip
  var_table <- data.table::copy(S7::prop(object, "var_table"))
  var_table[, names(res) := res]
  S7::prop(object, "var_table") <- var_table

  hvg <- switch(
    hvg_params$method,
    "vst" = order(res$var_std, decreasing = TRUE)[1:hvg_no],
    "dispersion" = order(res$dispersion, decreasing = TRUE)[1:hvg_no],
    "meanvarbin" = order(res$dispersion_scaled, decreasing = TRUE)[1:hvg_no],
    stop("Unknown HVG method: ", hvg_params$method)
  )

  object <- set_hvg(object, hvg = hvg)

  return(object)
}

#### without changing object state ---------------------------------------------

#' @method get_hvg_data_sc MetaCells
#'
#' @export
S7::method(get_hvg_data_sc, MetaCells) <- function(
  object,
  cell_ids = NULL,
  hvg_no = 3000L,
  hvg_params = params_sc_hvg(),
  streaming = NULL,
  .verbose = TRUE
) {
  checkmate::assertTRUE(S7::S7_inherits(object, MetaCells))
  checkmate::qassert(cell_ids, c("0", "S+"))
  checkmate::qassert(hvg_no, "I1")
  assertScHvg(hvg_params)
  checkmate::qassert(streaming, c("B1", "0"))
  checkmate::qassert(.verbose, c("B1", "I1[0,2]"))

  cell_indices <- if (is.null(cell_ids)) {
    NULL
  } else {
    get_cell_indices(object, cell_ids = cell_ids, rust_index = FALSE)
  }

  assay <- if (hvg_params$method == "vst") "raw" else "norm"

  count_list <- mc_counts_to_list(
    object = object,
    cell_indices = cell_indices,
    assay = assay
  )

  res <- with(
    hvg_params,
    rs_mc_hvg(
      sparse_data = count_list,
      hvg_method = method,
      loess_span = loess_span,
      binning = bin_method,
      n_bins = num_bin,
      clip_max = NULL,
      verbose = parse_verbosity(.verbose)
    )
  )

  var_table <- get_sc_var(object, cols = c("gene_idx", "gene_id"))

  build_hvg_table(
    var_table = var_table,
    res = res,
    hvg_no = hvg_no,
    hvg_method = hvg_params$method
  )
}


### pca ------------------------------------------------------------------------

# generic found in R/base_generics_sc.R

#' @method calculate_pca_sc MetaCells
#'
#' @export
#'
#' @importFrom zeallot %<-%
#' @importFrom magrittr %>%
S7::method(calculate_pca_sc, MetaCells) <- function(
  object,
  no_pcs,
  pca_params = params_sc_pca(),
  sparse_svd = FALSE,
  hvg = NULL,
  seed = 42L,
  .verbose = TRUE
) {
  checkmate::assertTRUE(S7::S7_inherits(object, MetaCells))
  checkmate::qassert(no_pcs, "I1")
  assertScPca(pca_params)
  checkmate::qassert(sparse_svd, "B1")
  checkmate::qassert(hvg, c("I+", "0"))
  checkmate::qassert(seed, "I1")
  checkmate::qassert(.verbose, c("B1", "I1[0,2]"))

  if ((length(get_hvg(object)) == 0) && is.null(hvg)) {
    warning(paste(
      "No HVGs identified in the object nor provided.",
      "Please run find_hvg_sc() or provide the indices of the HVG",
      "Returning object as is."
    ))
    return(object)
  }

  selected_hvg <- if (!is.null(hvg)) {
    if (.verbose) {
      message(
        paste(
          "HVGs provided.",
          "Will use these ones and set the internal HVG to the provided genes."
        )
      )
    }
    # the method here uses the R 1 indices
    object <- set_hvg(object, hvg)
    hvg
  } else {
    get_hvg(object)
  }

  count_list <- mc_counts_to_list(
    object = object,
    gene_indices = selected_hvg,
    assay = "norm"
  )

  clr_offsets <- if (pca_params$clr) {
    mc_get_clr_offsets(object)
  } else {
    NULL
  }

  zeallot::`%<-%`(
    c(pca_factors, pca_loadings, singular_values),
    rs_mc_pca(
      sparse_data = count_list,
      no_pcs = no_pcs,
      pca_params = pca_params,
      clr_offsets = clr_offsets,
      seed = seed,
      verbose = parse_verbosity(.verbose)
    )
  )

  object <- set_pca_factors(object, pca_factors)
  object <- set_pca_loadings(object, pca_loadings)
  object <- set_pca_singular_vals(object, singular_values[1:no_pcs])

  return(object)
}
