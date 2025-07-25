# parallelisation stuff --------------------------------------------------------

#' Helper function for core detection.
#'
#' @description
#' Identifies the number of cores/workers to use for parallel tasks. Will
#' default to half of the available cores, to a maximum of `abs_max_workers`.
#'
#' @param abs_max_workers Integer. Absolute maximum number of cores to use.
#' Defaults to `8L`.
#'
#' @returns Integer. Number of cores to use.
get_cores <- function(abs_max_workers = 8L) {
  checkmate::qassert(abs_max_workers, "I1")
  cores <- as.integer(parallel::detectCores() / 2)
  cores <- ifelse(cores >= abs_max_workers, abs_max_workers, cores)
  return(cores)
}

# user options -----------------------------------------------------------------

#' Helper function for user option selection
#'
#' @param options Character vector. The options the user has to choose from.
#'
#' @returns The chosen option
select_user_option <- function(options) {
  cat("Please select an option:\n")
  for (i in seq_along(options)) {
    cat(i, ": ", options[i], "\n")
  }

  while (TRUE) {
    cat("\nEnter the number of your choice: ")
    user_choice <- as.integer(readLines(n = 1))

    if (
      !is.na(user_choice) && user_choice > 0 && user_choice <= length(options)
    ) {
      return(options[user_choice])
    } else {
      cat("Invalid selection. Please try again.\n")
    }
  }
}

# gene set libraries -----------------------------------------------------------

## gene ontology ---------------------------------------------------------------

### loading data ---------------------------------------------------------------

#' Get the Gene Ontology data human
#'
#' @description
#' This function loads in gene ontology data stored in the package. This is
#' for humans only.
#'
#' @returns A list containing:
#' \itemize{
#'  \item go_info - data.table. The gene ontology identifier, name and namespace
#'  can be found in this one
#'  \item gene_ontology - data.table. The relationships between different
#'  gene ontology terms.
#'  \item go_to_genes - data.table. The gene ontology term to gene (ensembl id)
#'  relationships.
#' }
#'
#' @export
load_go_human_data <- function() {
  files_to_load <- c(
    "go_info" = "gene_ontology_info.parquet",
    "gene_ontology" = "go_relationship_info.parquet",
    "go_to_genes" = "go_to_genes.parquet"
  )

  files <- purrr::map(files_to_load, \(file_name) {
    path <- system.file(
      "extdata",
      file_name,
      package = "bixverse"
    )
    if (path == "") {
      stop(sprintf("The expected extdata was not found for %s", file_name))
    }
    df <- data.table::setDT(arrow::read_parquet(
      path
    ))
    df
  })

  return(files)
}

### processing -----------------------------------------------------------------

#' Process Gene Ontology data into the right format
#'
#' @description Helper function that takes in different files containing gene
#' ontology data and puts them together for various gene set enrichment methods
#' using the ontological information
#'
#' @param go_info data.table. Contains `go_id`, `go_name` and `namespace.`
#' @param go_genes data.table. Contains `go_id` and corresponding `ensembl_id`.
#' @param go_relationships data.table. Contains `parent`, `child` and
#' `relationship`
#'
#' @returns data.table ready for usage in [bixverse::gene_ontology_data()].
#'
#' @export
#' @import data.table
#' @importFrom magrittr `%>%`
process_go_data <- function(go_info, go_genes, go_relationships) {
  # scope
  . <- ind <- NULL

  # checks
  checkmate::assertDataTable(go_info)
  checkmate::assertDataTable(go_genes)
  checkmate::assertDataTable(go_relationships)
  checkmate::assertNames(
    names(go_info),
    must.include = c("go_id", "go_name", "namespace")
  )
  checkmate::assertNames(
    names(go_genes),
    must.include = c("go_id", "ensembl_id")
  )
  checkmate::assertNames(
    names(go_relationships),
    must.include = c("parent", "child", "relationship")
  )

  # function
  go_genes_summarised <- go_genes[, .(ensembl_id = list(ensembl_id)), .(go_id)]

  ancestry <- get_ontology_ancestry(parent_child_dt = go_relationships)

  ancestor_dt <- data.table::setDT(stack(ancestry$ancestors))[,
    ind := as.character(ind)
  ] %>%
    .[, .(ancestors = list(values)), .(ind)] %>%
    data.table::setnames(old = "ind", new = "go_id")

  go_levels <- get_go_levels(go_relationships)

  go_data_final <- Reduce(
    function(x, y) merge(x, y, by = "go_id"),
    list(
      go_info,
      go_genes_summarised,
      ancestor_dt,
      go_levels
    )
  ) %>%
    data.table::setorder(-depth)

  return(go_data_final)
}

#' Wrapper function to load and process the gene ontology data.
#'
#' @description
#' This function loads in gene ontology data stored in the package and processes
#' it into the format for [bixverse::gene_ontology_data()]. Wraps
#' [bixverse::load_go_human_data()] and [bixverse::process_go_data()] into one.
#'
#' @param filter_relationships Boolean. Shall the ontology be filtered to
#' only `"is_a"` and `"part_of"` relationships. Defaults to TRUE.
#' @param .verbose Boolean. Controls verbosity of the function.
#'
#' @returns A data.table
#'
#' @export
get_go_data_human <- function(filter_relationships = TRUE, .verbose = TRUE) {
  # checks
  checkmate::qassert(.verbose, "B1")

  if (.verbose) {
    message("Loading the data from the package.")
  }

  go_data <- load_go_human_data()

  go_info <- go_data$go_info
  go_genes <- go_data$go_to_genes
  go_relationships <- go_data$gene_ontology %>%
    setnames(
      old = c("from", "to"),
      new = c("parent", "child")
    )

  if (filter_relationships) {
    go_relationships <- go_relationships[relationship %in% c("is_a", "part_of")]
  }

  if (.verbose) {
    message("Processing data for the gene_ontology class.")
  }

  go_data_final <- process_go_data(
    go_info = go_info,
    go_genes = go_genes,
    go_relationships = go_relationships
  )

  return(go_data_final)
}

### helpers --------------------------------------------------------------------

#' Helper to get the Gene Ontology levels
#'
#' @description
#' Gets the ontology depth based on the edge dt
#'
#' @param edge_dt data.table. The gene ontology edge data (i.e., term
#' connections between the different terms)
#'
#' @return A data.table with the identifier and depth.
get_go_levels <- function(edge_dt) {
  # checks
  checkmate::assertDataTable(edge_dt)

  go_hierarchy_graph <- igraph::graph_from_data_frame(
    edge_dt,
    directed = TRUE
  )

  key_go_ids <- c("GO:0008150", "GO:0003674", "GO:0005575")

  distances <- vector(mode = "list", length = 3L)

  for (i in seq_along(key_go_ids)) {
    index <- which(igraph::V(go_hierarchy_graph)$name == key_go_ids[[i]])

    dfs_search <- igraph::dfs(
      go_hierarchy_graph,
      root = index,
      order.out = T,
      father = T,
      dist = T,
      unreachable = F,
      mode = "out"
    )$dist

    dfs_search[dfs_search < 0] <- NA

    distances[[i]] <- dfs_search
  }

  level <- matrixStats::rowMins(do.call(cbind, distances), na.rm = TRUE)

  level_dt <- data.table::data.table(
    go_id = names(level),
    depth = level
  )

  return(level_dt)
}

# sparse matrices --------------------------------------------------------------

#' Transform an upper triangle-stored matrix to a sparse one
#'
#' @param upper_triangle_vals Numerical vector. The values of the upper triangle
#' stored in a row major format.
#' @param shift Integer. Did you shift the diagonal up.
#' @param n Integer. Number of columns and rows of the symmetric matrix.
#'
#' @returns The sparse matrix.
#'
#' @export
upper_triangle_to_sparse <- function(upper_triangle_vals, shift, n) {
  # checks
  checkmate::qassert(upper_triangle_vals, "N+")
  checkmate::qassert(shift, "I1")
  checkmate::qassert(n, "I1")

  # functions
  data <- rs_upper_triangle_to_sparse(upper_triangle_vals, as.integer(shift), n)

  matrix = Matrix::sparseMatrix(
    i = data$row_indices + 1,
    p = data$col_ptr,
    x = unlist(data$data),
    dims = c(n, n)
  )

  return(matrix)
}

#' Helper function to transform the Rust-exported sparse matrices into R ones
#'
#' @param ls List. Needs to represent the (column) sparse data.
#'
#' @returns The sparseMatrix from the data.
#'
#' @export
sparse_list_to_mat <- function(ls) {
  # checks
  checkmate::assertList(ls, types = "numeric", names = "named")
  checkmate::assertNames(
    names(ls),
    must.include = c("data", "row_indices", "col_ptr", "ncol", "nrow")
  )

  # body
  sparse_mat <- Matrix::sparseMatrix(
    i = ls$row_indices + 1,
    p = ls$col_ptr,
    x = ls$data,
    dims = c(ls$nrow, ls$ncol)
  )

  return(sparse_mat)
}

# inflection points ------------------------------------------------------------

#' Identify the inflection point for elbow-like data
#'
#' @description
#' This function will identify the index of the inflection point of an arbitrary
#' series `x ~ y` via the biggest increase in the first derivative. Useful to
#' identify key points in elbow plots.
#'
#' @param x,y The x and y values.
#' @param span The span parameter for the loess function.
#'
#' @return A list containing:
#' \itemize{
#'   \item inflection_idx - Index of the inflection point
#'   \item gradient_change - Absolute change in the first derivative
#' }
get_inflection_point <- function(x, y, span = 0.5) {
  # Checks
  checkmate::assertNumeric(x, len = length(y))
  checkmate::assertNumeric(y, len = length(x))
  checkmate::qassert(span, "R+[0,1]")
  # Function body
  span <- max(0.1, min(1.0, span))
  fit <- tryCatch(
    {
      fit <- loess(
        y ~ x,
        span = span,
        family = "gaussian",
        degree = 2L,
        normalize = TRUE
      )
      list(fit = fit, warning = FALSE)
    },
    warning = function(w) {
      return(list(fit = NULL, warning = TRUE))
    }
  )

  # early return if the inflection point function is unhappy
  if (fit$warning) {
    return(list(
      inflection_idx = NULL,
      gradient_change = NULL
    ))
  }

  py <- predict(fit$fit, x)

  n <- length(x)
  gradient <- numeric(n)
  gradient[1] <- (py[2] - py[1]) / (x[2] - x[1])
  gradient[n] <- (py[n] - py[n - 1]) / (x[n] - x[n - 1])

  for (i in 2:(n - 1)) {
    gradient[i] <- (py[i + 1] - py[i - 1]) / (x[i + 1] - x[i - 1])
  }

  gradient_change <- diff(gradient)
  # One after the point with the biggest delta
  inflection_idx <- which.max(gradient_change) + 1

  return(
    list(inflection_idx = inflection_idx, gradient_change = gradient_change)
  )
}
