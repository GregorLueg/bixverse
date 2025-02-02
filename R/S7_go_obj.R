# S7 ----

gene_ontology_data <- S7::new_class(
  # Name
  name = "gene_ontology_data",
  # Properties, i.e., slots
  properties = list(
    go_info = S7::class_data.frame,
    go_to_genes = S7::class_list,
    ancestry = S7::class_list,
    levels = S7::class_list,
    min_genes = S7::class_integer
  ),

  #' Gene Ontology data
  #'
  #' @description
  #' This class is used to store the gene ontology information for usage in GSE
  #' elimination methods.
  #'
  #' - go_info: data.table. Contains the gene ontology identifiers and names.
  #' - go_to_genes: List. Contains the genes within each gene ontology term.
  #' - ancestry: List. Contains the ancestors for each gene ontology term.
  #' - levels: List. Which gene ontology terms sit at which level.
  #' - min_genes: Integer, the minimum genes in the gene ontology term to
  #' conduct the test.
  #'
  #' @param go_data_dt A data.table that contains the gene ontology information.
  #' This can for example be produced with `biomind_to_go_data()`.
  #' @param min_genes data.frame. Meta-data information in form of a data.frame.
  #'
  #' @return Returns the S7 object for further operations.
  #'
  #' @export
  constructor = function(go_data_dt, min_genes) {
    # Checks
    checkmate::assertDataTable(go_data_dt)
    checkmate::qassert(min_genes, "I1")
    go_data_dt <-
      go_data_dt[, `:=`(
        no_genes = purrr::map_dbl(ensembl_id, length),
        depth = sprintf("%02d", depth)
      )]

    go_data_dt <- go_data_dt[no_genes >= min_genes]

    go_info <- go_data_dt[, c("go_id", "go_name")]

    go_to_genes <- go_data_dt$ensembl_id
    names(go_to_genes) <- go_data_dt$go_id

    ancestry <- go_data_dt$ancestors
    names(ancestry) <- go_data_dt$go_id

    depth_df <- go_data_dt[, .(go_ids = list(go_id)), .(depth)]

    levels <- depth_df$go_ids
    names(levels) <- depth_df$depth

    # Finalise object
    S7::new_object(
      S7::S7_object(),
      go_info = go_info,
      go_to_genes = go_to_genes,
      ancestry = ancestry,
      levels = levels,
      min_genes = min_genes
    )
  }
)

# Raw data helpers ----

#' Wrapper to clean up the GO identifiers to ensembl identifer table.
#'
#' @param dt The data.table containing the pathway to gene edges from BioMind.
#'
#' @return data.table with go_id and a list_col of collapsed, associated ensembl
#' identifiers.
#'
#' @export
#'
#' @import data.table
#' @import foreach
#' @importFrom magrittr `%>%`
.get_go_genes_bm <- function(dt) {
  # Avoid binding erros
  `.` <- from <- to <- NULL
  # Checks
  checkmate::assertDataTable(dt)
  # Checks
  checkmate::assertDataTable(dt)
  # Function body
  df_final <- dt %>%
    setnames(., c(":START_ID", ":END_ID"), c("from", "to")) %>%
    .[, c("from", "to")] %>%
    .[from %like% "ENSG" & to %like% "GO:"] %>%
    .[, .(collapsed = list(unique(from))), to] %>%
    setnames(., c("to", "collapsed"), c("go_id", "ensembl_id"))

  df_final
}

#' Wrapper to get the ancestry and the ontology depth level.
#'
#' @param dt The data.table containing the ontology edges from BioMind.
#'
#' @return A list with two data.tables with the first containing the depth
#' information and the second the ancestry information
#'
#' @export
#'
#' @import data.table
#' @importFrom magrittr `%>%`
.get_go_depth_ancestry_bm <- function(dt) {
  # Bindings
  `.` <- from <- to <- i <- foreach <- `%do%` <- NULL
  # Checks
  checkmate::assertDataTable(dt)
  # Function body
  ontology <- dt %>%
    setnames(., c(":START_ID", ":END_ID"), c("from", "to")) %>%
    .[from %like% "GO:" &
      to %like% "GO:", c("from", "to")] %>%
    unique()

  # Create igraph
  go_ontology_igraph <-
    igraph::graph_from_data_frame(ontology)

  # Identify distance to key nodes
  key_nodes <- setNames(
    c("GO:0008150", "GO:0005575", "GO:0003674"),
    c(
      "biological_process",
      "cellular_component",
      "molecular_function"
    )
  )

  dist_to_nodes <- foreach(i = seq_along(key_nodes)) %do% {
    # Index
    index <-
      which(igraph::V(go_ontology_igraph)$name == key_nodes[[i]])
    dfs_search <- igraph::bfs(
      go_ontology_igraph,
      root = index,
      father = TRUE,
      dist = TRUE,
      unreachable = FALSE,
      mode = "out"
    )

    dfs_search
  }

  dist_to_nodes_combined <-
    do.call(cbind, purrr::map(dist_to_nodes, ~ .$dist)) %>%
    `colnames<-`(key_nodes)

  dist_to_nodes_combined <- dist_to_nodes_combined
  dist_to_nodes_combined[dist_to_nodes_combined < 0] <- NA

  dist_to_nodes <-
    matrixStats::rowMaxs(dist_to_nodes_combined, na.rm = TRUE)
  dist_to_nodes[is.infinite(dist_to_nodes)] <- 0

  level_df <- data.table(
    go_id = names(dist_to_nodes),
    depth = dist_to_nodes
  )

  # Get the ancestors
  ancestor_dt <- go_ontology_igraph %>%
    igraph::ego(
      order = igraph::vcount(go_ontology_igraph),
      mode = "out"
    ) %>%
    setNames(igraph::V(go_ontology_igraph)$name) %>%
    Map(f = names) %>%
    stack() %>%
    rev() %>%
    setNames(names(ontology)) %>%
    as.data.table() %>%
    .[from != to] %>%
    dplyr::mutate_all(as.character) %>%
    .[, .(ancestors = list(from)), .(to)] %>%
    dplyr::rename(go_id = to)

  list(level_df, ancestor_dt)
}


#' Generates the GO data from the BioMind pre-processed data.
#'
#' The function needs to be pointed to the processed BioMind data and will
#' automatically ingest and process the gene ontology data for usage for GSE
#' with elimination method.
#'
#' @param path_to_biomind_processed The path to the processed BioMind data with
#' the .parquet files in there.
#' @param .verbose
#'
#' @return A data.table with all of the data ready for ingestion into the
#' S7_GO_obj.
#'
#' @export
#'
#' @import data.table
#' @importFrom magrittr `%>%`
#' @importFrom zeallot `%<-%`
biomind_to_go_data <- function(path_to_biomind_processed, .verbose = TRUE) {
  # Bindings
  go_ancestry <- go_depth <- `%<-%` <- node_id <- node_name <- depth <- NULL
  # Check that expected files exist
  checkmate::qassert(path_to_biomind_processed, "S1")
  checkmate::qassert(.verbose, "B1")
  path_nodes <- file.path(
    path_to_biomind_processed,
    "nodes_OT_gene_ontology.parquet"
  )
  path_go_genes <-
    file.path(path_to_biomind_processed, "edges_OT_pathways.parquet")
  path_ontology <-
    file.path(path_to_biomind_processed, "edges_UP_ontologies.parquet")
  checkmate::assertFileExists(path_nodes)
  checkmate::assertFileExists(path_go_genes)
  checkmate::assertFileExists(path_ontology)
  # Function body
  if (.verbose) {
    print("1. Loading data in.")
  }
  go_nodes <- arrow::read_parquet(path_nodes) %>%
    as.data.table()
  edges_go_genes <- arrow::read_parquet(path_go_genes) %>%
    as.data.table()
  edges_ontology <- arrow::read_parquet(path_ontology) %>%
    as.data.table()

  if (.verbose) {
    print("2. Processing gene ontology to genes.")
  }

  go_genes <- .get_go_genes_bm()(edges_go_genes)

  if (.verbose) {
    print("3. Processing gene ontological information.")
  }

  c(go_ancestry, go_depth) %<-% .get_go_depth_ancestry_bm(edges_ontology)

  df_list <- list(
    go_nodes[, c("node_id", "node_name")] %>%
      dplyr::rename(go_id = node_id, go_name = node_name),
    go_genes,
    go_ancestry,
    go_depth
  )

  final_res <- purrr::reduce(df_list,
    merge,
    by = "go_id"
  ) %>%
    setorder(-depth)

  final_res
}
