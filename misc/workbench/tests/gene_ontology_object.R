library(arrow)
library(data.table)
library(foreach)

path_to_biomind_processed = "~/Desktop/biomind_downloads/processed_data/"

list.files(path_to_biomind_processed)

go_nodes = arrow::read_parquet(file.path(path_to_biomind_processed, "nodes_OT_gene_ontology.parquet")) %>%
  as.data.table()

go_to_genes = arrow::read_parquet(
  file.path(path_to_biomind_processed, "edges_OT_pathways.parquet"),
  col_select = c(":START_ID", ":END_ID")
) %>% as.data.table() %>%
  setnames(., c(":START_ID", ":END_ID"), c("from", "to"))

head(go_to_genes)

ontology = arrow::read_parquet(
  file.path(path_to_biomind_processed, "edges_UP_ontologies.parquet")
) %>% as.data.table() %>%
  setnames(., c(":START_ID", ":END_ID"), c("from", "to"))

go_to_genes_clean = go_to_genes[from %like% 'ENSG' &
                                  to %like% 'GO:'] %>%
  .[, .(collapsed = list(unique(from))), to] %>%
  setnames(., c("to", "collapsed"), c("go_id", "ensembl_id"))

ontology_clean = ontology[from %like% 'GO:' &
           to %like% 'GO:', c("from", "to")]

go_ontology_igraph = igraph::graph_from_data_frame(ontology_clean)

key_nodes = setNames(
  c("GO:0008150", "GO:0005575", "GO:0003674"),
  c(
    "biological_process",
    "cellular_component",
    "molecular_function"
  )
)

dist_to_nodes = foreach(i = seq_along(key_nodes)) %do% {
  # Index
  index = which(igraph::V(go_ontology_igraph)$name == key_nodes[[i]])
  dfs_search = igraph::bfs(
    go_ontology_igraph,
    root = index,
    father = T,
    dist = T,
    unreachable = F,
    mode = "out"
  )
}

dist_to_nodes_combined = do.call(cbind, purrr::map(dist_to_nodes, ~ .$dist )) %>%
  `colnames<-`(key_nodes)

dist_to_nodes_combined = dist_to_nodes_combined
dist_to_nodes_combined[dist_to_nodes_combined < 0] = NA

dist_to_nodes = matrixStats::rowMaxs(dist_to_nodes_combined, na.rm = T)
dist_to_nodes[is.infinite(dist_to_nodes)] = 0

level_df = data.table(
  go_id = names(dist_to_nodes),
  depth = dist_to_nodes
)


ancestor_DT = go_ontology_igraph %>%
  igraph::ego(order = igraph::vcount(go_ontology_igraph),
              mode = "out") %>%
  setNames(igraph::V(go_ontology_igraph)$name) %>%
  Map(f = names) %>%
  stack() %>%
  rev() %>%
  setNames(names(ontology_clean)) %>%
  as.data.table() %>%
  .[from != to] %>%
  dplyr::mutate_all(as.character) %>%
  .[, .(ancestors = list(from)), .(to)] %>%
  dplyr::rename(go_id = to)




go_data_final = purrr::reduce(
  list(
    go_nodes[, c('nodeID', 'node_name')] %>% dplyr::rename(go_id = nodeID, go_name = node_name),
    go_to_genes_clean,
    ancestor_DT,
    level_df
  ),
  merge,
  by = 'go_id'
) %>%
  setorder(-depth)
