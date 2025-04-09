
library(Matrix)
library(tictoc)
library(igraph)
library(reshape2)
library(polars)
library(data.table)
library(dplyr)
library(igraph)
library(purrr)
library(furrr)

data_path = "/Users/liesbeth/Datascience/data/processed/OpenTargets_platform"

ontology = pl$read_parquet(file.path(data_path,"OT_edges_disease_hierarchy.parquet"))$
  filter(pl$col("source:string") == "efo")$
  rename(":END_ID" = "to",
         ":START_ID" = "from",
         "relation_type:string" = "relation")$
  select("from", "to", "relation")
nodes = pl$read_parquet(file.path(data_path,"OT_nodes_diseases_phenotypes.parquet"))$
  filter(pl$col("source:string") == "efo")$
  rename("id:ID" = "id",
         "name:string" = "name")$
  select("id", "name")
total_size = nrow(nodes)
max_ic = -log2(1/total_size)

ontology_ancestors = ontology$filter(pl$col("relation") == "ancestor_of")

calculate_information_content <- function(ontology_ancestors){
  checkmate::assert(class(ontology_ancestors) == "RPolarsDataFrame")
  checkmate::assert("from" %in% ontology_ancestors$columns)

  information_content = ontology_ancestors$
    group_by("from")$agg(pl$len()$alias("nb_desc"))$
    with_columns(pl$lit(total_size)$alias("total_terms"),
                 pl$lit(max_ic)$alias("max_ic"))$
    with_columns((-((pl$col("nb_desc")$cast(pl$Float32))$div(pl$col("total_terms")))$log(2))$alias("ic"))$
    with_columns(pl$col("ic")$div(pl$col("max_ic"))$alias("norm_ic"))$
    rename("from" = "id")
  ## add nodes without descendants
  onto_no_descendants = nodes$filter(!pl$col("id")$is_in(information_content$select("id")$to_series()))$
    with_columns(pl$lit(0)$cast(pl$UInt32)$alias("nb_desc"),
                 pl$lit(total_size)$alias("total_terms"),
                 pl$lit(max_ic)$alias("max_ic"),
                 pl$lit(max_ic)$alias("ic"),
                 pl$lit(1)$cast(pl$Float64)$alias("norm_ic"))$
    drop("name")
  information_content = pl$concat(information_content, onto_no_descendants, how = "diagonal")
  information_content
}

information_content = calculate_information_content(ontology_ancestors)

library(bixverse)
terms = as.data.table(nodes)$id
ancestor_list = as.data.table(ontology_ancestors)
ancestor_list = split(ancestor_list$from, ancestor_list$to)
ic_list = split(as.data.table(information_content)$ic, as.data.table(information_content)$id)
tic()
bxv_sim <- rs_onto_similarity_both(terms = terms,
                   ancestor_list = ancestor_list,
                   ic_list = ic_list,
                   max_ic = max(unlist(ic_list)))
toc()

bxv_sim_mat <- bixverse:::upper_triangular_cor_mat$new(
  cor_coef = bxv_sim$resnik_sim,
  features = bxv_sim$terms,
  shift = 1L
)$get_cor_matrix()

## some functions for the R function
get_index <- function(name, graph){
  match(name, V(graph)$name)
}
LCA = function(distance, n1, n2){
  d = rowSums(distance[, c(n1, n2)])
  d = d[!is.infinite(d)]
  if(length(d) == 0){
    NA
  }else{
    names(d[which.min(d)])
  }
}
get_distances <- function(ontology){
  graph = igraph::graph_from_data_frame(ontology )
  distance = igraph::distances(graph, V(graph), mode="out")
}
get_similarity <- function(information_content, ontology, term1, term2){
  checkmate::assert(all(c(term1, term2) %in% information_content$select("id")$to_series()$to_list()))
  ## calculate LCA
  distance <- get_distances(ontology)
  lca = LCA(distance, n1, n2)

  resnik = information_content$filter(pl$col("id") == names(lca))$select("ic", "norm_ic")$to_data_frame()
  lin = (2*resnik$ic)/(sum(information_content$filter(pl$col("id")$is_in(c(term1, term2)))$select("ic")$to_data_frame()))

  return(data.table(term1 = term1, term2 = term2, resnik = resnik, lin.ic = lin))
}
# run it in R only
ontology = ontology$filter(pl$col("relation") == "parent_of") %>% as.data.table()
nodes = nodes %>% as.data.table()

tic()
d = get_distances(ontology)
md = reshape2::melt(d)

seq.int <- seq(1, nrow(md), 2000000)
results = list()
n = 1
for(i in 1:length(seq.int)){
  start = seq.int[i]
  end = seq.int[i+1]-1
  if(is.na(end)){
    end = nrow(md)
  }
  tmp = md[start:end,] %>%
    as.data.table() %>%
    .[!is.infinite(value) & Var1 != Var2]
  results[[n]] <- tmp
  n = n + 1
}
results <- do.call(rbind, results)

## so now we have all disease pairs that connect somewhere in the tree
plan(multisession, workers = 5)
options(future.globals.maxSize = 6 * 1024^3)

lca_all = future_map2(results$Var1, results$Var2,
                      ~{
                        LCA(d, .x, .y)
                      })
results <- results[, lca := unlist(lca_all)]

results_pl <- pl$DataFrame(results)$
  with_columns(pl$col("Var1")$cast(pl$String)$alias("node1"), pl$col("Var2")$cast(pl$String)$alias("node2"))$
  join(information_content$select("id", "ic", "norm_ic")$rename("ic" = "resnik", "norm_ic" = "norm_resnik"),
       left_on = "lca", right_on = "id")$
  join(information_content$select("id", "ic")$rename("ic" = "ic_var1"), left_on = "node1", right_on = "id")$
  join(information_content$select("id", "ic")$rename("ic" = "ic_var2"), left_on = "node2", right_on = "id")$
  with_columns(((2*pl$col("resnik"))$div(pl$col("ic_var1") + pl$col("ic_var2")))$alias("lin"))$
  select("node1", "node2", "lca", "resnik", "norm_resnik", "lin")
toc()


n1 = "MONDO_0008648"
n2 = "Orphanet_7"
lca = LCA(d, n1, n2)
lca
information_content$filter(pl$col("id") == lca)

