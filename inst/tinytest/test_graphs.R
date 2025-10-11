# constrained page-rank tests --------------------------------------------------

## data ------------------------------------------------------------------------

### vector form ----------------------------------------------------------------

node_names <- c(
  "receptor_a",
  "receptor_b",
  "kinase_c",
  "kinase_d",
  "tf_e",
  "tf_f",
  "target_gene_g",
  "target_gene_h",
  "target_gene_i"
)

node_types <- c(
  "receptor",
  "receptor",
  "kinase",
  "kinase",
  "tf",
  "tf",
  "target_gene",
  "target_gene",
  "target_gene"
)

from <- c(
  "receptor_a",
  "receptor_a",
  "receptor_b",
  "kinase_c",
  "kinase_d",
  "tf_e",
  "tf_e",
  "tf_f",
  "tf_f",
  "tf_f"
)

to <- c(
  "kinase_c",
  "kinase_d",
  "kinase_d",
  "tf_e",
  "tf_f",
  "receptor_a",
  "target_gene_g",
  "receptor_b",
  "target_gene_h",
  "target_gene_i"
)

edge_type <- c(
  rep("activation", 3),
  rep("phosphorylation", 2),
  rep("tf_activation", 5)
)

edge_weight <- rep(1, 10)

personalisation_vec_1 <- c(1, rep(0, 8))
personalisation_vec_2 <- c(0, 1, rep(0, 7))

### igraph form ----------------------------------------------------------------

node_df <- data.frame(
  name = node_names,
  type = node_types
)

edge_df <- data.frame(
  from = from,
  to = to,
  weight = edge_weight,
  type = edge_type
)

g <- igraph::graph_from_data_frame(
  d = edge_df,
  vertices = node_df,
  directed = TRUE
)

### expected results -----------------------------------------------------------

expected_res_no_constraints_receptor_1 <- c(
  0.30167290,
  0.03881697,
  0.12821904,
  0.16136227,
  0.10897813,
  0.13700106,
  0.04631570,
  0.03881697,
  0.03881697
)

expected_res_no_constraints_receptor_2 <- c(
  0,
  0.33536714,
  0,
  0.28495089,
  0,
  0.24220825,
  0,
  0.06873686,
  0.06873686
)

names(expected_res_no_constraints_receptor_1) <- names(
  expected_res_no_constraints_receptor_2
) <- node_names

## tests -----------------------------------------------------------------------

### no constraints -------------------------------------------------------------

#### rust ----------------------------------------------------------------------

res_no_constraints_receptor_1 <- rs_constrained_page_rank(
  node_names = node_names,
  node_types = node_types,
  from = from,
  to = to,
  weights = edge_weight,
  edge_type = edge_type,
  personalised = personalisation_vec_1,
  sink_nodes = NULL,
  sink_edges = NULL
)

res_no_constraints_receptor_2 <- rs_constrained_page_rank(
  node_names = node_names,
  node_types = node_types,
  from = from,
  to = to,
  weights = edge_weight,
  edge_type = edge_type,
  personalised = personalisation_vec_2,
  sink_nodes = NULL,
  sink_edges = NULL
)

# only equivalence
expect_equivalent(
  current = res_no_constraints_receptor_1,
  target = expected_res_no_constraints_receptor_1,
  tolerance = 1e-7,
  info = "No constraints - diffusion from receptor 1"
)

expect_equivalent(
  current = res_no_constraints_receptor_2,
  target = expected_res_no_constraints_receptor_2,
  tolerance = 1e-7,
  info = "No constraints - diffusion from receptor 2"
)

#### igraph --------------------------------------------------------------------

res_no_constraints_receptor_1 <- constrained_page_rank(
  graph = g,
  personalisation_vector = personalisation_vec_1
)

res_no_constraints_receptor_2 <- constrained_page_rank(
  graph = g,
  personalisation_vector = personalisation_vec_2
)

expect_equal(
  current = res_no_constraints_receptor_1,
  target = expected_res_no_constraints_receptor_1,
  tolerance = 1e-7,
  info = "No constraints - diffusion from receptor 1 (igraph)"
)

expect_equal(
  current = res_no_constraints_receptor_2,
  target = expected_res_no_constraints_receptor_2,
  tolerance = 1e-7,
  info = "No constraints - diffusion from receptor 2 (igraph)"
)

#### list variant --------------------------------------------------------------

personalisation_list <- list(
  first = personalisation_vec_1,
  second = personalisation_vec_2
)

list_results <- constrained_page_rank_ls(
  graph = g,
  personalisation_list = personalisation_list
)

expect_equal(
  current = names(personalisation_list),
  target = names(list_results),
  info = "List version "
)

expect_equal(
  current = list_results$first,
  target = expected_res_no_constraints_receptor_1,
  tolerance = 1e-7,
  info = "No constraints - diffusion from receptor 1 (igraph - list)"
)

expect_equal(
  current = list_results$second,
  target = expected_res_no_constraints_receptor_2,
  tolerance = 1e-7,
  info = "No constraints - diffusion from receptor 2 (igraph - list)"
)

#### errors --------------------------------------------------------------------

g.2 <- data.table::copy(g)
g.2 <- igraph::delete_vertex_attr(g.2, "type")

g.3 <- data.table::copy(g)
g.3 <- igraph::delete_edge_attr(g.3, "type")


expect_error(
  current = constrained_page_rank(
    graph = g.2,
    personalisation_vector = personalisation_vec_1
  ),
  info = "Igraph error version 1"
)

expect_error(
  current = constrained_page_rank(
    graph = g.3,
    personalisation_vector = personalisation_vec_1
  ),
  info = "Igraph error version 2"
)

### edge constrained -----------------------------------------------------------

#### rust ----------------------------------------------------------------------

# different versions of constraints
res_edge_constraint_receptor_1.1 <- rs_constrained_page_rank(
  node_names = node_names,
  node_types = node_types,
  from = from,
  to = to,
  weights = edge_weight,
  edge_type = edge_type,
  personalised = personalisation_vec_1,
  sink_nodes = NULL,
  sink_edges = c("activation")
)

res_edge_constraint_receptor_1.2 <- rs_constrained_page_rank(
  node_names = node_names,
  node_types = node_types,
  from = from,
  to = to,
  weights = edge_weight,
  edge_type = edge_type,
  personalised = personalisation_vec_1,
  sink_nodes = NULL,
  sink_edges = c("phosphorylation")
)


res_edge_constraint_receptor_2.1 <- rs_constrained_page_rank(
  node_names = node_names,
  node_types = node_types,
  from = from,
  to = to,
  weights = edge_weight,
  edge_type = edge_type,
  personalised = personalisation_vec_2,
  sink_nodes = NULL,
  sink_edges = c("activation")
)

res_edge_constraint_receptor_2.2 <- rs_constrained_page_rank(
  node_names = node_names,
  node_types = node_types,
  from = from,
  to = to,
  weights = edge_weight,
  edge_type = edge_type,
  personalised = personalisation_vec_2,
  sink_nodes = NULL,
  sink_edges = c("phosphorylation")
)

# 6 should be 0, i.e., not reachable
expect_true(
  current = sum(res_edge_constraint_receptor_1.1 == 0) == 6,
  info = "Edge constraints version 1 - diffusion from receptor 1"
)

# 4 should be 0, i.e., not reachable
expect_true(
  current = sum(res_edge_constraint_receptor_1.2 == 0) == 4,
  info = "Edge constraints version 2 - diffusion from receptor 1"
)

# 7 should be 0, i.e., not reachable
expect_true(
  current = sum(res_edge_constraint_receptor_2.1 == 0) == 7,
  info = "Edge constraints version 1 - diffusion from receptor 2"
)

# 6 should be 0, i.e., not reachable
expect_true(
  current = sum(res_edge_constraint_receptor_2.2 == 0) == 6,
  info = "Edge constraints version 2 - diffusion from receptor 2"
)

### node constraint ------------------------------------------------------------

# this should be VERY close to res_edge_constraint_receptor_1.2
# due to the implementation it won't be perfect
res_node_constraint_receptor_1 <- rs_constrained_page_rank(
  node_names = node_names,
  node_types = node_types,
  from = from,
  to = to,
  weights = edge_weight,
  edge_type = edge_type,
  personalised = personalisation_vec_1,
  sink_nodes = c("tf"),
  sink_edges = NULL
)

# this should be VERY close to res_edge_constraint_receptor_2.1
# due to the implementation it won't be perfect
res_node_constraint_receptor_2 <- rs_constrained_page_rank(
  node_names = node_names,
  node_types = node_types,
  from = from,
  to = to,
  weights = edge_weight,
  edge_type = edge_type,
  personalised = personalisation_vec_2,
  sink_nodes = c("kinase"),
  sink_edges = NULL
)

#
res_node_constraint_receptor_3 <- rs_constrained_page_rank(
  node_names = node_names,
  node_types = node_types,
  from = from,
  to = to,
  weights = edge_weight,
  edge_type = edge_type,
  personalised = personalisation_vec_2,
  sink_nodes = c("kinase"),
  sink_edges = c("activation")
)

expect_true(
  current = cor(
    res_node_constraint_receptor_1,
    res_edge_constraint_receptor_1.2
  ) >=
    0.99,
  info = "Node constraint - diffusion from receptor 1"
)

expect_true(
  current = cor(
    res_node_constraint_receptor_2,
    res_edge_constraint_receptor_2.1
  ) >=
    0.99,
  info = "Node constraint - diffusion from receptor 2"
)

expect_true(
  current = cor(
    res_node_constraint_receptor_2,
    res_node_constraint_receptor_3
  ) >=
    0.99,
  info = "Constraint comparison"
)

# label propogation ------------------------------------------------------------

## synthetic data --------------------------------------------------------------

# fmt: skip
edge_list <- c(
  1, 2,   # Sample 1 (A) connects to 2
  1, 5,   # Sample 1 (A) connects to 5
  2, 3,   # Sample 2 connects to 3 (B)
  2, 4,   # Sample 2 connects to 4
  3, 4,   # Sample 3 (B) connects to 4
  3, 6,   # Sample 3 (B) connects to 6
  4, 7,   # Sample 4 connects to 7
  5, 6,   # Sample 5 connects to 6
  5, 9,   # Sample 5 connects to 9
  6, 7,   # Sample 6 connects to 7
  7, 8,   # Sample 7 connects to 8 (C)
  8, 9,   # Sample 8 (C) connects to 9
  8, 10,  # Sample 8 (C) connects to 10
  9, 10   # Sample 9 connects to 10
)
edge_list <- as.integer(edge_list)

labels <- c("A", NA, "B", NA, NA, NA, NA, "C", NA, NA)

propagated_labels <- knn_graph_label_propagation(
  edge_list = edge_list,
  labels = labels
)

expect_equal(
  current = propagated_labels$final_labels,
  target = c("A", "B", "B", "B", "A", "B", "C", "C", "C", "C"),
  info = "graph label propagation working as expected"
)

# igraph comparisons -----------------------------------------------------------

if (!requireNamespace("igraph", quietly = TRUE)) {
  exit_file("BiocNeighbors, bluster and cluster not available")
}

## data ------------------------------------------------------------------------

edge_dt <- data.table::data.table(
  from = c("a", "b", "c", "d", "d"),
  to = c("b", "c", "d", "a", "e")
)

edge_dt_weighted <- data.table::data.table(
  from = c("a", "b", "c", "d", "d"),
  to = c("b", "c", "d", "a", "e"),
  weight = c(1, 1, 0.5, 0.4, 0.25)
)

unique_nodes <- unique(c(edge_dt$from, edge_dt$to))

personalised_v1 <- c(1, 0, 0, 0, 0)
personalised_v2 <- rep(1, 5) / 5

## tests -----------------------------------------------------------------------

# graphs
g_undir <- igraph::graph_from_data_frame(edge_dt, directed = FALSE)
g_dir <- igraph::graph_from_data_frame(edge_dt, directed = TRUE)
g_weighted <- igraph::graph_from_data_frame(
  edge_dt_weighted,
  directed = FALSE
)

# version 1 - igraph
igraph_res_undir_v1 <- igraph::page_rank(
  graph = g_undir,
  personalized = personalised_v1
)$vector

igraph_res_dir_v1 <- igraph::page_rank(
  graph = g_dir,
  personalized = personalised_v1
)$vector

igraph_res_weighted_v1 <- igraph::page_rank(
  graph = g_weighted,
  personalized = personalised_v1
)$vector

# version 2 - igraph
igraph_res_undir_v2 <- igraph::page_rank(
  graph = g_undir,
  personalized = personalised_v2
)$vector

igraph_res_dir_v2 <- igraph::page_rank(
  graph = g_dir,
  personalized = personalised_v2
)$vector

igraph_res_weighted_v2 <- igraph::page_rank(
  graph = g_weighted,
  personalized = personalised_v2
)$vector

# version 1 - rust
rs_res_undir_v1 <- rs_page_rank(
  node_names = unique_nodes,
  from = edge_dt$from,
  to = edge_dt$to,
  weights = NULL,
  personalised = personalised_v1,
  undirected = TRUE
)

rs_res_dir_v1 <- rs_page_rank(
  node_names = unique_nodes,
  from = edge_dt$from,
  to = edge_dt$to,
  weights = NULL,
  personalised = personalised_v1,
  undirected = FALSE
)

rs_res_weighted_v1 <- rs_page_rank(
  node_names = unique_nodes,
  from = edge_dt_weighted$from,
  to = edge_dt_weighted$to,
  weights = edge_dt_weighted$weight,
  personalised = personalised_v1,
  undirected = TRUE
)

# version 2 - rust
rs_res_undir_v2 <- rs_page_rank(
  node_names = unique_nodes,
  from = edge_dt$from,
  to = edge_dt$to,
  weights = NULL,
  personalised = personalised_v2,
  undirected = TRUE
)

rs_res_dir_v2 <- rs_page_rank(
  node_names = unique_nodes,
  from = edge_dt$from,
  to = edge_dt$to,
  weights = NULL,
  personalised = personalised_v2,
  undirected = FALSE
)

rs_res_weighted_v2 <- rs_page_rank(
  node_names = unique_nodes,
  from = edge_dt_weighted$from,
  to = edge_dt_weighted$to,
  weights = edge_dt_weighted$weight,
  personalised = personalised_v2,
  undirected = TRUE
)

# version 1
cor_undir_v1 <- cor(igraph_res_undir_v1, rs_res_undir_v1)
cor_dir_v1 <- cor(igraph_res_dir_v1, rs_res_dir_v1)
cor_weighted_v1 <- cor(igraph_res_weighted_v1, rs_res_weighted_v1)

expect_true(
  cor_undir_v1 > 0.99,
  info = "Rust personalsied Page Rank implementation undirected network (v1)."
)

expect_true(
  cor_dir_v1 > 0.99,
  info = "Rust personalsied Page Rank implementation directed network (v1)."
)

expect_true(
  cor_weighted_v1 > 0.99,
  info = "Rust personalsied Page Rank implementation weighted network (v1)."
)

# version 2
cor_undir_v2 <- cor(igraph_res_undir_v2, rs_res_undir_v2)
cor_dir_v2 <- cor(igraph_res_dir_v2, rs_res_dir_v2)
cor_weighted_v2 <- cor(igraph_res_weighted_v2, rs_res_weighted_v2)

expect_true(
  cor_undir_v2 > 0.99,
  info = "Rust personalsied Page Rank implementation undirected network (v2)."
)

expect_true(
  cor_dir_v2 > 0.99,
  info = "Rust personalsied Page Rank implementation directed network (v2)."
)

expect_true(
  cor_weighted_v2 > 0.99,
  info = "Rust personalsied Page Rank implementation weighted network (v2)."
)
