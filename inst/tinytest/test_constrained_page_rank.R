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
