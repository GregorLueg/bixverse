# graph diffusion tests --------------------------------------------------------

suppressPackageStartupMessages(library(data.table))
library(magrittr)

## data ------------------------------------------------------------------------

### artifical data -------------------------------------------------------------

set.seed(123)

barabasi_graph <- igraph::sample_pa(15)

node_names <- sprintf("node_%i", 1:15)

edge_data <- data.table::setDT(igraph::as_data_frame(barabasi_graph))[, `:=`(
  from = node_names[from],
  to = node_names[to]
)]

genes_1 <- c("node_1", "node_3", "node_10")
genes_2 <- c("node_2", "node_4", "node_6", "node_15")

diffusion_vector_1 <- rep(1, length(genes_1)) %>% `names<-`(genes_1)
diffusion_vector_2 <- rep(1, length(genes_2)) %>% `names<-`(genes_2)

### expected data --------------------------------------------------------------

#### single diffusions ---------------------------------------------------------

# manually calculated
res_single_diff <- c(
  0.233706474,
  0.058360543,
  0.156027594,
  0.051164216,
  0.049662626,
  0.013617091,
  0.005787264,
  0.177805709,
  0.066311728,
  0.100378284,
  0.008697917,
  0.013617091,
  0.005787264,
  0.008697917,
  0.050378284
) %>%
  `names<-`(node_names)

res_single_diff_z <- c(
  1.03753387,
  -0.77109687,
  1.46038021,
  -1.67922749,
  0.03132234,
  -1.17452575,
  -0.92969543,
  1.79531831,
  0.87785097,
  2.36333586,
  -0.96761295,
  -1.01965007,
  -0.82964529,
  -0.92708780,
  0.59262532
) %>%
  `names<-`(node_names)

expected_clusters_v1 <- c(rep("cluster_1", 5), rep("cluster_2", 3))
expected_seed_genes_v1 <- c(rep(2, 5), rep(1, 3))
expected_node_ids_v1 <- c(
  "node_1",
  "node_2",
  "node_3",
  "node_4",
  "node_9",
  "node_10",
  "node_15",
  "node_8"
)

expected_clusters_v2 <- rep("cluster_1", 2)
expected_seed_genes_v2 <- rep(1, 2)
expected_node_ids_v2 <- c("node_10", "node_8")

#### tied diffusion ------------------------------------------------------------

res_tied_diff <- c(
  0.103336196,
  0.058360543,
  0.034377991,
  0.051164216,
  0.021958942,
  0.013617091,
  0.005787264,
  0.103859695,
  0.014610646,
  0.029426913,
  0.008697917,
  0.013617091,
  0.005787264,
  0.008697917,
  0.050378284
) %>%
  `names<-`(node_names)

res_tied_diff_z <- c(
  0.14818944,
  -0.46793517,
  -0.50552495,
  -1.69536239,
  -0.09524208,
  -1.15013896,
  -1.07308734,
  2.21181170,
  -0.48982606,
  1.00214642,
  -1.33395314,
  -1.12305640,
  -1.03930116,
  -1.29465227,
  3.31025175
) %>%
  `names<-`(node_names)

tied_expected_node_ids_v1 <- c(
  "node_1",
  "node_2",
  "node_3",
  "node_4",
  "node_5",
  "node_10",
  "node_15",
  "node_8"
)

tied_expected_clusters_v2 <- c(rep("cluster_1", 3), rep("cluster_2", 2))
tied_expected_node_ids_v2 <- c(
  "node_11",
  "node_14",
  "node_4",
  "node_15",
  "node_8"
)
tied_expected_seed_nodes_1 <- rep(0, 5)
tied_expected_seed_nodes_2 <- rep(1, 5)

## tests -----------------------------------------------------------------------

### single diffusion -----------------------------------------------------------

test_class <- network_diffusions(edge_data, weighted = FALSE, directed = FALSE)

test_class <- diffuse_seed_nodes(test_class, diffusion_vector_1, "max")

expect_warning(
  current = community_detection(
    test_class,
    community_params = params_community_detection(
      min_seed_nodes = 0L,
      min_nodes = 2L,
      threshold_type = "pval_based",
      pval_threshold = 0.1
    )
  ),
  info = "diffusion class - warning when permutations were not calculated"
)

test_class <- permute_seed_nodes(
  test_class,
  perm_iters = 100L,
  .verbose = FALSE
)

single_diff <- S7::props(test_class, "diffusion_res")[[1]][node_names]
single_diff_z_scores <- S7::props(test_class, "diffusion_perm")[[1]][node_names]

expect_equal(
  current = single_diff,
  target = res_single_diff,
  info = paste(
    "single diffusion - diffusion vector"
  ),
  tolerance = 10e-6
)

expect_equal(
  current = single_diff_z_scores,
  target = res_single_diff_z,
  info = paste(
    "single diffusion - diffusion vector (Z scores)"
  ),
  tolerance = 10e-6
)

### community detection single diffusion ---------------------------------------

#### version 1 -----------------------------------------------------------------

test_class <- community_detection(
  test_class,
  community_params = params_community_detection(
    min_seed_nodes = 0L,
    min_nodes = 2L
  )
)

single_com_res_v1 <- get_results(test_class) %>% setorder(cluster_id, node_id)

expect_equal(
  current = single_com_res_v1$cluster_id,
  target = expected_clusters_v1,
  info = paste(
    "single diffusion (v1) - communities, cluster identifiers"
  )
)

expect_equal(
  current = single_com_res_v1$node_id,
  target = expected_node_ids_v1,
  info = paste(
    "single diffusion (v1) - communities, node identifiers"
  )
)

expect_equal(
  current = single_com_res_v1$seed_nodes_no,
  target = expected_seed_genes_v1,
  info = paste(
    "single diffusion (v1) - communities, seed nodes"
  )
)

expect_true(
  current = sum(single_com_res_v1$seed_node) == 3,
  info = paste(
    "single diffusion (v1) - communites, no seed nodes"
  )
)

#### version 2 -----------------------------------------------------------------

test_class <- community_detection(
  test_class,
  community_params = params_community_detection(
    min_seed_nodes = 0L,
    min_nodes = 2L,
    threshold_type = "pval_based",
    pval_threshold = 0.1
  )
)

single_com_res_v2 <- get_results(test_class) %>% setorder(cluster_id, node_id)

expect_equal(
  current = single_com_res_v2$cluster_id,
  target = expected_clusters_v2,
  info = paste(
    "single diffusion (v2) - communities, cluster identifiers"
  )
)

expect_equal(
  current = single_com_res_v2$node_id,
  target = expected_node_ids_v2,
  info = paste(
    "single diffusion (v2) - communities, node identifiers"
  )
)

expect_equal(
  current = single_com_res_v2$seed_nodes_no,
  target = expected_seed_genes_v2,
  info = paste(
    "single diffusion (v2) - communities, seed nodes"
  )
)

expect_true(
  current = sum(single_com_res_v2$seed_node) == 1,
  info = paste(
    "single diffusion (v2) - communites, no seed nodes"
  )
)

### tied diffusion -------------------------------------------------------------

test_class <- network_diffusions(edge_data, weighted = FALSE, directed = FALSE)

test_class <- tied_diffusion(
  test_class,
  diffusion_vector_1 = diffusion_vector_1,
  diffusion_vector_2 = diffusion_vector_2,
  summarisation = "max",
  score_aggregation = "min"
)

test_class <- permute_seed_nodes(
  test_class,
  perm_iters = 100L,
  .verbose = FALSE
)

tied_diff <- S7::props(test_class, "diffusion_res")[[1]][node_names]
tied_diff_z_scores <- S7::props(test_class, "diffusion_perm")[[1]][node_names]

expect_equal(
  current = tied_diff,
  target = res_tied_diff,
  info = paste(
    "tied diffusion - diffusion vector"
  ),
  tolerance = 10e-6
)

expect_equal(
  current = tied_diff_z_scores,
  target = res_tied_diff_z,
  info = paste(
    "tied diffusion - z score vector"
  ),
  tolerance = 10e-6
)

### community detection tied diffusion -----------------------------------------

#### version 1 -----------------------------------------------------------------

test_class <- community_detection(
  test_class,
  community_params = params_community_detection(
    min_seed_nodes = 0L,
    min_nodes = 2L
  )
)

tied_com_res_v1 <- get_results(test_class) %>% setorder(cluster_id, node_id)

expect_equal(
  current = tied_com_res_v1$cluster_id,
  target = expected_clusters_v1,
  info = paste(
    "tied diffusion (v1) - communities, cluster identifiers"
  )
)

expect_equal(
  current = tied_com_res_v1$node_id,
  target = tied_expected_node_ids_v1,
  info = paste(
    "tied diffusion (v1) - communities, node identifiers"
  )
)

expect_equal(
  current = tied_com_res_v1$seed_nodes_1,
  target = expected_seed_genes_v1,
  info = paste(
    "tied diffusion (v1) - communities, seed nodes set 1"
  )
)

expect_equal(
  current = tied_com_res_v1$seed_nodes_2,
  target = expected_seed_genes_v1,
  info = paste(
    "tied diffusion (v1) - communities, seed nodes set 2"
  )
)

expect_true(
  current = sum(tied_com_res_v1$seed_node_a) == 3,
  info = paste(
    "tied diffusion (v1) - communities, no seed nodes set 1"
  )
)

expect_true(
  current = sum(tied_com_res_v1$seed_node_b) == 3,
  info = paste(
    "tied diffusion (v1) - communities, no seed nodes set 2"
  )
)

#### version 2 -----------------------------------------------------------------

test_class <- community_detection(
  test_class,
  community_params = params_community_detection(
    min_seed_nodes = 0L,
    min_nodes = 2L,
    threshold_type = "pval_based",
    pval_threshold = 0.1
  )
)

tied_com_res_v2 <- get_results(test_class) %>% setorder(cluster_id, node_id)


expect_equal(
  current = tied_com_res_v2$cluster_id,
  target = tied_expected_clusters_v2,
  info = paste(
    "tied diffusion (v2) - communities, cluster identifiers"
  )
)

expect_equal(
  current = tied_com_res_v2$node_id,
  target = tied_expected_node_ids_v2,
  info = paste(
    "tied diffusion (v2) - communities, node identifiers"
  )
)

expect_equal(
  current = tied_com_res_v2$seed_nodes_1,
  target = tied_expected_seed_nodes_1,
  info = paste(
    "tied diffusion (v2) - communities, seed nodes set 1"
  )
)

expect_equal(
  current = tied_com_res_v2$seed_nodes_2,
  target = tied_expected_seed_nodes_2,
  info = paste(
    "tied diffusion (v2) - communities, seed nodes set 2"
  )
)

expect_true(
  current = sum(tied_com_res_v2$seed_node_a) == 0,
  info = paste(
    "tied diffusion (v2) - communities, no seed nodes set 1"
  )
)

expect_true(
  current = sum(tied_com_res_v2$seed_node_b) == 2,
  info = paste(
    "tied diffusion (v2) - communities, no seed nodes set 2"
  )
)
