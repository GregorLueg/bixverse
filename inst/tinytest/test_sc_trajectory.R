# sc trajectory ----------------------------------------------------------------

source("helper_sc.R", local = TRUE)

test_temp_dir <- sc_test_dir("trajectory")

## fixtures --------------------------------------------------------------------

# these mirror the fixtures the algorithm is calibrated against in bixverse-rs.
# genuinely noisy on purpose: a clean lattice makes every neighbourhood
# identical, which is not what the adaptive bandwidths are built for.

linear_manifold <- function(n) {
  t <- (seq_len(n) - 1) / (n - 1)
  m <- cbind(t * 10, runif(n, -0.3, 0.3), runif(n, -0.3, 0.3))
  rownames(m) <- sprintf("cell_%04i", seq_len(n))
  m
}

# the arms step away from the axis immediately. letting them start at y = 0 puts
# them within the noise of each other at the branch point and the kNN bridges
# the two. too fast and they detach from the trunk, which leaves the diffusion
# kernel with a degenerate eigenspace.
y_manifold <- function(trunk, arm, divergence = 1.4) {
  t_trunk <- (seq_len(trunk) - 1) / trunk * 5
  t_arm <- seq_len(arm) / arm * 5
  m <- rbind(
    cbind(t_trunk, runif(trunk, -0.15, 0.15)),
    cbind(5 + t_arm, t_arm * divergence + runif(arm, -0.15, 0.15)),
    cbind(5 + t_arm, -t_arm * divergence + runif(arm, -0.15, 0.15))
  )
  rownames(m) <- sprintf("cell_%04i", seq_len(nrow(m)))
  m
}

knn_of <- function(m) {
  generate_sc_knn(
    data = m,
    neighbours_params = params_sc_neighbours(
      knn = list(k = 15L, ann_dist = "euclidean", knn_method = "exhaustive")
    ),
    seed = 42L,
    .verbose = FALSE
  )
}

# n_eigs is pinned. the eigengap heuristic keeps six or more on fixtures this
# small, and the higher harmonics oscillate across the full scaled range, which
# inflates the arc length until the pruning stops cutting anything.
trajectory_params <- function(...) {
  params_sc_palantir(
    knn = 24L,
    num_waypoints = 250L,
    n_eigs = 3L,
    knn_params = list(knn_method = "exhaustive"),
    ...
  )
}

# tests ------------------------------------------------------------------------

## palantir on a linear trajectory ---------------------------------------------

n_linear <- 120L

set.seed(6L)
linear_data <- linear_manifold(n_linear)
linear_knn <- knn_of(linear_data)

expect_equal(
  current = linear_knn$dist_metric,
  target = "euclidean",
  info = "sc trajectory - the kNN graph carries squared euclidean distances"
)

linear_res <- rs_palantir(
  knn_data = linear_knn,
  palantir_params = trajectory_params(),
  early_cell = 0L,
  terminal_states = NULL,
  seed = 42L,
  verbose = 0L
)

expect_true(
  current = all(linear_res$pseudotime >= 0 & linear_res$pseudotime <= 1),
  info = "sc trajectory - pseudotime is min-max scaled to [0, 1]"
)

# the start cell is not pinned to 0, the refinement can put the minimum on a
# neighbour, but it should still sit at the very beginning of the trajectory
expect_true(
  current = linear_res$pseudotime[linear_res$start_cell + 1L] < 0.05,
  info = "sc trajectory - the start cell sits at the start of the trajectory"
)

expect_true(
  current = cor(linear_res$pseudotime, seq_len(n_linear)) > 0.95,
  info = "sc trajectory - pseudotime tracks the linear manifold"
)

expect_equal(
  current = length(linear_res$terminal_states),
  target = 1L,
  info = "sc trajectory - a linear trajectory has a single terminal state"
)

expect_equal(
  current = linear_res$repair_edges,
  target = 0L,
  info = "sc trajectory - the linear kNN graph needs no connectivity repair"
)

expect_equal(
  current = linear_res$waypoints[1],
  target = linear_res$start_cell,
  info = "sc trajectory - the start cell is the first waypoint"
)

### reproducibility ------------------------------------------------------------

linear_repeat <- rs_palantir(
  knn_data = linear_knn,
  palantir_params = trajectory_params(),
  early_cell = 0L,
  terminal_states = NULL,
  seed = 42L,
  verbose = 0L
)

expect_equal(
  current = linear_repeat$pseudotime,
  target = linear_res$pseudotime,
  info = "sc trajectory - the same seed gives identical pseudotime"
)

## palantir on a branching trajectory ------------------------------------------

trunk_n <- 200L
arm_n <- 200L

set.seed(6L)
y_data <- y_manifold(trunk_n, arm_n)
y_knn <- knn_of(y_data)

# terminal states supplied, so this exercises the Markov chain, the absorbing
# solve, the projection and the entropy without leaning on the detection
# heuristic, which a fixture this small does not determine well enough.
arm_tips <- c(trunk_n + arm_n - 1L, trunk_n + 2L * arm_n - 1L)

y_res <- rs_palantir(
  knn_data = y_knn,
  palantir_params = trajectory_params(),
  early_cell = 0L,
  terminal_states = arm_tips,
  seed = 42L,
  verbose = 0L
)

expect_equal(
  current = y_res$terminal_states,
  target = arm_tips,
  info = "sc trajectory - the provided terminal states are respected"
)

expect_equal(
  current = y_res$repair_edges,
  target = 0L,
  info = "sc trajectory - the arms stay attached to the trunk"
)

expect_equal(
  current = dim(y_res$branch_probs),
  target = c(trunk_n + 2L * arm_n, 2L),
  info = "sc trajectory - fate probabilities are cells x terminal states"
)

expect_true(
  current = all(y_res$branch_probs >= 0 & y_res$branch_probs <= 1),
  info = "sc trajectory - fate probabilities are in [0, 1]"
)

expect_true(
  current = y_res$branch_probs[trunk_n + arm_n, 1] > 0.75 &&
    y_res$branch_probs[trunk_n + 2L * arm_n, 2] > 0.75,
  info = "sc trajectory - each arm tip commits to its own fate"
)

# the root has both fates reachable, so its entropy should approach ln 2
expect_true(
  current = abs(y_res$entropy[1] - log(2)) < 0.05,
  info = "sc trajectory - the root sits close to maximum fate entropy"
)

expect_true(
  current = mean(y_res$entropy[c(
    (trunk_n + arm_n - 19L):(trunk_n + arm_n),
    (trunk_n + 2L * arm_n - 19L):(trunk_n + 2L * arm_n)
  )]) <
    mean(y_res$entropy[1:20]),
  info = "sc trajectory - entropy drops as cells commit to a fate"
)

## paga ------------------------------------------------------------------------

### a chain of clusters --------------------------------------------------------

chain_labels <- factor(
  rep(c("A", "B", "C"), each = n_linear / 3L),
  levels = c("A", "B", "C", "unused")
)

chain_paga <- rs_paga(
  knn_mat = get_knn_mat(linear_knn),
  partitions = as.integer(chain_labels) - 1L,
  n_partitions = nlevels(chain_labels)
)

chain_conn <- as.matrix(sparse_list_to_mat(chain_paga$connectivities))
chain_tree <- as.matrix(sparse_list_to_mat(chain_paga$connectivities_tree))

expect_equal(
  current = dim(chain_conn),
  target = c(4L, 4L),
  info = "sc trajectory - PAGA keeps empty factor levels"
)

expect_equal(
  current = chain_paga$sizes,
  target = c(40L, 40L, 40L, 0L),
  info = "sc trajectory - PAGA cluster sizes match the clustering"
)

expect_true(
  current = all(diag(chain_conn) == 0),
  info = "sc trajectory - the abstracted graph has a zero diagonal"
)

expect_true(
  current = isTRUE(all.equal(
    chain_conn,
    t(chain_conn),
    check.attributes = FALSE
  )),
  info = "sc trajectory - the abstracted graph is symmetric"
)

expect_true(
  current = chain_conn[1, 2] > 0 &&
    chain_conn[2, 3] > 0 &&
    chain_conn[1, 3] == 0,
  info = "sc trajectory - only adjacent clusters of a chain connect"
)

expect_equal(
  current = sum(chain_tree > 0) / 2,
  target = 2,
  info = "sc trajectory - the spanning forest of a chain has n - 1 edges"
)

### a branching topology -------------------------------------------------------

y_labels <- factor(
  c(rep("trunk", trunk_n), rep("arm_1", arm_n), rep("arm_2", arm_n)),
  levels = c("trunk", "arm_1", "arm_2")
)

y_paga <- rs_paga(
  knn_mat = get_knn_mat(y_knn),
  partitions = as.integer(y_labels) - 1L,
  n_partitions = nlevels(y_labels)
)

y_tree <- as.matrix(sparse_list_to_mat(y_paga$connectivities_tree))
dimnames(y_tree) <- list(levels(y_labels), levels(y_labels))

expect_equal(
  current = sum(y_tree > 0) / 2,
  target = 2,
  info = "sc trajectory - the spanning forest of a Y has two edges"
)

expect_true(
  current = y_tree["trunk", "arm_1"] > 0 &&
    y_tree["trunk", "arm_2"] > 0 &&
    y_tree["arm_1", "arm_2"] == 0,
  info = "sc trajectory - the spanning forest recovers the Y topology"
)

## parameter wrappers ----------------------------------------------------------

palantir_defaults <- params_sc_palantir()

expect_equal(
  current = palantir_defaults[c(
    "n_dcs",
    "n_eigs",
    "knn",
    "num_waypoints",
    "scale_components",
    "use_early_cell_as_start",
    "max_iterations",
    "branch_prob_threshold"
  )],
  target = list(
    n_dcs = 10L,
    n_eigs = NULL,
    knn = 30L,
    num_waypoints = 1200L,
    scale_components = TRUE,
    use_early_cell_as_start = TRUE,
    max_iterations = 25L,
    branch_prob_threshold = 0.01
  ),
  info = "sc trajectory - Palantir defaults mirror the Rust defaults"
)

expect_true(
  current = isTRUE(bixverse:::checkScPalantirParams(palantir_defaults)),
  info = "sc trajectory - the default Palantir params pass the check"
)

expect_error(
  current = params_sc_palantir(knn = 3L),
  info = "sc trajectory - knn below the Rust minimum is rejected"
)

expect_error(
  current = params_sc_palantir(max_iterations = 1L),
  info = "sc trajectory - max_iterations below the Rust minimum is rejected"
)

expect_error(
  current = params_sc_palantir(branch_prob_threshold = 1.5),
  info = "sc trajectory - an out-of-range branch threshold is rejected"
)

broken_params <- palantir_defaults
broken_params$knn <- 2L

expect_false(
  current = isTRUE(bixverse:::checkScPalantirParams(broken_params)),
  info = "sc trajectory - the check catches a tampered parameter list"
)

## result classes --------------------------------------------------------------

palantir_res <- new_palantir_res(
  rs_res = y_res,
  used_cells = rownames(y_data),
  modality = "rna"
)

expect_true(
  current = checkmate::testClass(palantir_res, "PalantirRes"),
  info = "sc trajectory - the Palantir result carries its class"
)

expect_equal(
  current = palantir_res$start_cell,
  target = rownames(y_data)[y_res$start_cell + 1L],
  info = "sc trajectory - the start cell is mapped back to a cell name"
)

expect_equal(
  current = palantir_res$terminal_states,
  target = rownames(y_data)[arm_tips + 1L],
  info = "sc trajectory - terminal states are mapped back to cell names"
)

expect_equal(
  current = colnames(palantir_res$branch_probs),
  target = palantir_res$terminal_states,
  info = "sc trajectory - fate probability columns name their terminal state"
)

expect_equal(
  current = palantir_res$pseudotime$cell_id,
  target = rownames(y_data),
  info = "sc trajectory - the pseudotime table is ordered as the kNN graph"
)

paga_res <- new_paga_res(
  rs_res = chain_paga,
  cluster_levels = levels(chain_labels),
  cluster_col = "test_clusters",
  modality = "rna"
)

expect_true(
  current = checkmate::testClass(paga_res, "PagaRes"),
  info = "sc trajectory - the PAGA result carries its class"
)

expect_equal(
  current = rownames(paga_res$connectivities),
  target = levels(chain_labels),
  info = "sc trajectory - the abstracted graph is named by cluster"
)

expect_equal(
  current = paga_res$sizes$n_cells,
  target = c(40L, 40L, 40L, 0L),
  info = "sc trajectory - the cluster sizes survive the wrapping"
)

## methods on a SingleCells object ---------------------------------------------

set.seed(123L)

single_cell_test_data <- sc_test_fixture()

sc_qc_param <- sc_test_qc_params(single_cell_test_data, target_size = 1000)

sc_object <- sc_test_object(
  test_temp_dir,
  single_cell_test_data,
  sc_qc_param = sc_qc_param
)

### early returns and guards ---------------------------------------------------

expect_warning(
  current = run_palantir_sc(sc_object, early_cell = "whatever"),
  info = "sc trajectory - Palantir warns when no kNN graph is present"
)

expect_warning(
  current = run_paga_sc(sc_object, cluster_col = "whatever"),
  info = "sc trajectory - PAGA warns when no kNN graph is present"
)

expect_error(
  current = run_paga_sc(sc_object, cluster_col = "whatever", modality = "adt"),
  info = "sc trajectory - non-RNA modalities need a multi-modal object"
)

### full run -------------------------------------------------------------------

sc_object <- find_hvg_sc(sc_object, hvg_no = 30L, .verbose = FALSE)
sc_object <- calculate_pca_sc(sc_object, no_pcs = 15L, .verbose = FALSE)
sc_object <- find_neighbours_sc(
  object = sc_object,
  neighbours_params = params_sc_neighbours(
    knn = list(k = 15L, ann_dist = "euclidean")
  ),
  .verbose = FALSE
)
sc_object <- find_clusters_sc(sc_object, res = 1.0, name = "clusters")

sc_knn <- get_knn_obj(sc_object)

paga_sc_res <- run_paga_sc(
  sc_object,
  cluster_col = "clusters",
  .verbose = FALSE
)

expect_true(
  current = checkmate::testClass(paga_sc_res, "PagaRes"),
  info = "sc trajectory - run_paga_sc returns a PagaRes"
)

expect_equal(
  current = sum(paga_sc_res$sizes$n_cells),
  target = sc_knn$n,
  info = "sc trajectory - every cell in the kNN graph lands in a cluster"
)

expect_equal(
  current = paga_sc_res$sizes$cluster,
  target = levels(as.factor(unlist(sc_object[["clusters"]]))),
  info = "sc trajectory - the cluster levels come back as names"
)

### paga plot data -------------------------------------------------------------

# "pca" is a virtual embedding name, so this needs no manifold
paga_pd <- extract_paga_plot_data(
  sc_object,
  paga_sc_res,
  embedding = "pca"
)

expect_equal(
  current = names(paga_pd),
  target = c("nodes", "edges"),
  info = "sc paga plot - the extractor returns nodes and edges"
)

expect_equal(
  current = attr(paga_pd, "embedding"),
  target = "pca",
  info = "sc paga plot - the embedding rides along as an attribute"
)

expect_equal(
  current = as.character(paga_pd$nodes$cluster),
  target = paga_sc_res$sizes$cluster,
  info = "sc paga plot - one node per cluster, in graph order"
)

expect_equal(
  current = paga_pd$nodes$n_cells,
  target = paga_sc_res$sizes$n_cells,
  info = "sc paga plot - the node sizes come from the PAGA result"
)

# the node has to sit where its cells sit
paga_embd <- extract_embedding_data(
  sc_object,
  embedding = "pca",
  obs_cols = "clusters"
)
first_cluster <- as.character(paga_pd$nodes$cluster)[1]
manual_centroid <- stats::median(
  paga_embd$dim_1[as.character(paga_embd$clusters) == first_cluster]
)

expect_equal(
  current = paga_pd$nodes$dim_1[1],
  target = manual_centroid,
  info = "sc paga plot - a node sits at the median of its cells"
)

expect_false(
  current = isTRUE(all.equal(
    extract_paga_plot_data(
      sc_object,
      paga_sc_res,
      embedding = "pca",
      centroid = "mean"
    )$nodes$dim_1,
    paga_pd$nodes$dim_1
  )),
  info = "sc paga plot - the mean and the median give different positions"
)

# both graphs are stored symmetrically, so every edge must appear exactly once
edge_keys <- purrr::map2_chr(
  paga_pd$edges$from,
  paga_pd$edges$to,
  \(a, b) paste(sort(c(a, b)), collapse = "|")
)

expect_equal(
  current = length(unique(edge_keys)),
  target = nrow(paga_pd$edges),
  info = "sc paga plot - each undirected edge is emitted once"
)

expect_true(
  current = all(paga_pd$edges$weight >= 0.01),
  info = "sc paga plot - the default threshold is applied"
)

expect_true(
  current = nrow(
    extract_paga_plot_data(
      sc_object,
      paga_sc_res,
      embedding = "pca",
      threshold = 0
    )$edges
  ) >=
    nrow(paga_pd$edges),
  info = "sc paga plot - a lower threshold keeps at least as many edges"
)

# the edge coordinates are what a segment layer draws, so they have to agree
# with the nodes they claim to join
node_idx <- match(paga_pd$edges$from, as.character(paga_pd$nodes$cluster))

expect_equal(
  current = paga_pd$edges$x,
  target = paga_pd$nodes$dim_1[node_idx],
  info = "sc paga plot - the edge start matches its source node"
)

paga_tree_pd <- extract_paga_plot_data(
  sc_object,
  paga_sc_res,
  embedding = "pca",
  tree_only = TRUE
)

expect_true(
  current = nrow(paga_tree_pd$edges) <= nrow(paga_pd$nodes) - 1L,
  info = "sc paga plot - the spanning forest holds at most n - 1 edges"
)

### node statistics ------------------------------------------------------------

sc_object[["paga_stat"]] <- seq_len(sc_knn$n) / sc_knn$n

paga_stat_pd <- extract_paga_plot_data(
  sc_object,
  paga_sc_res,
  embedding = "pca",
  node_stat_col = "paga_stat"
)

expect_true(
  current = "stat" %in% names(paga_stat_pd$nodes),
  info = "sc paga plot - a node statistic column is added on request"
)

expect_equal(
  current = paga_stat_pd$nodes$stat[1],
  target = stats::median(
    (seq_len(sc_knn$n) / sc_knn$n)[
      as.character(paga_embd$clusters) == first_cluster
    ]
  ),
  info = "sc paga plot - the node statistic uses the same summary as the position"
)

### guards ---------------------------------------------------------------------

expect_error(
  current = extract_paga_plot_data(
    sc_object,
    paga_sc_res,
    embedding = "pca",
    cluster_col = "not_the_one"
  ),
  pattern = "PAGA was run on",
  info = "sc paga plot - a mismatched cluster column is refused"
)

# PAGA retains empty factor levels, and an empty cluster has no centroid
paga_empty <- paga_sc_res
n_clusters <- nrow(paga_empty$connectivities)

pad_graph <- function(mat) {
  padded <- matrix(
    0,
    nrow = n_clusters + 1L,
    ncol = n_clusters + 1L,
    dimnames = rep(list(c(rownames(mat), "phantom")), 2L)
  )
  padded[seq_len(n_clusters), seq_len(n_clusters)] <- as.matrix(mat)
  Matrix::Matrix(padded, sparse = TRUE)
}

paga_empty$connectivities <- pad_graph(paga_empty$connectivities)
paga_empty$connectivities_tree <- pad_graph(paga_empty$connectivities_tree)
paga_empty$sizes <- rbind(
  paga_empty$sizes,
  data.table::data.table(cluster = "phantom", n_cells = 0L)
)

paga_empty_pd <- extract_paga_plot_data(
  sc_object,
  paga_empty,
  embedding = "pca"
)

expect_false(
  current = "phantom" %in% as.character(paga_empty_pd$nodes$cluster),
  info = "sc paga plot - a cluster holding no cells is dropped"
)

expect_true(
  current = all(is.finite(c(
    paga_empty_pd$nodes$dim_1,
    paga_empty_pd$nodes$dim_2
  ))),
  info = "sc paga plot - no NaN positions survive the empty cluster drop"
)

# a kNN graph over a subset of the cells no longer aligns with the obs table
sc_object_short <- set_knn(
  sc_object,
  new_sc_knn(
    knn_data = list(
      indices = sc_knn$indices[1:100, ],
      dist = sc_knn$dist[1:100, ],
      dist_metric = sc_knn$dist_metric
    ),
    used_cells = sc_knn$used_cells[1:100]
  )
)

expect_error(
  current = run_paga_sc(
    sc_object_short,
    cluster_col = "clusters",
    .verbose = FALSE
  ),
  pattern = "entries but the kNN graph covers",
  info = "sc trajectory - PAGA refuses a cluster column of the wrong length"
)

expect_error(
  current = run_palantir_sc(
    sc_object,
    early_cell = "not_a_cell",
    .verbose = FALSE
  ),
  pattern = "not found in the kNN graph",
  info = "sc trajectory - Palantir refuses an unknown early cell"
)

palantir_sc_res <- run_palantir_sc(
  sc_object,
  early_cell = sc_knn$used_cells[1],
  palantir_params = params_sc_palantir(
    knn = 24L,
    num_waypoints = 200L,
    n_eigs = 3L,
    use_early_cell_as_start = TRUE,
    knn_params = list(knn_method = "exhaustive")
  ),
  .verbose = FALSE
)

expect_true(
  current = checkmate::testClass(palantir_sc_res, "PalantirRes"),
  info = "sc trajectory - run_palantir_sc returns a PalantirRes"
)

expect_equal(
  current = palantir_sc_res$start_cell,
  target = sc_knn$used_cells[1],
  info = "sc trajectory - the requested early cell is used as start cell"
)

expect_equal(
  current = nrow(palantir_sc_res$pseudotime),
  target = sc_knn$n,
  info = "sc trajectory - one pseudotime value per cell in the kNN graph"
)

expect_true(
  current = all(palantir_sc_res$terminal_states %in% sc_knn$used_cells),
  info = "sc trajectory - terminal states are reported as cell names"
)

### terminal state labels ------------------------------------------------------

labelled_states <- palantir_sc_res$terminal_states
names(labelled_states) <- sprintf("fate_%i", seq_along(labelled_states))

palantir_labelled <- run_palantir_sc(
  sc_object,
  early_cell = sc_knn$used_cells[1],
  terminal_states = labelled_states,
  palantir_params = params_sc_palantir(
    knn = 24L,
    num_waypoints = 200L,
    n_eigs = 3L,
    use_early_cell_as_start = TRUE,
    knn_params = list(knn_method = "exhaustive")
  ),
  .verbose = FALSE
)

expect_equal(
  current = sort(colnames(palantir_labelled$branch_probs)),
  target = sort(names(labelled_states)),
  info = "sc trajectory - named terminal states label the fate columns"
)

expect_equal(
  current = unname(palantir_labelled$terminal_states),
  target = sort(unname(labelled_states)),
  info = "sc trajectory - the labelled run keeps the cell names as values"
)

expect_null(
  current = names(palantir_sc_res$terminal_states),
  info = "sc trajectory - detected terminal states carry no labels"
)

# magic ------------------------------------------------------------------------

## parameter wrappers ----------------------------------------------------------

magic_defaults <- params_sc_magic()

expect_equal(
  current = magic_defaults,
  target = list(
    n_steps = 3L,
    clip_threshold = 0.01,
    gene_batch_size = 1000L,
    layer = "norm",
    allow_large = FALSE
  ),
  info = "sc magic - the defaults mirror the Rust defaults"
)

expect_true(
  current = isTRUE(bixverse:::checkScMagicParams(magic_defaults)),
  info = "sc magic - the default params pass the check"
)

expect_error(
  current = params_sc_magic(n_steps = -1L),
  info = "sc magic - a negative step count is rejected"
)

expect_true(
  current = identical(params_sc_magic(n_steps = 0L)$n_steps, 0L),
  info = "sc magic - zero steps is legal and means no imputation"
)

expect_error(
  current = params_sc_magic(layer = "counts"),
  info = "sc magic - an unknown layer is rejected"
)

broken_magic <- magic_defaults
broken_magic$layer <- "logcounts"

expect_false(
  current = isTRUE(bixverse:::checkScMagicParams(broken_magic)),
  info = "sc magic - the check catches a tampered layer"
)

## running it ------------------------------------------------------------------

magic_genes <- get_gene_names(sc_object)[1:5]

sc_object <- run_magic_sc(
  sc_object,
  features = magic_genes,
  .verbose = FALSE
)

magic_layer <- get_magic(sc_object)

expect_true(
  current = checkmate::testClass(magic_layer, "ScMagic"),
  info = "sc magic - run_magic_sc writes an ScMagic layer"
)

expect_equal(
  current = dim(magic_layer$data),
  target = c(sc_knn$n, length(magic_genes)),
  info = "sc magic - the layer is cells x requested genes"
)

expect_equal(
  current = dimnames(magic_layer$data),
  target = list(get_cell_names(sc_object, filtered = TRUE), magic_genes),
  info = "sc magic - the layer is named by cell and gene"
)

expect_true(
  current = all(magic_layer$data >= 0),
  info = "sc magic - smoothing non-negative counts stays non-negative"
)

### provenance -----------------------------------------------------------------

magic_status <- get_sc_cache_status(sc_object)
magic_row <- magic_status[magic_status$artefact == "magic", ]

expect_equal(
  current = nrow(magic_row),
  target = 1L,
  info = "sc magic - the layer shows up in the cache status"
)

expect_false(
  current = magic_row$stale,
  info = "sc magic - a freshly written layer is not stale"
)

expect_equal(
  current = magic_row$from[[1]],
  target = magic_status$id[magic_status$artefact == "knn"],
  info = "sc magic - the layer records the kNN graph it came from"
)

sc_object_restale <- find_neighbours_sc(
  object = sc_object,
  neighbours_params = params_sc_neighbours(
    knn = list(k = 15L, ann_dist = "euclidean")
  ),
  .verbose = FALSE
)

restale_status <- get_sc_cache_status(sc_object_restale)

expect_true(
  current = restale_status$stale[restale_status$artefact == "magic"],
  info = "sc magic - re-running the neighbours marks the layer stale"
)

### the un-imputed path --------------------------------------------------------

# `n_steps = 0` is the same code path with the operator applied zero times, so
# it has to hand back exactly what the plain extractor reads off the store
sc_object_raw <- run_magic_sc(
  sc_object,
  features = magic_genes,
  magic_params = params_sc_magic(n_steps = 0L, clip_threshold = 0),
  .verbose = FALSE
)

norm_dt <- extract_gene_expression(sc_object, features = magic_genes)

expect_equal(
  current = unname(get_magic(sc_object_raw)$data),
  target = unname(as.matrix(norm_dt[, magic_genes, with = FALSE])),
  info = "sc magic - zero steps reproduces the un-imputed counts"
)

### the layer argument ---------------------------------------------------------

magic_dt <- extract_gene_expression(
  sc_object,
  features = magic_genes,
  layer = "magic"
)

expect_equal(
  current = magic_dt$cell_id,
  target = norm_dt$cell_id,
  info = "sc magic - the imputed layer keeps the cell order of the counts"
)

expect_false(
  current = isTRUE(all.equal(magic_dt$gene_001, norm_dt$gene_001)),
  info = "sc magic - three diffusion steps actually change the values"
)

expect_error(
  current = extract_gene_expression(
    sc_object,
    features = magic_genes,
    layer = "raw"
  ),
  info = "sc magic - an unknown layer is rejected by the extractor"
)

sc_object_no_magic <- remove_magic(sc_object)

expect_error(
  current = extract_gene_expression(
    sc_object_no_magic,
    features = magic_genes,
    layer = "magic"
  ),
  pattern = "No imputed layer found",
  info = "sc magic - reading a layer that was never written errors"
)

# gene trends ------------------------------------------------------------------

## parameter wrappers ----------------------------------------------------------

branch_defaults <- params_sc_branch_selection()

expect_equal(
  current = branch_defaults,
  target = list(q = 0.01, eps = 0.01, resolution = 500L),
  info = "sc gene trends - the branch selection defaults mirror the reference"
)

trend_defaults <- params_sc_gene_trends()

expect_equal(
  current = trend_defaults,
  target = list(
    resolution = 500L,
    weighting = "hard_mask",
    length_scale = 1.0,
    sigma = 1.0,
    jitter = 1e-6,
    max_jitter_retries = 3L,
    chunk_size = 2048L
  ),
  info = "sc gene trends - the trend defaults mirror the Rust defaults"
)

expect_true(
  current = isTRUE(bixverse:::checkScGeneTrendParams(trend_defaults)),
  info = "sc gene trends - the default trend params pass the check"
)

expect_error(
  current = params_sc_gene_trends(length_scale = 0),
  info = "sc gene trends - a zero length scale is rejected"
)

expect_error(
  current = params_sc_gene_trends(weighting = "gam"),
  info = "sc gene trends - an unknown weighting is rejected"
)

expect_error(
  current = params_sc_branch_selection(q = 1.5),
  info = "sc gene trends - an out-of-range quantile is rejected"
)

broken_trends <- trend_defaults
broken_trends$sigma <- -1

expect_false(
  current = isTRUE(bixverse:::checkScGeneTrendParams(broken_trends)),
  info = "sc gene trends - the check catches a tampered sigma"
)

## running it ------------------------------------------------------------------

trend_resolution <- 50L

trends_res <- run_gene_trends_sc(
  sc_object,
  palantir_res = palantir_sc_res,
  features = magic_genes,
  gene_trend_params = params_sc_gene_trends(resolution = trend_resolution),
  .verbose = FALSE
)

expect_true(
  current = checkmate::testClass(trends_res, "GeneTrendsRes"),
  info = "sc gene trends - run_gene_trends_sc returns a GeneTrendsRes"
)

expect_equal(
  current = trends_res$branches,
  target = colnames(palantir_sc_res$branch_probs),
  info = "sc gene trends - the branches follow the fate probability columns"
)

expect_equal(
  current = nrow(trends_res$trends),
  target = as.integer(trend_resolution) *
    length(magic_genes) *
    length(trends_res$branches),
  info = "sc gene trends - one row per grid point, gene and branch"
)

expect_equal(
  current = levels(trends_res$trends$gene),
  target = magic_genes,
  info = "sc gene trends - the gene factor follows the requested order"
)

expect_true(
  current = all(is.finite(trends_res$trends$expression)),
  info = "sc gene trends - the fitted values are finite"
)

expect_true(
  current = all(purrr::map_lgl(
    trends_res$branch_cells,
    \(cells) all(cells %in% sc_knn$used_cells)
  )),
  info = "sc gene trends - the branch cells come back as cell names"
)

expect_equal(
  current = purrr::map_int(trends_res$branch_cells, length),
  target = trends_res$run_info$n_cells,
  info = "sc gene trends - the reported cell counts match the selections"
)

### against the imputed layer --------------------------------------------------

trends_magic <- run_gene_trends_sc(
  sc_object,
  palantir_res = palantir_sc_res,
  use_magic = TRUE,
  gene_trend_params = params_sc_gene_trends(resolution = trend_resolution),
  .verbose = FALSE
)

expect_equal(
  current = levels(trends_magic$trends$gene),
  target = magic_genes,
  info = "sc gene trends - a NULL feature list falls back to the whole layer"
)

expect_true(
  current = trends_magic$params$use_magic,
  info = "sc gene trends - the result records which source it fitted"
)

### weighting and smoothing ----------------------------------------------------

# a different Rust path: every cell enters every branch weighted by its fate
# probability rather than by a hard mask
trends_weighted <- run_gene_trends_sc(
  sc_object,
  palantir_res = palantir_sc_res,
  features = magic_genes,
  gene_trend_params = params_sc_gene_trends(
    resolution = trend_resolution,
    weighting = "fate_probability"
  ),
  .verbose = FALSE
)

expect_true(
  current = all(is.finite(trends_weighted$trends$expression)),
  info = "sc gene trends - the fate probability weighting fits finite values"
)

trends_tight <- run_gene_trends_sc(
  sc_object,
  palantir_res = palantir_sc_res,
  features = magic_genes,
  gene_trend_params = params_sc_gene_trends(
    resolution = trend_resolution,
    length_scale = 0.1
  ),
  .verbose = FALSE
)

expect_false(
  current = isTRUE(all.equal(
    trends_tight$trends$expression,
    trends_res$trends$expression
  )),
  info = "sc gene trends - a shorter length scale changes the posterior"
)

### guards ---------------------------------------------------------------------

expect_error(
  current = run_gene_trends_sc(
    sc_object,
    palantir_res = palantir_sc_res,
    use_magic = TRUE,
    features = c(magic_genes, "not_a_gene"),
    .verbose = FALSE
  ),
  pattern = "not in the imputed layer",
  info = "sc gene trends - genes outside the imputed layer are named"
)

expect_error(
  current = run_gene_trends_sc(
    sc_object,
    palantir_res = palantir_sc_res,
    .verbose = FALSE
  ),
  pattern = "only optional when",
  info = "sc gene trends - a NULL feature list needs the imputed layer"
)

palantir_renamed <- palantir_sc_res
palantir_renamed$pseudotime$cell_id <- sprintf(
  "other_%i",
  seq_len(nrow(palantir_renamed$pseudotime))
)

expect_error(
  current = run_gene_trends_sc(
    sc_object,
    palantir_res = palantir_renamed,
    features = magic_genes,
    .verbose = FALSE
  ),
  pattern = "missing from the Palantir result",
  info = "sc gene trends - a mismatched cell set is caught by name"
)

# clean up ---------------------------------------------------------------------

sc_test_cleanup(test_temp_dir)
