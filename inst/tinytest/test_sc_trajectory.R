# sc trajectory ----------------------------------------------------------------

test_temp_dir <- file.path(
  tempdir(),
  "trajectory"
)

dir.create(test_temp_dir, recursive = TRUE, showWarnings = FALSE)
stopifnot("Test directory does not exist" = dir.exists(test_temp_dir))

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

single_cell_test_data <- generate_single_cell_test_data()

sc_qc_param <- params_sc_min_quality(
  min_unique_genes = 45L,
  min_lib_size = 300L,
  min_cells = 500L,
  target_size = 1000
)

sc_object <- SingleCells(dir_data = test_temp_dir)

sc_object <- load_r_data(
  object = sc_object,
  counts = single_cell_test_data$counts,
  obs = single_cell_test_data$obs,
  var = single_cell_test_data$var,
  sc_qc_param = sc_qc_param,
  streaming = 0L,
  .verbose = FALSE
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
