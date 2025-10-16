use petgraph::graph::{NodeIndex, UnGraph};
use petgraph::visit::EdgeRef;
use rand::prelude::*;
use rustc_hash::{FxHashMap, FxHashSet};

/////////////
// Helpers //
/////////////

/// Track community state efficiently
struct ClusterState {
    communities: Vec<usize>,
    comm_degree_sums: FxHashMap<usize, f32>,
}

impl ClusterState {
    fn new(n: usize, degrees: &[f32]) -> Self {
        let communities: Vec<usize> = (0..n).collect();
        let mut comm_degree_sums = FxHashMap::default();

        for node in 0..n {
            *comm_degree_sums.entry(node).or_insert(0.0) += degrees[node];
        }

        Self {
            communities,
            comm_degree_sums,
        }
    }

    fn move_node(&mut self, node: usize, from_comm: usize, to_comm: usize, degree: f32) {
        self.communities[node] = to_comm;

        *self.comm_degree_sums.get_mut(&from_comm).unwrap() -= degree;
        *self.comm_degree_sums.entry(to_comm).or_insert(0.0) += degree;
    }

    #[inline]
    fn get_degree_sum(&self, comm: usize) -> f32 {
        *self.comm_degree_sums.get(&comm).unwrap_or(&0.0)
    }

    #[inline]
    fn get_community(&self, node: usize) -> usize {
        self.communities[node]
    }
}

/////////////
// Louvain //
/////////////

/// Louvain community detection
///
/// ### Params
///
/// * `graph` - Undirected pet graph
/// * `resolution` - Resolution parameter for the Louvain clustering
/// * `seed` - Seed for reproducibility purposes
///
/// ### Returns
///
/// Vector of communitiies
pub fn louvain_clustering(graph: &UnGraph<(), f32>, resolution: f32, seed: usize) -> Vec<usize> {
    let n = graph.node_count();
    let mut rng = StdRng::seed_from_u64(seed as u64);
    let m = graph.edge_count() as f32;

    let degrees: Vec<f32> = (0..n)
        .map(|i| graph.edges(NodeIndex::new(i)).count() as f32)
        .collect();

    let mut state = ClusterState::new(n, &degrees);
    let mut improved = true;
    let mut iteration = 0;

    // Reuse allocation
    let mut neighbour_comms: FxHashMap<usize, f32> = FxHashMap::default();

    while improved && iteration < 100 {
        improved = false;
        let mut node_order: Vec<usize> = (0..n).collect();
        node_order.shuffle(&mut rng);

        for &node in &node_order {
            let current_comm = state.get_community(node);
            let k_i = degrees[node];

            // Calculate neighbouring communities
            neighbour_comms.clear();
            for edge in graph.edges(NodeIndex::new(node)) {
                let neighbour = edge.target().index();
                if neighbour != node {
                    let comm = state.get_community(neighbour);
                    *neighbour_comms.entry(comm).or_insert(0.0) += 1.0;
                }
            }

            // Find best community
            let mut best_comm = current_comm;
            let mut best_delta = 0.0f32;

            for (&comm, &weight) in &neighbour_comms {
                if comm == current_comm {
                    continue;
                }

                // O(1) lookup instead of O(n) scan!
                let sigma_tot = state.get_degree_sum(comm);
                let delta = (weight - resolution * k_i * sigma_tot / (2.0 * m)) / (2.0 * m);

                if delta > best_delta {
                    best_delta = delta;
                    best_comm = comm;
                }
            }

            if best_comm != current_comm && best_delta > 1e-10 {
                state.move_node(node, current_comm, best_comm, k_i);
                improved = true;
            }
        }

        iteration += 1;
    }

    // Relabel communities to be consecutive
    let unique_comms: FxHashSet<usize> = state.communities.iter().copied().collect();
    let mut comm_map: FxHashMap<usize, usize> = FxHashMap::default();
    for (idx, &comm) in unique_comms.iter().enumerate() {
        comm_map.insert(comm, idx);
    }

    state.communities.iter().map(|&c| comm_map[&c]).collect()
}
