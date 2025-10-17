use petgraph::graph::UnGraph;
use petgraph::visit::EdgeRef;
use rand::prelude::*;

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
    let res_over_two_m = resolution / (2.0 * m);

    let mut adj: Vec<Vec<(usize, f32)>> = vec![Vec::new(); n];
    let mut degrees = vec![0.0f32; n];
    for edge in graph.edge_references() {
        let a = edge.source().index();
        let b = edge.target().index();
        let w = *edge.weight();
        adj[a].push((b, w));
        adj[b].push((a, w));
        degrees[a] += w;
        degrees[b] += w;
    }

    let mut communities: Vec<usize> = (0..n).collect();
    let mut comm_degree_sums = degrees.clone();
    let mut improved = true;
    let mut iteration = 0;
    let mut neighbour_weights: Vec<f32> = vec![0.0; n];
    let mut touched: Vec<usize> = Vec::with_capacity(n);

    while improved && iteration < 100 {
        improved = false;
        let mut node_order: Vec<usize> = (0..n).collect();
        node_order.shuffle(&mut rng);

        for &node in &node_order {
            let current_comm = unsafe { *communities.get_unchecked(node) };
            let k_i = unsafe { *degrees.get_unchecked(node) };

            touched.clear();

            for &(neighbour, weight) in unsafe { adj.get_unchecked(node) } {
                let comm = unsafe { *communities.get_unchecked(neighbour) };
                if unsafe { *neighbour_weights.get_unchecked(comm) } == 0.0 {
                    touched.push(comm);
                }
                unsafe { *neighbour_weights.get_unchecked_mut(comm) += weight };
            }

            let mut best_comm = current_comm;
            let mut best_delta = 0.0f32;

            for &comm in &touched {
                if comm != current_comm {
                    let weight = unsafe { *neighbour_weights.get_unchecked(comm) };
                    let sigma_tot = unsafe { *comm_degree_sums.get_unchecked(comm) };
                    let delta = weight - res_over_two_m * k_i * sigma_tot;

                    if delta > best_delta {
                        best_delta = delta;
                        best_comm = comm;
                    }
                }
            }

            for &comm in &touched {
                unsafe { *neighbour_weights.get_unchecked_mut(comm) = 0.0 };
            }

            if best_comm != current_comm && best_delta > 1e-10 {
                unsafe {
                    *communities.get_unchecked_mut(node) = best_comm;
                    *comm_degree_sums.get_unchecked_mut(current_comm) -= k_i;
                    *comm_degree_sums.get_unchecked_mut(best_comm) += k_i;
                }
                improved = true;
            }
        }
        iteration += 1;
    }

    let mut seen = vec![false; n];
    let mut label = 0;
    let mut comm_map = vec![0; n];
    for &c in &communities {
        if !seen[c] {
            seen[c] = true;
            comm_map[c] = label;
            label += 1;
        }
    }
    communities.iter().map(|&c| comm_map[c]).collect()
}
