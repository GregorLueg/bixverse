use petgraph::graph::UnGraph;
use petgraph::visit::EdgeRef;
use rand::prelude::*;

use crate::core::graph::graph_structures::*;

////////////////////////
// Louvain - PetGraph //
////////////////////////

/// Louvain community detection
///
/// This version works on undirected PetGraphs
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
#[allow(dead_code)]
pub fn louvain_petgraph(
    graph: &UnGraph<(), f32>,
    resolution: f32,
    iters: usize,
    seed: usize,
) -> Vec<usize> {
    let n = graph.node_count();
    if n == 0 {
        return Vec::new();
    }

    let mut rng = StdRng::seed_from_u64(seed as u64);
    let m: f32 = graph.edge_references().map(|e| *e.weight()).sum::<f32>() / 2.0;
    let res_over_two_m = resolution / (2.0 * m);

    let mut adj: Vec<Vec<(u32, f32)>> = vec![Vec::new(); n];
    let mut degrees = vec![0.0f32; n];

    for edge in graph.edge_references() {
        let (a, b) = (edge.source().index(), edge.target().index());
        let w = *edge.weight();
        adj[a].push((b as u32, w));
        adj[b].push((a as u32, w));
        degrees[a] += w;
        degrees[b] += w;
    }

    let mut communities: Vec<u32> = (0..n as u32).collect();
    let mut comm_degree_sums = degrees.clone();
    let mut neighbour_weights: Vec<f32> = vec![0.0; n];
    let mut comm_active = vec![false; n];
    let mut active_comms = Vec::with_capacity(256);
    let mut node_order: Vec<u32> = (0..n as u32).collect();

    for _iteration in 0..iters {
        let mut move_count = 0;
        node_order.shuffle(&mut rng);

        for &node in &node_order {
            let node_idx = node as usize;
            let current_comm = communities[node_idx] as usize;
            let k_i = degrees[node_idx];
            let k_i_scaled = k_i * res_over_two_m;

            for &(neighbour, weight) in &adj[node_idx] {
                let comm = communities[neighbour as usize] as usize;
                if !comm_active[comm] {
                    comm_active[comm] = true;
                    active_comms.push(comm);
                }
                neighbour_weights[comm] += weight;
            }

            let mut best_comm = current_comm;
            let mut best_delta = 0.0f32;

            for &comm in &active_comms {
                if comm != current_comm {
                    let delta = neighbour_weights[comm] - k_i_scaled * comm_degree_sums[comm];
                    if delta > best_delta {
                        best_delta = delta;
                        best_comm = comm;
                    }
                }
            }

            for &comm in &active_comms {
                neighbour_weights[comm] = 0.0;
                comm_active[comm] = false;
            }
            active_comms.clear();

            if best_comm != current_comm && best_delta > 1e-10 {
                communities[node_idx] = best_comm as u32;
                comm_degree_sums[current_comm] -= k_i;
                comm_degree_sums[best_comm] += k_i;
                move_count += 1;
            }
        }

        if move_count == 0 {
            break;
        }
    }

    let mut comm_map = vec![u32::MAX; n];
    let mut label = 0u32;
    communities.iter_mut().for_each(|c| {
        let idx = *c as usize;
        if comm_map[idx] == u32::MAX {
            comm_map[idx] = label;
            label += 1;
        }
        *c = comm_map[idx];
    });

    communities.iter().map(|&c| c as usize).collect()
}

///////////////////////////
// Louvain - SparseGraph //
///////////////////////////

/// Louvain version for
pub fn louvain_sparse_graph(
    graph: &SparseGraph,
    resolution: f32,
    max_iter: usize,
    seed: usize,
) -> Vec<usize> {
    let n = graph.get_node_number();
    if n == 0 {
        return Vec::new();
    }

    let mut rng = StdRng::seed_from_u64(seed as u64);

    let m: f32 = (0..n)
        .map(|i| {
            graph
                .get_neighbours(i)
                .1
                .iter()
                .map(|w| w.to_f32())
                .sum::<f32>()
        })
        .sum::<f32>()
        / 2.0;

    let res_over_two_m = resolution / (2.0 * m);

    let mut degrees = vec![0.0f32; n];
    for i in 0..n {
        degrees[i] = graph.get_neighbours(i).1.iter().map(|w| w.to_f32()).sum();
    }

    let mut communities: Vec<u32> = (0..n as u32).collect();
    let mut comm_degree_sums = degrees.clone();
    let mut neighbour_weights = vec![0.0f32; n];
    let mut comm_active = vec![false; n];
    let mut active_comms = Vec::with_capacity(256);
    let mut node_order: Vec<u32> = (0..n as u32).collect();

    for _ in 0..max_iter {
        let mut move_count = 0;
        node_order.shuffle(&mut rng);

        for &node in &node_order {
            let node_idx = node as usize;
            let current_comm = communities[node_idx] as usize;
            let k_i = degrees[node_idx];
            let k_i_scaled = k_i * res_over_two_m;

            let (neighbours, weights) = graph.get_neighbours(node_idx);

            for (&neighbour, &weight) in neighbours.iter().zip(weights.iter()) {
                let comm = communities[neighbour] as usize;
                if !comm_active[comm] {
                    comm_active[comm] = true;
                    active_comms.push(comm);
                }
                neighbour_weights[comm] += weight.to_f32();
            }

            let mut best_comm = current_comm;
            let mut best_delta = 0.0f32;

            for &comm in &active_comms {
                if comm != current_comm {
                    let delta = neighbour_weights[comm] - k_i_scaled * comm_degree_sums[comm];
                    if delta > best_delta {
                        best_delta = delta;
                        best_comm = comm;
                    }
                }
            }

            for &comm in &active_comms {
                neighbour_weights[comm] = 0.0;
                comm_active[comm] = false;
            }
            active_comms.clear();

            if best_comm != current_comm && best_delta > 1e-10 {
                communities[node_idx] = best_comm as u32;
                comm_degree_sums[current_comm] -= k_i;
                comm_degree_sums[best_comm] += k_i;
                move_count += 1;
            }
        }

        if move_count == 0 {
            break;
        }
    }

    let mut comm_map = vec![u32::MAX; n];
    let mut label = 0u32;
    for c in &mut communities {
        let idx = *c as usize;
        if comm_map[idx] == u32::MAX {
            comm_map[idx] = label;
            label += 1;
        }
        *c = comm_map[idx];
    }

    communities.iter().map(|&c| c as usize).collect()
}
