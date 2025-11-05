use petgraph::graph::UnGraph;
use petgraph::visit::EdgeRef;
use rand::prelude::*;
use rayon::prelude::*;
use rustc_hash::FxHashMap;
use std::cmp::Ordering;
use std::collections::BinaryHeap;
use std::time::Instant;

use crate::core::graph::graph_structures::*;

/////////////
// Louvain //
/////////////

//////////////////////
// Petgraph version //
//////////////////////

/// Louvain community detection
///
/// This version works on undirected PetGraphs
///
/// ### Params
///
/// * `graph` - Undirected pet graph
/// * `resolution` - Resolution parameter for the Louvain clustering
/// * `iter` - Numbers of iterations for the algorithm
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

/////////////////
// SparseGraph //
/////////////////

/// Louvain community detection
///
/// This version works on sparse graphs
///
/// ### Params
///
/// * `graph` - Undirected sparse graph
/// * `resolution` - Resolution parameter for the Louvain clustering
/// * `iter` - Numbers of iterations for the algorithm
/// * `seed` - Seed for reproducibility purposes
///
/// ### Returns
///
/// Vector of communitiies
pub fn louvain_sparse_graph(
    graph: &SparseGraph,
    resolution: f32,
    max_iter: usize,
    seed: usize,
) -> Vec<usize> {
    assert!(
        !graph.is_directed(),
        "Louvain does not work for directed graphs!"
    );

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

//////////////
// Walktrap //
//////////////

#[derive(Clone, Copy)]
pub enum Linkage {
    /// Average Linkage
    Average,
    /// Complete linkage
    Complete,
}

/// Parse the linkage distance
///
/// ### Params
///
/// * `s` - The string to parse
///
/// ### Return
///
/// The Option of the `LinkageDist`
pub fn parse_linkage_distance(s: &str) -> Option<Linkage> {
    match s.to_lowercase().as_str() {
        "average" => Some(Linkage::Average),
        "complete" => Some(Linkage::Complete),
        _ => None,
    }
}

/// Wrapper for f32 to use in BinaryHeap (min-heap)
#[derive(Clone, Copy)]
struct OrderedFloat(f32);

impl PartialEq for OrderedFloat {
    /// Partial equality implementation
    ///
    /// ### Params
    ///
    /// * `other` - Other self
    ///
    /// ### Returns
    ///
    /// Boolean if same
    fn eq(&self, other: &Self) -> bool {
        self.0 == other.0
    }
}

impl Eq for OrderedFloat {}

impl PartialOrd for OrderedFloat {
    /// Partial ordering implementation
    ///
    /// ### Params
    ///
    /// * `other` - Other self
    ///
    /// ### Returns
    ///
    /// Option of `Ordering`
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for OrderedFloat {
    //// Compare implementation
    ///
    /// ### Params
    ///
    /// * `other` - Other self
    ///
    /// ### Returns
    ///
    /// `Ordering`
    fn cmp(&self, other: &Self) -> Ordering {
        other.0.partial_cmp(&self.0).unwrap_or(Ordering::Equal)
    }
}

/////////////
// Helpers //
/////////////

/// Compute transition probability matrix
///
/// ### Params
///
/// * `graph` - SparseGraph structure
///
/// ### Returns
///
/// The transition probabilities from a given node to the others with their
/// weights
fn compute_transition_matrix(graph: &SparseGraph) -> Vec<Vec<(usize, f32)>> {
    let n = graph.get_node_number();

    (0..n)
        .into_par_iter()
        .map(|i| {
            let (neighbours, weights) = graph.get_neighbours(i);
            let degree: f32 = weights.iter().map(|w| w.to_f32()).sum();

            if degree > 0.0 {
                neighbours
                    .iter()
                    .zip(weights.iter())
                    .map(|(&j, &w)| (j, w.to_f32() / degree))
                    .collect()
            } else {
                Vec::new()
            }
        })
        .collect()
}

/// Compute distance based on random walk distributions
///
/// ### Params
///
/// * `all_probs` - Probability vectors for all nodes
/// * `i` - Index of node i
/// * `j` - Index of node j
/// * `degrees` - Precomputed degrees for all nodes
///
/// ### Returns
///
/// Distance
#[inline(always)]
fn walk_distance(all_probs: &[Vec<f32>], i: usize, j: usize, degrees: &[f32]) -> f32 {
    let n = all_probs[0].len();
    let mut dist_sq = 0.0;

    for k in 0..n {
        let deg = degrees[k];
        if deg > 1e-10 {
            let diff = all_probs[i][k] - all_probs[j][k];
            dist_sq += (diff * diff) / deg;
        }
    }

    dist_sq.sqrt()
}

/// Compute random walk distances between all pairs within neighbourhood
///
/// ### Params
///
/// * `graph` - The SparseGraph
/// * `transition_probs` - Pre-calculated transition probabilities
/// * `walk_length` - Length of the walk
///
/// ### Returns
///
/// HashMap with the walk distances between the different points
fn compute_walk_distances(
    graph: &SparseGraph,
    transition_probs: &[Vec<(usize, f32)>],
    walk_length: usize,
) -> FxHashMap<(usize, usize), f32> {
    let n = graph.get_node_number();

    // Precompute degrees once
    let degrees: Vec<f32> = (0..n)
        .map(|k| {
            let (_, weights) = graph.get_neighbours(k);
            weights.iter().map(|w| w.to_f32()).sum()
        })
        .collect();

    // Compute all probability distributions
    let all_probs: Vec<Vec<f32>> = (0..n)
        .into_par_iter()
        .map(|i| {
            let mut probs = vec![0.0f32; n];
            let mut new_probs = vec![0.0f32; n];
            probs[i] = 1.0;

            for _ in 0..walk_length {
                new_probs.fill(0.0);
                for node in 0..n {
                    let p = probs[node];
                    if p > 1e-10 {
                        for &(neighbour, trans_p) in &transition_probs[node] {
                            new_probs[neighbour] += p * trans_p;
                        }
                    }
                }
                std::mem::swap(&mut probs, &mut new_probs);
            }
            probs
        })
        .collect();

    // Compute pairwise distances (upper triangle only)
    let estimated_capacity = (n * (n - 1)) / 2;
    let all_distances: Vec<_> = (0..n)
        .into_par_iter()
        .flat_map(|i| {
            (i + 1..n)
                .filter_map(|j| {
                    let dist = walk_distance(&all_probs, i, j, &degrees);
                    if dist < 1e6 {
                        Some(((i, j), dist))
                    } else {
                        None
                    }
                })
                .collect::<Vec<_>>()
        })
        .collect();

    let mut result = FxHashMap::with_capacity_and_hasher(estimated_capacity, Default::default());
    result.extend(all_distances);
    result
}

/// Compute distance between two communities
///
/// ### Params
///
/// * `node_distances` - HashMap with all of the distances
/// * `comm_1` - Indices of members in community 1
/// * `comm_2` - Indices of members in community 2
/// * `linkage` - Linkage type
///
/// ### Returns
///
/// Distance between the members of the two communities
fn compute_community_distance(
    node_distances: &FxHashMap<(usize, usize), f32>,
    comm_1: &[usize],
    comm_2: &[usize],
    linkage: Linkage,
) -> f32 {
    match linkage {
        Linkage::Average => {
            let mut sum = 0.0;
            let mut count = 0;

            for &i in comm_1 {
                for &j in comm_2 {
                    let key = if i < j { (i, j) } else { (j, i) };
                    if let Some(&d) = node_distances.get(&key) {
                        sum += d;
                        count += 1;
                    }
                }
            }

            if count > 0 {
                sum / count as f32
            } else {
                1e10
            }
        }
        Linkage::Complete => {
            let mut max_dist: f32 = 0.0;

            for &i in comm_1 {
                for &j in comm_2 {
                    let key = if i < j { (i, j) } else { (j, i) };
                    if let Some(&d) = node_distances.get(&key) {
                        max_dist = max_dist.max(d);
                    }
                }
            }

            if max_dist > 0.0 {
                max_dist
            } else {
                1e10
            }
        }
    }
}

/// WalkTrap community detection
///
/// ### Params
///
/// * `graph` - The SparseGraph.
/// * `walk_length` - The walk length for the random walkers.
/// * `num_clusters` - Number of communities to return.
/// * `linkage_dist` - The type of Linkage distance to use. One of `"average"`
///   or `"complete"`.
/// * `verbose` - Controls verbosity of the function.
///
/// ### Returns
///
/// The community memberships
pub fn walktrap_sparse_graph(
    graph: &SparseGraph,
    walk_length: usize,
    num_clusters: usize,
    linkage_dist: &str,
    verbose: bool,
) -> Vec<usize> {
    let n = graph.get_node_number();
    if n == 0 {
        return Vec::new();
    }
    if n <= num_clusters {
        return (0..n).collect();
    }

    let walktrap_start = Instant::now();

    // transition probabilities
    let start_transition_probs = Instant::now();

    let transition_probs = compute_transition_matrix(graph);

    let end_transition_probs = start_transition_probs.elapsed();

    if verbose {
        println!(
            "Calculated transition probabilities: {:.2?}",
            end_transition_probs
        );
    }

    // distance calculations
    let start_distance_calc = Instant::now();

    let rw_distances = compute_walk_distances(graph, &transition_probs, walk_length);

    let end_distance_calc = start_distance_calc.elapsed();

    if verbose {
        println!(
            "Calculated Random Walk distances: {:.2?}",
            end_distance_calc
        );
    }

    let linkage_dist = parse_linkage_distance(linkage_dist).unwrap_or(Linkage::Complete);

    let start_linkage = Instant::now();

    let mut community_map: Vec<usize> = (0..n).collect();
    let mut community_sizes: Vec<usize> = vec![1; n];
    let mut active_communities: Vec<bool> = vec![true; n];
    let mut community_generation: Vec<usize> = vec![0; n];
    let mut num_active = n;

    let mut distances: FxHashMap<(usize, usize), f32> =
        FxHashMap::with_capacity_and_hasher(rw_distances.len(), Default::default());

    let mut merge_queue: BinaryHeap<(OrderedFloat, usize, usize, usize, usize)> =
        BinaryHeap::with_capacity(rw_distances.len());

    for (&(i, j), &dist) in &rw_distances {
        if dist.is_finite() {
            distances.insert((i, j), dist);
            merge_queue.push((
                OrderedFloat(dist),
                i,
                j,
                community_generation[i],
                community_generation[j],
            ));
        }
    }

    let mut community_members: Vec<Vec<usize>> = (0..n).map(|i| vec![i]).collect();

    while num_active > num_clusters {
        let (dist, mut c1, mut c2) = loop {
            if let Some((d, a, b, gen_a, gen_b)) = merge_queue.pop() {
                if active_communities[a]
                    && active_communities[b]
                    && community_generation[a] == gen_a
                    && community_generation[b] == gen_b
                {
                    break (d.0, a, b);
                }
            } else {
                break (f32::INFINITY, 0, 0);
            }
        };

        if !dist.is_finite() {
            break;
        }

        if c1 > c2 {
            std::mem::swap(&mut c1, &mut c2);
        }

        let new_size = community_sizes[c1] + community_sizes[c2];

        let mut new_members = std::mem::take(&mut community_members[c1]);
        new_members.extend_from_slice(&community_members[c2]);

        active_communities[c1] = false;
        active_communities[c2] = false;
        num_active -= 1;

        community_sizes.push(new_size);
        active_communities.push(true);
        community_members.push(new_members);
        community_generation.push(0);

        let new_comm = community_sizes.len() - 1;

        for &member in &community_members[new_comm] {
            community_map[member] = new_comm;
        }

        let active_comms: Vec<usize> = (0..active_communities.len())
            .filter(|&i| active_communities[i] && i != new_comm)
            .collect();

        let new_distances: Vec<_> = active_comms
            .par_iter()
            .filter_map(|&other| {
                let dist_new = compute_community_distance(
                    &rw_distances,
                    &community_members[new_comm],
                    &community_members[other],
                    linkage_dist,
                );

                if dist_new.is_finite() {
                    let key = if new_comm < other {
                        (new_comm, other)
                    } else {
                        (other, new_comm)
                    };
                    Some((key, dist_new, new_comm, other))
                } else {
                    None
                }
            })
            .collect();

        for (key, dist_new, comm1, comm2) in new_distances {
            distances.insert(key, dist_new);
            merge_queue.push((
                OrderedFloat(dist_new),
                comm1,
                comm2,
                community_generation[comm1],
                community_generation[comm2],
            ));
        }
    }

    let end_linkage = start_linkage.elapsed();

    if verbose {
        println!("Finished merging the communities: {:.2?}", end_linkage);
    }

    let mut final_labels = vec![0; n];
    let mut label_map: FxHashMap<usize, usize> =
        FxHashMap::with_capacity_and_hasher(num_clusters, Default::default());
    let mut next_label = 0;

    for node in 0..n {
        let comm = community_map[node];
        let label = *label_map.entry(comm).or_insert_with(|| {
            let l = next_label;
            next_label += 1;
            l
        });
        final_labels[node] = label;
    }

    let walktrap_end = walktrap_start.elapsed();

    if verbose {
        println!(
            "Finished WalkTrap communitie detection: {:.2?}",
            walktrap_end
        );
    }

    final_labels
}
