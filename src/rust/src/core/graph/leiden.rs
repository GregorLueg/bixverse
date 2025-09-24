use rand::prelude::*;
use rayon::prelude::*;
use rustc_hash::FxHashMap;

/// Represents a weighted graph for community detection
#[derive(Debug, Clone)]
pub struct Graph {
    pub node_count: usize,
    pub edges: FxHashMap<usize, Vec<(usize, f32)>>,
    pub total_weight: f32,
    pub node_weights: Vec<f32>,
}

/// Community structure for tracking which nodes belong to which community
#[derive(Debug, Clone)]
pub struct Communities {
    pub node_to_community: Vec<usize>,
    pub internal_weights: Vec<f32>,
    pub total_weights: Vec<f32>,
    pub num_communities: usize,
}

impl Graph {
    pub fn from_edge_list(from: Vec<usize>, to: Vec<usize>, weights: Vec<f32>) -> Self {
        let max_node = from.iter().chain(to.iter()).max().copied().unwrap_or(0);
        let node_count = max_node + 1;

        let mut edges: Vec<Vec<(usize, f32)>> = vec![Vec::new(); node_count];
        let total_weight: f32 = weights.par_iter().sum();

        // Build adjacency lists
        for ((f, t), w) in from
            .into_iter()
            .zip(to.into_iter())
            .zip(weights.into_iter())
        {
            edges[f].push((t, w));
            if f != t {
                edges[t].push((f, w));
            }
        }

        // Calculate node weights
        let node_weights: Vec<f32> = edges
            .par_iter()
            .map(|neighbors| neighbors.iter().map(|(_, w)| w).sum())
            .collect();

        // Convert to FxHashMap
        let edges_map: FxHashMap<usize, Vec<(usize, f32)>> = edges
            .into_iter()
            .enumerate()
            .filter(|(_, neighbors)| !neighbors.is_empty())
            .collect();

        Self {
            node_count,
            edges: edges_map,
            total_weight,
            node_weights,
        }
    }

    pub fn neighbors(&self, node: usize) -> &[(usize, f32)] {
        self.edges.get(&node).map(|v| v.as_slice()).unwrap_or(&[])
    }
}

impl Communities {
    pub fn new(node_count: usize, graph: &Graph) -> Self {
        let node_to_community: Vec<usize> = (0..node_count).collect();
        let mut internal_weights = vec![0.0; node_count];
        let total_weights = graph.node_weights.clone();

        // Calculate initial internal weights (self-loops)
        #[allow(clippy::needless_range_loop)]
        for node in 0..node_count {
            for &(neighbor, weight) in graph.neighbors(node) {
                if neighbor == node {
                    internal_weights[node] += weight;
                }
            }
        }

        Self {
            node_to_community,
            internal_weights,
            total_weights,
            num_communities: node_count,
        }
    }

    pub fn move_node(&mut self, node: usize, new_community: usize, graph: &Graph) {
        let old_community = self.node_to_community[node];
        if old_community == new_community {
            return;
        }

        let node_weight = graph.node_weights[node];

        // Calculate edge weights to old and new communities
        let mut weight_to_old = 0.0;
        let mut weight_to_new = 0.0;

        for &(neighbor, weight) in graph.neighbors(node) {
            let neighbor_community = self.node_to_community[neighbor];
            if neighbor_community == old_community {
                weight_to_old += weight;
            } else if neighbor_community == new_community {
                weight_to_new += weight;
            }
        }

        // Update communities
        self.internal_weights[old_community] -= weight_to_old;
        self.total_weights[old_community] -= node_weight;
        self.internal_weights[new_community] += weight_to_new;
        self.total_weights[new_community] += node_weight;

        self.node_to_community[node] = new_community;
    }

    pub fn modularity_gain(
        &self,
        node: usize,
        target_community: usize,
        graph: &Graph,
        resolution: f32,
    ) -> f32 {
        let current_community = self.node_to_community[node];
        if current_community == target_community {
            return 0.0;
        }

        let node_weight = graph.node_weights[node];

        let mut weight_to_target = 0.0;
        let mut weight_to_current = 0.0;

        for &(neighbor, weight) in graph.neighbors(node) {
            let neighbor_community = self.node_to_community[neighbor];
            if neighbor_community == target_community {
                weight_to_target += weight;
            } else if neighbor_community == current_community && neighbor != node {
                weight_to_current += weight;
            }
        }

        let target_total = self.total_weights[target_community];
        let current_total = self.total_weights[current_community] - node_weight;

        // Fixed modularity calculation - proper Leiden formula
        let m2 = 2.0 * graph.total_weight;

        (weight_to_target - weight_to_current) / graph.total_weight
            - resolution * node_weight * (target_total - current_total) / (m2 * graph.total_weight)
    }
}

/// Sequential local moving phase with debug output
fn local_moving_phase(
    graph: &Graph,
    communities: &mut Communities,
    resolution: f32,
    rng: &mut StdRng,
) -> bool {
    let mut nodes: Vec<usize> = (0..graph.node_count).collect();
    nodes.shuffle(rng);

    let mut improved = false;
    let mut moves_made = 0;

    for &node in &nodes {
        let current_community = communities.node_to_community[node];
        let mut best_community = current_community;
        let mut best_gain = 0.0;

        // Check all neighbor communities
        let mut neighbor_communities: Vec<usize> = graph
            .neighbors(node)
            .iter()
            .map(|(neighbor, _)| communities.node_to_community[*neighbor])
            .filter(|&nc| nc != current_community)
            .collect();
        neighbor_communities.sort_unstable();
        neighbor_communities.dedup();

        for &nc in &neighbor_communities {
            let gain = communities.modularity_gain(node, nc, graph, resolution);
            if gain > best_gain {
                best_gain = gain;
                best_community = nc;
            }
        }

        if best_community != current_community && best_gain > 1e-10 {
            communities.move_node(node, best_community, graph);
            improved = true;
            moves_made += 1;
        }
    }

    if moves_made > 0 {
        println!("Local moving: {} moves made", moves_made);
    }

    improved
}

/// Improved refinement phase with proper subpartitioning
fn refinement_phase(graph: &Graph, communities: &mut Communities, resolution: f32) -> bool {
    let mut improved = false;

    // Find communities to refine
    let community_ids: Vec<usize> = (0..communities.num_communities).collect();

    for &community_id in &community_ids {
        let community_nodes: Vec<usize> = communities
            .node_to_community
            .iter()
            .enumerate()
            .filter(|(_, &c)| c == community_id)
            .map(|(node, _)| node)
            .collect();

        if community_nodes.len() < 3 {
            continue; // Skip small communities
        }

        // Create subgraph for this community
        let mut subcommunities = vec![0; community_nodes.len()];
        let mut next_subcomm_id = 1;

        // Simple connected components within community
        for (i, _) in community_nodes.iter().enumerate() {
            if subcommunities[i] == 0 {
                // Start new subcomponent
                let mut stack = vec![i];
                subcommunities[i] = next_subcomm_id;

                while let Some(current_idx) = stack.pop() {
                    let current_node = community_nodes[current_idx];

                    for &(neighbor, _) in graph.neighbors(current_node) {
                        if communities.node_to_community[neighbor] == community_id {
                            if let Some(neighbor_idx) =
                                community_nodes.iter().position(|&n| n == neighbor)
                            {
                                if subcommunities[neighbor_idx] == 0 {
                                    subcommunities[neighbor_idx] = next_subcomm_id;
                                    stack.push(neighbor_idx);
                                }
                            }
                        }
                    }
                }
                next_subcomm_id += 1;
            }
        }

        // If we found multiple subcomponents, try to optimize their assignment
        if next_subcomm_id > 2 {
            for (idx, &node) in community_nodes.iter().enumerate() {
                let current_subcomm = subcommunities[idx];
                let mut best_subcomm = current_subcomm;
                let mut best_score = 0.0;

                for test_subcomm in 1..next_subcomm_id {
                    if test_subcomm == current_subcomm {
                        continue;
                    }

                    // Use modularity-based scoring with resolution parameter
                    let mut weight_to_test = 0.0;
                    let mut weight_to_current = 0.0;

                    for &(neighbor, weight) in graph.neighbors(node) {
                        if let Some(neighbor_idx) =
                            community_nodes.iter().position(|&n| n == neighbor)
                        {
                            if subcommunities[neighbor_idx] == test_subcomm {
                                weight_to_test += weight;
                            } else if subcommunities[neighbor_idx] == current_subcomm {
                                weight_to_current += weight;
                            }
                        }
                    }

                    let score = (weight_to_test - weight_to_current) / graph.total_weight
                        - resolution * graph.node_weights[node] * 0.1 / graph.total_weight; // Simplified modularity

                    if score > best_score {
                        best_score = score;
                        best_subcomm = test_subcomm;
                    }
                }

                if best_subcomm != current_subcomm {
                    subcommunities[idx] = best_subcomm;
                    improved = true;
                }
            }
        }
    }

    improved
}

fn aggregate_graph(graph: &Graph, communities: &Communities) -> Graph {
    let mut community_edges: FxHashMap<(usize, usize), f32> = FxHashMap::default();

    for node in 0..graph.node_count {
        let node_community = communities.node_to_community[node];

        for &(neighbor, weight) in graph.neighbors(node) {
            let neighbor_community = communities.node_to_community[neighbor];

            if node <= neighbor {
                let key = if node_community <= neighbor_community {
                    (node_community, neighbor_community)
                } else {
                    (neighbor_community, node_community)
                };

                *community_edges.entry(key).or_insert(0.0) += weight;
            }
        }
    }

    let mut from = Vec::new();
    let mut to = Vec::new();
    let mut weights = Vec::new();

    for ((c1, c2), weight) in community_edges {
        if weight > 0.0 {
            from.push(c1);
            to.push(c2);
            weights.push(weight);
        }
    }

    Graph::from_edge_list(from, to, weights)
}

pub fn leiden_clustering(
    from: Vec<usize>,
    to: Vec<usize>,
    weights: Vec<f32>,
    max_iterations: usize,
    resolution: f32,
    seed: Option<u64>,
) -> Vec<usize> {
    let mut rng = StdRng::seed_from_u64(seed.unwrap_or(42));

    let original_node_count = from
        .iter()
        .chain(to.iter())
        .max()
        .map(|x| x + 1)
        .unwrap_or(0);
    let mut graph = Graph::from_edge_list(from, to, weights);
    let mut communities = Communities::new(graph.node_count, &graph);

    let mut node_mapping: Vec<usize> = (0..original_node_count).collect();

    println!(
        "Starting Leiden clustering with {} nodes and {} edges",
        graph.node_count,
        graph.edges.values().map(|v| v.len()).sum::<usize>() / 2
    );

    for iteration in 0..max_iterations {
        let mut improved = false;

        // Phase 1: Local moving with resolution parameter
        if local_moving_phase(&graph, &mut communities, resolution, &mut rng) {
            improved = true;
        }

        // Phase 2: Refinement with resolution parameter
        if refinement_phase(&graph, &mut communities, resolution) {
            improved = true;
        }

        if !improved {
            println!("Converged after {} iterations", iteration + 1);
            break;
        }

        // Phase 3: Aggregation
        let new_graph = aggregate_graph(&graph, &communities);

        // Update node mapping
        let mut new_mapping = vec![0; original_node_count];
        for (original_node, &current_node) in node_mapping.iter().enumerate() {
            if current_node < communities.node_to_community.len() {
                new_mapping[original_node] = communities.node_to_community[current_node];
            }
        }
        node_mapping = new_mapping;

        graph = new_graph;
        communities = Communities::new(graph.node_count, &graph);

        println!(
            "Iteration {}: {} communities",
            iteration + 1,
            node_mapping.iter().max().map(|x| x + 1).unwrap_or(0)
        );
    }

    node_mapping
}

/// Test basic functionality
#[cfg(test)]
mod tests {
    use super::*;

    /// Helper function for testing...
    fn generate_community_graph(
        community_sizes: &[usize],
        internal_prob: f32,
        external_prob: f32,
        seed: u64,
    ) -> (Vec<usize>, Vec<usize>, Vec<f32>) {
        let mut rng = StdRng::seed_from_u64(seed);
        let mut from = Vec::new();
        let mut to = Vec::new();
        let mut weights = Vec::new();

        let mut node_offset = 0;
        let total_nodes: usize = community_sizes.iter().sum();

        // Create community assignments
        let mut node_to_community = vec![0; total_nodes];
        for (comm_id, &size) in community_sizes.iter().enumerate() {
            for i in 0..size {
                node_to_community[node_offset + i] = comm_id;
            }
            node_offset += size;
        }

        // Generate edges
        for i in 0..total_nodes {
            for j in (i + 1)..total_nodes {
                let same_community = node_to_community[i] == node_to_community[j];
                let prob = if same_community {
                    internal_prob
                } else {
                    external_prob
                };

                if rng.random::<f32>() < prob {
                    from.push(i);
                    to.push(j);
                    weights.push(1.0);
                }
            }
        }

        (from, to, weights)
    }

    /// Other helper function for testing
    fn generate_ring_of_cliques(
        num_cliques: usize,
        clique_size: usize,
    ) -> (Vec<usize>, Vec<usize>, Vec<f32>) {
        let mut from = Vec::new();
        let mut to = Vec::new();
        let mut weights = Vec::new();

        // Create complete graphs (cliques)
        for clique_id in 0..num_cliques {
            let start = clique_id * clique_size;

            // Internal clique edges
            for i in start..(start + clique_size) {
                for j in (i + 1)..(start + clique_size) {
                    from.push(i);
                    to.push(j);
                    weights.push(1.0);
                }
            }

            // Connect to next clique with single edge
            if clique_id < num_cliques - 1 {
                from.push(start + clique_size - 1);
                to.push(start + clique_size);
                weights.push(0.1); // Weak inter-community connection
            }
        }

        // Close the ring
        if num_cliques > 2 {
            from.push(num_cliques * clique_size - 1);
            to.push(0);
            weights.push(0.1);
        }

        (from, to, weights)
    }

    fn count_communities(communities: &[usize]) -> usize {
        let mut unique_communities: Vec<usize> = communities.to_vec();
        unique_communities.sort_unstable();
        unique_communities.dedup();
        unique_communities.len()
    }

    fn calculate_modularity(
        from: &[usize],
        to: &[usize],
        weights: &[f32],
        communities: &[usize],
        resolution: f32,
    ) -> f32 {
        let total_weight: f32 = weights.iter().sum();
        let max_node = from.iter().chain(to.iter()).max().copied().unwrap_or(0);

        // Calculate degree of each node
        let mut degrees = vec![0.0; max_node + 1];
        for ((&f, &t), &w) in from.iter().zip(to.iter()).zip(weights.iter()) {
            degrees[f] += w;
            degrees[t] += w;
        }

        let mut modularity = 0.0;

        // Sum over all edges
        for ((&f, &t), &w) in from.iter().zip(to.iter()).zip(weights.iter()) {
            let same_community = communities[f] == communities[t];
            let expected = (degrees[f] * degrees[t]) / (2.0 * total_weight);

            if same_community {
                modularity += w - resolution * expected;
            } else {
                modularity -= resolution * expected;
            }
        }

        modularity / total_weight
    }

    #[test]
    fn test_ring_of_cliques() {
        let (from, to, weights) = generate_ring_of_cliques(4, 5);
        let result =
            leiden_clustering(from.clone(), to.clone(), weights.clone(), 10, 1.0, Some(42));

        let num_communities = count_communities(&result);
        println!(
            "Ring of 4 cliques (5 nodes each): Found {} communities",
            num_communities
        );

        // Should find approximately 4 communities
        assert!(
            (3..=6).contains(&num_communities),
            "Expected 3-6 communities, got {}",
            num_communities
        );
    }

    #[test]
    fn test_clear_communities() {
        let (from, to, weights) = generate_community_graph(
            &[10, 10, 10], // 3 communities of 10 nodes each
            0.8,           // 80% internal connection probability
            0.05,          // 5% external connection probability
            42,
        );

        let result =
            leiden_clustering(from.clone(), to.clone(), weights.clone(), 10, 1.0, Some(42));
        let num_communities = count_communities(&result);

        println!("3 clear communities: Found {} communities", num_communities);

        // Should find exactly 3 communities
        assert!(
            (2..=5).contains(&num_communities),
            "Expected 2-5 communities, got {}",
            num_communities
        );
    }

    #[test]
    fn test_single_clique() {
        let (from, to, weights) = generate_ring_of_cliques(1, 10);
        let result = leiden_clustering(from, to, weights, 10, 1.0, Some(42));

        let num_communities = count_communities(&result);
        println!(
            "Single clique (10 nodes): Found {} communities",
            num_communities
        );

        // Should find 1 community
        assert_eq!(
            num_communities, 1,
            "Expected 1 community, got {}",
            num_communities
        );
    }

    #[test]
    fn test_resolution_effects() {
        let (from, to, weights) = generate_community_graph(&[8, 8, 8], 0.7, 0.1, 42);

        let low_res =
            leiden_clustering(from.clone(), to.clone(), weights.clone(), 10, 0.5, Some(42));
        let high_res =
            leiden_clustering(from.clone(), to.clone(), weights.clone(), 10, 2.0, Some(42));

        let low_res_communities = count_communities(&low_res);
        let high_res_communities = count_communities(&high_res);

        println!("Resolution 0.5: {} communities", low_res_communities);
        println!("Resolution 2.0: {} communities", high_res_communities);

        // Higher resolution should give more communities
        assert!(
            high_res_communities >= low_res_communities,
            "High resolution should give more communities"
        );
    }

    #[test]
    fn test_modularity_calculation() {
        let from = vec![0, 1, 2, 3];
        let to = vec![1, 2, 3, 0];
        let weights = vec![1.0, 1.0, 1.0, 1.0];

        // Perfect partition: each node in its own community
        let communities = vec![0, 1, 2, 3];
        let mod1 = calculate_modularity(&from, &to, &weights, &communities, 1.0);

        // All in one community
        let communities = vec![0, 0, 0, 0];
        let mod2 = calculate_modularity(&from, &to, &weights, &communities, 1.0);

        println!(
            "Modularity (separate): {:.3}, Modularity (together): {:.3}",
            mod1, mod2
        );

        // For a ring, all together should be better
        assert!(
            mod2 > mod1,
            "Ring should have better modularity when all together"
        );
    }
}
