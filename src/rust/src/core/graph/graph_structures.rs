#![allow(dead_code)]

use crate::core::data::sparse_structures::*;
use half::f16;

/////////////////////////////
// Sparse graph structures //
/////////////////////////////

/// Structure representation of a sparse graph
///
/// ### Fields
///
/// * `adjacency` - Sparse CSR representation of the graph with `f16` as
///   weights.
/// * `num_nodes` - Number of nodes represented in the graph.
#[derive(Clone, Debug)]
pub struct SparseGraph {
    /// Adjacency matrix in CSR (symmetric for undirected graphs)
    adjacency: CompressedSparseData<f16>,
    num_nodes: usize,
    directed: bool,
}

impl SparseGraph {
    pub fn new(num_nodes: usize, adjacency: CompressedSparseData<f16>, directed: bool) -> Self {
        Self {
            adjacency,
            num_nodes,
            directed,
        }
    }

    /// Helper function to get the neighbours and weights
    ///
    /// ### Params
    ///
    /// * `node` - Index of the node for which to get the neighbours
    ///
    /// ### Return
    ///
    /// Tuple of `(neighbour_indices, edge_weights)`
    #[inline]
    pub fn get_neighbours(&self, node: usize) -> (&[usize], &[f16]) {
        let start = self.adjacency.indptr[node];
        let end = self.adjacency.indptr[node + 1];
        (
            &self.adjacency.indices[start..end],
            &self.adjacency.data[start..end],
        )
    }

    /// Get the node degree
    ///
    /// ### Params
    ///
    /// * `node` - Index of the node to get the node degree from
    ///
    /// ### Return
    ///
    /// The node degree for this node.
    #[inline]
    pub fn get_node_degree(&self, node: usize) -> usize {
        self.adjacency.indptr[node + 1] - self.adjacency.indptr[node]
    }

    /// Get total weight
    ///
    /// ### Params
    ///
    /// * `undirected` - Is the graph undirected or not.
    ///
    /// ### Returns
    ///
    /// The total edge weight
    pub fn total_weight(&self) -> f32 {
        let mut total = 0.0_f32;
        for i in 0..self.num_nodes {
            let (_, weights) = self.get_neighbours(i);
            total += weights.iter().map(|w| w.to_f32()).sum::<f32>();
        }
        if !self.directed {
            total * 0.5
        } else {
            total
        }
    }

    /// Expose the number of nodes
    ///
    /// ### Returns
    ///
    /// The number of nodes in the graph
    pub fn get_node_number(&self) -> usize {
        self.num_nodes
    }

    /// Expose if graph is directed
    ///
    /// ### Returns
    ///
    /// Boolean indicating if graph is directed
    pub fn is_directed(&self) -> bool {
        self.directed
    }
}

//////////////////////////////////
// Hierachical graph structures //
//////////////////////////////////

/// Represents one level in the hierarchy
///
/// ### Fields
///
/// * `graph` - The sparse graph
/// * `node_map` - Maps fine node -> coarse node
/// * `coarse_weights` - Number of fine nodes that mapped to each coarse node
pub struct CoarseLevel {
    graph: SparseGraph,
    node_map: Vec<usize>,
    coarse_weights: Vec<usize>,
}

/// Multi-level graph hierarchy
///
/// ### Fields
///
/// * `finest` - Original finest level graph
/// * `levels` - Coarse levels
pub struct GraphHierarchy {
    finest: SparseGraph,
    levels: Vec<CoarseLevel>,
}

impl GraphHierarchy {
    /// Build a Graph hierarchy based on an initial graph
    ///
    /// This one assumes undirected graphs!!!
    ///
    /// ### Params
    ///
    /// * `graph` - The initial SparseGraph
    /// * `reduction_ratio` - By how much to reduce per level the number of
    ///   nodes. Typically 0.5.
    /// * `max_levels` - Maximum depth in terms of coarsion.
    ///
    /// ### Returns
    ///
    /// Initialised self.
    pub fn build(graph: SparseGraph, reduction_ratio: f32, max_levels: usize) -> Self {
        assert!(
            !graph.is_directed(),
            "This is implemented for undirected graphs!"
        );
        let mut levels = Vec::new();
        let mut current_graph = graph.clone();

        for _ in 0..max_levels {
            let coarse_level = Self::coarsen_level(&current_graph);

            let reduction = coarse_level.graph.num_nodes as f32 / current_graph.num_nodes as f32;

            if reduction > reduction_ratio || coarse_level.graph.num_nodes < 100 {
                break;
            }

            current_graph = coarse_level.graph.clone();
            levels.push(coarse_level);
        }

        Self {
            finest: graph,
            levels,
        }
    }

    /// Coarsen one level using heavy-edge matching
    ///
    /// Generates a CoarseLevel of a given graph via heavy edge matching.
    ///
    /// ### Params
    ///
    /// * `graph` -
    fn coarsen_level(graph: &SparseGraph) -> CoarseLevel {
        let matching = Self::heavy_edge_matching(graph);
        let (coarse_graph, node_map, coarse_weights) = Self::contract_graph(graph, &matching);

        CoarseLevel {
            graph: coarse_graph,
            node_map,
            coarse_weights,
        }
    }

    /// Match nodes via heavy edges
    ///
    /// ### Params
    ///
    /// * `graph` - SparseGraph for which to identify the matching nodes
    ///
    /// ### Returns
    ///
    /// A vector of nodes to combine.
    fn heavy_edge_matching(graph: &SparseGraph) -> Vec<Option<usize>> {
        let n = graph.num_nodes;
        let mut matching = vec![None; n];
        let mut matched = vec![false; n];

        for i in 0..n {
            if matched[i] {
                continue;
            }

            let (neighbour, weights) = graph.get_neighbours(i);

            // Find heaviest unmatched neighbour
            let mut best_neighbour = None;
            let mut best_weight = f16::from_f32(0.0);

            for (&j, &w) in neighbour.iter().zip(weights.iter()) {
                if !matched[j] && j != i && w > best_weight {
                    best_weight = w;
                    best_neighbour = Some(j)
                }
            }

            if let Some(j) = best_neighbour {
                matching[i] = Some(j);
                matching[j] = Some(i);
                matched[i] = true;
                matched[j] = true;
            }
        }

        matching
    }

    fn contract_graph(
        graph: &SparseGraph,
        matching: &[Option<usize>],
    ) -> (SparseGraph, Vec<usize>, Vec<usize>) {
        let n = graph.num_nodes;
        let mut node_map = vec![0; n];
        let mut coarse_id = 0;

        // assign coarse node IDs
        for i in 0..n {
            if node_map[i] == 0 || i == 0 {
                match matching[i] {
                    Some(j) if j > i => {
                        node_map[i] = coarse_id;
                        node_map[j] = coarse_id;
                        coarse_id += 1;
                    }
                    None => {
                        // i unchanced
                        node_map[i] = coarse_id;
                        coarse_id += 1;
                    }
                    _ => {} // already assigned
                }
            }
        }

        let num_coarse = coarse_id;
        let mut coarse_weights = vec![0; num_coarse];

        for i in 0..n {
            coarse_weights[node_map[i]] += 1;
        }

        // Build coarse graph edges
        let mut rows = Vec::new();
        let mut cols = Vec::new();
        let mut vals = Vec::new();

        for i in 0..n {
            let ci = node_map[i];
            let (neighbours, weights) = graph.get_neighbours(i);

            for (&j, &w) in neighbours.iter().zip(weights.iter()) {
                let cj = node_map[j];
                if ci != cj {
                    rows.push(ci);
                    cols.push(cj);
                    vals.push(w.to_f32());
                }
            }
        }

        // Sort to group duplicates
        let mut edges: Vec<(usize, usize, f32)> = rows
            .into_iter()
            .zip(cols)
            .zip(vals)
            .map(|((r, c), v)| (r, c, v))
            .collect();
        edges.sort_unstable_by_key(|(r, c, _)| (*r, *c));

        // Sum duplicates
        let mut deduped_rows = Vec::new();
        let mut deduped_cols = Vec::new();
        let mut deduped_vals = Vec::new();

        if !edges.is_empty() {
            let (mut curr_r, mut curr_c, mut curr_v) = edges[0];

            for &(r, c, v) in &edges[1..] {
                if r == curr_r && c == curr_c {
                    curr_v += v;
                } else {
                    deduped_rows.push(curr_r);
                    deduped_cols.push(curr_c);
                    deduped_vals.push(f16::from_f32(curr_v));
                    (curr_r, curr_c, curr_v) = (r, c, v);
                }
            }
            deduped_rows.push(curr_r);
            deduped_cols.push(curr_c);
            deduped_vals.push(f16::from_f32(curr_v));
        }

        let coarse_csr = coo_to_csr(
            &deduped_rows,
            &deduped_cols,
            &deduped_vals,
            (num_coarse, num_coarse),
        );

        let coarse_graph = SparseGraph::new(num_coarse, coarse_csr, false);

        (coarse_graph, node_map, coarse_weights)
    }

    /// Project coarse solution to fine level
    ///
    /// ### Params
    ///
    /// * `level` - Which coarse level (0 = first coarsening)
    /// * `coarse_labels` - Labels at coarse level
    pub fn prolong(&self, level: usize, coarse_labels: &[usize]) -> Vec<usize> {
        if level >= self.levels.len() {
            return coarse_labels.to_vec();
        }

        let node_map = &self.levels[level].node_map;
        node_map.iter().map(|&ci| coarse_labels[ci]).collect()
    }

    /// Get coarsest graph
    pub fn coarsest(&self) -> &SparseGraph {
        self.levels.last().map(|l| &l.graph).unwrap_or(&self.finest)
    }
}
