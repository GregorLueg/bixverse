use faer::{Mat, MatRef};
use half::f16;
use petgraph::graph::{Graph, NodeIndex, UnGraph};
use rayon::prelude::*;
use std::cmp::Ordering;
use std::collections::BinaryHeap;

use crate::assert_symmetric_mat;
use crate::core::data::sparse_structures::*;
use crate::core::graph::graph_structures::SparseGraph;

/// Enum for the approximate nearest neighbour search
#[derive(Clone, Debug, Copy, PartialEq)]
pub enum AnnDist {
    /// Euclidean distance
    Euclidean,
    /// Cosine distance
    Cosine,
}

/// Parsing the approximate nearest neighbour distance
///
/// ### Params
///
/// * `s` - The string that defines the tied summarisation type
///
/// ### Results
///
/// The `AnnDist` defining the distance metric to use for the approximate
/// neighbour search.
pub fn parse_ann_dist(s: &str) -> Option<AnnDist> {
    match s.to_lowercase().as_str() {
        "euclidean" => Some(AnnDist::Euclidean),
        "cosine" => Some(AnnDist::Cosine),
        _ => None,
    }
}

/// kNN symmetrisation method
#[derive(Clone, Copy)]
pub enum KnnSymmetrisation {
    /// Only intersecting nearest neigbhbours will be considered
    Intersection,
    /// The union of nearest neighbours will be considered
    Union,
}

/// Helper function to parse the SEACell graph generation
///
/// ### Params
///
/// * `s` - Type of graph to build
///
/// ### Returns
///
/// Option of the SeaCellGraphGen
pub fn parse_knn_symmetrisation(s: &str) -> Option<KnnSymmetrisation> {
    match s.to_lowercase().as_str() {
        "intersection" => Some(KnnSymmetrisation::Intersection),
        "union" => Some(KnnSymmetrisation::Union),
        _ => None,
    }
}

///////////////////////
// KNN and Laplacian //
///////////////////////

/// Helper struct for KNN with heap
///
/// ### Fields
///
/// * `index` - Index position of that neighbours
/// * `similiarity` - Similarity value for that Neighbour
#[derive(Debug)]
pub struct SimilarityItem {
    pub index: usize,
    pub similarity: f64,
}

/// Equality Trait for SimilarityItem
impl Eq for SimilarityItem {}

/// `PartialEq` trait for `SimilarityItem`
///
/// Check equality between two `SimilarityItem` items in terms of similarity
///
/// ### Returns
///
/// `true` if they are the same.
impl PartialEq for SimilarityItem {
    fn eq(&self, other: &Self) -> bool {
        self.similarity == other.similarity
    }
}

/// Ord trait `SimilarityItem`
///
/// How to order `SimilarityItem`
impl Ord for SimilarityItem {
    fn cmp(&self, other: &Self) -> Ordering {
        // Reverse for min-heap (we want to keep highest similarities)
        other
            .similarity
            .partial_cmp(&self.similarity)
            .unwrap_or(Ordering::Equal)
    }
}

/// PartialOrd trait `SimilarityItem`
///
/// How to order `SimilarityItem`
impl PartialOrd for SimilarityItem {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

/// Generate a KNN graph adjacency matrix from a similarity matrix
///
/// ### Params
///
/// * `similarities` - The symmetric similarity matrix.
/// * `k` - Number of neighbours to take
///
/// ### Returns
///
/// The KNN adjacency matrix
pub fn get_knn_graph_adj(similarities: &MatRef<f64>, k: usize) -> Mat<f64> {
    assert_symmetric_mat!(similarities);

    let n = similarities.nrows();
    let mut adjacency: Mat<f64> = Mat::zeros(n, n);

    // Parallelize across rows
    let rows: Vec<Vec<(usize, f64)>> = (0..n)
        .into_par_iter()
        .map(|i| {
            let mut heap = BinaryHeap::with_capacity(k + 1);

            // Use min-heap to keep top-k similarities
            for j in 0..n {
                if i != j {
                    let sim = similarities[(i, j)];
                    heap.push(SimilarityItem {
                        index: j,
                        similarity: sim,
                    });

                    if heap.len() > k {
                        heap.pop(); // Remove smallest
                    }
                }
            }

            heap.into_iter()
                .map(|item| (item.index, item.similarity))
                .collect()
        })
        .collect();

    // Fill adjacency matrix
    for (i, neighbors) in rows.iter().enumerate() {
        for &(j, sim) in neighbors {
            adjacency[(i, j)] = sim;
        }
    }

    // Symmetrize in parallel
    for i in 0..n {
        for j in i + 1..n {
            let val = (adjacency[(i, j)] + adjacency[(j, i)]) / 2.0;
            adjacency[(i, j)] = val;
            adjacency[(j, i)] = val;
        }
    }

    adjacency
}

/// Generate a Laplacian matrix from an adjacency matrix
///
/// ### Params
///
/// * `adjacency` - The symmetric adjacency matrix.
///
/// ### Returns
///
/// The Laplacian matrix
pub fn adjacency_to_laplacian(adjacency: &MatRef<f64>, normalise: bool) -> Mat<f64> {
    assert_symmetric_mat!(adjacency);
    let n = adjacency.nrows();

    // Compute degrees
    let degrees: Vec<f64> = (0..n)
        .map(|i| adjacency.row(i).iter().sum::<f64>())
        .collect();

    if !normalise {
        // Unnormalised: L = D - A
        let mut laplacian = adjacency.cloned();
        for i in 0..n {
            laplacian[(i, i)] = degrees[i] - adjacency[(i, i)];
            for j in 0..n {
                if i != j {
                    laplacian[(i, j)] = -adjacency[(i, j)];
                }
            }
        }
        laplacian
    } else {
        // Normalised: L_sym = I - D^(-1/2) * A * D^(-1/2)
        let inv_sqrt_d: Vec<f64> = degrees
            .iter()
            .map(|&d| if d > 1e-10 { 1.0 / d.sqrt() } else { 0.0 })
            .collect();

        let mut laplacian = Mat::zeros(n, n);
        for i in 0..n {
            laplacian[(i, i)] = 1.0;
            for j in 0..n {
                laplacian[(i, j)] -= inv_sqrt_d[i] * adjacency[(i, j)] * inv_sqrt_d[j];
            }
        }
        laplacian
    }
}

///////////////////////////
// KNN label propagation //
///////////////////////////

/// Structure to store KNN graphs and do label propagation
///
/// ### Fields
///
/// * `offsets` - Stores the offsets for node i's neighbours
/// * `neighbours`- Flat array with neighbour indices
/// * `weights` - Normalised edge weights
#[derive(Debug, Clone)]
pub struct KnnLabPropGraph {
    pub offsets: Vec<usize>,
    pub neighbours: Vec<usize>,
    pub weights: Vec<f32>,
}

impl KnnLabPropGraph {
    /// Generate the KnnLabelPropGraph
    ///
    /// ### Params
    ///
    /// * `edges` - edge list in form of [node_1, node_2, node_3, ...] which
    ///   indicates alternating pairs (node_1, node_2), etc in terms of edges.
    /// * `n_nodes` - Number of nodes in the graph
    ///
    /// ### Returns
    ///
    /// Self with the data stored in the structure.
    pub fn from_edge_list(edges: &[usize], n_nodes: usize) -> Self {
        // generate an adjaceny matrix for normalisation; could be faer
        let mut adj: Vec<Vec<(usize, f32)>> = vec![vec![]; n_nodes];

        for chunk in edges.chunks(2) {
            let (u, v) = (chunk[0], chunk[1]);
            adj[u].push((v, 1_f32));
            adj[v].push((u, 1_f32));
        }

        for neighbours in &mut adj {
            let sum = neighbours.len() as f32;
            for (_, w) in neighbours {
                *w /= sum;
            }
        }

        // conversion to CSR for better cache locality and look-ups
        let mut offsets: Vec<usize> = vec![0];
        let mut neighbours: Vec<usize> = Vec::with_capacity(edges.len());
        let mut weights: Vec<f32> = Vec::with_capacity(edges.len());

        for node_neighbours in adj {
            for (neighbour, weight) in node_neighbours {
                neighbours.push(neighbour);
                weights.push(weight);
            }
            offsets.push(neighbours.len())
        }

        Self {
            offsets,
            neighbours,
            weights,
        }
    }

    /// Label spreading algorithm
    ///
    /// Function will spread the labels (one hot encoding for categorical data)
    /// over the graph. The input needs to be of structure:
    ///
    /// class 1 -> `[1.0, 0.0, 0.0, 0.0]`
    ///
    /// class 2 -> `[0.0, 1.0, 0.0, 0.0]`
    ///
    /// unlabelled -> `[0.0, 0.0, 0.0, 0.0]`
    ///
    /// ### Params
    ///
    /// * `labels` - One-hot encoded group membership. All zeroes == unlabelled.
    /// * `mask` - Boolean indicating which samples are unlabelled.
    /// * `alpha` - Controls the spreading. Usually between 0.9 to 0.95. Larger
    ///   values goes further labelling, smaller values are more conversative.
    ///
    /// ### Returns
    ///
    /// A Vec<Vec<f32>> with the probabilities of a given group being of that
    /// class.
    pub fn label_spreading(
        &self,
        labels: &[Vec<f32>],
        mask: &[bool],
        alpha: f32,
        iterations: usize,
        tolerance: f32,
    ) -> Vec<Vec<f32>> {
        let n = labels.len();
        let num_classes = labels[0].len();
        let mut y = labels.to_vec();
        let mut y_new = vec![vec![0.0; num_classes]; n];

        for _ in 0..iterations {
            y_new.par_iter_mut().enumerate().for_each(|(node, y_dist)| {
                let start = self.offsets[node];
                let end = self.offsets[node + 1];

                y_dist.fill(0.0);
                for i in start..end {
                    let neighbor_dist = &y[self.neighbours[i]];
                    for c in 0..num_classes {
                        y_dist[c] += self.weights[i] * neighbor_dist[c];
                    }
                }
            });

            let max_change = y
                .par_iter_mut()
                .enumerate()
                .map(|(i, y_dist)| {
                    let mut max_diff = 0.0_f32;

                    if mask[i] {
                        // unlabeled - pure propagation

                        for c in 0..num_classes {
                            max_diff = max_diff.max((y_new[i][c] - y_dist[c]).abs());
                            y_dist[c] = y_new[i][c];
                        }
                    } else {
                        // labeled - anchor to original

                        for c in 0..num_classes {
                            let new_val = alpha * labels[i][c] + (1.0 - alpha) * y_new[i][c];
                            max_diff = max_diff.max((new_val - y_dist[c]).abs());
                            y_dist[c] = new_val;
                        }
                    }

                    max_diff
                })
                .reduce(|| 0.0, f32::max);

            if max_change < tolerance {
                break;
            }
        }

        y
    }
}

/////////////////////
// kNN to PetGraph //
/////////////////////

/// Convert kNN indices to undirected petgraph
///
/// ### Params
///
/// * `knn` - kNN indices. Excludes self.
///
/// ### Returns
///
/// Undirected pet graph, i.e., `UnGraph`
#[allow(dead_code)]
pub fn knn_to_pet_graph(knn: &[Vec<usize>]) -> UnGraph<(), f32> {
    let n_nodes = knn.len();
    let mut graph = Graph::new_undirected();
    let nodes: Vec<NodeIndex> = (0..n_nodes).map(|_| graph.add_node(())).collect();

    // Collect and normalise edges
    let mut edges: Vec<(usize, usize)> = Vec::new();
    for (i, neighbours) in knn.iter().enumerate() {
        for &j in neighbours {
            let edge = if i < j { (i, j) } else { (j, i) };
            edges.push(edge);
        }
    }

    // Sort and deduplicate
    edges.sort_unstable();
    edges.dedup();

    // Add to graph
    for (i, j) in edges {
        graph.add_edge(nodes[i], nodes[j], 1.0);
    }

    graph
}

/////////////////////////
// kNN to sparse graph //
/////////////////////////

/// Convert kNN indices to undirected SparseGraph
///
/// ### Params
///
/// * `knn` - kNN indices. Excludes self.
///
/// ### Returns
///
/// Undirected sparse graph with symmetric CSR representation
pub fn knn_to_sparse_graph(knn: &[Vec<usize>]) -> SparseGraph {
    let n_nodes = knn.len();
    let mut edges = Vec::new();

    // Collect all edges in both directions
    for (i, neighbours) in knn.iter().enumerate() {
        for &j in neighbours {
            edges.push((i, j));
            edges.push((j, i)); // Symmetric
        }
    }

    // Sort to group duplicates
    edges.sort_unstable();

    // Deduplicate and sum weights
    let mut rows = Vec::new();
    let mut cols = Vec::new();
    let mut vals = Vec::new();

    if !edges.is_empty() {
        let (mut curr_r, mut curr_c) = edges[0];
        let mut weight = 1.0f32;

        for &(r, c) in &edges[1..] {
            if r == curr_r && c == curr_c {
                weight += 1.0;
            } else {
                rows.push(curr_r);
                cols.push(curr_c);
                vals.push(f16::from_f32(weight));
                (curr_r, curr_c) = (r, c);
                weight = 1.0;
            }
        }
        rows.push(curr_r);
        cols.push(curr_c);
        vals.push(f16::from_f32(weight));
    }

    let adjacency = coo_to_csr(&rows, &cols, &vals, (n_nodes, n_nodes));

    SparseGraph::new(n_nodes, adjacency, false)
}
