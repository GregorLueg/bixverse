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
