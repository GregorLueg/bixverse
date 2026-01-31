use crate::core::data::sparse_structures::*;
use crate::core::graph::graph_structures::SparseGraph;
use half::f16;

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
