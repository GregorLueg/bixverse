use faer::MatRef;
use instant_distance::{Builder, Point as DistancePoint, Search};
use rayon::prelude::*;

use crate::core::methods::coremo::intersection_size_sorted;

////////////////
// Structures //
////////////////

#[derive(Clone, Debug)]
struct Point(Vec<f32>);

impl DistancePoint for Point {
    fn distance(&self, other: &Self) -> f32 {
        self.0
            .iter()
            .zip(other.0.iter())
            .map(|(a, b)| (a - b).powi(2))
            .sum::<f32>()
            .sqrt()
    }
}

/////////////
// Helpers //
/////////////

/// Transform a kNN graph to an edge list
///
/// ### Params
///
/// * `knn_graph` - kNN graph in form of a slice of `Vec<usize>`
///
/// ### Returns
///
/// A tuple of `(Vec<from>, Vec<to>)`
pub fn knn_to_edge_list(knn_graph: &[Vec<usize>]) -> (Vec<usize>, Vec<usize>) {
    let total_edges: usize = knn_graph.par_iter().map(|neighbors| neighbors.len()).sum();

    let mut from = Vec::with_capacity(total_edges);
    let mut to = Vec::with_capacity(total_edges);

    for (i, neighbors) in knn_graph.iter().enumerate() {
        for &neighbor in neighbors {
            from.push(i);
            to.push(neighbor);
        }
    }

    (from, to)
}

/// Transform an sNN graph into an edge list
///
/// ### Params
///
/// * `snn_edges` - sNN edges based on a Vec of tuples that indicate
///   `(from, to, weight)`
///
/// ### Returns
///
/// A tuple of `(Vec<from>, Vec<to>, Vec<weight>)`
pub fn snn_to_edge_list(snn_edges: Vec<(usize, usize, f32)>) -> (Vec<usize>, Vec<usize>, Vec<f32>) {
    let mut from = Vec::with_capacity(snn_edges.len());
    let mut to = Vec::with_capacity(snn_edges.len());
    let mut weights = Vec::with_capacity(snn_edges.len());

    for (f, t, w) in snn_edges {
        from.push(f);
        to.push(t);
        weights.push(w);
    }

    (from, to, weights)
}

////////////////////
// Main functions //
////////////////////

/// Get the kNN graph based on some embedding
///
/// This function generates the kNN graph based on approximate nearest neighbour
/// search based on the HNSW algorithm (hierarchical navigable small world).
///
/// ### Params
///
/// * `mat` - Matrix in which rows represent the samples and columns the
///   respective embeddings for that sample
/// * `no_neighbours` - Number of neighbours for the KNN graph
/// * `seed` - Seed for the HNSW algorithm
///
/// ### Returns
///
/// The k-nearest neighbours based on the HNSW algorithm
pub fn generate_knn(mat: MatRef<f32>, no_neighbours: usize, seed: usize) -> Vec<Vec<usize>> {
    let n_samples = mat.nrows();

    let points: Vec<Point> = (0..n_samples)
        .into_par_iter()
        .map(|i| Point(mat.row(i).iter().cloned().collect()))
        .collect();

    let map = Builder::default()
        .seed(seed as u64)
        .build(points, (0..n_samples).collect::<Vec<_>>());

    (0..n_samples)
        .into_par_iter()
        .map(|i| {
            let query_point = Point(mat.row(i).iter().cloned().collect());
            let mut search = Search::default();
            let mut nearest_neighbours: Vec<usize> = map
                .search(&query_point, &mut search)
                .skip(1)
                .take(no_neighbours)
                .map(|item| *item.value)
                .collect();
            nearest_neighbours.sort_unstable();
            nearest_neighbours
        })
        .collect()
}

/// Generate an sNN graph based on the kNN graph
///
/// ### Params
///
/// * `knn_graph` - K-nearest neighbours data
/// * `no_neighbours` - Number of neighbours in the kNN graph
/// * `pruning` - Below which Jaccard similarity to prune the edge
///
/// ### Returns
///
/// A a vector of tuples `(index_cell_i, index_cell_j, weight)`
pub fn generate_snn(
    knn_graph: &[Vec<usize>],
    no_neighbours: usize,
    pruning: f32,
) -> Vec<(usize, usize, f32)> {
    let n_f32 = no_neighbours as f32;
    (0..knn_graph.len())
        .into_par_iter()
        .flat_map(|i| {
            ((i + 1)..knn_graph.len())
                .into_par_iter()
                .map(move |j| {
                    let weight = intersection_size_sorted(&knn_graph[i], &knn_graph[j]) as f32;
                    (i, j, weight / n_f32)
                })
                .filter(|(_, _, weight)| *weight > pruning)
        })
        .collect()
}
