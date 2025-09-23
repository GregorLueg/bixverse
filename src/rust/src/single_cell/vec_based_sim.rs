use faer::MatRef;
use instant_distance::{Builder, Point as DistancePoint, Search};
use rayon::prelude::*;

use crate::core::graph::knn::AnnoyIndex;
use crate::core::methods::coremo::intersection_size_sorted;

///////////
// Enums //
///////////

////////////////
// Structures //
////////////////

#[derive(Clone, Debug)]
struct Point(Vec<f32>);

impl DistancePoint for Point {
    /// Distance function. This is Euclidean distance without squaring for
    /// speed gains. Does not change the rank order in KNN generation.
    fn distance(&self, other: &Self) -> f32 {
        let mut sum = 0.0f32;

        for i in 0..self.0.len() {
            let diff = self.0[i] - other.0[i];
            sum += diff * diff;
        }
        sum
    }
}

/////////////
// Helpers //
/////////////

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
pub fn generate_knn_hnsw(mat: MatRef<f32>, no_neighbours: usize, seed: usize) -> Vec<Vec<usize>> {
    let n_samples = mat.nrows();
    let points: Vec<Point> = (0..n_samples)
        .into_par_iter()
        .map(|i| Point(mat.row(i).iter().cloned().collect()))
        .collect();

    let map = Builder::default()
        .seed(seed as u64)
        .build(points.clone(), (0..n_samples).collect::<Vec<_>>());

    let res: Vec<Vec<usize>> = points
        .par_iter()
        .enumerate()
        .map(|(_i, point)| {
            // Fix warning: prefix with _
            let mut search = Search::default();
            let mut nearest_neighbours: Vec<usize> = map
                .search(point, &mut search)
                .take(no_neighbours + 1)
                .map(|item| *item.value)
                .collect();

            nearest_neighbours.retain(|&x| x != _i);
            nearest_neighbours.truncate(no_neighbours);
            nearest_neighbours.sort_unstable();
            nearest_neighbours
        })
        .collect();
    res
}

pub fn generate_knn_annoy(mat: MatRef<f32>, no_neighbours: usize, seed: usize) -> Vec<Vec<usize>> {
    let index = AnnoyIndex::new(mat, 100, seed); // More trees

    (0..mat.nrows())
        .into_par_iter()
        .map(|i| {
            let query_vec: Vec<f32> = mat.row(i).iter().cloned().collect();
            let search_k = Some(no_neighbours * 200); // Higher search budget
            let mut neighbors = index.query(&query_vec, no_neighbours + 1, search_k);
            neighbors.retain(|&x| x != i);
            neighbors.truncate(no_neighbours);
            neighbors.sort_unstable();
            neighbors
        })
        .collect()
}

/// Generate an sNN graph based on the kNN graph
///
/// ### Params
///
/// * `knn_graph` - K-nearest neighbours data
/// * `no_neighbours` - Number of neighbours in the kNN graph
/// * `pruning` - Below which Jaccard similarity to prune the edge. In this case
///   the weight is set to `0`.
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
            knn_graph[i]
                .iter()
                .filter_map(move |&j| {
                    if i < j {
                        let weight =
                            intersection_size_sorted(&knn_graph[i], &knn_graph[j]) as f32 / n_f32;
                        Some((i, j, if weight >= pruning { weight } else { 0.0 }))
                    } else {
                        None
                    }
                })
                .collect::<Vec<_>>()
        })
        .collect()
}

///////////
// Tests //
///////////

#[cfg(test)]
mod tests {
    use super::*;
    use faer::Mat;
    use std::collections::HashSet;

    fn create_clustered_data() -> Mat<f32> {
        // Create 3 well-separated clusters
        let mut data = Vec::new();

        // Cluster 1: around (0, 0)
        for i in 0..50 {
            data.push(vec![
                (i as f32 * 0.1) % 2.0 - 1.0,
                (i as f32 * 0.15) % 2.0 - 1.0,
            ]);
        }

        // Cluster 2: around (10, 10)
        for i in 0..50 {
            data.push(vec![
                10.0 + (i as f32 * 0.1) % 2.0 - 1.0,
                10.0 + (i as f32 * 0.15) % 2.0 - 1.0,
            ]);
        }

        // Cluster 3: around (-10, 10)
        for i in 0..50 {
            data.push(vec![
                -10.0 + (i as f32 * 0.1) % 2.0 - 1.0,
                10.0 + (i as f32 * 0.15) % 2.0 - 1.0,
            ]);
        }

        Mat::from_fn(150, 2, |i, j| data[i][j])
    }

    fn jaccard_similarity(a: &[usize], b: &[usize]) -> f32 {
        let set_a: HashSet<_> = a.iter().collect();
        let set_b: HashSet<_> = b.iter().collect();
        let intersection = set_a.intersection(&set_b).count();
        let union = set_a.union(&set_b).count();
        intersection as f32 / union as f32
    }

    #[test]
    fn test_no_self_neighbors() {
        let data = create_clustered_data();
        let hnsw_result = generate_knn_hnsw(data.as_ref(), 5, 42);
        let annoy_result = generate_knn_annoy(data.as_ref(), 5, 42);

        for (i, neighbors) in hnsw_result.iter().enumerate() {
            assert!(
                !neighbors.contains(&i),
                "HNSW: Node {} found itself in neighbors",
                i
            );
        }

        for (i, neighbors) in annoy_result.iter().enumerate() {
            assert!(
                !neighbors.contains(&i),
                "Annoy: Node {} found itself in neighbors",
                i
            );
        }
    }

    #[test]
    fn test_neighbor_count() {
        let data = create_clustered_data();
        let k = 10;
        let hnsw_result = generate_knn_hnsw(data.as_ref(), k, 42);
        let annoy_result = generate_knn_annoy(data.as_ref(), k, 42);

        for neighbors in &hnsw_result {
            assert_eq!(neighbors.len(), k, "HNSW: Wrong number of neighbors");
        }

        for neighbors in &annoy_result {
            assert_eq!(neighbors.len(), k, "Annoy: Wrong number of neighbors");
        }
    }

    #[test]
    fn test_cluster_structure() {
        let data = create_clustered_data();
        let hnsw_result = generate_knn_hnsw(data.as_ref(), 5, 42);
        let annoy_result = generate_knn_annoy(data.as_ref(), 5, 42);

        // Test that nodes in same cluster find each other
        // Node 0 (cluster 1) should have neighbors mostly in range [0, 49]
        let hnsw_cluster1_neighbors: HashSet<_> =
            hnsw_result[0].iter().filter(|&&n| n < 50).collect();
        let annoy_cluster1_neighbors: HashSet<_> =
            annoy_result[0].iter().filter(|&&n| n < 50).collect();

        assert!(
            hnsw_cluster1_neighbors.len() >= 3,
            "HNSW: Node 0 should find cluster neighbors"
        );
        assert!(
            annoy_cluster1_neighbors.len() >= 3,
            "Annoy: Node 0 should find cluster neighbors"
        );
    }

    #[test]
    fn test_deterministic_results() {
        let data = create_clustered_data();

        let hnsw1 = generate_knn_hnsw(data.as_ref(), 5, 42);
        let hnsw2 = generate_knn_hnsw(data.as_ref(), 5, 42);
        assert_eq!(
            hnsw1, hnsw2,
            "HNSW results should be deterministic with same seed"
        );

        let annoy1 = generate_knn_annoy(data.as_ref(), 5, 42);
        let annoy2 = generate_knn_annoy(data.as_ref(), 5, 42);
        assert_eq!(
            annoy1, annoy2,
            "Annoy results should be deterministic with same seed"
        );
    }

    #[test]
    fn test_algorithm_overlap() {
        let data = create_clustered_data();
        let hnsw_result = generate_knn_hnsw(data.as_ref(), 10, 42);
        let annoy_result = generate_knn_annoy(data.as_ref(), 10, 42);

        let mut total_jaccard = 0.0;
        let mut valid_comparisons = 0;

        for i in 0..hnsw_result.len() {
            let jaccard = jaccard_similarity(&hnsw_result[i], &annoy_result[i]);
            if jaccard > 0.0 {
                total_jaccard += jaccard;
                valid_comparisons += 1;
            }
        }

        let avg_jaccard = total_jaccard / valid_comparisons as f32;
        println!("Average Jaccard similarity: {:.3}", avg_jaccard);

        // With well-separated clusters, both should find similar patterns
        assert!(
            avg_jaccard > 0.2,
            "Algorithms should have reasonable overlap (got {:.3})",
            avg_jaccard
        );
    }

    #[test]
    fn test_snn_weights_nonzero() {
        let data = create_clustered_data();
        let hnsw_knn = generate_knn_hnsw(data.as_ref(), 10, 42);
        let annoy_knn = generate_knn_annoy(data.as_ref(), 10, 42);

        let hnsw_snn = generate_snn(&hnsw_knn, 10, 0.1);
        let annoy_snn = generate_snn(&annoy_knn, 10, 0.1);

        let hnsw_nonzero = hnsw_snn.iter().filter(|(_, _, w)| *w > 0.0).count();
        let annoy_nonzero = annoy_snn.iter().filter(|(_, _, w)| *w > 0.0).count();

        println!("HNSW non-zero sNN edges: {}", hnsw_nonzero);
        println!("Annoy non-zero sNN edges: {}", annoy_nonzero);

        assert!(hnsw_nonzero > 0, "HNSW should produce non-zero sNN weights");
        assert!(
            annoy_nonzero > 0,
            "Annoy should produce non-zero sNN weights"
        );
    }

    #[test]
    fn test_sorted_neighbors() {
        let data = create_clustered_data();
        let hnsw_result = generate_knn_hnsw(data.as_ref(), 5, 42);
        let annoy_result = generate_knn_annoy(data.as_ref(), 5, 42);

        for neighbors in &hnsw_result {
            let mut sorted = neighbors.clone();
            sorted.sort_unstable();
            assert_eq!(*neighbors, sorted, "HNSW neighbors should be sorted");
        }

        for neighbors in &annoy_result {
            let mut sorted = neighbors.clone();
            sorted.sort_unstable();
            assert_eq!(*neighbors, sorted, "Annoy neighbors should be sorted");
        }
    }
}
