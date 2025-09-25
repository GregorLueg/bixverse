use faer::MatRef;
use instant_distance::{Builder, Point as DistancePoint, Search};
use rayon::prelude::*;
use rustc_hash::{FxHashMap, FxHashSet};
use std::time::Instant;

use crate::core::graph::annoy::AnnoyIndex;

///////////
// Enums //
///////////

/// SNN similarity method
#[derive(Clone, Copy)]
pub enum SnnSimilarityMethod {
    /// This will calculate the Jaccard similarity as weight
    Intersection,
    /// This will calculate the Rank version as a weight
    Rank,
}

/// Enum for the different methods
pub enum KnnSearch {
    /// Annoy-based
    Annoy,
    /// Hierarchical Navigable Small World
    Hnsw,
}

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

/// Helper function to get the KNN method
///
/// ### Params
///
/// * `s` - Type of KNN algorithm to use
///
/// ### Returns
///
/// Option of the HvgMethod (some not yet implemented)
pub fn get_knn_method(s: &str) -> Option<KnnSearch> {
    match s.to_lowercase().as_str() {
        "annoy" => Some(KnnSearch::Annoy),
        "hnsw" => Some(KnnSearch::Hnsw),
        _ => None,
    }
}

/// Helper function to get the type of sNN similarity
///
/// ### Params
///
/// * `s` - Type of SNN similarity to use
///
/// ### Returns
///
/// Option of the SnnSimilarityMethod
pub fn get_snn_similiarity_method(s: &str) -> Option<SnnSimilarityMethod> {
    match s.to_lowercase().as_str() {
        "jaccard" => Some(SnnSimilarityMethod::Intersection),
        "rank" => Some(SnnSimilarityMethod::Rank),
        _ => None,
    }
}

////////////////////
// Main functions //
////////////////////

/// Get the kNN graph based on HNSW
///
/// This function generates the kNN graph via an approximate nearest neighbour
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
        .map(|(i, point)| {
            let mut search = Search::default();
            let mut nearest_neighbours: Vec<usize> = map
                .search(point, &mut search)
                .take(no_neighbours + 1)
                .map(|item| *item.value)
                .collect();

            nearest_neighbours.retain(|&x| x != i);
            nearest_neighbours.truncate(no_neighbours);
            nearest_neighbours
        })
        .collect();

    res
}

/// Get the kNN graph based on Annoy
///
/// This function generates the kNN graph based via an approximate nearest
/// neighbour search based on the Annoy algorithm (or a version thereof).
///
/// ### Params
///
/// * `mat` - Matrix in which rows represent the samples and columns the
///   respective embeddings for that sample
/// * `no_neighbours` - Number of neighbours for the KNN graph.
/// * `n_trees` - Number of trees to use for the search.
/// * `search_budget` - Search budget per given query.
/// * `seed` - Seed for the HNSW algorithm
///
/// ### Returns
///
/// The k-nearest neighbours based on the HNSW algorithm
pub fn generate_knn_annoy(
    mat: MatRef<f32>,
    no_neighbours: usize,
    n_trees: usize,
    search_budget: usize,
    seed: usize,
) -> Vec<Vec<usize>> {
    let start_index = Instant::now();

    let index = AnnoyIndex::new(mat, n_trees, seed);

    let end_index = start_index.elapsed();

    println!("Alloy index generation : {:.2?}", end_index);

    let res: Vec<Vec<usize>> = (0..mat.nrows())
        .into_par_iter()
        .map(|i| {
            let query_vec: Vec<f32> = mat.row(i).iter().cloned().collect();
            let search_k = Some(no_neighbours * search_budget);
            let mut neighbors = index.query(&query_vec, no_neighbours + 1, search_k);
            neighbors.retain(|&x| x != i);
            neighbors.truncate(no_neighbours);
            neighbors
        })
        .collect();

    res
}

/// Generate an sNN graph based on the kNN graph (full)
///
/// This version will compare all cells against all cells and generate an edge
/// if any neighbours are shared. This yields way denser graphs and is the
/// approach taken in the `bluster` R package to generate the sNN.
///
/// ### Params
///
/// * `knn_graph` - K-nearest neighbours data as a flat vector in column-major.
/// * `no_neighbours` - Number of neighbours in the kNN graph
/// * `pruning` - Below which Jaccard similarity to prune the edge. In this case
///   the weight is set to `0`.
/// * `method` - Which similarity method to use
///
/// ### Returns
///
/// A tuple with `(<edges>, <weights>)`. The edges are stored in a way that the
/// the first edge points goes from the first element to the second, the second
/// edge from the third to the fourth, etc.
pub fn generate_snn_full(
    flat_knn: &[usize],
    k: usize,
    n_samples: usize,
    pruning: f32,
    method: SnnSimilarityMethod,
) -> (Vec<usize>, Vec<f32>) {
    let mut reverse_mappings: Vec<Vec<(usize, usize)>> = vec![Vec::new(); n_samples];

    for i in 0..n_samples {
        reverse_mappings[i].push((i, 0));

        for neighbor_idx in 0..k {
            let neighbor = flat_knn[neighbor_idx * n_samples + i];
            reverse_mappings[neighbor].push((i, neighbor_idx + 1));
        }
    }

    let results: Vec<(usize, usize, f32)> = (0..n_samples)
        .into_par_iter()
        .flat_map(|j| {
            let mut scores = vec![0.0f32; n_samples];
            let mut added = Vec::new();

            for i in 0..=k {
                let cur_neighbor = if i == 0 {
                    j
                } else {
                    flat_knn[(i - 1) * n_samples + j]
                };

                for &(othernode, other_rank) in &reverse_mappings[cur_neighbor] {
                    if othernode < j {
                        match method {
                            SnnSimilarityMethod::Rank => {
                                let combined_rank = (i + other_rank) as f32;
                                if scores[othernode] == 0.0 {
                                    scores[othernode] = combined_rank;
                                    added.push(othernode);
                                } else if combined_rank < scores[othernode] {
                                    scores[othernode] = combined_rank;
                                }
                            }
                            SnnSimilarityMethod::Intersection => {
                                if scores[othernode] == 0.0 {
                                    added.push(othernode);
                                }
                                scores[othernode] += 1.0;
                            }
                        }
                    }
                }
            }

            added
                .into_iter()
                .filter_map(|othernode| {
                    let weight = match method {
                        SnnSimilarityMethod::Rank => {
                            let preliminary = k as f32 - scores[othernode] / 2.0;
                            let raw_weight = preliminary.max(1e-6);
                            raw_weight / k as f32
                        }
                        SnnSimilarityMethod::Intersection => {
                            scores[othernode] / (2.0 * (k as f32 + 1.0) - scores[othernode])
                        }
                    };

                    if weight >= pruning {
                        Some((j, othernode, weight))
                    } else {
                        None
                    }
                })
                .collect::<Vec<_>>()
        })
        .collect();

    let mut edges = Vec::with_capacity(results.len() * 2);
    let mut weights = Vec::with_capacity(results.len());

    for (i, j, weight) in results {
        edges.push(i);
        edges.push(j);
        weights.push(weight);
    }

    (edges, weights)
}

/// Generate an sNN graph based on the kNN graph (limited)
///
/// This version will only compare cells to the neighbouring cells and
/// deduplicate edges in taking the maximum weight between two given cells.
///
/// ### Params
///
/// * `knn_graph` - K-nearest neighbours data as a flat vector in column-major.
/// * `no_neighbours` - Number of neighbours in the kNN graph
/// * `pruning` - Below which Jaccard similarity to prune the edge. In this case
///   the weight is set to `0`.
/// * `method` - Which similarity method to use
///
/// ### Returns
///
/// A tuple with `(<edges>, <weights>)`. The edges are stored in a way that the
/// the first edge points goes from the first element to the second, the second
/// edge from the third to the fourth, etc.
pub fn generate_snn_limited(
    flat_knn: &[usize],
    k: usize,
    n_samples: usize,
    pruning: f32,
    method: SnnSimilarityMethod,
) -> (Vec<usize>, Vec<f32>) {
    // We need to use a hashmap to store unique edges (smaller_idx, larger_idx)
    // -> weight
    let edge_map: FxHashMap<(usize, usize), f32> = (0..n_samples)
        .into_par_iter()
        .flat_map(|i| {
            let mut edges = Vec::new();

            // only consider edges to this cell's k nearest neighbors
            for neighbor_idx in 0..k {
                let j = flat_knn[neighbor_idx * n_samples + i];

                // Calculate sNN similarity between cell i and its neighbor j
                let weight = match method {
                    SnnSimilarityMethod::Intersection => {
                        // Get neighbors of both cells
                        let neighbors_i: FxHashSet<usize> = (0..k)
                            .map(|idx| flat_knn[idx * n_samples + i])
                            .chain(std::iter::once(i)) // include self
                            .collect();

                        let neighbors_j: FxHashSet<usize> = (0..k)
                            .map(|idx| flat_knn[idx * n_samples + j])
                            .chain(std::iter::once(j)) // include self
                            .collect();

                        let intersection_count =
                            neighbors_i.intersection(&neighbors_j).count() as f32;
                        intersection_count / (2.0 * (k as f32 + 1.0) - intersection_count)
                        // Jaccard
                    }
                    SnnSimilarityMethod::Rank => {
                        // build ranks i
                        let mut ranks_i = FxHashMap::default();
                        ranks_i.insert(i, 0); // self at rank 0
                        for (rank, neighbor) in
                            (0..k).map(|idx| flat_knn[idx * n_samples + i]).enumerate()
                        {
                            ranks_i.insert(neighbor, rank + 1);
                        }

                        // build ranks j
                        let mut ranks_j = FxHashMap::default();
                        ranks_j.insert(j, 0); // self at rank 0
                        for (rank, neighbor) in
                            (0..k).map(|idx| flat_knn[idx * n_samples + j]).enumerate()
                        {
                            ranks_j.insert(neighbor, rank + 1);
                        }

                        // find minimum combined rank of shared neighbors
                        let min_combined_rank = ranks_i
                            .keys()
                            .filter(|&neighbor| ranks_j.contains_key(neighbor))
                            .map(|neighbor| ranks_i[neighbor] + ranks_j[neighbor])
                            .min()
                            .unwrap_or(2 * k)
                            as f32;

                        let preliminary = k as f32 - min_combined_rank / 2.0;
                        let raw_weight = preliminary.max(1e-6);

                        raw_weight / k as f32
                    }
                };

                if weight >= pruning {
                    // Store edge with smaller index first to ensure uniqueness
                    let edge_key = if i < j { (i, j) } else { (j, i) };
                    edges.push((edge_key, weight));
                }
            }

            edges
        })
        .collect::<Vec<_>>()
        .into_iter()
        .fold(FxHashMap::default(), |mut acc, (edge_key, weight)| {
            // Keep the maximum weight if we see the same edge multiple times
            acc.entry(edge_key)
                .and_modify(|existing_weight| {
                    if weight > *existing_weight {
                        *existing_weight = weight;
                    }
                })
                .or_insert(weight);
            acc
        });

    let mut edges = Vec::with_capacity(edge_map.len() * 2);
    let mut weights = Vec::with_capacity(edge_map.len());

    for ((i, j), weight) in edge_map {
        edges.push(i);
        edges.push(j);
        weights.push(weight);
    }

    (edges, weights)
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
        let annoy_result = generate_knn_annoy(data.as_ref(), 5, 100, 200, 42);

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
        let annoy_result = generate_knn_annoy(data.as_ref(), k, 100, 200, 42);

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
        let annoy_result = generate_knn_annoy(data.as_ref(), 5, 100, 200, 42);

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

        let annoy1 = generate_knn_annoy(data.as_ref(), 5, 100, 200, 42);
        let annoy2 = generate_knn_annoy(data.as_ref(), 5, 100, 200, 42);
        assert_eq!(
            annoy1, annoy2,
            "Annoy results should be deterministic with same seed"
        );
    }

    #[test]
    fn test_algorithm_overlap() {
        let data = create_clustered_data();
        let hnsw_result = generate_knn_hnsw(data.as_ref(), 10, 42);
        let annoy_result = generate_knn_annoy(data.as_ref(), 10, 100, 200, 42);

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
    fn test_sorted_neighbors() {
        let data = create_clustered_data();
        let hnsw_result = generate_knn_hnsw(data.as_ref(), 5, 42);
        let annoy_result = generate_knn_annoy(data.as_ref(), 5, 100, 200, 42);

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
