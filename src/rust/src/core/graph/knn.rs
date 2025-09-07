use faer::{Mat, MatRef};
use rayon::prelude::*;
use std::cmp::Ordering;
use std::collections::BinaryHeap;

use crate::assert_symmetric_mat;

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
struct SimilarityItem {
    index: usize,
    similarity: f64,
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
pub fn adjacency_to_laplacian(adjacency: &MatRef<f64>) -> Mat<f64> {
    assert_symmetric_mat!(adjacency);

    let n = adjacency.nrows();

    let mut laplacian = adjacency.cloned();

    for i in 0..n {
        let degree = adjacency.row(i).iter().sum::<f64>();
        laplacian[(i, i)] = degree - adjacency[(i, i)];

        for j in 0..n {
            if i != j {
                laplacian[(i, j)] = -adjacency[(i, j)];
            }
        }
    }

    laplacian
}
