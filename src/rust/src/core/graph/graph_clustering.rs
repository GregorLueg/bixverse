use faer::{Mat, MatRef};
use rand::prelude::*;

use crate::core::graph::knn::*;

/////////////
// Helpers //
/////////////

/// Helper function to get K-means clustering
///
/// ### Params
///
/// * `data` - The data on which to apply k-means clustering
/// * `k` - Number of clusters
/// * `max_iters` - Maximum number of iterations
/// * `seed` - Random seed for reproducibility
///
/// ### Return
///
/// Vector with usizes, indicating cluster membership
fn kmeans(data: &MatRef<f64>, k: usize, max_iters: usize, seed: usize) -> Vec<usize> {
    let n = data.nrows();
    let d = data.ncols();
    let mut labels = vec![0; n];
    let mut centroids = Mat::zeros(k, d);

    let mut rng = StdRng::seed_from_u64(seed as u64);
    centroids
        .as_mut()
        .row_mut(0)
        .copy_from(data.row(rng.random_range(0..n)));

    for i in 1..k {
        let mut distances = vec![f64::INFINITY; n];
        for j in 0..n {
            for c in 0..i {
                let dist: f64 = (0..d)
                    .map(|dim| (data[(j, dim)] - centroids[(c, dim)]).powi(2))
                    .sum::<f64>();
                distances[j] = distances[j].min(dist);
            }
        }
        let idx = distances
            .iter()
            .enumerate()
            .max_by(|a, b| a.1.partial_cmp(b.1).unwrap())
            .unwrap()
            .0;
        centroids.as_mut().row_mut(i).copy_from(data.row(idx));
    }

    // Lloyd's algorithm
    for _ in 0..max_iters {
        // Assign to nearest centroid
        for i in 0..n {
            let mut min_dist = f64::INFINITY;
            for j in 0..k {
                let dist: f64 = (0..d)
                    .map(|dim| (data[(i, dim)] - centroids[(j, dim)]).powi(2))
                    .sum();
                if dist < min_dist {
                    min_dist = dist;
                    labels[i] = j;
                }
            }
        }

        // Update centroids
        let mut counts = vec![0; k];
        let mut new_centroids = Mat::zeros(k, d);
        for i in 0..n {
            counts[labels[i]] += 1;
            for j in 0..d {
                new_centroids[(labels[i], j)] += data[(i, j)];
            }
        }
        for i in 0..k {
            if counts[i] > 0 {
                for j in 0..d {
                    new_centroids[(i, j)] /= counts[i] as f64;
                }
            }
        }
        centroids = new_centroids;
    }

    labels
}

////////////////////
// Main functions //
////////////////////

/// Spectral clustering
///
/// ### Params
///
/// * `similarities` - The matrix of similarities
/// * `k_neighbours` - Number of neighbours to consider
/// * `n_cluster` - Number of clusters to detect.
/// * `max_iters` - Maximum iterations for the k-means clustering.
/// * `seed` - For reproducibility purposes in the centroid initialisation
///
/// ### Returns
///
/// Vector with usizes, indicating cluster membership
pub fn spectral_clustering(
    similarities: &MatRef<f64>,
    k_neighbours: usize,
    n_clusters: usize,
    max_iters: usize,
    seed: usize,
) -> Vec<usize> {
    let adjacency = get_knn_graph_adj(similarities, k_neighbours);

    let laplacian = adjacency_to_laplacian(&adjacency.as_ref(), true);

    let eigendecomp = laplacian.eigen().unwrap();
    let eigenvalues = eigendecomp.S().column_vector();
    let eigenvectors = eigendecomp.U();

    let mut indices: Vec<usize> = (0..eigenvalues.nrows()).collect();
    indices.sort_by(|&a, &b| eigenvalues[a].re.partial_cmp(&eigenvalues[b].re).unwrap());

    let mut features = Mat::zeros(similarities.nrows(), n_clusters);
    for i in 0..similarities.nrows() {
        for j in 0..n_clusters {
            features[(i, j)] = eigenvectors[(i, indices[j])].re;
        }
    }

    for i in 0..features.nrows() {
        let norm: f64 = (0..n_clusters)
            .map(|j| features[(i, j)].powi(2))
            .sum::<f64>()
            .sqrt();
        if norm > 1e-10 {
            for j in 0..n_clusters {
                features[(i, j)] /= norm;
            }
        }
    }

    kmeans(&features.as_ref(), n_clusters, max_iters, seed)
}
