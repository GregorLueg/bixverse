use faer::MatRef;
use rayon::prelude::*;
use rustc_hash::FxHashMap;

///////////////////////
// Cluster stability //
///////////////////////

/// Helper functions to calculate the intersection of sorted usize vectors
///
/// ### Params
///
/// * `a` - (Sorted) slice of usize
/// * `b` - (Sorted) slice of usize
///
/// ### Returns
///
/// Intersection between the two sorted slices.
fn intersection_size_sorted(a: &[usize], b: &[usize]) -> usize {
    let mut count = 0;
    let mut i = 0;
    let mut j = 0;

    while i < a.len() && j < b.len() {
        match a[i].cmp(&b[j]) {
            std::cmp::Ordering::Equal => {
                count += 1;
                i += 1;
                j += 1;
            }
            std::cmp::Ordering::Less => i += 1,
            std::cmp::Ordering::Greater => j += 1,
        }
    }

    count
}

/// Function that assesses the cluster stability.
///
/// ### Params
///
/// * `cluster_matrix` - A matrix with the columns representing the bootstraps,
///   the rows the features and the values which cluster the feature belongs to.
///
/// ### Returns
///
/// Returns tuple of (average Jaccard similarities, standard deviations of the
/// Jaccard similarities).
pub fn cluster_stability(cluster_matrix: &MatRef<i32>) -> Vec<(f64, f64)> {
    let n_features = cluster_matrix.nrows();
    let n_iter = cluster_matrix.ncols();

    // Pre-compute cluster membership maps for all bootstraps
    let bootstrap_cluster_maps: Vec<FxHashMap<i32, Vec<usize>>> = (0..n_iter)
        .into_par_iter()
        .map(|boot_idx| {
            let mut clusters_map: FxHashMap<i32, Vec<usize>> = FxHashMap::default();
            for feature_idx in 0..n_features {
                let cluster_id = cluster_matrix[(feature_idx, boot_idx)];
                clusters_map
                    .entry(cluster_id)
                    .or_default()
                    .push(feature_idx);
            }
            clusters_map
        })
        .collect();

    // Process features in parallel with optimized memory usage
    (0..n_features)
        .into_par_iter()
        .map(|feature_idx| {
            let n_pairs = (n_iter * (n_iter - 1)) / 2;
            let mut jaccard_scores = Vec::with_capacity(n_pairs);

            for i in 0..(n_iter - 1) {
                for j in (i + 1)..n_iter {
                    let cluster_i = cluster_matrix[(feature_idx, i)];
                    let cluster_j = cluster_matrix[(feature_idx, j)];

                    let members_i = bootstrap_cluster_maps[i]
                        .get(&cluster_i)
                        .map(|v| v.as_slice())
                        .unwrap_or(&[]);

                    let members_j = bootstrap_cluster_maps[j]
                        .get(&cluster_j)
                        .map(|v| v.as_slice())
                        .unwrap_or(&[]);

                    let intersection_size = intersection_size_sorted(members_i, members_j);
                    let union_size = members_i.len() + members_j.len() - intersection_size;

                    let jaccard = if union_size == 0 {
                        0.0
                    } else {
                        intersection_size as f64 / union_size as f64
                    };
                    jaccard_scores.push(jaccard);
                }
            }

            let mean_jaccard = jaccard_scores.iter().sum::<f64>() / jaccard_scores.len() as f64;
            let variance = jaccard_scores
                .iter()
                .map(|x| (x - mean_jaccard).powi(2))
                .sum::<f64>()
                / jaccard_scores.len() as f64;
            let std_jaccard = variance.sqrt();

            (mean_jaccard, std_jaccard)
        })
        .collect()
}
