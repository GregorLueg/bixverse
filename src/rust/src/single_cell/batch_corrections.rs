use rayon::prelude::*;
use rustc_hash::FxHashMap;
use statrs::distribution::ChiSquared;
use statrs::distribution::ContinuousCDF;

// use crate::core::base::stats::hypergeom_pval;

//////////
// kBET //
//////////

/// Calculate kBET-based mixing scores on kNN data
///
/// ### Params
///
/// * `knn_data` - KNN data. Outer vector represents the cells, while the inner
///   vector represents
/// * `batches` - Vector indicating the batches.
/// * `pval_threshold`- Value below which we consider the lack of mixing
///   signficant.
///
/// ### Return
///
/// Boolean vector indicating which cells have not the mixing across batches
/// as one would expect from background distributions.
pub fn kbet(knn_data: &Vec<Vec<usize>>, batches: &Vec<usize>, pval_threshold: f64) -> Vec<bool> {
    let mut batch_counts = FxHashMap::default();
    for &batch in batches {
        *batch_counts.entry(batch).or_insert(0) += 1;
    }
    let total = batches.len() as f64;
    let batch_ids: Vec<usize> = batch_counts.keys().copied().collect();
    let dof = (batch_ids.len() - 1) as f64;

    knn_data
        .par_iter()
        .map(|neighbours| {
            let k = neighbours.len() as f64;
            let mut neighbours_count = FxHashMap::default();
            for &neighbour_idx in neighbours {
                *neighbours_count.entry(batches[neighbour_idx]).or_insert(0) += 1;
            }

            // Chi-square test: Σ (observed - expected)² / expected
            let mut chi_square = 0.0;
            for &batch_id in &batch_ids {
                let expected = k * (batch_counts[&batch_id] as f64 / total);
                let observed = *neighbours_count.get(&batch_id).unwrap_or(&0) as f64;
                chi_square += (observed - expected).powi(2) / expected;
            }

            // Compute p-value from chi-square distribution
            let p_value = 1.0 - ChiSquared::new(dof).unwrap().cdf(chi_square);

            // Low p-value = poor mixing (reject null hypothesis of good mixing)
            p_value <= pval_threshold
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_kbet_no_mixing() {
        // Perfect batch separation with k=10 neighbors
        let knn_data = vec![
            vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9], // batch 0: all batch 0 neighbors
            vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
            vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
            vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
            vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
            vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
            vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
            vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
            vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
            vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
            vec![10, 11, 12, 13, 14, 15, 16, 17, 18, 19], // batch 1: all batch 1 neighbors
            vec![10, 11, 12, 13, 14, 15, 16, 17, 18, 19],
            vec![10, 11, 12, 13, 14, 15, 16, 17, 18, 19],
            vec![10, 11, 12, 13, 14, 15, 16, 17, 18, 19],
            vec![10, 11, 12, 13, 14, 15, 16, 17, 18, 19],
            vec![10, 11, 12, 13, 14, 15, 16, 17, 18, 19],
            vec![10, 11, 12, 13, 14, 15, 16, 17, 18, 19],
            vec![10, 11, 12, 13, 14, 15, 16, 17, 18, 19],
            vec![10, 11, 12, 13, 14, 15, 16, 17, 18, 19],
            vec![10, 11, 12, 13, 14, 15, 16, 17, 18, 19],
        ];
        let batches = vec![0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1];

        let result = kbet(&knn_data, &batches, 0.05);
        let poor_mixing_count = result.iter().filter(|&&x| x).count();

        println!(
            "No mixing: {}/{} cells have poor mixing",
            poor_mixing_count,
            result.len()
        );
        assert!(poor_mixing_count > 15); // Should be almost all
    }

    #[test]
    fn test_kbet_perfect_mixing() {
        // Perfect 50/50 mixing with k=10
        let knn_data = vec![
            vec![0, 1, 2, 3, 4, 10, 11, 12, 13, 14], // 5 from each batch
            vec![0, 1, 2, 3, 4, 10, 11, 12, 13, 14],
            vec![0, 1, 2, 3, 4, 10, 11, 12, 13, 14],
            vec![0, 1, 2, 3, 4, 10, 11, 12, 13, 14],
            vec![0, 1, 2, 3, 4, 10, 11, 12, 13, 14],
            vec![0, 1, 2, 3, 4, 10, 11, 12, 13, 14],
            vec![0, 1, 2, 3, 4, 10, 11, 12, 13, 14],
            vec![0, 1, 2, 3, 4, 10, 11, 12, 13, 14],
            vec![0, 1, 2, 3, 4, 10, 11, 12, 13, 14],
            vec![0, 1, 2, 3, 4, 10, 11, 12, 13, 14],
            vec![0, 1, 2, 3, 4, 10, 11, 12, 13, 14],
            vec![0, 1, 2, 3, 4, 10, 11, 12, 13, 14],
            vec![0, 1, 2, 3, 4, 10, 11, 12, 13, 14],
            vec![0, 1, 2, 3, 4, 10, 11, 12, 13, 14],
            vec![0, 1, 2, 3, 4, 10, 11, 12, 13, 14],
            vec![0, 1, 2, 3, 4, 10, 11, 12, 13, 14],
            vec![0, 1, 2, 3, 4, 10, 11, 12, 13, 14],
            vec![0, 1, 2, 3, 4, 10, 11, 12, 13, 14],
            vec![0, 1, 2, 3, 4, 10, 11, 12, 13, 14],
            vec![0, 1, 2, 3, 4, 10, 11, 12, 13, 14],
        ];
        let batches = vec![0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1];

        let result = kbet(&knn_data, &batches, 0.05);
        let poor_mixing_count = result.iter().filter(|&&x| x).count();

        println!(
            "Perfect mixing: {}/{} cells have poor mixing",
            poor_mixing_count,
            result.len()
        );
        assert!(poor_mixing_count < 3); // Should be almost none
    }

    #[test]
    fn test_kbet_partial_mixing() {
        // 90/10 mixing - much more skewed
        let knn_data = vec![
            vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 10], // 9 batch 0, 1 batch 1
            vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 10],
            vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 10],
            vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 10],
            vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 10],
            vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 10],
            vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 10],
            vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 10],
            vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 10],
            vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 10],
            vec![10, 11, 12, 13, 14, 15, 16, 17, 18, 0], // 9 batch 1, 1 batch 0
            vec![10, 11, 12, 13, 14, 15, 16, 17, 18, 0],
            vec![10, 11, 12, 13, 14, 15, 16, 17, 18, 0],
            vec![10, 11, 12, 13, 14, 15, 16, 17, 18, 0],
            vec![10, 11, 12, 13, 14, 15, 16, 17, 18, 0],
            vec![10, 11, 12, 13, 14, 15, 16, 17, 18, 0],
            vec![10, 11, 12, 13, 14, 15, 16, 17, 18, 0],
            vec![10, 11, 12, 13, 14, 15, 16, 17, 18, 0],
            vec![10, 11, 12, 13, 14, 15, 16, 17, 18, 0],
            vec![10, 11, 12, 13, 14, 15, 16, 17, 18, 0],
        ];
        let batches = vec![0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1];

        let result = kbet(&knn_data, &batches, 0.05);
        let poor_mixing_count = result.iter().filter(|&&x| x).count();

        println!(
            "Partial mixing: {}/{} cells have poor mixing",
            poor_mixing_count,
            result.len()
        );
        assert!(poor_mixing_count > 10); // 90/10 split should show significant poor mixing
    }
}
