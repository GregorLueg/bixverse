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
///
/// ### Return
///
/// Numerical vector indicating with the p-values from the ChiSquare test
pub fn kbet(knn_data: &Vec<Vec<usize>>, batches: &Vec<usize>) -> Vec<f64> {
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
            1.0 - ChiSquared::new(dof).unwrap().cdf(chi_square)
        })
        .collect()
}
