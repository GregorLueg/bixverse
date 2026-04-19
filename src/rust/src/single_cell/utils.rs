//! Helpers for the single cell part.

use bixverse_rs::single_cell::sc_processing::hvg::HvgDispersionRes;
use extendr_api::*;

/// Helper function to flatten the dispersion results to a List
///
/// ### Params
///
/// * `results` - A vector of `HvgDispersionRes`.
///
/// ### Returns
///
/// The results list
pub fn flatten_dispersion_batches(results: Vec<HvgDispersionRes>) -> List {
    let n_genes = results[0].mean.len();
    let total_len = n_genes * results.len();
    let mut mean_flat = Vec::with_capacity(total_len);
    let mut disp_flat = Vec::with_capacity(total_len);
    let mut disp_scaled_flat = Vec::with_capacity(total_len);
    let mut bin_flat = Vec::with_capacity(total_len);
    let mut batch_idx = Vec::with_capacity(total_len);
    let mut gene_idx = Vec::with_capacity(total_len);

    for (batch, res) in results.into_iter().enumerate() {
        mean_flat.extend(res.mean);
        disp_flat.extend(res.dispersion);
        disp_scaled_flat.extend(res.dispersion_scaled);
        bin_flat.extend(res.bin);
        batch_idx.extend(vec![batch as i32; n_genes]);
        gene_idx.extend(0..n_genes as i32);
    }

    list!(
        mean = mean_flat,
        dispersion = disp_flat,
        dispersion_scaled = disp_scaled_flat,
        bin = bin_flat,
        batch = batch_idx,
        gene_idx = gene_idx
    )
}
