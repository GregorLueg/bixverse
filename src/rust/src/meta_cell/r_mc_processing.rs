//! Functions to apply on meta cells. These are in-memory versions of the
//! variants used for single cells due to the massive compression achieved
//! with meta cells.

use bixverse_rs::prelude::*;
use bixverse_rs::single_cell::mc_processing::hvg_pca::*;
use bixverse_rs::single_cell::sc_processing::hvg::*;
use bixverse_rs::utils::r_rust_interface::list_to_sparse_matrix;
use extendr_api::*;

use crate::meta_cell::utils::cast_compressed_sparse_data;

/////////////
// extendR //
/////////////

extendr_module! {
    mod r_mc_processing;
    // hvg and pca
    fn rs_mc_hvg;
}

///////////////////////////
// Highly variable genes //
///////////////////////////

/// Meta cells highly variable genes
///
/// @param sparse_data A named list that needs to have `data`, `indptr`,
/// `indices`, `nrow`, `ncol` and `format`.
/// @param hvg_method String. Which HVG detection method to use. Options
/// are `c("vst", "meanvarbin", "dispersion")`.
/// @param loess_span Numeric. The span parameter for the loess function
/// (only used for `"vst"`).
/// @param clip_max Optional clipping number. Defaults to `sqrt(no_cells)` if
/// not provided (only used for `"vst"`).
/// @param binning String. The binning strategy for the `meanvarbin` and
/// `dispersion` methods. One of `c("equal_width", "equal_frequency")`.
/// @param n_bins Integer. Number of bins for the `meanvarbin` and
/// `dispersion` methods.
///
/// @return A list with the HVG statistics. If `hvg_method == "vst"`:
/// \itemize{
///   \item mean - The average expression of the gene.
///   \item var - The variance of the gene.
///   \item var_exp - The expected variance of the gene.
///   \item var_std - The standardised variance of the gene.
/// }
/// For `"meanvarbin"` and `"dispersion"`:
/// \itemize{
///   \item mean - The average expression of the gene.
///   \item dispersion - The dispersion of the gene.
///   \item dispersion_scaled - The scaled dispersion per bin per gene.
///   \item bin - The bin of the gene.
/// }
///
/// @export
#[extendr]
fn rs_mc_hvg(
    sparse_data: List,
    hvg_method: &str,
    loess_span: f64,
    binning: String,
    n_bins: usize,
    clip_max: Option<f32>,
) -> List {
    let sparse: CompressedSparseData2<f64, f64> = list_to_sparse_matrix(sparse_data, true);
    let sparse = cast_compressed_sparse_data(sparse);

    let hvg_type = parse_hvg_method(hvg_method)
        .ok_or_else(|| format!("Invalid HVG method: {}", hvg_method))
        .unwrap();

    match hvg_type {
        HvgMethod::Vst => {
            let res = get_hvg_vst_from_sparse(&sparse, loess_span as f32, clip_max);
            list!(
                mean = res.mean,
                var = res.var,
                var_exp = res.var_exp,
                var_std = res.var_std
            )
        }
        HvgMethod::MeanVarBin => {
            let res = get_hvg_mvb_from_sparse(&sparse, &binning, n_bins);
            list!(
                mean = res.mean,
                dispersion = res.dispersion,
                dispersion_scaled = res.dispersion_scaled,
                bin = res.bin
            )
        }
        HvgMethod::Dispersion => {
            let res = get_hvg_dispersion_from_sparse(&sparse, &binning, n_bins);
            list!(
                mean = res.mean,
                dispersion = res.dispersion,
                dispersion_scaled = res.dispersion_scaled,
                bin = res.bin
            )
        }
    }
}
