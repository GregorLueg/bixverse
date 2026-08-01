//! Rust <> R interface for spatially variable gene detection.

use bixverse_rs::prelude::*;
use bixverse_rs::spatial::sp_processing::morans_i::{MoranIRes, MoransI};
use bixverse_rs::spatial::sp_processing::sparkx::{SparkX, SparkXParams, SparkXRes};
use bixverse_rs::spatial::sp_processing::svg_utils::SpatialSvgParams;
use extendr_api::*;

use crate::spatial::utils::{adjacency_from_r_lists, coords_from_r_matrix, row_major_to_r_matrix};

////////////////////
// extendr Module //
////////////////////

extendr_module! {
    mod r_sp_processing;
    fn rs_sp_morans_i;
    fn rs_sp_sparkx;
}

/////////////////
// Moran's I //
/////////////////

/// Detect spatially variable genes with Moran's I
///
/// @description
/// `r lifecycle::badge("experimental")`
/// Per-gene Moran's I against a precomputed spatial graph, tested under the
/// Cliff-Ord normality null. Expectation and variance depend only on the spot
/// count and the graph weights, so both are scalars across genes and get
/// computed once. Per-gene work is a single sparse quadratic form.
///
/// Index spaces differ between the two spot arguments and mixing them up is
/// the classic failure here. `spots_to_keep` is **global** (0-based positions
/// in the count store), while `neighbours` and `weights` are **local**
/// (0-based, positional against `spots_to_keep`). The output of
/// [bixverse::rs_sp_build_graph()] is already local, so the two compose.
///
/// @param f_path_gene String. Path to the `counts_genes.bin` file.
/// @param spots_to_keep Integer. Global 0-based spot indices making up this
/// sample.
/// @param neighbours List of integer vectors. Local 0-based neighbour indices,
/// one element per spot. Must be as long as `spots_to_keep`.
/// @param weights List of numeric vectors. Edge weights, aligned element for
/// element with `neighbours`.
/// @param gene_indices Integer. On-disk 0-based gene indices to test.
/// @param svg_params List. With element `assay`, one of `c("raw", "norm")`.
/// @param streaming Logical. If `TRUE`, gene chunks are loaded in batches of
/// 1000 rather than all at once.
/// @param verbose Integer. `0L` - quiet; `1L` - normal verbosity; `2L` -
/// detailed verbosity.
///
/// @return A list with the following elements
/// \itemize{
///   \item gene_idx - Integer. On-disk gene indices that were tested.
///   \item i - Numeric. Moran's I per gene. `NaN` for constant genes.
///   \item e_i - Numeric scalar. Expected I under the null.
///   \item var_i - Numeric scalar. Variance of I under the null.
///   \item z - Numeric. Standardised statistic per gene.
///   \item pval - Numeric. Two-sided p-value per gene.
///   \item fdr - Numeric. BH-adjusted p-value per gene.
/// }
///
/// @references
/// Cliff & Ord, Spatial Autocorrelation, 1973
///
/// @export
///
/// @keywords internal
#[extendr]
#[allow(clippy::too_many_arguments)]
fn rs_sp_morans_i(
    f_path_gene: String,
    spots_to_keep: Vec<i32>,
    neighbours: List,
    weights: List,
    gene_indices: Vec<i32>,
    svg_params: List,
    streaming: bool,
    verbose: usize,
) -> extendr_api::Result<List> {
    let spots_to_keep: Vec<usize> = spots_to_keep.r_int_convert();
    let gene_indices: Vec<usize> = gene_indices.r_int_convert();
    let params = SpatialSvgParams::from_r_list(svg_params)?;

    let (neighbours, mut weights) = adjacency_from_r_lists(&neighbours, &weights)?;

    let reader = ParallelSparseReader::new(&f_path_gene).to_extendr()?;

    let morans =
        MoransI::new(&reader, &spots_to_keep, &neighbours, &mut weights, params).to_extendr()?;

    let res: MoranIRes = if streaming {
        morans
            .compute_all_genes_streaming(&gene_indices, verbose)
            .to_extendr()?
    } else {
        morans
            .compute_all_genes(&gene_indices, verbose)
            .to_extendr()?
    };

    Ok(list!(
        gene_idx = res.gene_idx.iter().map(|&x| x as i32).collect::<Vec<i32>>(),
        i = res.i,
        e_i = res.e_i,
        var_i = res.var_i,
        z = res.z,
        pval = res.pval,
        fdr = res.fdr
    ))
}

////////////
// SPARKX //
////////////

/// Detect spatially variable genes with SPARK-X
///
/// @description
/// `r lifecycle::badge("experimental")`
/// Kernel-based non-parametric test for spatial expression patterns. Each
/// kernel gets a low-rank factorisation of the centred kernel matrix, the
/// per-gene statistic comes out of a projection against that factor, and the
/// per-kernel p-values are combined with the Cauchy test.
///
/// An empty kernel bank triggers the paper-default five Gaussian plus five
/// cosine bandwidths, derived from quantiles of the pairwise coordinate
/// distance distribution on a sub-sample.
///
/// @param f_path_gene String. Path to the `counts_genes.bin` file.
/// @param spots_to_keep Integer. Global 0-based spot indices making up this
/// sample.
/// @param coords Numeric matrix. Two columns, `x` then `y`, in full-resolution
/// pixels. One row per spot, aligned with `spots_to_keep`.
/// @param gene_indices Integer. On-disk 0-based gene indices to test.
/// @param sparkx_params List. With the following elements:
/// \itemize{
///   \item kernels - List or `NULL`. Each element a list with `kernel`
///   (`"gaussian"` or `"cosine"`) and `bandwidth`. `NULL` or empty uses the
///   defaults derived from the coordinates.
///   \item n_landmarks - Integer. Nystroem landmarks per Gaussian kernel.
///   \item bandwidth_subsample - Integer. Sub-sample size for the default
///   bandwidth estimation.
/// }
/// @param seed Integer. Seed for landmark selection and bandwidth
/// sub-sampling.
/// @param streaming Logical. If `TRUE`, gene chunks are loaded in batches of
/// 1000 rather than all at once.
/// @param verbose Integer. `0L` - quiet; `1L` - normal verbosity; `2L` -
/// detailed verbosity.
///
/// @return A list with the following elements
/// \itemize{
///   \item gene_idx - Integer. On-disk gene indices that were tested.
///   \item kernels - Character. Kernel labels, one per column of the two
///   matrices below.
///   \item stat_per_kernel - Numeric matrix, genes x kernels.
///   \item pval_per_kernel - Numeric matrix, genes x kernels.
///   \item pval_combined - Numeric. Cauchy-combined p-value per gene.
///   \item fdr - Numeric. BH-adjusted `pval_combined`.
/// }
///
/// @references
/// Zhu, Sun & Zhou, Genome Biology, 2021
///
/// @export
///
/// @keywords internal
#[extendr]
#[allow(clippy::too_many_arguments)]
fn rs_sp_sparkx(
    f_path_gene: String,
    spots_to_keep: Vec<i32>,
    coords: RMatrix<f64>,
    gene_indices: Vec<i32>,
    sparkx_params: List,
    seed: usize,
    streaming: bool,
    verbose: usize,
) -> extendr_api::Result<List> {
    let spots_to_keep: Vec<usize> = spots_to_keep.r_int_convert();
    let gene_indices: Vec<usize> = gene_indices.r_int_convert();
    let coordinates = coords_from_r_matrix(&coords)?;
    let params = SparkXParams::from_r_list(sparkx_params)?;

    let reader = ParallelSparseReader::new(&f_path_gene).to_extendr()?;

    let sparkx = SparkX::new(&reader, &spots_to_keep, &coordinates, params, seed).to_extendr()?;

    let res: SparkXRes = if streaming {
        sparkx
            .compute_all_genes_streaming(&gene_indices, verbose)
            .to_extendr()?
    } else {
        sparkx
            .compute_all_genes(&gene_indices, verbose)
            .to_extendr()?
    };

    let n_genes = res.gene_idx.len();

    Ok(list!(
        gene_idx = res.gene_idx.iter().map(|&x| x as i32).collect::<Vec<i32>>(),
        kernels = res.kernels,
        stat_per_kernel = row_major_to_r_matrix(&res.stat_per_kernel, n_genes, res.n_kernels),
        pval_per_kernel = row_major_to_r_matrix(&res.pval_per_kernel, n_genes, res.n_kernels),
        pval_combined = res.pval_combined,
        fdr = res.fdr
    ))
}
