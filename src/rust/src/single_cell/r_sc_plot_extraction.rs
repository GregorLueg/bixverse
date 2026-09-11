use extendr_api::*;

use bixverse_rs::prelude::*;
use bixverse_rs::single_cell::sc_data::plotting::*;

/////////////
// extendR //
/////////////

extendr_module! {
    mod r_sc_plot_extraction;
    fn rs_extract_counts_plots;
    fn rs_extract_several_genes_plots;
    fn rs_extract_grouped_gene_stats;
}

//////////////////////////
// Plot data extractors //
//////////////////////////

/// Helper to extract single cell counts as a dense vector for plotting
///
/// @description
/// `r lifecycle::badge("experimental")`
/// Extract dense counts for this given gene.
///
/// @param f_path String. Path to the `counts_genes.bin` file.
/// @param cell_indices Integer positions (0-indexed!) that defines the cells
/// to keep.
/// @param gene_index Integer. Gene index position to return (0-indexed!).
/// @param norm Boolean. Shall normalised counts be returned.
/// @param scale Boolean. Shall the normalised counts be z-scored across the
/// selected cells.
/// @param clip Optional float. Clips the Z-scores to `[-clip, clip]`. Only
/// used if `scale = TRUE`.
///
/// @returns Numerical vector with one expression value per cell in
/// `cell_indices`.
///
/// @export
///
/// @keywords internal
#[extendr]
fn rs_extract_counts_plots(
    f_path: &str,
    cell_indices: &[i32],
    gene_index: usize,
    norm: bool,
    scale: bool,
    clip: Option<f32>,
) -> Result<Vec<f64>> {
    let cell_indices = cell_indices.r_int_convert();
    let reader = ParallelSparseReader::new(f_path).to_extendr()?;

    let counts = if norm {
        extract_norm_counts(&reader, &cell_indices, gene_index, scale, clip).to_extendr()?
    } else {
        let raw_counts = extract_raw_counts(&reader, &cell_indices, gene_index).to_extendr()?;
        raw_counts.iter().map(|x| *x as f32).collect()
    };

    Ok(counts.r_float_convert())
}

/// Helper to extract single cell counts for several genes
///
/// @description
/// `r lifecycle::badge("experimental")`
/// Extract the normalised single cell counts of several genes at once.
///
/// @param f_path String. Path to the `counts_genes.bin` file.
/// @param cell_indices Integer positions (0-indexed!) that defines the cells
/// to keep.
/// @param gene_indices Integer vector. Gene index positions to return
/// (0-indexed!).
/// @param scale Boolean. Shall the normalised counts be z-scored per gene
/// across the selected cells.
/// @param clip Optional float. Clips the Z-scores to `[-clip, clip]`. Only
/// used if `scale = TRUE`.
///
/// @returns A list of numerical vectors, one per gene in `gene_indices`, each
/// with one normalised value per cell in `cell_indices`.
///
/// @export
///
/// @keywords internal
#[extendr]
fn rs_extract_several_genes_plots(
    f_path: &str,
    cell_indices: &[i32],
    gene_indices: &[i32],
    scale: bool,
    clip: Option<f32>,
) -> Result<List> {
    let cell_indices = cell_indices.r_int_convert();
    let gene_indices = gene_indices.r_int_convert();
    let reader = ParallelSparseReader::new(f_path).to_extendr()?;

    let all_counts = extract_norm_counts_multi(&reader, &cell_indices, &gene_indices, scale, clip)
        .to_extendr()?;

    let mut res = List::new(all_counts.len());

    for i in 0..res.len() {
        let vec_i = all_counts[i].clone().r_float_convert();
        res.set_elt(i, Robj::from(vec_i))?;
    }

    Ok(res)
}

/// Calculates the gene statistics for a set of cell groups and genes
///
/// @description
/// `r lifecycle::badge("experimental")`
/// Helper function to extract data for dot plots and/or heatmaps.
///
/// @param f_path String. Path to the `counts_genes.bin` file.
/// @param cell_indices Integer positions (0-indexed!) that defines the cells
/// to keep.
/// @param gene_indices Integer vector. Gene index positions to return
/// (0-indexed!).
/// @param group_ids Integer vector. Group of each cell in `cell_indices`, as
/// an index into `group_levels` (0-indexed!). Same length as `cell_indices`.
/// @param group_levels Character vector. The group labels.
///
/// @returns A list with the following elements:
/// \itemize{
///   \item grp_label - The group labels, i.e. `group_levels`.
///   \item mean_exp - Mean normalised expression per gene and group over all
///   cells of the group (zeros included), row-major (genes x groups).
///   \item perc_exp - Fraction (`[0, 1]`) of cells in the group with a
///   non-zero count, row-major (genes x groups).
/// }
///
/// @export
///
/// @keywords internal
#[extendr]
fn rs_extract_grouped_gene_stats(
    f_path: &str,
    cell_indices: &[i32],
    gene_indices: &[i32],
    group_ids: &[i32],
    group_levels: Vec<String>,
) -> Result<List> {
    let cell_indices = cell_indices.r_int_convert();
    let gene_indices = gene_indices.r_int_convert();
    let group_ids = group_ids.r_int_convert();
    let reader = ParallelSparseReader::new(f_path).to_extendr()?;

    let gene_res: GroupedGeneStats = extract_grouped_gene_stats(
        &reader,
        &cell_indices,
        &gene_indices,
        &group_ids,
        &group_levels,
    )
    .to_extendr()?;

    Ok(list!(
        grp_label = gene_res.group_labels,
        mean_exp = gene_res.mean_expression.r_float_convert(),
        perc_exp = gene_res.pct_expressed.r_float_convert(),
    ))
}
