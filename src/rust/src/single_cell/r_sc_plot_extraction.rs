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
/// @param f_path String. Path to the `counts_genes.bin` file.
/// @param cell_indices Integer positions (0-indexed!) that defines the cells
/// to keep.
/// @param gene_index Integer. Gene index position to return (0-indexed!).
/// @param norm Boolean. Shall normalised counts be returned.
/// @param scale Boolean. Shall the normalised counts be scaled.
/// @param clip Optional float. Clipping for the Z-scores if scale is set to
/// `TRUE`
///
/// @returns The dense vector of expression values for this gene.
///
/// @export
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

    let counts = if norm {
        let raw_counts = extract_raw_counts(f_path, &cell_indices, gene_index).to_extendr()?;
        raw_counts.iter().map(|x| *x as f32).collect()
    } else {
        extract_norm_counts(f_path, &cell_indices, gene_index, scale, clip).to_extendr()?
    };

    Ok(counts.r_float_convert())
}

/// Helper to extract single cell counts for several genes
///
/// @param f_path String. Path to the `counts_genes.bin` file.
/// @param cell_indices Integer positions (0-indexed!) that defines the cells
/// to keep.
/// @param gene_indices Integer. Gene index position to return (0-indexed!).
/// @param scale Boolean. Shall the normalised counts be scaled.
/// @param clip Optional float. Clipping for the Z-scores if scale is set to
/// `TRUE`
///
/// @returns A list of dense vectors of the normalised counts.
///
/// @export
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

    let all_counts = extract_norm_counts_multi(f_path, &cell_indices, &gene_indices, scale, clip)
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
/// Helper function to extract data for dot plots and/or heatmaps.
///
/// @param f_path String. Path to the `counts_genes.bin` file.
/// @param cell_indices Integer positions (0-indexed!) that defines the cells
/// to keep.
/// @param gene_indices Integer. Gene index position to return (0-indexed!).
/// @param group_ids Integer. The levels of the data. (0-indexed!)
/// @param group_levels String. Name of the factors.
///
/// @returns A list with the following elements:
/// \itemize{
///   \item grp_label - The label of that group
///   \item mean_exp - Vector of mean expression values in row major (genes x
///   n_levels)
///   \item perc_exp - Vector of proportions of cells with expression in row
///   major (genes x n_levels)
/// }
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

    let gene_res: GroupedGeneStats = extract_grouped_gene_stats(
        f_path,
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
