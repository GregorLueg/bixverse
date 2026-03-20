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
) -> Vec<f64> {
    let cell_indices = cell_indices.r_int_convert();

    let counts = if norm {
        let raw_counts = extract_raw_counts(f_path, &cell_indices, gene_index);
        raw_counts.iter().map(|x| *x as f32).collect()
    } else {
        extract_norm_counts(f_path, &cell_indices, gene_index, scale, clip)
    };

    counts.r_float_convert()
}

/// Helper to extract single cell counts for several genes
///
/// @param f_path String. Path to the `counts_genes.bin` file.
/// @param cell_indices Integer positions (0-indexed!) that defines the cells
/// to keep.
/// @param gene_index Integer. Gene index position to return (0-indexed!).
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
) -> extendr_api::Result<List> {
    let cell_indices = cell_indices.r_int_convert();
    let gene_indices = gene_indices.r_int_convert();

    let all_counts = extract_norm_counts_multi(f_path, &cell_indices, &gene_indices, scale, clip);

    let mut res = List::new(cell_indices.len());

    for i in 0..res.len() {
        let vec_i = all_counts[i].clone().r_float_convert();
        res.set_elt(i, Robj::from(vec_i))?;
    }

    Ok(res)
}

fn rs_extract_grouped_gene_stats(
    f_path: &str,
    cell_indices: &[i32],
    gene_indices: &[i32],
    group_ids: Vec<String>,
) {
}
