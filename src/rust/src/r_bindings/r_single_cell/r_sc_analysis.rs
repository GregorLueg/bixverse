use extendr_api::*;
use faer::Mat;

use crate::single_cell::dge_aucs::*;
use crate::utils::r_rust_interface::faer_to_r_matrix;

/// Calculate DGEs between cells based on Mann Whitney stats
///
/// @description
/// The function will take two sets of cell indices and calculate the
/// differential gene expression based on the Mann Whitney test between the
/// two groups.
///
/// @param f_path String. Path to the `counts_cells.bin` file.
/// @param cell_indices_1 Integer. Index positions (0-indexed) of the cells
/// of group 1.
/// @param cell_indices_2 Integer. Index positions (0-indexed) of the cells
/// of group 2.
/// @param min_prop Minimum proportion of expression in at least one of the
/// two groups to be tested.
/// @alternative String. One of `c("twosided", "greater", "less")`. Null
/// hypothesis.
/// @param verbose Boolean. Controls verbosity of the function.
///
/// @return A list with the following elements
/// \itemize{
///   \item lfc - Log fold changes between the two groups.
///   \item prop1 - Proportion of cells expressing the gene in group 1.
///   \item prop2 - Proportion of cells expressing the gene in group 2.
///   \item z_scores - Z-scores based on the Mann Whitney statistic.
///   \item p_values - P-values of the Mann Whitney statistic.
///   \item fdr - False discovery rate after BH adjustment
///   \item genes_to_keep - Boolean indicating which genes were tested.
/// }
///
/// @export
#[extendr]
fn rs_calculate_dge_mann_whitney(
    f_path: String,
    cell_indices_1: &[i32],
    cell_indices_2: &[i32],
    min_prop: f64,
    alternative: String,
    verbose: bool,
) -> extendr_api::Result<List> {
    let cell_indices_1 = cell_indices_1
        .iter()
        .map(|x| *x as usize)
        .collect::<Vec<usize>>();
    let cell_indices_2 = cell_indices_2
        .iter()
        .map(|x| *x as usize)
        .collect::<Vec<usize>>();

    let dge_results: DgeMannWhitneyRes = calculate_dge_grps_mann_whitney(
        &f_path,
        &cell_indices_1,
        &cell_indices_2,
        min_prop as f32,
        &alternative,
        verbose,
    )?;

    Ok(list!(
        lfc = dge_results
            .lfc
            .iter()
            .map(|x| *x as f64)
            .collect::<Vec<f64>>(),
        prop1 = dge_results
            .prop1
            .iter()
            .map(|x| *x as f64)
            .collect::<Vec<f64>>(),
        prop2 = dge_results
            .prop2
            .iter()
            .map(|x| *x as f64)
            .collect::<Vec<f64>>(),
        z_scores = dge_results.z_scores,
        p_values = dge_results.p_vals,
        fdr = dge_results.fdr,
        genes_to_keep = dge_results.genes_to_keep
    ))
}

#[extendr]
fn rs_aucell(
    f_path: String,
    gs_list: List,
    auc_type: &str,
    streaming: bool,
    verbose: bool,
) -> extendr_api::Result<extendr_api::RArray<f64, [usize; 2]>> {
    let mut gs_indices: Vec<Vec<usize>> = Vec::with_capacity(gs_list.len());

    for i in 0..gs_list.len() {
        let r_obj = gs_list.elt(i).unwrap();
        let int = r_obj
            .as_integer_vector()
            .unwrap()
            .iter()
            .map(|x| *x as usize)
            .collect::<Vec<usize>>();
        gs_indices.push(int);
    }

    let res = if streaming {
        calculate_aucell_streaming(&f_path, &gs_indices, auc_type, verbose)?
    } else {
        calculate_aucell(&f_path, &gs_indices, auc_type, verbose)?
    };

    let auc_mat = Mat::from_fn(res.len(), res[0].len(), |i, j| res[i][j] as f64);

    Ok(faer_to_r_matrix(auc_mat.as_ref()))
}

extendr_module! {
    mod r_sc_analysis;
    fn rs_calculate_dge_mann_whitney;
    fn rs_aucell;
}
