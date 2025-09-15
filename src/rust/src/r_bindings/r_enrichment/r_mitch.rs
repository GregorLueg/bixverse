use extendr_api::prelude::*;
use rayon::prelude::*;

use crate::core::base::stats::calc_fdr;
use crate::core::enrichment::mitch::*;
use crate::utils::general::flatten_vector;
use crate::utils::r_rust_interface::*;

/// Calculate mitch enrichment leveraging Rust under the hood
///
/// @param x Numerical matrix. Each column represents on the contrasts you
/// wish to test for and the rows represent the gene statistics per contrast.
/// @param pathway_list Named list. Each element represents one of the pathways
/// to test for.
/// @param min_size Integer. Minimum size of gene the gene set to be tested for.
///
/// @return A list with the following elements:
///  \itemize{
///     \item pathway_names - The name of the pathway.
///     \item pathway_sizes The size of the pathway.
///     \item manova_pvals - The p-value of the MANOVA test.
///     \item anova_pvals The p-values of the ANOVA test on top of the MANOVA
///     results. Total length = `ncol(x)` * number of pathways.
///     \item scores - The scores for each pathway set, contrast. Same length
///     as `anova_pvals`.
///     \item s_dist - Calculated distances from the hypotenuse.
///     \item sd - SDs of the scores.
/// }
///
/// @export
#[extendr]
fn rs_mitch_calc(
    x: RMatrix<f64>,
    pathway_list: List,
    min_size: usize,
) -> extendr_api::Result<List> {
    let data = r_matrix_to_faer(&x);
    let ranked_data = mitch_rank(&data);
    let r_row_names = x.get_rownames().unwrap();
    let mut row_names = Vec::with_capacity(r_row_names.len());
    for e in r_row_names.iter() {
        row_names.push(e.to_string());
    }

    let (pathway_names, pathway_indices) =
        prepare_mitch_pathways(&row_names, pathway_list, min_size)?;

    let mitch_res: Vec<MitchResult<'_>> = pathway_names
        .par_iter()
        .zip(pathway_indices.par_iter())
        .map(|(p_name, p_idx)| process_mitch_pathway(ranked_data.as_ref(), p_name, p_idx))
        .collect();

    let mut pathway_names = Vec::with_capacity(mitch_res.len());
    let mut pathway_sizes = Vec::with_capacity(mitch_res.len());
    let mut manova_pvals = Vec::with_capacity(mitch_res.len());
    let mut anova_pvals = Vec::with_capacity(mitch_res.len());
    let mut scores = Vec::with_capacity(mitch_res.len());
    let mut s_dist = Vec::with_capacity(mitch_res.len());
    let mut sd = Vec::with_capacity(mitch_res.len());

    for res in mitch_res {
        pathway_names.push(res.pathway_name.to_string());
        pathway_sizes.push(res.pathway_size);
        manova_pvals.push(res.manova_pval);
        anova_pvals.push(res.anova_pvals);
        scores.push(res.scores);
        s_dist.push(res.s_dist);
        sd.push(res.mysd);
    }

    let anova_pvals = flatten_vector(anova_pvals);
    let scores = flatten_vector(scores);

    Ok(list!(
        pathway_names = pathway_names,
        pathway_sizes = pathway_sizes,
        manova_pval = manova_pvals.clone(),
        manova_fdr = calc_fdr(&manova_pvals),
        anova_pvals = anova_pvals,
        scores = scores,
        s_dist = s_dist,
        sd = sd
    ))
}

extendr_module! {
    mod r_mitch;
    fn rs_mitch_calc;
}
