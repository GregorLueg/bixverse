use extendr_api::prelude::*;

use rayon::prelude::*;
use rustc_hash::FxHashMap;

use crate::core::enrichment::gsea::*;
use crate::utils::r_rust_interface::r_named_vec_data;

//////////////////////
// Helper functions //
//////////////////////

/// Calculates the traditional GSEA enrichment score
///
/// @param stats Named numerical vector. Needs to be sorted. The gene level statistics.
/// @param pathway_r String vector. The genes in the pathway.
///
/// @return The enrichment score
///
/// @export
#[extendr]
fn rs_calc_es(stats: Robj, pathway_r: Vec<String>) -> extendr_api::Result<f64> {
    let vec_data = r_named_vec_data(stats)?;

    let index: Vec<usize> = vec_data
        .0
        .iter()
        .enumerate()
        .filter_map(|(a, b)| if pathway_r.contains(b) { Some(a) } else { None })
        .collect();

    Ok(calculate_es(&vec_data.1, &index))
}

/// Helper function to rapidly retrieve the indices of the gene set members
///
/// @param gene_universe Character Vector. The genes represented in the gene universe.
/// @param pathway_list List. A named list with each element containing the genes for this
/// pathway.
///
/// @return Returns a list with the index positions of the gene set genes in the gene universe.
/// Importantly, these are indexed to R's 1-indexing!
///
/// @export
#[extendr]
fn rs_get_gs_indices(gene_universe: Vec<String>, pathway_list: List) -> extendr_api::Result<List> {
    // HashMap for fast look ups
    let gene_map: FxHashMap<&str, usize> = gene_universe
        .iter()
        .enumerate()
        .map(|(i, gene)| (gene.as_str(), i))
        .collect();

    let mut result_list = List::new(pathway_list.len());
    if let Some(names) = pathway_list.names() {
        result_list.set_names(names)?;
    }

    for i in 0..pathway_list.len() {
        let element = pathway_list.elt(i)?;
        if let Some(internal_vals) = element.as_string_vector() {
            let mut indices = Vec::with_capacity(internal_vals.len());

            for gene in &internal_vals {
                if let Some(&idx) = gene_map.get(gene.as_str()) {
                    indices.push((idx as i32) + 1);
                }
            }

            indices.sort_unstable();
            result_list.set_elt(i, Robj::from(indices))?;
        } else {
            result_list.set_elt(i, element)?;
        }
    }

    Ok(result_list)
}

////////////////////
// Core functions //
////////////////////

/// Rust implementation of the fgsea::calcGseaStat() function
///
/// @param stats Numeric vector. The gene level statistic. Needs to
/// sorted in descending nature.
/// @param gs_idx Integer vector. The indices of the gene set genes.
/// @param gsea_param Float. The GSEA parameter. Usually defaults to 1.0.
/// @param return_leading_edge Boolean. Return the leading edge indices.
/// @param return_all_extremes Boolean. Shall the extreme values be returned
/// for plotting.
///
/// @return List with the following elements
/// \itemize{
///     \item gene_stat Enrichment score for that gene set
///     \item leading_edge Indicies of the leading edge genes.
///     \item top Top values of the curve.
///     \item bottom Bottom values of the curve.
/// }
///
/// @export
#[extendr]
fn rs_calc_gsea_stats(
    stats: &[f64],
    gs_idx: &[i32],
    gsea_param: f64,
    return_leading_edge: bool,
    return_all_extremes: bool,
) -> List {
    let res: GseaStats = calc_gsea_stats(
        stats,
        gs_idx,
        gsea_param,
        return_leading_edge,
        return_all_extremes,
        true,
    );

    list!(
        es = res.es,
        leading_edge = res.leading_edge,
        top = res.top,
        bottom = res.bottom
    )
}

/// Helper function to generate traditional GSEA-based permutations
///
/// @param stats Numeric vector. The gene level statistic. Needs to
/// sorted in descending nature.
/// @param pathway_scores Numeric vector. The enrichment scores for the
/// pathways
/// @param pathway_sizes Integer vector. The sizes of the pathways.
/// @param iters Integer. Number of permutations.
/// @param seed Integer For reproducibility purposes
///
/// @return List with the following elements
/// \itemize{
///     \item es Enrichment scores for the gene sets
///     \item nes Normalised enrichment scores for the gene sets
///     \item pvals The calculated p-values.
///     \item n_more_extreme Number of times the enrichment score was
///     bigger or smaller than the permutation (pending sign).
///     \item size Pathway size.
/// }
///
/// @export
#[extendr]
fn rs_calc_gsea_stat_traditional_batch(
    stats: &[f64],
    pathway_scores: &[f64],
    pathway_sizes: &[i32],
    iters: usize,
    seed: u64,
) -> extendr_api::Result<List> {
    let pathway_sizes: Vec<usize> = pathway_sizes.iter().map(|x| *x as usize).collect();

    let batch_res: GseaBatchResults =
        calc_gsea_stat_traditional_batch(stats, pathway_scores, &pathway_sizes, iters, seed);

    let gsea_res: GseaResults<'_> =
        calculate_nes_es_pval(pathway_scores, &pathway_sizes, &batch_res);

    Ok(list!(
        es = gsea_res.es,
        nes = gsea_res.nes,
        pvals = gsea_res.pvals,
        n_more_extreme = gsea_res.n_more_extreme,
        size = gsea_res.size
    ))
}

/// Helper function to generate fgsea simple-based permutations
///
/// @param stats Numeric vector. The gene level statistic. Needs to
/// sorted in descending nature.
/// @param pathway_scores Numeric vector. The enrichment scores for the
/// pathways
/// @param pathway_sizes Integer vector. The sizes of the pathways.
/// @param iters Integer. Number of permutations.
/// @param gsea_param Float. The Gene Set Enrichment parameter.
/// @param return_add_stats Boolean. Returns additional statistics
/// necessary for the multi-level calculations.
/// @param seed Integer. For reproducibility purposes
///
/// @return List with the following elements
/// \itemize{
///     \item es Enrichment scores for the gene sets
///     \item nes Normalised enrichment scores for the gene sets
///     \item pvals The calculated p-values.
///     \item n_more_extreme Number of times the enrichment score was
///     bigger or smaller than the permutation (pending sign).
///     \item size Pathway size.
/// }
///
/// If `return_add_stats` is set to true, there is additional elements in the
/// list:
/// \itemize{
///     \item le_zero Number of times the permutation was less than zero.
///     \item ge_zero Number of times the permutation was greater than zero.
/// }
///
/// @export
#[extendr]
fn rs_calc_gsea_stat_cumulative_batch(
    stats: &[f64],
    pathway_scores: &[f64],
    pathway_sizes: &[i32],
    iters: usize,
    gsea_param: f64,
    return_add_stats: bool,
    seed: u64,
) -> extendr_api::Result<List> {
    // Convert indices from R - keep as 1-based for algorithm
    let pathway_sizes: Vec<usize> = pathway_sizes.iter().map(|x| *x as usize).collect();

    let batch_res: GseaBatchResults = calc_gsea_stat_cumulative_batch(
        stats,
        pathway_scores,
        &pathway_sizes,
        iters,
        gsea_param,
        seed,
    );

    let gsea_res: GseaResults<'_> =
        calculate_nes_es_pval(pathway_scores, &pathway_sizes, &batch_res);

    let res = if return_add_stats {
        list!(
            es = gsea_res.es,
            nes = gsea_res.nes,
            pvals = gsea_res.pvals,
            n_more_extreme = gsea_res.n_more_extreme,
            le_zero = gsea_res.le_zero,
            ge_zero = gsea_res.ge_zero,
            size = gsea_res.size
        )
    } else {
        list!(
            es = gsea_res.es,
            nes = gsea_res.nes,
            pvals = gsea_res.pvals,
            n_more_extreme = gsea_res.n_more_extreme,
            size = gsea_res.size
        )
    };

    Ok(res)
}

/// Calculates p-values for pre-processed data
///
/// @param stats Named numerical vector. Needs to be sorted. The gene level statistics.
/// @param es Numerical vector. The enrichment scores of the pathways of that specific size
/// @param pathway_size Integer. The size of the pathways to test.
/// @param sample_size Integer. The size of the random gene sets to test against.
/// @param seed Integer. Random seed.
/// @param eps Float. Boundary for calculating the p-value.
/// @param sign Boolean. Used for the only positive or only negative score version.
///
/// @return List with the following elements:
/// \itemize{
///     \item pvals The pvalues.
///     \item is_cp_ge_half Flag indicating if conditional probability is ≥ 0.5. Indicates
///     overesimation of the p-values.
/// }
///
/// @export
#[extendr]
fn rs_calc_multi_level(
    stats: Robj,
    es: &[f64],
    pathway_size: &[i32],
    sample_size: usize,
    seed: u64,
    eps: f64,
    sign: bool,
) -> extendr_api::Result<List> {
    let (_, ranks) = r_named_vec_data(stats)?;
    let pathway_size: Vec<usize> = pathway_size
        .iter()
        .map(|&x| x.try_into().unwrap_or(0))
        .collect();

    let chunk_size = std::cmp::max(1, es.len() / (rayon::current_num_threads() * 4));

    let res: Vec<(f64, bool)> = es
        .par_iter()
        .zip(pathway_size.par_iter())
        .chunks(chunk_size)
        .flat_map(|chunk| {
            let mut local_results = Vec::with_capacity(chunk.len());

            for (es_i, size_i) in chunk {
                let res_i =
                    fgsea_multilevel_helper(*es_i, &ranks, *size_i, sample_size, seed, eps, sign);
                local_results.push(res_i);
            }

            local_results
        })
        .collect();

    let mut pvals: Vec<f64> = Vec::with_capacity(res.len());
    let mut is_cp_ge_half: Vec<bool> = Vec::with_capacity(res.len());

    for (pval, cp_flag) in res {
        pvals.push(pval);
        is_cp_ge_half.push(cp_flag);
    }

    Ok(list!(pvals = pvals, is_cp_ge_half = is_cp_ge_half))
}

/// Calculates the simple and multi error for fgsea multi level
///
/// @param n_more_extreme Integer vector. The number of times the ES was larger than the
/// permutations.
/// @param nperm Integer. Number of permutations.
/// @param sample_size Integer. Number of samples.
///
/// @return List with the following elements:
/// \itemize{
///     \item simple_err Vector of simple errors.
///     \item multi_err Vector of multi errors.
/// }
///
/// @export
#[extendr]
fn rs_simple_and_multi_err(n_more_extreme: &[i32], nperm: usize, sample_size: usize) -> List {
    // Conversion needed
    let n_more_extreme: Vec<usize> = n_more_extreme.iter().map(|x| *x as usize).collect();

    let res: MultiLevelErrRes = calc_simple_and_multi_error(&n_more_extreme, nperm, sample_size);

    list!(simple_err = res.0, multi_err = res.1)
}

extendr_module! {
    mod r_gsea;
    fn rs_calc_es;
    fn rs_get_gs_indices;
    fn rs_calc_gsea_stats;
    fn rs_calc_gsea_stat_cumulative_batch;
    fn rs_calc_gsea_stat_traditional_batch;
    fn rs_calc_multi_level;
    fn rs_simple_and_multi_err;
}
