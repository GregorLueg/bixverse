use extendr_api::prelude::*;

use std::collections::HashMap;

use crate::helpers_fgsea::*;
use crate::utils_r_rust::{r_list_to_str_vec, r_named_vec_data};
use crate::utils_rust::{array_max, array_min, cumsum};

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
///
/// @export
#[extendr]
fn rs_get_gs_indices(gene_universe: Vec<String>, pathway_list: List) -> extendr_api::Result<List> {
    // HashMap for fast look ups
    let gene_map: HashMap<&str, usize> = gene_universe
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

/// Traditional Gene Set Enrichment analysis (in Rust)
///
/// @description This function serves as an internal control. It implements the
/// gene set enrichment analysis in the traditional way without the convex approximations
/// used in fgsea implementations. TODO: should also return leading edge genes.
///
/// @param stats Named numerical vector. Needs to be sorted. The gene level statistics.
/// @param iters Integer. Number of permutations to test for.
/// @param pathway_list List. A named list with each element containing the genes for this
/// pathway.
/// @param seed Integer. For reproducibility purposes
///
/// @return List with the following elements
/// \itemize{
///     \item es Enrichment scores for the gene sets
///     \item nes Normalised enrichment scores for the gene sets
///     \item pvals The calculated p-values.
/// }
///
/// @export
#[extendr]
fn rs_gsea_traditional(
    stats: Robj,
    iters: usize,
    pathway_list: List,
    seed: u64,
) -> extendr_api::Result<List> {
    let vec_data = r_named_vec_data(stats)?;
    // Transform to  Vec of Vecs
    let gene_sets = r_list_to_str_vec(pathway_list)?;

    let results = calc_gsea_traditional(&vec_data.1, &vec_data.0, gene_sets, iters, seed);

    let mut es_vec = Vec::with_capacity(results.len());
    let mut nes_vec = Vec::with_capacity(results.len());
    let mut pvals_vec = Vec::with_capacity(results.len());

    for (es, nes, pval) in results {
        es_vec.push(es);
        nes_vec.push(nes);
        pvals_vec.push(pval);
    }

    Ok(list!(es = es_vec, nes = nes_vec, pvals = pvals_vec))
}

/// Rust implementation of the fgsea::calcGseaStat() function
///
/// @param stats Numeric vector. The gene level statistic. Needs to
/// sorted in descending nature.
/// @param gs_idx Integer vector. The indices of the gene set genes.
/// @param gsea_param Float. The GSEA parameter. Usually defaults to 1.0.
/// @param return_leading_edge Boolean. Return the leading edge indices.
///
/// @return List with the following elements
/// \itemize{
///     \item gene_stat Enrichment score for that gene set
///     \item leading_edge Indicies of the leading edge genes.
/// }
///
/// @export
#[extendr]
fn rs_calc_gsea_stats(
    stats: &[f64],
    gs_idx: &[i32],
    gsea_param: f64,
    return_leading_edge: bool,
) -> List {
    let n = stats.len();
    let m = gs_idx.len();
    let r_adj: Vec<f64> = gs_idx
        .iter()
        .map(|i| stats[(*i - 1) as usize].abs().powf(gsea_param))
        .collect();
    let nr: f64 = r_adj.iter().sum();
    let r_cum_sum: Vec<f64> = if nr == 0.0 {
        gs_idx
            .iter()
            .enumerate()
            .map(|(i, _)| i as f64 / r_adj.len() as f64)
            .collect()
    } else {
        cumsum(&r_adj).iter().map(|x| x / nr).collect()
    };
    let top_tmp: Vec<f64> = gs_idx
        .iter()
        .enumerate()
        .map(|(i, x)| (*x as f64 - (i + 1) as f64) / (n as f64 - m as f64))
        .collect();
    let tops: Vec<f64> = r_cum_sum
        .iter()
        .zip(top_tmp.iter())
        .map(|(x1, x2)| x1 - x2)
        .collect();
    let bottoms: Vec<f64> = if nr == 0.0 {
        tops.iter().map(|x| x - (1.0 / m as f64)).collect()
    } else {
        tops.iter()
            .zip(r_adj.iter())
            .map(|(top, adj)| top - (adj / nr))
            .collect()
    };
    let max_p = array_max(&tops);
    let min_p = array_min(&bottoms);

    let gene_stat = if max_p == -min_p {
        0.0
    } else if max_p > -min_p {
        max_p
    } else {
        min_p
    };
    let leading_edge = if return_leading_edge {
        if max_p > -min_p {
            let max_idx = bottoms
                .iter()
                .enumerate()
                .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
                .map(|(idx, _)| idx)
                .unwrap_or(0);

            gs_idx.iter().take(max_idx + 1).cloned().collect()
        } else if max_p < -min_p {
            let min_idx = bottoms
                .iter()
                .enumerate()
                .min_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
                .map(|(idx, _)| idx)
                .unwrap_or(0);
            gs_idx.iter().skip(min_idx).cloned().rev().collect()
        } else {
            Vec::new()
        }
    } else {
        Vec::new()
    };

    list!(es = gene_stat, leading_edge = leading_edge)
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
/// @param seed Integer For
///
/// @return List with the following elements
/// \itemize{
///     \item le_es ...
///     \item ge_es ...
///     \item le_zero ...
///     \item ge_zero ...
///     \item le_zero_sum ...
///     \item ge_zero_sum ...
/// }
///
/// @export
#[extendr]
pub fn rs_calc_gsea_stat_cumulative_batch(
    stats: &[f64],
    pathway_scores: &[f64],
    pathway_sizes: &[i32],
    iters: usize,
    gsea_param: f64,
    seed: u64,
) -> extendr_api::Result<List> {
    // Convert indices from R - keep as 1-based for algorithm
    let pathway_sizes: Vec<usize> = pathway_sizes.iter().map(|x| *x as usize).collect();

    let res = calc_gsea_stat_cumulative_batch(
        stats,
        pathway_scores,
        &pathway_sizes,
        iters,
        gsea_param,
        seed,
    )?;

    Ok(list!(
        le_es = res.le_es,
        ge_es = res.ge_es,
        le_zero = res.le_zero,
        ge_zero = res.ge_zero,
        le_zero_sum = res.le_zero_sum,
        ge_zero_sum = res.ge_zero_sum
    ))
}

extendr_module! {
    mod fun_fgsea;
    fn rs_get_gs_indices;
    fn rs_calc_es;
    fn rs_gsea_traditional;
    fn rs_calc_gsea_stats;
    fn rs_calc_gsea_stat_cumulative_batch;
}
