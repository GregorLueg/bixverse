use extendr_api::prelude::*;

use std::collections::HashMap;

use crate::helpers_fgsea::*;
use crate::helpers_geom_elim::*;
use crate::utils_r_rust::r_named_vec_data;
use crate::utils_rust::flatten_vector;

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
    let res = calc_gsea_stats(stats, gs_idx, gsea_param, return_leading_edge, true);

    list!(es = res.0, leading_edge = res.1)
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
/// @param seed Integer For reproducibility purposes
///
/// @return List with the following elements
/// \itemize{
///     \item es The enrichment scores for the pathway
///     \item nes The normalised enrichment scores for the pathway
///     \item pvals The p-values for this pathway based on permutation
///     testing
///     \item size The pathway sizes.
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
    )?;

    let gsea_res: GseaResults<'_> =
        calculate_nes_es_pval(pathway_scores, &pathway_sizes, &batch_res);

    Ok(list!(
        es = gsea_res.es,
        nes = gsea_res.nes,
        pvals = gsea_res.pvals,
        size = gsea_res.size
    ))
}

/// Run fgsea simple method for gene ontology with elimination method
///
/// @export
#[allow(clippy::too_many_arguments)]
#[extendr]
fn rs_geom_elim_fgsea(
    stats: Robj,
    levels: Vec<String>,
    go_obj: Robj,
    gsea_param: f64,
    elim_threshold: f64,
    min_size: usize,
    max_size: usize,
    iters: usize,
    seed: u64,
    debug: bool,
) -> extendr_api::Result<List> {
    let vec_data = r_named_vec_data(stats)?;

    let (go_to_gene, ancestors_map, levels_map) = prepare_go_data(go_obj)?;

    let mut go_obj = GeneOntology::new(go_to_gene, &ancestors_map, &levels_map);

    let n = vec_data.1.len();

    let shared_perm =
        create_perm_es_simple(&vec_data.1, gsea_param, iters, max_size, n, seed, true);

    let go_shared_perm = GeneOntologyRandomPerm::new(&shared_perm);

    if debug {
        println!("Shared permutations generated");
    }

    let mut go_ids: Vec<Vec<String>> = Vec::with_capacity(levels.len());
    let mut es: Vec<Vec<f64>> = Vec::with_capacity(levels.len());
    let mut nes: Vec<Vec<Option<f64>>> = Vec::with_capacity(levels.len());
    let mut size: Vec<Vec<usize>> = Vec::with_capacity(levels.len());
    let mut pvals: Vec<Vec<f64>> = Vec::with_capacity(levels.len());
    let mut leading_edge: Vec<Vec<Vec<i32>>> = Vec::with_capacity(levels.len());

    for level in levels {
        if debug {
            println!("I am processing level: {:?}?", level);
        }

        let level_res = process_ontology_level_fgsea_simple(
            &vec_data.1,
            &vec_data.0,
            gsea_param,
            &level,
            &mut go_obj,
            &go_shared_perm,
            min_size,
            max_size,
            elim_threshold,
            debug,
        )?;

        go_ids.push(level_res.go_ids);
        es.push(level_res.es);
        nes.push(level_res.nes);
        size.push(level_res.size);
        pvals.push(level_res.pvals);
        leading_edge.push(level_res.leading_edge);
    }

    let go_ids: Vec<String> = flatten_vector(go_ids);
    let es: Vec<f64> = flatten_vector(es);
    let nes: Vec<Option<f64>> = flatten_vector(nes);
    let size: Vec<usize> = flatten_vector(size);
    let pvals: Vec<f64> = flatten_vector(pvals);
    let leading_edge: Vec<Vec<i32>> = flatten_vector(leading_edge);

    let mut leading_edge_list = List::new(leading_edge.len());

    for (i, x) in leading_edge.iter().enumerate() {
        leading_edge_list.set_elt(i, Robj::from(x.clone()))?;
    }

    Ok(list!(
        go_id = go_ids,
        es = es,
        nes = nes,
        size = size,
        pvals = pvals,
        leading_edge = leading_edge_list
    ))
}

extendr_module! {
    mod fun_fgsea;
    fn rs_calc_es;
    fn rs_get_gs_indices;
    fn rs_calc_gsea_stats;
    fn rs_calc_gsea_stat_cumulative_batch;
    fn rs_calc_gsea_stat_traditional_batch;
    fn rs_geom_elim_fgsea;
}
