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
        n_more_extreme = gsea_res.n_more_extreme,
        size = gsea_res.size
    ))
}

/// Calculates p-values for pre-processed data
///
/// @param es Numerical vector. The enrichment scores of the pathways of that specific size
/// @param stats Named numerical vector. Needs to be sorted. The gene level statistics.
/// @param pathway_size Integer. The size of the pathways to test.
/// @param sample_size Integer. The size of the random gene sets to test against.
/// @param seed Integer. Random seed.
/// @param eps Float. Boundary for calculating the p-value.
/// @param sign Boolean. Bit unclear what this is supposed to do. Original documentation says
/// `This option will be used in future implementations.`, but is used in the function.
///
/// @return List with the following elements:
/// \itemize{
///     \item pvals The pvalues.
///     \item is_cp_ge_half Flag indicating if conditional probability is ≥0.5. Indicates
///     overesimation of the p-values.
/// }
///
/// @export
#[extendr]
fn rs_calc_multi_level(
    es: &[f64],
    stats: Robj,
    pathway_size: usize,
    sample_size: usize,
    seed: u64,
    eps: f64,
    sign: bool,
) -> extendr_api::Result<List> {
    let (_, ranks) = r_named_vec_data(stats)?;

    let res: GseaMultiLevelresults =
        fgsea_multilevel_helper(es, &ranks, pathway_size, sample_size, seed, eps, sign);

    Ok(list!(pvals = res.pvals, is_cp_ge_half = res.is_cp_ge_half))
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

/////////////////////////
// Core elim functions //
/////////////////////////

/// Run fgsea simple method for gene ontology with elimination method
///
/// @param stats Named numerical vector. Needs to be sorted. The gene level statistics.
/// @param levels A character vector representing the levels to iterate through.
/// The order will be the one the iterations are happening in.
/// @param go_obj The gene_ontology_data S7 class. See [bixverse::gene_ontology_data()].
/// @param gsea_param Float. The GSEA parameter. Usually defaults to 1.0.
/// @param elim_threshold p-value below which the elimination procedure shall be
/// applied to the ancestors.
/// @param min_size Minimum size of the gene ontology term for testing.
/// @param max_size Maximum size of the gene ontology term for testing. Setting this
/// parameter to large values will slow the function down.
/// @param iters Integer. Number of random permutations for the fgsea simple method
/// to use
/// @param seed Integer. For reproducibility purposes.
/// @param debug Boolean that will provide additional console information for
/// debugging purposes.
///
/// @return List with the following elements
/// \itemize{
///     \item go_ids The name of the tested gene ontology identifer.
///     \item es The enrichment scores for the pathway
///     \item nes The normalised enrichment scores for the pathway
///     \item size The pathway sizes (after elimination!).
///     \item pvals The p-values for this pathway based on permutation
///     testing
///     \item n_more_extreme Number of times the enrichment score was
///     bigger or smaller than the permutation (pending sign).
///     \item leading_edge A list of the index positions of the leading edge
///     genes for this given GO term.
/// }
///
/// @export
#[allow(clippy::too_many_arguments)]
#[extendr]
fn rs_geom_elim_fgsea_simple(
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

    let stat_name_indices: HashMap<&String, usize> = vec_data
        .0
        .iter()
        .enumerate()
        .map(|(i, name)| (name, i))
        .collect();

    if debug {
        println!("Shared permutations generated");
    }

    let mut go_ids: Vec<Vec<String>> = Vec::with_capacity(levels.len());
    let mut es: Vec<Vec<f64>> = Vec::with_capacity(levels.len());
    let mut nes: Vec<Vec<Option<f64>>> = Vec::with_capacity(levels.len());
    let mut size: Vec<Vec<usize>> = Vec::with_capacity(levels.len());
    let mut pvals: Vec<Vec<f64>> = Vec::with_capacity(levels.len());
    let mut n_more_extreme: Vec<Vec<usize>> = Vec::with_capacity(levels.len());
    let mut leading_edge: Vec<Vec<Vec<i32>>> = Vec::with_capacity(levels.len());

    for level in levels {
        if debug {
            println!("I am processing level: {:?}?", level);
        }

        let level_res = process_ontology_level_fgsea_simple(
            &vec_data.1,
            &stat_name_indices,
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
        n_more_extreme.push(level_res.n_more_extreme);
        leading_edge.push(level_res.leading_edge);
    }

    let go_ids: Vec<String> = flatten_vector(go_ids);
    let es: Vec<f64> = flatten_vector(es);
    let nes: Vec<Option<f64>> = flatten_vector(nes);
    let size: Vec<usize> = flatten_vector(size);
    let pvals: Vec<f64> = flatten_vector(pvals);
    let n_more_extreme: Vec<usize> = flatten_vector(n_more_extreme);
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
        n_more_extreme = n_more_extreme,
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
    fn rs_calc_multi_level;
    fn rs_geom_elim_fgsea_simple;
    fn rs_simple_and_multi_err;
}
