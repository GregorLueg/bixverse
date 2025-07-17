use extendr_api::prelude::*;

use rayon::prelude::*;
use rustc_hash::FxBuildHasher;
use rustc_hash::FxHashSet;

use crate::helpers::geom_elim::*;
use crate::helpers::hypergeom::*;
use crate::utils::general::flatten_vector;
use crate::utils::r_rust_interface::{r_list_to_hash_vec, r_list_to_str_vec};
use crate::utils::utils_stats::calc_fdr;

/// Run a single hypergeometric test.
///
/// @description Given a set of target genes, this is a Rust implementation of
/// an hypergeometric test testing for overenrichment of the target genes in the
/// gene sets. WARNING! Incorrect use can cause kernel crashes. Wrapper around
/// the Rust functions with type checks are provided in the package.
///
/// @param target_genes String vector. Represents the target gene set.
/// @param gene_sets List. Contains the strings that represent the gene sets to
/// test against.
/// @param gene_universe String vector. The features representing the gene universe
/// from which the target genes and gene sets are sampled from.
/// @param min_overlap Optional integer. Shall a filter be applied on the minimum of
/// overlappign genes.
/// @param fdr_threshold Optional float. Shall a filter be applied for the maximum
/// tolerated FDR.
///
/// @return A list containing:
///  \itemize{
///   \item pvals - The p-values from the hypergeometric test
///   \item odds_ratios - The calculated odds ratios
///   \item hits - The size of the overlap
///   \item gene_set_lengths - The length of the gene sets.
///   \item fdr - The FDR calculated across the gene sets.
///   \item to_keep - Indices of the gene sets that passed (optional) thresholds.
/// }
///
/// @export
#[extendr]
fn rs_hypergeom_test(
    target_genes: Vec<String>,
    gene_sets: List,
    gene_universe: Vec<String>,
    min_overlap: Option<usize>,
    fdr_threshold: Option<f64>,
) -> extendr_api::Result<List> {
    let gene_sets = r_list_to_hash_vec(gene_sets)?;

    let mut target_genes_set =
        FxHashSet::with_capacity_and_hasher(target_genes.len(), FxBuildHasher);
    for s in target_genes {
        target_genes_set.insert(s);
    }

    let res: HypergeomResult = hypergeom_helper(&target_genes_set, &gene_sets, &gene_universe);
    let (filtered_res, to_keep) = filter_gse_results(res, min_overlap, fdr_threshold);

    Ok(list!(
        pvals = filtered_res.pval,
        odds_ratios = filtered_res.odds_ratio,
        hits = filtered_res.hits,
        gene_set_lengths = filtered_res.gs_length,
        fdr = filtered_res.fdr,
        to_keep = to_keep,
    ))
}

/// Run a hypergeometric test over a list of target genes
///
/// @description Given a list of target gene sets, this function will test for
/// each of the individual target genes the hypergeoemetric enrichment against
/// the specified gene sets. WARNING! Incorrect use can cause kernel crashes.
/// Wrapper around the Rust functions with type checks are provided in the
/// package.
///
/// @param target_genes_list A character vector representing the target gene set.
/// @param gene_sets A list of strings that represent the gene sets to test
/// against.
/// @param gene_universe A character vector representing the gene universe from
/// which the target genes and gene sets are sampled from.
///
/// @return A list containing:
///  \itemize{
///   \item pvals - The p-values from the hypergeometric test.
///   \item fdr - The FDRs for each target gene calculated across all gene sets.
///   \item odds ratios - The calculated odds ratios
///   \item hits - The size of the overlap between the target gene set and individual
///   gene sets.
///   \item gene_set_lengths - The length of the gene sets.
///   \item to_keep - indices
/// }
///
/// @export
#[extendr]
fn rs_hypergeom_test_list(
    target_genes_list: List,
    gene_sets: List,
    gene_universe: Vec<String>,
    min_overlap: Option<usize>,
    fdr_threshold: Option<f64>,
) -> extendr_api::Result<List> {
    let gene_sets = r_list_to_hash_vec(gene_sets)?;
    let target_genes_list = r_list_to_hash_vec(target_genes_list)?;

    let (hypergeom_results, indices_vectors): (Vec<HypergeomResult>, Vec<Vec<usize>>) =
        target_genes_list
            .par_iter()
            .map(|x_i| {
                let res_i: HypergeomResult = hypergeom_helper(x_i, &gene_sets, &gene_universe);
                filter_gse_results(res_i, min_overlap, fdr_threshold)
            })
            .unzip();

    let tests_passed: Vec<usize> = indices_vectors.iter().map(|x| x.len()).collect();

    let mut pvals = Vec::with_capacity(hypergeom_results.len());
    let mut fdr = Vec::with_capacity(hypergeom_results.len());
    let mut odds_ratio = Vec::with_capacity(hypergeom_results.len());
    let mut hits = Vec::with_capacity(hypergeom_results.len());
    let mut gs_lengths = Vec::with_capacity(hypergeom_results.len());

    for res in hypergeom_results {
        pvals.push(res.pval);
        fdr.push(res.fdr);
        odds_ratio.push(res.odds_ratio);
        hits.push(res.hits);
        gs_lengths.push(res.gs_length);
    }

    let pvals = flatten_vector(pvals);
    let fdr = flatten_vector(fdr);
    let odds_ratio = flatten_vector(odds_ratio);
    let hits = flatten_vector(hits);
    let gs_lengths = flatten_vector(gs_lengths);
    let to_keep = flatten_vector(indices_vectors);

    Ok(list!(
        pvals = pvals,
        fdr = fdr,
        odds_ratios = odds_ratio,
        hits = hits,
        gene_set_lengths = gs_lengths,
        to_keep = to_keep,
        tests_passed = tests_passed
    ))
}

/// Run hypergeometric enrichment over the gene ontology
///
/// @description This function implements a Rust version of the gene ontology
/// enrichment with elimination: the starting point are the leafs of the
/// ontology and hypergeometric tests will first conducted there. Should the
/// hypergeometric test p-value be below a certain threshold, the genes of that
/// gene ontology term will be removed from all ancestors. WARNING! Incorrect
/// use can cause kernel crashes. Wrapper around the Rust functions with type
/// checks are provided in the package.
///
/// @param target_genes A character vector representing the target gene set.
/// @param levels A character vector representing the levels to iterate through.
/// The order will be the one the iterations are happening in.
/// @param go_obj The gene_ontology_data S7 class. See [bixverse::gene_ontology_data()].
/// @param gene_universe_length The length of the gene universe.
/// @param min_genes number of minimum genes for the gene ontology term to be
/// tested.
/// @param elim_threshold p-value below which the elimination procedure shall be
/// applied to the ancestors.
/// @param min_overlap Optional minimum overlap threshold.
/// @param fdr_threshold Optional fdr threshold.
///
/// @return A list containing:
///  \itemize{
///   \item go_ids - The gene ontology identifier.
///   \item pvals - The calculated odds ratios.
///   \item odds_ratios - The calculated odds ratios.
///   \item overlap - The size of the overlap.
///   \item gene_set_lengths - The length of the gene sets.
/// }
///
/// @export
#[allow(clippy::too_many_arguments)]
#[extendr]
fn rs_gse_geom_elim(
    target_genes: Vec<String>,
    levels: Vec<String>,
    go_obj: Robj,
    gene_universe_length: usize,
    min_genes: usize,
    elim_threshold: f64,
    min_overlap: Option<usize>,
    fdr_threshold: Option<f64>,
) -> extendr_api::Result<List> {
    let (go_to_gene, ancestors_map, levels_map) = prepare_go_data(go_obj)?;

    let mut go_obj = GeneOntology::new(go_to_gene, &ancestors_map, &levels_map);

    let mut go_ids: Vec<Vec<String>> = Vec::with_capacity(levels.len());
    let mut pvals: Vec<Vec<f64>> = Vec::with_capacity(levels.len());
    let mut odds_ratios: Vec<Vec<f64>> = Vec::with_capacity(levels.len());
    let mut hits: Vec<Vec<usize>> = Vec::with_capacity(levels.len());
    let mut gene_set_lengths: Vec<Vec<usize>> = Vec::with_capacity(levels.len());

    let target_set: FxHashSet<String> = target_genes.iter().cloned().collect();

    for level in levels.iter() {
        let level_res: GoElimLevelResults = process_ontology_level(
            &target_set,
            level,
            &mut go_obj,
            min_genes,
            gene_universe_length,
            elim_threshold,
        );
        go_ids.push(level_res.go_ids);
        pvals.push(level_res.pvals);
        odds_ratios.push(level_res.odds_ratios);
        hits.push(level_res.hits);
        gene_set_lengths.push(level_res.gene_set_lengths);
    }

    let go_ids: Vec<_> = flatten_vector(go_ids);
    let pvals: Vec<_> = flatten_vector(pvals);
    let odds_ratios: Vec<_> = flatten_vector(odds_ratios);
    let hits: Vec<_> = flatten_vector(hits);
    let gene_set_lengths: Vec<_> = flatten_vector(gene_set_lengths);
    let fdr = calc_fdr(&pvals);

    let filtered_res: GoElimFinalResults = filter_go_res(
        go_ids,
        pvals,
        fdr,
        odds_ratios,
        hits,
        gene_set_lengths,
        min_overlap,
        fdr_threshold,
    );

    Ok(list!(
        go_ids = filtered_res.go_ids,
        pvals = filtered_res.pvals,
        fdr = filtered_res.fdr,
        odds_ratios = filtered_res.odds_ratios,
        hits = filtered_res.hits,
        gene_set_lengths = filtered_res.gs_length
    ))
}

/// Run hypergeometric enrichment a list of target genes over the gene ontology
///
/// This function implements a Rust version of the gene ontology enrichment with
/// elimination: the starting point are the leafs of the ontology and
/// hypergeometric tests will first conducted there. Should the hypergeometric
/// test p-value be below a certain threshold, the genes of that gene ontology
/// term will be removed from all ancestors. This function is designed to
/// leverage Rust-based threading for parallel processing of a list of target
/// genes. WARNING! Incorrect use can cause kernel crashes. Wrapper around the
/// Rust functions with type checks are provided in the package.
///
/// @param target_genes_list A list of target genes against which to run the
/// method.
/// @param levels A character vector representing the levels to iterate through.
/// The order will be the one the iterations are happening in.
/// @param go_obj The gene_ontology_data S7 class. See [bixverse::gene_ontology_data()].
/// @param gene_universe_length The length of the gene universe.
/// @param min_genes number of minimum genes for the gene ontology term to be
/// tested.
/// @param elim_threshold p-value below which the elimination procedure shall
/// be applied to the ancestors.
/// @param min_overlap Optional minimum overlap threshold.
/// @param fdr_threshold Optional fdr threshold.
///
/// @return A list containing:
///  \itemize{
///   \item go_ids - The gene ontology identifier.
///   \item pvals - The calculated odds ratios.
///   \item fdrs - The calculated fdrs.
///   \item odds_ratios - The calculated odds ratios.
///   \item overlap - The size of the overlap.
///   \item gene_set_lengths - The length of the gene sets.
///   \item no_test - The number of tests for that target set that passed the
///   thresholds.
/// }
///
/// @export
#[allow(clippy::too_many_arguments)]
#[extendr]
fn rs_gse_geom_elim_list(
    target_genes_list: List,
    levels: Vec<String>,
    go_obj: Robj,
    gene_universe_length: usize,
    min_genes: usize,
    elim_threshold: f64,
    min_overlap: Option<usize>,
    fdr_threshold: Option<f64>,
) -> extendr_api::Result<List> {
    // Prepare various variables
    let target_genes_list = r_list_to_str_vec(target_genes_list)?;

    // Prepare the data
    let go_data = prepare_go_data(go_obj)?;

    let (geom_res, test_passed): (Vec<GoElimFinalResults>, Vec<usize>) = target_genes_list
        .par_iter()
        .map(|targets| {
            // Create necessary mutables
            let mut go_obj = GeneOntology::new(go_data.0.clone(), &go_data.1, &go_data.2);

            let mut go_ids: Vec<Vec<String>> = Vec::with_capacity(levels.len());
            let mut pvals: Vec<Vec<f64>> = Vec::with_capacity(levels.len());
            let mut odds_ratios: Vec<Vec<f64>> = Vec::with_capacity(levels.len());
            let mut hits: Vec<Vec<usize>> = Vec::with_capacity(levels.len());
            let mut gene_set_lengths: Vec<Vec<usize>> = Vec::with_capacity(levels.len());

            let target_set: FxHashSet<String> = targets.iter().cloned().collect();

            // Iterate over the levels
            for level in levels.iter() {
                let level_res: GoElimLevelResults = process_ontology_level(
                    &target_set,
                    level,
                    &mut go_obj,
                    min_genes,
                    gene_universe_length,
                    elim_threshold,
                );

                go_ids.push(level_res.go_ids);
                pvals.push(level_res.pvals);
                odds_ratios.push(level_res.odds_ratios);
                hits.push(level_res.hits);
                gene_set_lengths.push(level_res.gene_set_lengths);
            }

            // Flatten the vectors
            let go_ids: Vec<_> = flatten_vector(go_ids);
            let pvals: Vec<_> = flatten_vector(pvals);
            let odds_ratios: Vec<_> = flatten_vector(odds_ratios);
            let hits: Vec<_> = flatten_vector(hits);
            let gene_set_lengths: Vec<_> = flatten_vector(gene_set_lengths);
            let fdr = calc_fdr(&pvals);

            let filtered_res: GoElimFinalResults = filter_go_res(
                go_ids,
                pvals,
                fdr,
                odds_ratios,
                hits,
                gene_set_lengths,
                min_overlap,
                fdr_threshold,
            );

            let passed_tests = filtered_res.go_ids.len();

            (filtered_res, passed_tests)
        })
        .unzip();

    let mut go_ids_final: Vec<Vec<_>> = Vec::with_capacity(geom_res.len());
    let mut pvals_final = Vec::with_capacity(geom_res.len());
    let mut fdrs_final = Vec::with_capacity(geom_res.len());
    let mut odds_ratios_final = Vec::with_capacity(geom_res.len());
    let mut hits_final = Vec::with_capacity(geom_res.len());
    let mut gene_set_lengths_final = Vec::with_capacity(geom_res.len());

    for res in geom_res {
        go_ids_final.push(res.go_ids);
        pvals_final.push(res.pvals);
        fdrs_final.push(res.fdr);
        odds_ratios_final.push(res.odds_ratios);
        hits_final.push(res.hits);
        gene_set_lengths_final.push(res.gs_length);
    }

    let go_ids_final: Vec<_> = flatten_vector(go_ids_final);
    let pvals_final: Vec<_> = flatten_vector(pvals_final);
    let fdrs_final: Vec<_> = flatten_vector(fdrs_final);
    let odds_ratios_final: Vec<_> = flatten_vector(odds_ratios_final);
    let hits_final: Vec<_> = flatten_vector(hits_final);
    let gene_set_lengths_final: Vec<_> = flatten_vector(gene_set_lengths_final);

    Ok(list!(
        go_ids = go_ids_final,
        pvals = pvals_final,
        fdr = fdrs_final,
        odds_ratios = odds_ratios_final,
        hits = hits_final,
        gene_set_lengths = gene_set_lengths_final,
        no_test = test_passed
    ))
}

extendr_module! {
    mod r_hypergeom;
    fn rs_hypergeom_test;
    fn rs_hypergeom_test_list;
    fn rs_gse_geom_elim;
    fn rs_gse_geom_elim_list;
}
