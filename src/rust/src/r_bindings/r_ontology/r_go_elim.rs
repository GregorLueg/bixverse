use extendr_api::prelude::*;
use rayon::prelude::*;
use rustc_hash::{FxHashMap, FxHashSet};

use crate::core::enrichment::gsea::*;
use crate::core::ontology::go_elim::*;
use crate::utils::general::flatten_vector;
use crate::utils::r_rust_interface::*;

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

    let mut go_obj: GeneOntology<'_> = GeneOntology::new(go_to_gene, &ancestors_map, &levels_map);

    let mut go_res: Vec<GoElimLevelResults> = Vec::with_capacity(levels.len());

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
        go_res.push(level_res)
    }

    let filtered_res: GoElimFinalResults = finalise_go_res(&go_res, min_overlap, fdr_threshold);

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
            let mut go_obj: GeneOntology<'_> =
                GeneOntology::new(go_data.0.clone(), &go_data.1, &go_data.2);

            let mut go_res: Vec<GoElimLevelResults> = Vec::with_capacity(levels.len());

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

                go_res.push(level_res);
            }

            let filtered_res: GoElimFinalResults =
                finalise_go_res(&go_res, min_overlap, fdr_threshold);

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

/// Run fgsea simple method for gene ontology with elimination method
///
/// @param stats Named numerical vector. Needs to be sorted. The gene level statistics.
/// @param levels A character vector representing the levels to iterate through.
/// The order will be the one the iterations are happening in.
/// @param go_obj The gene_ontology_data S7 class. See [bixverse::gene_ontology_data()].
/// @param gsea_params List. The GSEA parameters, see [bixverse::params_gsea()]
/// wrapper function. This function generates a list containing:
/// \itemize{
///     \item min_size - Integer. Minimum size for the gene sets.
///     \item max_size - Integer. Maximum size for the gene sets.
///     \item gsea_param - Float. The GSEA parameter. Defaults to `1.0`.
///     \item sample_size - Integer. Number of samples to iterate through for the
///     multi-level implementation of fgsea.
///     \item eps - Float. Boundary for calculating the p-value. Used for the multi-
///     level implementation of fgsea.
/// }
/// @param elim_threshold p-value below which the elimination procedure shall be
/// applied to the ancestors.
/// @param iters Integer. Number of random permutations for the fgsea simple method
/// to use
/// @param seed Integer. For reproducibility purposes.
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
///     \item le_zero Number of times the permutation was less than zero.
///     \item ge_zero Number of times the permutation was greater than zero.
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
    gsea_params: List,
    elim_threshold: f64,
    iters: usize,
    seed: u64,
) -> extendr_api::Result<List> {
    let vec_data = r_named_vec_data(stats)?;

    let (go_to_gene, ancestors_map, levels_map) = prepare_go_data(go_obj)?;

    let gsea_params = prepare_gsea_params(gsea_params);

    let mut go_obj = GeneOntology::new(go_to_gene, &ancestors_map, &levels_map);

    let n = vec_data.1.len();

    let shared_perm = create_perm_es_simple(
        &vec_data.1,
        gsea_params.gsea_param,
        iters,
        gsea_params.max_size,
        n,
        seed,
        true,
    );

    let go_shared_perm = GeneOntologyRandomPerm::new(&shared_perm);

    let stat_name_indices: FxHashMap<&String, usize> = vec_data
        .0
        .iter()
        .enumerate()
        .map(|(i, name)| (name, i))
        .collect();

    let mut go_ids: Vec<Vec<String>> = Vec::with_capacity(levels.len());
    let mut es = Vec::with_capacity(levels.len());
    let mut nes: Vec<Vec<Option<f64>>> = Vec::with_capacity(levels.len());
    let mut size: Vec<Vec<usize>> = Vec::with_capacity(levels.len());
    let mut pvals: Vec<Vec<f64>> = Vec::with_capacity(levels.len());
    let mut n_more_extreme: Vec<Vec<usize>> = Vec::with_capacity(levels.len());
    let mut ge_zero: Vec<Vec<usize>> = Vec::with_capacity(levels.len());
    let mut le_zero: Vec<Vec<usize>> = Vec::with_capacity(levels.len());
    let mut leading_edge: Vec<Vec<Vec<i32>>> = Vec::with_capacity(levels.len());

    for level in levels {
        let level_res = process_ontology_level_fgsea_simple(
            &vec_data.1,
            &stat_name_indices,
            &level,
            &mut go_obj,
            &go_shared_perm,
            &gsea_params,
            elim_threshold,
        )?;

        go_ids.push(level_res.go_ids);
        es.push(level_res.es);
        nes.push(level_res.nes);
        size.push(level_res.size);
        pvals.push(level_res.pvals);
        n_more_extreme.push(level_res.n_more_extreme);
        ge_zero.push(level_res.ge_zero);
        le_zero.push(level_res.le_zero);
        leading_edge.push(level_res.leading_edge);
    }

    let go_ids: Vec<String> = flatten_vector(go_ids);
    let es: Vec<f64> = flatten_vector(es);
    let nes: Vec<Option<f64>> = flatten_vector(nes);
    let size: Vec<usize> = flatten_vector(size);
    let pvals: Vec<f64> = flatten_vector(pvals);
    let n_more_extreme: Vec<usize> = flatten_vector(n_more_extreme);
    let ge_zero: Vec<usize> = flatten_vector(ge_zero);
    let le_zero: Vec<usize> = flatten_vector(le_zero);
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
        ge_zero = ge_zero,
        le_zero = le_zero,
        leading_edge = leading_edge_list
    ))
}

extendr_module! {
    mod r_go_elim;
    fn rs_gse_geom_elim;
    fn rs_gse_geom_elim_list;
    fn rs_geom_elim_fgsea_simple;
}
