use extendr_api::prelude::*;

use rayon::prelude::*;
use rustc_hash::FxBuildHasher;
use rustc_hash::FxHashSet;

use crate::core::enrichment::oea::*;
use crate::utils::general::flatten_vector;
use crate::utils::r_rust_interface::r_list_to_hash_vec;

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
/// @param min_overlap Optional integer. Shall a filter be applied on the minimum of
/// overlappign genes.
/// @param fdr_threshold Optional float. Shall a filter be applied for the maximum
/// tolerated FDR.
///
/// @return A list containing:
///  \itemize{
///   \item pvals - The p-values from the hypergeometric test.
///   \item fdr - The FDRs for each target gene calculated across all gene sets.
///   \item odds ratios - The calculated odds ratios
///   \item hits - The size of the overlap between the target gene set and individual
///   gene sets.
///   \item gene_set_lengths - The length of the gene sets.
///   \item to_keep - Indices of the tests that passed.
///   \item tests_passed - How many tests passed the filter criteria for that target
///   set.
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

extendr_module! {
    mod r_oea;
    fn rs_hypergeom_test;
    fn rs_hypergeom_test_list;
}
