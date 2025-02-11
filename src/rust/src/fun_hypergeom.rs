use extendr_api::prelude::*;
use rayon::prelude::*; 

use crate::helpers_hypergeom::*;
use crate::helpers_geom_elim::*;
use crate::utils_r_rust::r_list_to_str_vec;
use crate::utils_rust::flatten_vector;


/// A type alias that can be returned by par_iter() functions.
type GoElimLevelResultsIter = (Vec<String>, Vec<f64>, Vec<f64>, Vec<u64>, Vec<u64>);

/// Run a single hypergeometric test.
/// 
/// Given a set of target genes, this is a Rust implementation of an hypergeometric test testing for overenrichment
/// of the target genes in the gene sets.
/// 
/// @param target_genes A character vector representing the target gene set.
/// @param gene_sets A list of strings that represent the gene sets to test against.
/// @param gene_universe A character vector representing the gene universe from which the target genes
/// and gene sets are sampled from.
/// 
/// @returns A list with the following elements: 
/// - pvals, the p-values from the hypergeometric test 
/// - odds ratios, the calculated odds ratios
/// - overlap, the size of the overlap,
/// - gene_set_lengths, the length of the gene sets.
/// 
/// @export
#[extendr]
fn rs_hypergeom_test(
  target_genes: Vec<String>,
  gene_sets: List,
  gene_universe: Vec<String>
) -> List {
  let gene_sets = r_list_to_str_vec(gene_sets);

  let res: HypergeomResult = hypergeom_helper(
    &target_genes,
    &gene_sets,
    &gene_universe
  );

  list!(
    pvals = res.0, 
    odds_ratios = res.1,
    hits = res.2,
    gene_set_lengths = res.3
  )
}

/// Run a hypergeometric test over a list of target genes
/// 
/// Given a list of target gene sets, this function will test for each of the individual 
/// 
/// @param target_genes A character vector representing the target gene set.
/// @param gene_sets A list of strings that represent the gene sets to test against.
/// @param gene_universe A character vector representing the gene universe from which the target genes
/// and gene sets are sampled from.
/// 
/// @returns A list with the following elements: 
/// - pvals, the p-values from the hypergeometric test 
/// - odds ratios, the calculated odds ratios
/// - overlap, the size of the overlap,
/// - gene_set_lengths, the length of the gene sets.
/// 
/// @export
#[extendr]
fn rs_hypergeom_test_list(
  target_genes: List,
  gene_sets: List,
  gene_universe: Vec<String>
) -> List {
  let gene_sets = r_list_to_str_vec(gene_sets);
  let target_genes = r_list_to_str_vec(target_genes);

  let res: Vec<HypergeomResult> = target_genes
    .par_iter()
    .map(|x_i| {
      let res_i: HypergeomResult = hypergeom_helper(
        x_i,
        &gene_sets,
        &gene_universe
      );
      res_i
    })
    .collect();

  let mut pvals = Vec::new();
  let mut odds_ratios = Vec::new();
  let mut hits = Vec::new();
  let mut gene_set_lengths = Vec::new();

  for (pval, odds_ratio, hit, gene_set_length) in res {
    pvals.push(pval);
    odds_ratios.push(odds_ratio);
    hits.push(hit);
    gene_set_lengths.push(gene_set_length);
  }

  let pvals: Vec<_> = flatten_vector(pvals);
  let odds_ratios: Vec<_> = flatten_vector(odds_ratios);
  let hits: Vec<_> = flatten_vector(hits);
  let gene_set_lengths: Vec<_> = flatten_vector(gene_set_lengths);
  
  list!(
    pvals = pvals, 
    odds_ratios = odds_ratios,
    hits = hits,
    gene_set_lengths = gene_set_lengths
  )
}

/// Run hypergeometric enrichment over the gene ontology
/// 
/// This function implements a Rust version of the gene ontology enrichment with elimination:
/// the starting point are the leaves of the ontology and hypergeometric tests will first conducted there.
/// Should the hypergeometric test p-value be below a certain threshold, the genes of that gene ontology
/// term will be removed from all ancestors.
/// 
/// @param target_genes A character vector representing the target gene set.
/// @param go_to_genes A named list with the gene identifers as elements and gene ontology identifiers as 
/// names.
/// @param ancestors A named list with the go identifiers of all ancestors as elements and the gene ontology
/// identifiers as names.
/// @param levels A named list with the go identifiers of that ontology level as elements and the level name
/// as names. IMPORTANT! This list needs to be ordered in the right way!
/// @param gene_universe_length The length of the gene universe.
/// @param min_genes number of minimum genes for the gene ontology term to be tested.
/// @param elim_threshold p-value below which the elimination procedure shall be applied to the ancestors.
/// @param debug boolean that will provide additional console information for debugging purposes.
/// 
/// @returns A list with the following elements: 
/// - go_ids, the gene ontology identifier
/// - pvals, the p-values from the hypergeometric test 
/// - odds ratios, the calculated odds ratios
/// - overlap, the size of the overlap,
/// - gene_set_lengths, the length of the gene sets.
/// 
/// @export
#[extendr]
fn rs_gse_geom_elim(
  target_genes: Vec<String>,
  go_to_genes: List,
  ancestors: List,
  levels: List,
  gene_universe_length: u64,
  min_genes: i64,
  elim_threshold: f64,
  debug: bool,
) -> List {
  let level_ids: Vec<String> = levels
    .clone()
    .iter()
    .map(|(n, _)| {
      n.to_string()
    })
    .collect();

  let go_data = prepare_go_data(
    go_to_genes,
    ancestors,
    levels,
  );

  let mut go_obj = GeneOntology {
    go_to_gene: go_data.0,
    ancestors: go_data.1,
    levels: go_data.2,
  };
  
  let mut go_ids: Vec<Vec<String>> = Vec::new();
  let mut pvals: Vec<Vec<f64>> = Vec::new();
  let mut odds_ratios: Vec<Vec<f64>> = Vec::new();
  let mut hits: Vec<Vec<u64>> = Vec::new();
  let mut gene_set_lengths: Vec<Vec<u64>> = Vec::new();

  for level in level_ids.iter() {
    let level_res = process_ontology_level(
      target_genes.clone(),
      level,
      &mut go_obj,
      min_genes,
      gene_universe_length,
      elim_threshold,
      debug
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

  list!(
    go_ids = go_ids,
    pvals = pvals, 
    odds_ratios = odds_ratios,
    hits = hits,
    gene_set_lengths = gene_set_lengths
  )
}


/// Run hypergeometric enrichment a list of target genes over the gene ontology
/// 
/// This function implements a Rust version of the gene ontology enrichment with elimination:
/// the starting point are the leaves of the ontology and hypergeometric tests will first conducted there.
/// Should the hypergeometric test p-value be below a certain threshold, the genes of that gene ontology
/// term will be removed from all ancestors. This function is designed to leverage Rust-based threading
/// for parallel processing of a list of target genes.
/// 
/// @param target_genes_list A list of target genes against which to run the method.
/// @param go_to_genes A named list with the gene identifers as elements and gene ontology identifiers as 
/// names.
/// @param ancestors A named list with the go identifiers of all ancestors as elements and the gene ontology
/// identifiers as names.
/// @param levels A named list with the go identifiers of that ontology level as elements and the level name
/// as names. IMPORTANT! This list needs to be ordered in the right way!
/// @param gene_universe_length The length of the gene universe.
/// @param min_genes number of minimum genes for the gene ontology term to be tested.
/// @param elim_threshold: p-value below which the elimination procedure shall be applied to the ancestors.
/// @param debug boolean that will provide additional console information for debugging purposes.
/// 
/// @returns A list with the following elements: 
/// - go_ids, the gene ontology identifier
/// - pvals, the p-values from the hypergeometric test 
/// - odds ratios, the calculated odds ratios
/// - overlap, the size of the overlap,
/// - gene_set_lengths, the length of the gene sets.
/// - no_test, the number of tests that were conducted against target_gene_list. First element indicates
/// how many values belong to the first target_genes set in the list, etc.
/// 
/// @export
#[extendr]
fn rs_gse_geom_elim_list(
  target_genes_list: List,
  go_to_genes: List,
  ancestors: List,
  levels: List,
  gene_universe_length: u64,
  min_genes: i64,
  elim_threshold: f64,
  debug: bool,
) -> List {
  // Prepare various variables
  let target_genes_list = r_list_to_str_vec(target_genes_list);
  let level_ids: Vec<String> = levels
    .clone()
    .iter()
    .map(|(n, _)| {
      n.to_string()
    })
    .collect();
  // Prepare the data
  let go_data = prepare_go_data(
    go_to_genes,
    ancestors,
    levels,
  );

  let res: Vec<GoElimLevelResultsIter> = target_genes_list
    .par_iter()
    .map(|targets| {
      // Create necessary mutables
      let mut go_obj = GeneOntology {
        go_to_gene: go_data.0.clone(),
        ancestors: go_data.1.clone(),
        levels: go_data.2.clone(),
      };

      let mut go_ids: Vec<Vec<String>> = Vec::new();
      let mut pvals: Vec<Vec<f64>> = Vec::new();
      let mut odds_ratios: Vec<Vec<f64>> = Vec::new();
      let mut hits: Vec<Vec<u64>> = Vec::new();
      let mut gene_set_lengths: Vec<Vec<u64>> = Vec::new();

      // Iterate over the levels
      for level in level_ids.iter() {
        let level_res = process_ontology_level(
          targets.clone(),
          level,
          &mut go_obj,
          min_genes,
          gene_universe_length,
          elim_threshold,
          debug
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
    
      (go_ids, pvals, odds_ratios, hits, gene_set_lengths)
    })
    .collect();

  let mut go_ids_final: Vec<Vec<_>> = Vec::new();
  let mut pvals_final = Vec::new();
  let mut odds_ratios_final = Vec::new();
  let mut hits_final = Vec::new();
  let mut gene_set_lengths_final = Vec::new();
  let mut no_tests = Vec::new();

  for (
    go_ids, 
    pval, 
    odds_ratio, 
    hit, 
    gene_set_length
  ) in res {
    no_tests.push(go_ids.len());
    go_ids_final.push(go_ids);
    pvals_final.push(pval);
    odds_ratios_final.push(odds_ratio);
    hits_final.push(hit);
    gene_set_lengths_final.push(gene_set_length);
  }

  let go_ids_final: Vec<_>= flatten_vector(go_ids_final);
  let pvals_final: Vec<_> = flatten_vector(pvals_final);
  let odds_ratios_final: Vec<_> = flatten_vector(odds_ratios_final);
  let hits_final: Vec<_> = flatten_vector(hits_final);
  let gene_set_lengths_final: Vec<_> = flatten_vector(gene_set_lengths_final);

  list!(
    go_ids = go_ids_final,
    pvals = pvals_final, 
    odds_ratios = odds_ratios_final,
    hits = hits_final,
    gene_set_lengths = gene_set_lengths_final,
    no_test = no_tests
  )
}

extendr_module! {
    mod fun_hypergeom;
    fn rs_hypergeom_test;
    fn rs_hypergeom_test_list;
    fn rs_gse_geom_elim;
    fn rs_gse_geom_elim_list;
}
