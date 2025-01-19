use extendr_api::prelude::*;
use rayon::prelude::*;
use crate::hypergeom_helpers::*;
use crate::geom_elim_helpers::*;
use crate::utils_r_rust::r_list_to_str_vec;


/// Run a single hypergeometric test.
/// 
/// Given a set of target genes, this is a Rust implementation of an hypergeometric test testing for overenrichment
/// of the target genes in the gene sets.
/// 
/// @param target_genes: A character vector representing the target gene set.
/// @param gene_sets: A list of strings that represent the gene sets to test against.
/// @param gene_universe: A character vector representing the gene universe from which the target genes
/// and gene sets are sampled from.
/// 
/// @returns A list with the following elements: pvals, odds ratios, overlap and the length of the gene set.
/// 
/// @export
#[extendr]
fn rs_hypergeom_test(
  target_genes: Vec<String>,
  gene_sets: List,
  gene_universe: Vec<String>
) -> List {
  let gene_sets = r_list_to_str_vec(gene_sets);

  let res = hypergeom_helper(
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
/// @param target_genes: A character vector representing the target gene set.
/// @param gene_sets: A list of strings that represent the gene sets to test against.
/// @param gene_universe: A character vector representing the gene universe from which the target genes
/// and gene sets are sampled from.
/// 
/// @returns A list with the following elements: pvals, odds ratios, overlap and the length of the gene set.
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

  let res: Vec<(Vec<f64>, Vec<f64>, Vec<u64>, Vec<u64>)> = target_genes
    .par_iter()
    .map(|x_i| {
      let res_i: (Vec<f64>, Vec<f64>, Vec<u64>, Vec<u64>) = hypergeom_helper(
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

  let pvals: Vec<_> = pvals
    .into_iter()
    .flatten()
    .collect();
  let odds_ratios: Vec<_> = odds_ratios
    .into_iter()
    .flatten()
    .collect();
  let hits: Vec<_> = hits
    .into_iter()
    .flatten()
    .collect();
  let gene_set_lengths: Vec<_> = gene_set_lengths
    .into_iter()
    .flatten()
    .collect();
  
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
/// Should the hypergeometric test
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
  // Get the levels
  let level_ids: Vec<String> = levels
    .clone()
    .iter()
    .map(|(n, _)| {
      n.to_string()
    })
    .collect();

  let mut go_obj = generate_go_structure(
    go_to_genes,
    ancestors,
    levels
  );
  
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

    go_ids.push(level_res.0);
    pvals.push(level_res.1);
    odds_ratios.push(level_res.2);
    hits.push(level_res.3);
    gene_set_lengths.push(level_res.4);
  }
  
  let go_ids: Vec<_> = go_ids
    .into_iter()
    .flatten()
    .collect();
  let pvals: Vec<_> = pvals
    .into_iter()
    .flatten()
    .collect();
  let odds_ratios: Vec<_> = odds_ratios
    .into_iter()
    .flatten()
    .collect();
  let hits: Vec<_> = hits
    .into_iter()
    .flatten()
    .collect();
  let gene_set_lengths: Vec<_> = gene_set_lengths
    .into_iter()
    .flatten()
    .collect();

  list!(
    go_ids = go_ids,
    pvals = pvals, 
    odds_ratios = odds_ratios,
    hits = hits,
    gene_set_lengths = gene_set_lengths
  )
}


extendr_module! {
    mod hypergeom_fun;
    fn rs_hypergeom_test;
    fn rs_hypergeom_test_list;
    fn rs_gse_geom_elim;
}