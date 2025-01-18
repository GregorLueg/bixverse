mod hypergeom_helpers;
mod utils_r_rust;

use extendr_api::prelude::*;
use rayon::prelude::*;
use hypergeom_helpers::*;
use utils_r_rust::r_list_to_str_vec;

/// Run a single hypergeometric test.
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

/// Run a single hypergeometric test.
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


extendr_module! {
    mod BIXverse;
    fn rs_hypergeom_test;
    fn rs_hypergeom_test_list;
}