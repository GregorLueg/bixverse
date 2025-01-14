use extendr_api::prelude::*;
use statrs::distribution::Hypergeometric;
use statrs::distribution::DiscreteCDF;
use std::collections::HashSet;


/// Calculate the p-value of a hypergeometric test.
/// @export
#[extendr]
fn hypergeom_pval(
  q: u64, 
  m: u64, 
  n: u64, 
  k: u64
) -> f64 {
  let pval = if q == 0 {
    // Special case of no hits. Due to being -1 here, the p-value returns as 0.
    1.0
  } else {
    let population = m + n;
    let successes = m;     
    let draws = k;  
    let dist: Hypergeometric = Hypergeometric::new(
      population, 
      successes, 
      draws
    )
    .unwrap();
    1.0 - dist.cdf(q - 1)
  };

  return pval
}

/// Calculate the odds ratio
/// @export
#[extendr]
fn hypergeom_odds_ratio(
  a1_b1: u64,
  a0_b1: u64,
  a1_b0: u64,
  a0_b0: u64,
) -> f64 {
  let odds_ratio = (a1_b1 as f64 / a0_b1 as f64 ) / (a1_b0 as f64 / a0_b0 as f64);

  return odds_ratio
}

/// Create an integer vector of overlaps
/// @export
#[extendr]
fn count_hits(
  gene_set_list: List,
  target_genes: Vec<String>
) -> Vec<u64>{
  let target_genes_str: Vec<&str> = target_genes
    .iter()
    .map(|s| &**s)
    .collect();
  let target_genes_hash: HashSet<_> = target_genes_str
    .into_iter()
    .collect(); 

  let hits = gene_set_list
    .into_iter()
    .map(|(_, xi)| {
      let x_vec = xi
        .as_str_vector()
        .unwrap();
      let x_hash: HashSet<_> = x_vec
        .into_iter()
        .collect();
      let intersection = x_hash
        .intersection(&target_genes_hash)
        .count() as u64;
      intersection
    })
    .collect();

  return hits
}

/// Run a single hypergeometric test.
/// @export
#[extendr]
fn hypergeom_test(
  target_genes: Vec<String>,
  gene_sets: List,
  gene_universe: Vec<String>
) -> List {
  let gene_universe_length = gene_universe
    .clone()
    .into_iter()
    .collect::<Vec<_>>().len() as u64;
  let trials = target_genes
    .clone()
    .into_iter()
    .collect::<Vec<_>>().len() as u64;
  let gene_set_lengths: Vec<u64> = gene_sets
    .clone()
    .into_iter()
    .map(|(_, xi)| {
      let len = xi
        .as_str_vector()
        .unwrap()
        .len() as u64;
      len
    })
    .collect();
  let hits = count_hits(gene_sets, target_genes);

  let pvals: Vec<f64> = hits
    .clone()
    .iter()
    .zip(gene_set_lengths.iter())
    .map(|(hit, gene_set_length)| {
      hypergeom_pval(
        *hit, 
        *gene_set_length, 
        gene_universe_length - *gene_set_length, 
        trials
      )
    })
    .collect();
  let odds_ratios: Vec<f64> = hits
    .clone()
    .iter()
    .zip(gene_set_lengths.iter())
    .map(|(hit, gene_set_length)| {
      hypergeom_odds_ratio(
        *hit,
        *gene_set_length - *hit,
        trials - *hit,
        gene_universe_length - *gene_set_length - trials + *hit,
      )
    })
    .collect();

  list!(
    pvals = Doubles::from_values(pvals), 
    odds_ratios = odds_ratios,
    hits = hits,
    gene_set_lengths = gene_set_lengths
  )
}

/// Run hypergeometric tests over a list of target gene sets.
/// @export
#[extendr]
fn hypergeom_test_list(
  target_gene_lists: List,
  gene_sets: List,
  gene_universe: Vec<String>,
) -> List {
  target_gene_lists
    .into_iter()
    .map(|(_, xi)| {
      // Convert xi to Vec<String>
      let target_genes_list_i: Vec<String> = xi
        .clone()
        .as_str_vector()
        .unwrap_or_default() // Handle invalid or empty elements
        .into_iter()
        .map(|s| s.to_string())
        .collect();

        // println!("You guessed: {:?}", target_genes);

      // Call hypergeom_test for each list
      hypergeom_test(
        target_genes_list_i,
        gene_sets.clone(),
        gene_universe.clone(),
      )
    })
    .collect()
}

// Macro to generate exports.
// This ensures exported functions are registered with R.
// See corresponding C code in entrypoint.c.
extendr_module! {
    mod BIXverse;
    fn hypergeom_pval;
    fn hypergeom_odds_ratio;
    fn count_hits;
    fn hypergeom_test;
    fn hypergeom_test_list;
}