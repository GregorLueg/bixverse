use statrs::distribution::Hypergeometric;
use statrs::distribution::DiscreteCDF;
use std::collections::HashSet;


/// Calculate the p-value of a hypergeometric test.
pub fn hypergeom_pval(
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

/// Calculate odds ratios
pub fn hypergeom_odds_ratio(
  a1_b1: u64,
  a0_b1: u64,
  a1_b0: u64,
  a0_b0: u64,
) -> f64 {
  let odds_ratio = (a1_b1 as f64 / a0_b1 as f64 ) / (a1_b0 as f64 / a0_b0 as f64);

  return odds_ratio
}

/// Count the number of hits for the hypergeometric tests
pub fn count_hits(
  gene_set_list: &Vec<Vec<String>>,
  target_genes: &Vec<String>
) -> Vec<u64> {
  let target_genes_hash: HashSet<_> = target_genes
    .into_iter()
    .collect();
  let hits: Vec<u64> = gene_set_list
    .into_iter()
    .map(|s| {
      let s_hash: HashSet<_> = s
        .into_iter()
        .collect();
      let intersection = s_hash
        .intersection(&target_genes_hash)
        .count() as u64;
      intersection
    })
    .collect();
  hits
}

/// Helper function to generate the 
pub fn hypergeom_helper(
  target_genes: &Vec<String>,
  gene_sets: &Vec<Vec<String>>,
  gene_universe: &Vec<String>
) -> (Vec<f64>, Vec<f64>, Vec<u64>, Vec<u64>){
  let gene_universe_length = gene_universe
    .into_iter()
    .collect::<Vec<_>>()
    .len() as u64;
  let trials = target_genes
    .into_iter()
    .collect::<Vec<_>>()
    .len() as u64;
  let gene_set_lengths = gene_sets
    .into_iter()
    .map(|s| {
      s.len() as u64
    })
    .collect::<Vec<u64>>();
  let hits = count_hits(gene_sets, target_genes);
  let pvals: Vec<f64> = hits
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

  (pvals, odds_ratios, hits, gene_set_lengths)
}