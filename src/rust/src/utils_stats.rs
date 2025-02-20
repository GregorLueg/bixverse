use rand::prelude::*;
use rand::seq::SliceRandom;
use std::collections::HashSet;
use statrs::distribution::{Normal, ContinuousCDF};
use crate::utils_rust::*;


/// Split a vector randomly into two chunks with one being [..x] and the other [x..]
pub fn split_vector_randomly(
  vec: Vec<f64>, 
  x: usize, 
  seed: u64
) -> (Vec<f64>, Vec<f64>) {
  let mut rng = StdRng::seed_from_u64(seed);
  let mut shuffled = vec.clone();
  shuffled.shuffle(&mut rng);

  let first_set = shuffled[..x].to_vec();
  let second_set = shuffled[x..].to_vec();

  (first_set, second_set)
}


/// Calculate the set similarity. Options are Jaccard (similarity_index = False)
/// or the similarity index calculation.
pub fn set_similarity(
  s_1: &HashSet<String>,
  s_2: &HashSet<String>,
  overlap_coefficient: bool
) -> f64{
  let i = s_1
    .intersection(s_2)
    .count() as u64;
  let u = if overlap_coefficient {
    std::cmp::min(s_1.len(), s_2.len()) as u64
  } else {
    s_1.union(s_2).count() as u64
  };
  i as f64 / u as f64
}

/// Calculate the Hedge's g effect size and its standard error
pub fn hedge_g_effect(
  grp_a: &[f64],
  grp_b: &[f64],
  small_sample_correction: bool
) -> (f64, f64) {
  let len_a = grp_a.len();
  let len_b = grp_b.len();
  let total_n = (len_a + len_b) as f64;
  let mean_a = array_f64_mean(grp_a);
  let mean_b = array_f64_mean(grp_b);
  let var_a = array_f64_var(grp_a);
  let var_b = array_f64_var(grp_b);
  // Calculate the Hedge's G effect
  let pooled_sd = (((len_a - 1) as f64 * var_a + (len_b - 1) as f64 * var_b) / ((len_a - 1) as f64 + (len_b - 1) as f64)).sqrt();
  let effect_size = (mean_a - mean_b) / pooled_sd;
  let effect_size = if small_sample_correction {
    let correction_factor = ((total_n - 3.0) / (total_n - 2.25)) * ((total_n - 2.0) / total_n).sqrt();
    correction_factor * effect_size
  } else {
    effect_size
  };
  let standard_error = (total_n / (len_a as f64 * len_b as f64) + (effect_size.powi(2) / 2.0 * total_n)).sqrt();
  (effect_size, standard_error)
}

/// Transform Z-scores into p-values (assuming normality).
pub fn z_scores_to_pval(
  z_scores: &[f64]
) -> Vec<f64> {
  let normal = Normal::new(0.0, 1.0).unwrap();
  z_scores
    .iter()
    .map(|&z| {
      let abs_z = z.abs();
      2.0 * (1.0 - normal.cdf(abs_z))
    })
    .collect()
}

// /// Calculate the Jaccard or Set similarity over a vector of HashSets.
// pub fn set_similarity_iter(
//   s_1: HashSet<String>,
//   other_s: Vec<HashSet<String>>,
//   similarity_index: bool,
// ) -> Vec<f64>{
//   let sims: Vec<f64>= other_s
//     .into_iter()
//     .map(|hash_i| {
//       set_similarity(hash_i, &s_1, similarity_index)
//     })
//     .collect();

//   sims  
// }