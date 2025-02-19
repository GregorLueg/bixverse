use rand::prelude::*;
use rand::seq::SliceRandom;
use std::collections::HashSet;
use statrs::distribution::{Normal, ContinuousCDF};

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


/// Calculate the set similarity. Options are Jaccard (similarity_index = False) or the similarity index calculation.
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

pub fn z_scores_to_pval(
  z_scores: &[f64]
) -> Vec<f64> {
    let normal = Normal::new(0.0, 1.0).unwrap();
    z_scores.iter()
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