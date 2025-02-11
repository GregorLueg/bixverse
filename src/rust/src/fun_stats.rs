use extendr_api::prelude::*;
use rand::prelude::*;
// use rand::seq::SliceRandom;
use rayon::prelude::*; 

use crate::utils_stats::split_vector_randomly;
// use std::collections::HashSet;
// use crate::utils_r_rust::r_list_to_str_vec;


// /// Set similarities
// /// 
// /// This function calculates the Jaccard or similarity index between a given 
// /// string vector and a list of other string vectors.
// /// 
// /// @param string The String vector against which to calculate the set similarities.
// /// @param string_list The list of character vectors for which to calculate the set similarities. 
// /// @param similarity_index Shall the similarity index instead of the Jaccard similarity be calculated.
// /// 
// /// @export
// #[extendr]
// fn rs_set_sim_list(
//   string: Vec<String>,
//   string_list: List,
//   similarity_index: bool,
// ) -> Vec<f64> {
//     let string_vec = r_list_to_str_vec(string_list);
//     let string: HashSet<_> = string
//       .iter()
//       .collect();
//     let values: Vec<(u64, u64)> = string_vec
//       .iter()
//       .map(|s| {
//         let s_hash: HashSet<_> = s
//         .iter()
//         .collect();
//         let i = s_hash.intersection(&string).count() as u64;
//         let u = if similarity_index {
//           std::cmp::min(s_hash.len(), string.len()) as u64
//         } else {
//           s_hash.union(&string).count() as u64
//         };
//         (i, u)
//       })
//       .collect();

//     let mut sim: Vec<f64> = Vec::new();

//     for (i, u) in values {
//       let sim_i: f64 = (i as f64) / (u as f64);
//       sim.push(sim_i)
//     }

//     sim
//   }

/// Fast AUC calculation
/// 
/// This function calculates rapidly AUCs based on an approximation.
/// 
/// @param pos_scores The scores of your hits.
/// @param neg_scores The scores of your non-hits.
/// @param iters Number of iterations to run the function for. Recommended size: 10,000.
/// @param random_seed Seed.
/// 
/// @return The AUC.
/// 
/// @export
#[extendr]
fn rs_fast_auc(
  pos_scores: Vec<f64>,
  neg_scores: Vec<f64>,
  iters: usize, 
  seed: u64
) -> f64 {
  let mut rng = StdRng::seed_from_u64(seed);
  let mut count = 0;

  for _ in 0..iters {
    let pos_sample = pos_scores.choose(&mut rng).unwrap();
    let neg_sample = neg_scores.choose(&mut rng).unwrap();
    if pos_sample > neg_sample {
      count += 1;
    }
  }

  count as f64 / iters as f64
}

/// Create random AUCs
/// 
/// This function creates a random set of AUCs based on a score vector and
/// a size of the positive set. This can be used for permutation-based estimation
/// of Z-scores and subsequently p-values.
/// 
/// @param score_vec The overall vector of scores.
/// @param size_pos The size of the hits represented in the score_vec.
/// @param random_iters Number of random AUCs to generate.
/// @param auc_iters Number of random iterations to approximate the AUCs.
/// @param seed Seed.
/// 
/// @return A vector of random AUCs based the score vector and size of the positive set.
/// 
/// @export
#[extendr]
fn rs_create_random_aucs(
  score_vec: Vec<f64>,
  size_pos: usize,
  random_iters: usize,
  auc_iters: usize,
  seed: u64,
) -> Vec<f64> {
  let iter_vec: Vec<usize> = (0..random_iters).collect();

  let random_aucs: Vec<_> = iter_vec
    .par_iter()
    .map(|x|{
      let scores = split_vector_randomly(
        score_vec.clone(), 
        size_pos, 
        *x as u64 + seed
      );

      rs_fast_auc(
        scores.0, 
        scores.1, 
        auc_iters, 
        *x as u64 + 1 + seed
      )
    })
    .collect();

  random_aucs
}


/// Calculate the OT harmonic sum
/// 
/// @param x The numeric vector (should be between 0 and 1) for which to 
/// calculate the harmonic sum
/// 
/// @return Returns the harmonic sum according to the OT calculation.
/// 
/// @export
#[extendr]
fn rs_ot_harmonic_sum(
  mut x: Vec<f64>
) -> f64 {
  x.sort_by(|a, b| b.partial_cmp(a).unwrap());

  let harmonic_sum: f64 = x.iter()
    .enumerate()
    .map(|(i, x)| {
      x / (i + 1).pow(2) as f64
    })
    .sum();

  let max_sum: f64 = vec![1; x.len()].into_iter()
    .enumerate()
    .map(|(i, x)|{
      x as f64 / (i + 1).pow(2) as f64
    })
    .sum();

  harmonic_sum / max_sum
}


extendr_module! {
    mod fun_stats;
    // fn rs_set_sim_list;
    fn rs_fast_auc;
    fn rs_create_random_aucs;
    fn rs_ot_harmonic_sum;
}