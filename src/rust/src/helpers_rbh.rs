use std::collections::{HashMap, HashSet};
use ndarray::*;

use crate::utils_stats::set_similarity;
use crate::utils_rust::*;

/// RBH containing module name 1, module name 2 and the similarity
pub type RbhTriplet = Vec<(String, String, f64)>;

/// Calculates the reciprocal best hits based on set similarities.
pub fn calculate_rbh_set(
  origin_modules: &HashMap<String, HashSet<String>>,
  target_modules: &HashMap<String, HashSet<String>>,
  similiarity_index: bool,
  min_similarity: f64,
  debug: bool
) -> RbhTriplet {
  let similarities: Vec<Vec<f64>> = origin_modules
    .clone()
    .into_values()
    .map(|module_i| {

      let res_i: Vec<f64> = target_modules
        .clone()
        .into_values()
        .map(|module_l| {
          set_similarity(
            &module_i, 
            &module_l,
            similiarity_index
         ) 
        })
        .collect();

      res_i

    })
    .collect();

  let names_origin: Vec<&String>= origin_modules
    .keys()
    .collect();

  let names_targets: Vec<&String> = target_modules
    .keys()
    .collect();

  let mat_data: Vec<f64> = flatten_vector(similarities);

  let max_sim = array_f64_max(&mat_data);

  if debug {
    println!(
      "The observed max_similarity was {}.\nThe values are {:?}", 
      max_sim,
      mat_data,
    );
  }
  

  let result = if max_sim < min_similarity {
    if debug {
      println!("No similarity passed the threshold.\n\n")
    }
    vec![("NA".to_string(), "NA".to_string(), -1.0)]
  } else {
    let sim_mat= Array::from_shape_vec(
      (names_origin.len(), names_targets.len()), 
      mat_data
    ).unwrap();

    let row_maxima: Vec<f64> = sim_mat
      .rows()
      .into_iter()
      .map(|x| {array_f64_max(&x.to_vec())})
      .collect();

    let col_maxima: Vec<f64> = sim_mat
      .columns()
      .into_iter()
      .map(|x| {array_f64_max(&x.to_vec())})
      .collect();

    let (rows, cols) = sim_mat.dim();

    let mut matching_pairs: Vec<(String, String, f64)> = Vec::new();

    for r in 0..rows {
      for c in 0..cols {
        let value = sim_mat[(r, c)];
        if value == row_maxima[r] && value == col_maxima[c] {
          matching_pairs.push((names_origin[r].to_string(), names_targets[c].to_string(), value));
        }
      }
    }

    if debug {
      println!(
        "A total of {} RBH pairs were identified.\n\n", 
        matching_pairs.len()
      );
    }

    if !matching_pairs.is_empty() {
      matching_pairs
    } else {
      vec![("NA".to_string(), "NA".to_string(), -1.0)]
    }
  };

  result
}