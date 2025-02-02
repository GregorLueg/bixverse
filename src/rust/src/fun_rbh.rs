use extendr_api::prelude::*;
use rayon::prelude::*; 

use crate::utils_r_rust::{NestedHashMap, r_nested_list_to_rust};
use crate::helpers_rbh::*;
use crate::utils_rust::flatten_vector;

/// Structure to store the RBH results.
pub struct RbhResult {
    pub origin: String,
    pub target: String,
    pub origin_modules: Vec<String>,
    pub target_modules: Vec<String>,
    pub similarities: Vec<f64>
}


/// @export
#[extendr]
fn rs_rbh_sets(
  module_list: List,
  similiarity_index: bool,
  min_similarity: f64,
  debug: bool
) -> List {
  let module_list: NestedHashMap = r_nested_list_to_rust(module_list);
  // Pull out all the keys
  let origins: Vec<String>  = module_list
    .clone()
    .into_keys()
    .collect();

  // Sliding window approach to generate tuples to iterate over in parallel with Rayon
  // Likely not super memmory efficient, but does the job and allows for subsequent parallelisation.
  let origins_split: Vec<(String, Vec<String>)> = origins
    .iter()
    .enumerate()
    .map(|(i, first)| {
      let rest: Vec<String> = origins.iter().skip(i + 1).map(|s| s.to_string()).collect();
      (first.to_string(), rest)
    })
    .take_while(|(_, rest)| !rest.is_empty())
    .collect();

  let rbh_results: Vec<Vec<RbhResult>> = origins_split
    .par_iter()
    .map(|(origin_module, target_modules)|{
      // Parallel iterator starts here
      let origin_module_data = module_list
        .get(origin_module)
        .unwrap();

      // Iterate over the remaining target modules
      let res: Vec<RbhResult> =   target_modules
        .iter()
        .map(|target| {
          let target_module_data = module_list
            .get(target)
            .unwrap();

          let rbh_res = calculate_rbh_set(
            origin_module_data,
            target_module_data,
            similiarity_index,
            min_similarity,
            debug
          );

          let mut origin_modules = Vec::new();
          let mut target_modules = Vec::new();
          let mut similarities = Vec::new();

          for (origin, target, similarity) in rbh_res {
            origin_modules.push(origin);
            target_modules.push(target);
            similarities.push(similarity)
          }

          RbhResult{
            origin: origin_module.to_string(),
            target: target.to_string(),
            origin_modules,
            target_modules,
            similarities
          }
        })
        .collect();

      res
    })
    .collect();

  let rbh_results_unnested: Vec<_> = rbh_results
    .iter()
    .flatten()
    .collect();

  let mut origin = Vec::new();
  let mut target = Vec::new();
  let mut comparisons = Vec::new();
  let mut origin_modules = Vec::new();
  let mut target_modules = Vec::new();
  let mut similarity = Vec::new();

  for module in rbh_results_unnested {
    origin.push(module.origin.clone());
    target.push(module.target.clone());
    origin_modules.push(module.origin_modules.clone());
    target_modules.push(module.target_modules.clone());
    similarity.push(module.similarities.clone());
    comparisons.push(module.similarities.len());
  }
  
  let origin_modules: Vec<_> = flatten_vector(origin_modules);
  let target_modules: Vec<_> = flatten_vector(target_modules);
  let similarity: Vec<_> = flatten_vector(similarity);


  list!(
    origin = origin,
    target = target,
    comparisons = comparisons,
    origin_modules = origin_modules,
    target_modules = target_modules,
    similarity = similarity
  )
}

extendr_module! {
    mod fun_rbh;
    fn rs_rbh_sets;
}