use extendr_api::prelude::*;
use std::collections::{HashMap, HashSet};

pub fn r_list_to_hashmap(
  r_list: List
) -> HashMap<String, Vec<String>> {
  let tuple_array: Vec<_> = r_list
    .into_iter()
    .map(|(n, s)| {
      let n = n.to_string();
      let s_vec = s
        .as_string_vector()
        .unwrap();
    (n, s_vec)
    })
    .collect();
  
    tuple_array.into_iter().collect()
}

pub fn r_list_to_hashmap_set(
  r_list: List,
) -> HashMap<String, HashSet<String>> {
  let tuple_array: Vec<(String, HashSet<String>)> = r_list
    .into_iter()
    .map(|(n, s)| {
      let s_vec = s.as_string_vector().unwrap();
      let s_hash: HashSet<_> = s_vec.into_iter().collect();
      (n.to_string(), s_hash)
    })
    .collect();

  tuple_array.into_iter().collect()
}


/// List of Strings as Robj to Vec<Vec<String>>
pub fn r_list_to_str_vec(
  r_list: List
) -> Vec<Vec<String>> {
  r_list
    .into_iter()
    .map(|(_, s)| {
        let s_vec = s
          .as_string_vector()
          .unwrap();
        s_vec
    })
    .collect()
}