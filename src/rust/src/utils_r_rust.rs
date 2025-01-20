use extendr_api::prelude::*;
use std::collections::{HashMap, HashSet};

/// Transforms a Robj List into a Hashmap
/// 
/// #### Arguments
/// * r_list: The R list to transform. The function is expecting Strings as content of the list
/// 
/// #### Returns
/// * A Hashmap with the names of the list as keys and the character vectors stored as String arrays.
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

/// Transforms a Robj List into a Hashmap with the values as Hashset
/// 
/// #### Arguments
/// * r_list: The R list to transform. The function is expecting Strings as content of the list
/// 
/// #### Returns
/// * A Hashmap with the names of the list as keys and the character vectors stored as HashSets.
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


/// Transforms a Robj List into an array of String arrays.
/// 
/// #### Arguments
/// * r_list: The R list to transform. The function is expecting Strings as content of the list
/// 
/// #### Returns
/// * An array of String arrays.
pub fn r_list_to_str_vec(
  r_list: List
) -> Vec<Vec<String>> {
  r_list
    .into_iter()
    .map(|(_, s)| {
        s
          .as_string_vector()
          .unwrap()
    })
    .collect()
}