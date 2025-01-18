use extendr_api::prelude::*;

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