use extendr_api::prelude::*;
use faer::Mat;
use std::collections::{HashMap, HashSet};

/// A double nested HashMap
pub type NestedHashMap = HashMap<String, HashMap<String, HashSet<String>>>;

/// Transforms a Robj List into a Hashmap
pub fn r_list_to_hashmap(r_list: List) -> extendr_api::Result<HashMap<String, Vec<String>>> {
    let tuple_array: Vec<(String, Vec<String>)> = r_list
        .into_iter()
        .map(|(n, s)| {
            let n = n.to_string();
            let s_vec = s.as_string_vector().ok_or_else(|| {
                Error::Other(format!(
                    "Failed to convert value for key '{}' to string vector",
                    n
                ))
            })?;
            Ok((n, s_vec))
        })
        .collect::<extendr_api::Result<_>>()?;
    Ok(tuple_array.into_iter().collect())
}

/// Transforms a Robj List into a Hashmap with the values as Hashset
pub fn r_list_to_hashmap_set(
    r_list: List,
) -> extendr_api::Result<HashMap<String, HashSet<String>>> {
    let tuple_array: Vec<(String, HashSet<String>)> = r_list
        .into_iter()
        .map(|(n, s)| {
            let s_vec = s.as_string_vector().ok_or_else(|| {
                Error::Other(format!(
                    "Failed to convert value for key '{}' to string vector",
                    n
                ))
            })?;
            let s_hash: HashSet<_> = s_vec.into_iter().collect();
            Ok((n.to_string(), s_hash))
        })
        .collect::<extendr_api::Result<_>>()?;

    Ok(tuple_array.into_iter().collect())
}

/// Transforms a Robj List into an array of String arrays.
pub fn r_list_to_str_vec(r_list: List) -> extendr_api::Result<Vec<Vec<String>>> {
    r_list
        .into_iter()
        .map(|(n, s)| {
            s.as_string_vector().ok_or_else(|| {
                Error::Other(format!(
                    "Failed to convert value to string vector at key '{}'",
                    n
                ))
            })
        })
        .collect()
}

/// Transforms a Robj nested list into a nested hashmap
pub fn r_nested_list_to_rust(r_nested_list: List) -> extendr_api::Result<NestedHashMap> {
    let outer_list: Vec<(String, HashMap<String, HashSet<String>>)> = r_nested_list
        .into_iter()
        .map(|(key, value)| {
            let inner_list = value.as_list().ok_or_else(|| {
                Error::Other(format!("Failed to convert value for key '{}' to list", key))
            })?;
            let value = r_list_to_hashmap_set(inner_list)?;
            Ok((key.to_string(), value))
        })
        .collect::<extendr_api::Result<_>>()?;

    Ok(outer_list.into_iter().collect())
}

/// Transform an R matrix to a Faer one
pub fn r_matrix_to_faer(x: &RMatrix<f64>) -> faer::Mat<f64> {
    let ncol = x.ncols();
    let nrow = x.nrows();
    let data = x.data();

    Mat::from_fn(nrow, ncol, |i, j| data[i + j * nrow])
}

/// Transform a faer into an R matrix
pub fn faer_to_r_matrix(x: faer::MatRef<f64>) -> extendr_api::RArray<f64, [usize; 2]> {
    let nrow = x.nrows();
    let ncol = x.ncols();
    RArray::new_matrix(nrow, ncol, |row, column| x[(row, column)])
}
