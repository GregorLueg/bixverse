use extendr_api::prelude::*;
use faer::MatRef;
use rustc_hash::{FxHashMap, FxHashSet};
use std::collections::BTreeMap;

///////////////////
// Maps and Sets //
///////////////////

/// Type alias for a double nested HashMap
pub type NestedHashMap = FxHashMap<String, FxHashMap<String, FxHashSet<String>>>;

/// Type alias for double nested BtreeMap
pub type NestedBtreeMap = BTreeMap<String, BTreeMap<String, FxHashSet<String>>>;

/// Transforms a Robj List into a Hashmap
pub fn r_list_to_hashmap(r_list: List) -> extendr_api::Result<FxHashMap<String, Vec<String>>> {
    let mut result = FxHashMap::default();

    for (n, s) in r_list {
        let s_vec = s.as_string_vector().ok_or_else(|| {
            Error::Other(format!(
                "Failed to convert value for key '{}' to string vector",
                n
            ))
        })?;
        result.insert(n.to_string(), s_vec);
    }

    Ok(result)
}

/// Transforms a Robj List into a Hashmap with the values as Hashset
pub fn r_list_to_hashmap_set(
    r_list: List,
) -> extendr_api::Result<FxHashMap<String, FxHashSet<String>>> {
    let mut result = FxHashMap::default();

    for (n, s) in r_list {
        let s_vec = s.as_string_vector().ok_or_else(|| {
            Error::Other(format!(
                "Failed to convert value for key '{}' to string vector",
                n
            ))
        })?;
        let mut s_hash = FxHashSet::default();
        for item in s_vec {
            s_hash.insert(item);
        }
        result.insert(n.to_string(), s_hash);
    }

    Ok(result)
}

// Transforms an Robj nested list into a nested hashmap
#[allow(dead_code)]
pub fn r_nested_list_to_rust(r_nested_list: List) -> extendr_api::Result<NestedHashMap> {
    let mut result = FxHashMap::default();
    for (n, obj) in r_nested_list {
        let inner_list = obj.as_list().ok_or_else(|| {
            Error::Other(format!("Failed to convert value for key '{}' to list", n))
        })?;
        let inner_hashmap = r_list_to_hashmap_set(inner_list)?;
        result.insert(n.to_string(), inner_hashmap);
    }
    Ok(result)
}

/// Transform a Robj List into a BTreeMap with the values as HashSet
/// Import where ordering of the values matters
pub fn r_list_to_btree_set(
    r_list: List,
) -> extendr_api::Result<BTreeMap<String, FxHashSet<String>>> {
    let mut result = BTreeMap::new();
    for (n, s) in r_list {
        let s_vec = s.as_string_vector().ok_or_else(|| {
            Error::Other(format!(
                "Failed to convert value for key '{}' to string vector",
                n
            ))
        })?;
        let mut s_hash = FxHashSet::default();
        for item in s_vec {
            s_hash.insert(item);
        }
        result.insert(n.to_string(), s_hash);
    }
    Ok(result)
}

/// Transform an Robj nested list into a nested Btreemap
pub fn r_nested_list_to_btree_nest(r_nested_list: List) -> extendr_api::Result<NestedBtreeMap> {
    let mut result = BTreeMap::new();

    for (n, obj) in r_nested_list {
        let inner_list = obj.as_list().ok_or_else(|| {
            Error::Other(format!("Failed to convert value for key '{}' to list", n))
        })?;
        let inner_tree = r_list_to_btree_set(inner_list)?;
        result.insert(n.to_string(), inner_tree);
    }

    Ok(result)
}

/////////////
// Vectors //
/////////////

// Error handling for named numeric conversion
#[derive(Debug)]
pub enum NamedVecError {
    NotNumeric,
    NoNames,
    MissingValues,
}

impl std::fmt::Display for NamedVecError {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            NamedVecError::NotNumeric => write!(f, "Input is not a numeric vector"),
            NamedVecError::NoNames => write!(f, "Vector has no names attribute"),
            NamedVecError::MissingValues => write!(f, "Vector contains missing values"),
        }
    }
}

impl std::error::Error for NamedVecError {}

impl From<NamedVecError> for extendr_api::Error {
    fn from(err: NamedVecError) -> Self {
        extendr_api::Error::Other(err.to_string())
    }
}

/// Type alias for named vectors
pub type NamedNumericVec = (Vec<String>, Vec<f64>);

/// Transforms a Robj List into an array of String arrays.
pub fn r_list_to_str_vec(r_list: List) -> extendr_api::Result<Vec<Vec<String>>> {
    let mut result = Vec::with_capacity(r_list.len());

    for (n, s) in r_list.into_iter() {
        let s_vec = s.as_string_vector().ok_or_else(|| {
            Error::Other(format!(
                "Failed to convert value to string vector at key '{}'",
                n
            ))
        })?;
        result.push(s_vec);
    }

    Ok(result)
}

/// Get the names and numeric values from a named R vector
pub fn r_named_vec_data(named_vec: Robj) -> extendr_api::Result<NamedNumericVec> {
    let values = named_vec
        .as_real_vector()
        .ok_or(NamedVecError::NotNumeric)?;

    let names_attr = named_vec.names().ok_or(NamedVecError::NoNames)?;

    let names: Vec<String> = names_attr.into_iter().map(|s| s.to_string()).collect();

    Ok((names, values))
}

//////////////
// Matrices //
//////////////

/// Transform an R matrix to a Faer one
pub fn r_matrix_to_faer(x: &RMatrix<f64>) -> faer::MatRef<'_, f64> {
    let ncol = x.ncols();
    let nrow = x.nrows();
    let data = x.data();

    MatRef::from_column_major_slice(data, nrow, ncol)
}

/// Transform an R integer matrix to a faer one
pub fn r_matrix_to_faer_i32(x: &RMatrix<i32>) -> faer::MatRef<'_, i32> {
    let ncol = x.ncols();
    let nrow = x.nrows();
    let data = x.data();

    MatRef::from_column_major_slice(data, nrow, ncol)
}

/// Transform a faer into an R matrix
pub fn faer_to_r_matrix(x: faer::MatRef<f64>) -> extendr_api::RArray<f64, [usize; 2]> {
    let nrow = x.nrows();
    let ncol = x.ncols();

    RArray::new_matrix(nrow, ncol, |row, column| x[(row, column)])
}
