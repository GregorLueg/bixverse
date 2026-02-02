use extendr_api::prelude::*;
use faer::{Mat, MatRef};
use rustc_hash::{FxBuildHasher, FxHashMap, FxHashSet};
use std::ops::{Add, Mul};

use crate::core::data::sparse_structures::*;
use crate::utils::traits::FaerRType;

//////////////////
// Type aliases //
//////////////////

/// Type alias for a double nested HashMap
pub type NestedHashMap = FxHashMap<String, FxHashMap<String, FxHashSet<String>>>;

////////////
// Errors //
////////////

/// Error handling for named numeric conversion
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

/// Transforms a Robj List into a Hashmap with HashSet values
///
/// This function assumes that the R list contains string vector!
///
/// ### Params
///
/// * `r_list` - R list that has names and contains string vectors.
///
/// ### Returns
///
/// A HashMap with as keys the names of the list and values as HashSets.
pub fn r_list_to_hashmap_set(
    r_list: List,
) -> extendr_api::Result<FxHashMap<String, FxHashSet<String>>> {
    let mut result = FxHashMap::with_capacity_and_hasher(r_list.len(), FxBuildHasher);

    for (n, s) in r_list {
        let s_vec = s.as_string_vector().ok_or_else(|| {
            Error::Other(format!(
                "Failed to convert value for key '{}' to string vector",
                n
            ))
        })?;
        let mut s_hash = FxHashSet::with_capacity_and_hasher(s_vec.len(), FxBuildHasher);
        for item in s_vec {
            s_hash.insert(item);
        }
        result.insert(n.to_string(), s_hash);
    }

    Ok(result)
}

/// Transforms an Robj nested list into a nested HashMap containing further HashMap
///
/// A helper that generates a nested HashMap from a nested R list.
///
/// ### Params
///
/// * `r_nested_list` - A named R list that contains named lists with String vectors.
///
/// ### Returns
///
/// Returns a `NestedHashMap`
#[allow(dead_code)]
pub fn r_nested_list_to_nested_hashmap(r_nested_list: List) -> extendr_api::Result<NestedHashMap> {
    let mut result = FxHashMap::with_capacity_and_hasher(r_nested_list.len(), FxBuildHasher);
    for (n, obj) in r_nested_list {
        let inner_list = obj.as_list().ok_or_else(|| {
            Error::Other(format!("Failed to convert value for key '{}' to list", n))
        })?;
        let inner_hashmap = r_list_to_hashmap_set(inner_list)?;
        result.insert(n.to_string(), inner_hashmap);
    }
    Ok(result)
}

/////////////
// Vectors //
/////////////

/// Transform an R matrix to a Faer one
///
/// ### Params
///
/// * `x` - The R matrix to transform into a faer MatRef (with `f64`)
///
/// ### Returns
///
/// The faer `MatRef` from the original R matrix.
pub fn r_matrix_to_faer<T>(x: &RMatrix<T>) -> MatRef<'_, T>
where
    T: Copy + Clone,
    extendr_api::Robj: for<'a> extendr_api::AsTypedSlice<'a, T>,
{
    let ncol = x.ncols();
    let nrow = x.nrows();
    let data = x.data();
    MatRef::from_column_major_slice(data, nrow, ncol)
}

/// Transform a faer into an R matrix
///
/// ### Params
///
/// * `x` - faer `MatRef` matrix to transform into an R matrix
///
/// ###
///
/// The R matrix based on the faer matrix.
pub fn faer_to_r_matrix<T: FaerRType>(
    x: faer::MatRef<T>,
) -> extendr_api::RArray<T::RType, [usize; 2]> {
    T::to_r_matrix(x)
}

/// Transform an R matrix into a f32 one
///
/// ### Params
///
/// * `x` - R matrix with f64.
///
/// ### Returns
///
/// A faer Mat with f32
pub fn r_matrix_to_faer_fp32(x: &RMatrix<f64>) -> Mat<f32> {
    let ncol = x.ncols();
    let nrow = x.nrows();
    let data = x.data();
    let data_fp32 = data.iter().map(|x| *x as f32).collect::<Vec<f32>>();
    Mat::from_fn(nrow, ncol, |i, j| data_fp32[i + j * nrow])
}

/// Transform a SparseColumnMatrix to an R list
///
/// ### Params
///
/// * `sparse` - SparseColumnMatrix structure
///
/// ### Returns
///
/// R list with the following slots:
/// * `data` - The values
/// * `row_indices` - The row indices
/// * `col_ptr` - The column pointers
/// * `ncol` - Number of columns
/// * `nrow` - Number of rows
pub fn sparse_data_to_list<T>(sparse: CompressedSparseData<T>) -> List
where
    T: Into<Robj> + Clone + Default + Into<f64> + Sync + Add + PartialEq + Mul,
{
    let data: Vec<f64> = sparse
        .data
        .into_iter()
        .map(|x| {
            let robj: Robj = x.into();
            robj.as_real().unwrap()
        })
        .collect();
    let indptr = sparse
        .indptr
        .iter()
        .map(|x| *x as i32)
        .collect::<Vec<i32>>();
    let indices = sparse
        .indices
        .iter()
        .map(|x| *x as i32)
        .collect::<Vec<i32>>();
    let nrow = sparse.shape.0;
    let ncol = sparse.shape.1;

    list!(
        data = data,
        indptr = indptr,
        indices = indices,
        nrow = nrow,
        ncol = ncol
    )
}

/// Transform an R list storing CSR/C data into CompressedSparseData
///
/// ### Params
///
/// * `r_list` - R list that has the following elements: `indptr`, `indices`,
///   `data`, `nrow`, `ncol` and `format`
///
/// ### Returns
///
/// The CompressedSparseData Rust object with the data
pub fn list_to_sparse_matrix<T>(r_list: List) -> CompressedSparseData<T>
where
    T: Clone + Default + TryFrom<Robj> + Into<u32>,
{
    let r_data = r_list.into_hashmap();

    let indptr: Vec<usize> = r_data
        .get("indptr")
        .and_then(|v| v.as_integer_slice())
        .unwrap()
        .iter()
        .map(|&x| x as usize)
        .collect();

    let indices: Vec<usize> = r_data
        .get("indices")
        .and_then(|v| v.as_integer_slice())
        .unwrap()
        .iter()
        .map(|&x| x as usize)
        .collect();

    let data: Vec<T> = r_data
        .get("data")
        .unwrap()
        .as_real_slice()
        .unwrap()
        .iter()
        .map(|&x| T::try_from(Robj::from(x)).ok().unwrap())
        .collect();

    let nrow = r_data.get("nrow").and_then(|v| v.as_integer()).unwrap() as usize;
    let ncol = r_data.get("ncol").and_then(|v| v.as_integer()).unwrap() as usize;

    let format_str = r_data.get("format").and_then(|v| v.as_str()).unwrap();
    let cs_type = match format_str {
        "csr" => CompressedSparseFormat::Csr,
        "csc" => CompressedSparseFormat::Csc,
        _ => panic!("Unknown format"),
    };

    CompressedSparseData {
        data,
        indices,
        indptr,
        cs_type,
        data_2: None,
        shape: (nrow, ncol),
    }
}
