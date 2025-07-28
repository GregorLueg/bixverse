use crate::helpers::structs_sparse::SparseColumnMatrix;
use extendr_api::prelude::*;
use faer::{Mat, MatRef};
use rustc_hash::{FxBuildHasher, FxHashMap, FxHashSet};
use std::collections::BTreeMap;

//////////////////
// Type aliases //
//////////////////

/// Type alias for a double nested HashMap
pub type NestedHashMap = FxHashMap<String, FxHashMap<String, FxHashSet<String>>>;

/// Type alias for double nested BtreeMap
pub type NestedBtreeMap = BTreeMap<String, BTreeMap<String, FxHashSet<String>>>;

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

///////////
// Lists //
///////////

/// Transforms a Robj List into a Hashmap
///
/// This function assumes that the R list contains string vector!
///
/// ### Params
///
/// * `r_list` - R list that has names and contains string vectors.
///
/// ### Returns
///
/// A HashMap with as keys the names of the list and values the string vectors.
pub fn r_list_to_hashmap(r_list: List) -> extendr_api::Result<FxHashMap<String, Vec<String>>> {
    let mut result = FxHashMap::with_capacity_and_hasher(r_list.len(), FxBuildHasher);

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

/// Transforms an R list to a vector of HashSets
///
/// ### Params
///
/// * `r_list` - A named R list that contains named lists with String vectors.
///
/// ### Returns
///
/// Returns a Vector of FxHashSets
pub fn r_list_to_hash_vec(r_list: List) -> extendr_api::Result<Vec<FxHashSet<String>>> {
    let mut res = Vec::with_capacity(r_list.len());
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
        res.push(s_hash)
    }

    Ok(res)
}

/// Transform a Robj List into a BTreeMap with the values as HashSet
///
/// Use where ordering of the values matters as the HashMaps have non-deterministic
/// ordering
///
/// ### Params
///
/// * `r_list` - R list that has names and contains string vectors.
///
/// ### Returns
///
/// A BTreeMap with as keys the names of the list and values as HashSets.
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
        let mut s_hash = FxHashSet::with_capacity_and_hasher(s_vec.len(), FxBuildHasher);
        for item in s_vec {
            s_hash.insert(item);
        }
        result.insert(n.to_string(), s_hash);
    }
    Ok(result)
}

/// Transform an Robj nested list into a nested BtreeMap
///
/// A helper that generates a nested BTreeMap from a nested R list.
///
/// ### Params
///
/// * `r_nested_list` - A named R list that contains named lists with String vectors.
///
/// ### Returns
///
/// Returns a `NestedBtreeMap`
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

/// Type alias for named numeric vectors
///
/// ### Fields
///
/// * `0` The names of the vector
/// * `1` The values of the vector
pub type NamedNumericVec = (Vec<String>, Vec<f64>);

/// Transforms a Robj List into an array of String arrays.
///
/// ### Params
///
/// * `r_list` - R list that has names and contains string vectors.
///
/// ### Returns
///
/// A vector of vectors with Strings
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
///
/// ### Params
///
/// * `named_vec` - Robj that represents a named numeric in R
///
/// ### Returns
///
/// The `NamedNumericVec` type alias.
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

/// Structure to store named matrices
///
/// ### Fields
///
/// * `col_names` - A BTreeMap representing the column names and the col indices
/// * `row_names` - A BTreeMap representing the row names and the row indices
/// * `values` - A faer matrix reference representing the matrix values.
#[derive(Clone, Debug)]
pub struct NamedMatrix<'a> {
    pub col_names: BTreeMap<String, usize>,
    pub row_names: BTreeMap<String, usize>,
    pub values: faer::MatRef<'a, f64>,
}

#[allow(dead_code)]
impl<'a> NamedMatrix<'a> {
    /// Generate a new matrix with the feature and sample names stored in the structure
    ///
    /// ### Paramx
    ///
    /// * `x` - An RMatrix.
    ///
    /// ### Returns
    ///
    /// Initialised structure
    pub fn new(x: &'a RMatrix<f64>) -> Self {
        let col_names: BTreeMap<String, usize> = x
            .get_colnames()
            .unwrap()
            .iter()
            .enumerate()
            .map(|(i, s)| (s.to_string(), i))
            .collect();
        let row_names: BTreeMap<String, usize> = x
            .get_rownames()
            .unwrap()
            .iter()
            .enumerate()
            .map(|(i, s)| (s.to_string(), i))
            .collect();
        let mat = r_matrix_to_faer(x);
        NamedMatrix {
            col_names,
            row_names,
            values: mat,
        }
    }

    /// Return a submatrix based on the row names and columns to select.
    ///
    /// If no rows or columns are specified, returns the full matrix. If empty slices are
    /// provided, returns None.
    ///
    /// ### Params
    ///
    /// * `rows_to_select` - Names of the row features
    /// * `cols_to_select` - Names of the column features
    ///
    /// ### Returns
    ///
    /// The faer Matrix with the values from the rows and columns.
    pub fn get_sub_mat(
        &self,
        rows_to_select: Option<&[&str]>,
        cols_to_select: Option<&[&str]>,
    ) -> Option<Mat<f64>> {
        // Determine which rows to select
        let row_indices: Vec<usize> = match rows_to_select {
            None => (0..self.values.nrows()).collect(), // All rows
            Some(rows) => {
                if rows.is_empty() {
                    return None;
                }
                let mut indices = Vec::new();
                for &row_name in rows {
                    if let Some(&index) = self.row_names.get(row_name) {
                        indices.push(index);
                    }
                }
                if indices.is_empty() {
                    return None;
                }
                indices
            }
        };

        // Determine which columns to select
        let col_indices: Vec<usize> = match cols_to_select {
            None => (0..self.values.ncols()).collect(), // All columns
            Some(cols) => {
                if cols.is_empty() {
                    return None;
                }
                let mut indices = Vec::new();
                for &col_name in cols {
                    if let Some(&index) = self.col_names.get(col_name) {
                        indices.push(index);
                    }
                }
                if indices.is_empty() {
                    return None;
                }
                indices
            }
        };

        // Create new matrix by copying values
        let mut result = Mat::<f64>::zeros(row_indices.len(), col_indices.len());
        for (new_row, &old_row) in row_indices.iter().enumerate() {
            for (new_col, &old_col) in col_indices.iter().enumerate() {
                result[(new_row, new_col)] = self.values[(old_row, old_col)];
            }
        }

        Some(result)
    }

    /// Convenience method to get the full matrix
    ///
    /// ### Returns
    ///
    /// Faer matrix representing the full data.
    pub fn get_full_mat(&self) -> Mat<f64> {
        self.get_sub_mat(None, None).unwrap()
    }

    /// Convenience method to get submatrix with only row selection
    ///
    /// ### Params
    ///
    /// * `rows_to_select` - Names of the row features
    ///
    /// ### Returns
    ///
    /// Faer matrix with the selected rows.
    pub fn get_rows(&self, rows_to_select: &[&str]) -> Option<Mat<f64>> {
        self.get_sub_mat(Some(rows_to_select), None)
    }

    /// Convenience method to get submatrix with only column selection
    ///
    /// ### Params
    ///
    /// * `cols_to_select` - Names of the column features
    ///
    /// ### Returns
    ///
    /// Faer matrix with the selected columns.
    pub fn get_cols(&self, cols_to_select: &[&str]) -> Option<Mat<f64>> {
        self.get_sub_mat(None, Some(cols_to_select))
    }

    /// Get column names as references (for temporary use within same scope)
    ///
    /// ### Returns
    ///
    /// The column names as a reference
    pub fn get_col_names_refs(&self) -> Vec<&String> {
        self.col_names.keys().collect()
    }

    /// Get row names as references (for temporary use within same scope)
    ///
    /// ### Returns
    ///
    /// The row names as a reference
    pub fn get_row_names_refs(&self) -> Vec<&String> {
        self.row_names.keys().collect()
    }

    /// Get the row indices based on names
    ///
    /// ### Params
    ///
    /// * `rows_to_select` - Names of the row features
    pub fn get_row_indices(&self, rows_to_select: Option<&[&str]>) -> Vec<usize> {
        let row_indices: Vec<usize> = match rows_to_select {
            None => (0..self.values.nrows()).collect(), // All rows
            Some(rows) => {
                if rows.is_empty() {
                    return vec![];
                }
                let mut indices = Vec::new();
                for &row_name in rows {
                    if let Some(&index) = self.row_names.get(row_name) {
                        indices.push(index);
                    }
                }
                if indices.is_empty() {
                    return vec![];
                }
                indices
            }
        };

        row_indices
    }

    /// Get the col indices based on names
    ///
    /// ### Params
    ///
    /// * `rows_to_select` - Names of the row features
    pub fn get_col_indices(&self, cols_to_select: Option<&[&str]>) -> Vec<usize> {
        let col_indices: Vec<usize> = match cols_to_select {
            None => (0..self.values.ncols()).collect(), // All columns
            Some(cols) => {
                if cols.is_empty() {
                    return vec![];
                }
                let mut indices = Vec::new();
                for &col_name in cols {
                    if let Some(&index) = self.col_names.get(col_name) {
                        indices.push(index);
                    }
                }
                if indices.is_empty() {
                    return vec![];
                }
                indices
            }
        };

        col_indices
    }
}

/// Transform an R matrix to a Faer one
///
/// ### Params
///
/// * `x` - The R matrix to transform into a faer MatRef (with `f64`)
///
/// ### Returns
///
/// The faer `MatRef` from the original R matrix.
pub fn r_matrix_to_faer(x: &RMatrix<f64>) -> faer::MatRef<f64> {
    let ncol = x.ncols();
    let nrow = x.nrows();
    let data = x.data();

    MatRef::from_column_major_slice(data, nrow, ncol)
}

/// Transform an R matrix to a Faer one (i32)
///
/// ### Params
///
/// * `x` - The R matrix to transform into a faer MatRef (with `i32`)
///
/// ### Returns
///
/// The faer `MatRef` from the original R matrix.
pub fn r_matrix_to_faer_i32(x: &RMatrix<i32>) -> faer::MatRef<i32> {
    let ncol = x.ncols();
    let nrow = x.nrows();
    let data = x.data();

    MatRef::from_column_major_slice(data, nrow, ncol)
}

/// Transform an R matrix into a nested vector of booleans
///
/// ### Params
///
/// * `x` - The R matrix to transform into a vector of vectors with booleans
///
/// ### Returns
///
/// The nested vector with the outer vector representing the columns.
pub fn r_matrix_to_vec_bool(x: &RMatrix<Rbool>) -> Vec<Vec<bool>> {
    let ncol = x.ncols();
    let nrow = x.nrows();
    let data = x.data();

    (0..ncol)
        .map(|j| (0..nrow).map(|i| data[i + j * nrow].to_bool()).collect())
        .collect()
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
pub fn faer_to_r_matrix(x: faer::MatRef<f64>) -> extendr_api::RArray<f64, [usize; 2]> {
    let nrow = x.nrows();
    let ncol = x.ncols();

    RArray::new_matrix(nrow, ncol, |row, column| x[(row, column)])
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
pub fn sparse_matrix_to_list<T>(sparse: SparseColumnMatrix<T>) -> List
where
    T: Into<Robj>,
{
    let data: Vec<f64> = sparse
        .data
        .into_iter()
        .map(|x| {
            let robj: Robj = x.into();
            robj.as_real().unwrap()
        })
        .collect();
    let row_indices: Vec<usize> = sparse.row_indices;
    let col_ptr: Vec<usize> = sparse.col_ptrs;
    let nrow = sparse.nrow;
    let ncol = sparse.ncol;

    list![
        data = data,
        row_indices = row_indices,
        col_ptr = col_ptr,
        ncol = ncol,
        nrow = nrow
    ]
}
