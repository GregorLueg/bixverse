use faer::traits::ComplexField;
use faer::{Mat, MatRef};
use rayon::prelude::*;
use rustc_hash::FxHashMap;
use std::cmp::PartialEq;
use std::marker::Sync;
use std::ops::{Add, AddAssign, DivAssign, Mul};

/////////////
// Helpers //
/////////////

/// Counts the zeroes in a given faer matrix
///
/// ### Params
///
/// * `mat` - The respective faer matrix
///
/// ### Returns
///
/// A tuple with the first being the total zeroes, the second the zeroes per
/// row and the last element being the column zeroes.
pub fn count_zeroes<T>(mat: &MatRef<T>) -> (usize, Vec<usize>, Vec<usize>)
where
    T: Default + PartialEq,
{
    let (nrow, ncol) = mat.shape();
    let mut total_zeroes = 0_usize;
    let mut row_zeroes = vec![0_usize; nrow];
    let mut col_zeroes = vec![0_usize; ncol];

    let zero = T::default();

    for j in 0..ncol {
        for i in 0..nrow {
            let val = unsafe { mat.get_unchecked(i, j) };
            if *val == zero {
                total_zeroes += 1;
                row_zeroes[i] += 1;
                col_zeroes[j] += 1;
            }
        }
    }

    (total_zeroes, row_zeroes, col_zeroes)
}

////////////////
// Structures //
////////////////

/// Structure for SparseColumnMatrices
///
/// ### Fields
///
/// * `data` - Vector with the data.
/// * `row_indices` - The row indices of the data.
/// * `col_ptrs` - The column pointers of the data.
/// * `ncol` - Original number of columns.
/// * `nrow` - Original number of rows.
#[derive(Debug, Clone)]
pub struct SparseColumnMatrix<T> {
    pub data: Vec<T>,
    pub row_indices: Vec<usize>,
    pub col_ptrs: Vec<usize>,
    pub ncol: usize,
    pub nrow: usize,
}

#[allow(dead_code)]
impl<T> SparseColumnMatrix<T>
where
    T: Clone + Default + PartialEq + ComplexField + From<f64>,
{
    /// Generate a new sparse column matrix from values pre-computed data
    ///
    /// ### Params
    ///
    /// * `data` - Slice of the data.
    /// * `row_indices` - Slice of the row indices of the data.
    /// * `col_ptrs` - Slice of the column pointers of the data.
    /// * `ncol` - Original number of columns.
    /// * `nrow` - Original number of rows.
    pub fn new(
        data: &[T],
        row_indices: &[usize],
        col_ptrs: &[usize],
        ncol: usize,
        nrow: usize,
    ) -> Self {
        Self {
            data: data.to_vec(),
            row_indices: row_indices.to_vec(),
            col_ptrs: col_ptrs.to_vec(),
            ncol,
            nrow,
        }
    }

    /// Convert a faer dense matrix to sparse column format
    ///
    /// ### Params
    ///
    /// * `dense` - The original dense matrix.
    pub fn from_dense_matrix(dense: MatRef<T>) -> Self {
        let ncol = dense.ncols();
        let nrow = dense.nrows();

        let mut values = Vec::new();
        let mut row_indices = Vec::new();
        let mut col_ptrs = Vec::with_capacity(ncol + 1);

        col_ptrs.push(0_usize);

        for col in 0..ncol {
            for row in 0..nrow {
                let value = dense.get(row, col).clone();
                if value != T::default() {
                    values.push(value);
                    row_indices.push(row);
                }
            }
            col_ptrs.push(values.len());
        }

        Self {
            data: values,
            row_indices,
            col_ptrs,
            ncol,
            nrow,
        }
    }

    /// To a dense faer matrix
    ///
    /// ### Returns
    ///
    /// Returns a dense faer matrix.
    pub fn to_dense_matrix(&self) -> Mat<T> {
        let mut dense = Mat::zeros(self.nrow, self.ncol);

        for col in 0..self.ncol {
            let start = self.col_ptrs[col];
            let end = self.col_ptrs[col + 1];

            for idx in start..end {
                let row = self.row_indices[idx];
                let value = &self.data[idx];
                *dense.get_mut(row, col) = value.clone();
            }
        }

        dense
    }

    /// Return the number of non-zero values
    ///
    /// ### Returns
    ///
    /// Return the total number of NNZ values in the data
    pub fn nnz(&self) -> usize {
        self.data.len()
    }

    /// Create a sparse matrix from a row major upper triangle stored value
    ///
    /// This is a helper function to transform potentially sparse symmetric matrices
    /// stored as upper-triangles into a sparse matrix format in Rust.
    ///
    /// ### Params
    ///
    /// * `upper_triangle` - Represents the values of the upper triangle in
    ///   row major formant
    /// * `n` - Original nrows and ncols.
    /// * `include_diagonal` - Are the diagonal values included.
    pub fn from_upper_triangle_sym(upper_triangle: &[T], n: usize, include_diagonal: bool) -> Self {
        let mut values = Vec::new();
        let mut row_indices = Vec::new();
        let mut col_ptrs = Vec::with_capacity(n + 1);

        col_ptrs.push(0);

        for col in 0..n {
            for row in 0..n {
                let value = if row == col && !include_diagonal {
                    T::from(1.0)
                } else if row < col {
                    // Upper triangle (row < col)
                    let offset = if include_diagonal {
                        row * n - row * (row + 1) / 2 + col
                    } else {
                        row * (n - 1) - row * (row + 1) / 2 + col - 1
                    };
                    upper_triangle[offset].clone()
                } else if row == col {
                    // Diagonal (when included)
                    let offset = row * n - row * (row + 1) / 2 + col;
                    upper_triangle[offset].clone()
                } else {
                    // Lower triangle (col < row) - symmetric
                    let offset = if include_diagonal {
                        col * n - col * (col + 1) / 2 + row
                    } else {
                        col * (n - 1) - col * (col + 1) / 2 + row - 1
                    };
                    upper_triangle[offset].clone()
                };

                if value != T::default() {
                    values.push(value);
                    row_indices.push(row);
                }
            }
            col_ptrs.push(values.len());
        }

        Self {
            data: values,
            row_indices,
            col_ptrs,
            ncol: n,
            nrow: n,
        }
    }
}

//////////////////////////////
// Sparse format conversion //
//////////////////////////////

/// Type to describe the CompressedSparseFormat
#[derive(Debug, Clone)]
pub enum CompressedSparseFormat {
    /// CSC-formatted data
    Csc,
    /// CSR-formatted data
    Csr,
}

#[allow(dead_code)]
impl CompressedSparseFormat {
    /// Returns boolean if it's CSC
    pub fn is_csc(&self) -> bool {
        matches!(self, CompressedSparseFormat::Csc)
    }
    /// Returns boolean if it's CSR
    pub fn is_csr(&self) -> bool {
        matches!(self, CompressedSparseFormat::Csr)
    }
}

/// Structure to store compressed sparse data of either type
///
/// ### Fields
///
/// * `data` - The values
/// * `indices` - The indices of the values
/// * `indptr` - The index pointers
/// * `cs_type` - Is the data stored in `Csr` or `Csc`.
/// * `data_2` - An optional second data layer
/// * `shape` - The shape of the underlying matrix
#[derive(Debug, Clone)]
pub struct CompressedSparseData<T, U = T>
where
    T: Clone + Default,
    U: Clone + Default,
{
    pub data: Vec<T>,
    pub indices: Vec<usize>,
    pub indptr: Vec<usize>,
    pub cs_type: CompressedSparseFormat,
    pub data_2: Option<Vec<U>>,
    pub shape: (usize, usize),
}

impl<T, U> CompressedSparseData<T, U>
where
    T: Clone + Default + Into<f64> + Sync + Add + PartialEq + Mul,
    U: Clone + Default + Sync + Add + PartialEq + Mul,
{
    /// Generate a nes CSC version of the matrix
    ///
    /// ### Params
    ///
    /// * `data` - The underlying data
    /// * `indices` - The index positions (in this case row indices)
    /// * `indptr` - The index pointer (in this case the column index pointers)
    /// * `data2` - An optional second layer
    #[allow(dead_code)]
    pub fn new_csc(
        data: &[T],
        indices: &[usize],
        indptr: &[usize],
        data2: Option<&[U]>,
        shape: (usize, usize),
    ) -> Self {
        Self {
            data: data.to_vec(),
            indices: indices.to_vec(),
            indptr: indptr.to_vec(), // Fixed: was using indices instead of indptr
            cs_type: CompressedSparseFormat::Csc,
            data_2: data2.map(|d| d.to_vec()),
            shape,
        }
    }

    /// Generate a nes CSR version of the matrix
    ///
    /// ### Params
    ///
    /// * `data` - The underlying data
    /// * `indices` - The index positions (in this case row indices)
    /// * `indptr` - The index pointer (in this case the column index pointers)
    /// * `data2` - An optional second layer
    pub fn new_csr(
        data: &[T],
        indices: &[usize],
        indptr: &[usize],
        data2: Option<&[U]>,
        shape: (usize, usize),
    ) -> Self {
        Self {
            data: data.to_vec(),
            indices: indices.to_vec(),
            indptr: indptr.to_vec(), // Fixed: was using indices instead of indptr
            cs_type: CompressedSparseFormat::Csr,
            data_2: data2.map(|d| d.to_vec()),
            shape,
        }
    }

    /// Transform from CSC to CSR or vice versa
    ///
    /// ### Returns
    ///
    /// The transformed/transposed version
    pub fn transform(&self) -> Self {
        match self.cs_type {
            CompressedSparseFormat::Csc => csc_to_csr(self),
            CompressedSparseFormat::Csr => csr_to_csc(self),
        }
    }

    /// Transpose and convert
    ///
    /// This is a helper to deal with the h5ad madness. Takes in for example
    /// a genes x cell CSR matrix from h5ad and transforms it into a cell x
    /// genes CSR matrix which bixverse expects. Same for CSC.
    ///
    /// ### Returns
    ///
    /// The transformed/transposed version
    pub fn transpose_and_convert(&self) -> Self {
        match self.cs_type {
            CompressedSparseFormat::Csr => {
                // convert first and then switch around
                let csc_version = csr_to_csc(self);
                CompressedSparseData {
                    data: csc_version.data,
                    indices: csc_version.indices,
                    indptr: csc_version.indptr,
                    cs_type: CompressedSparseFormat::Csr, // relabel as CSR
                    data_2: csc_version.data_2,
                    shape: (self.shape.1, self.shape.0), // swap dimensions
                }
            }
            CompressedSparseFormat::Csc => {
                // no conversion needed here! simple transpose is enough...
                CompressedSparseData {
                    data: self.data.clone(),
                    indices: self.indices.clone(),
                    indptr: self.indptr.clone(),
                    cs_type: CompressedSparseFormat::Csr,
                    data_2: self.data_2.clone(),
                    shape: (self.shape.1, self.shape.0),
                }
            }
        }
    }

    /// Transpose the matrix
    #[allow(dead_code)]
    pub fn transpose_from_h5ad(&self) -> Self {
        CompressedSparseData {
            data: self.data.clone(),
            indices: self.indices.clone(),
            indptr: self.indptr.clone(),
            cs_type: self.cs_type.clone(),
            data_2: self.data_2.clone(),
            shape: (self.shape.1, self.shape.0),
        }
    }

    /// Returns the shape of the matrix
    ///
    /// ### Returns
    ///
    /// A tuple of `(nrow, ncol)`
    pub fn shape(&self) -> (usize, usize) {
        self.shape
    }

    /// Returns the NNZ
    ///
    /// ### Returns
    ///
    /// The number of NNZ
    pub fn get_nnz(&self) -> usize {
        self.data.len()
    }

    /// Return the second layer
    ///
    /// If this does not exist, the function will panic
    ///
    /// ### Returns
    ///
    /// Vector of the second layer
    pub fn get_data2_unsafe(&self) -> Vec<U> {
        self.data_2.clone().unwrap()
    }
}

/// Transforms a CompressedSparseData that is CSC to CSR
///
/// ### Params
///
/// * `sparse_data` - The CompressedSparseData you want to transform
pub fn csc_to_csr<T, U>(sparse_data: &CompressedSparseData<T, U>) -> CompressedSparseData<T, U>
where
    T: Clone + Default + Into<f64> + Sync + Add + PartialEq + Mul,
    U: Clone + Default + Sync + Add + PartialEq + Mul,
{
    let (nrow, _) = sparse_data.shape();
    let nnz = sparse_data.get_nnz();
    let mut row_ptr = vec![0; nrow + 1];

    for &r in &sparse_data.indices {
        row_ptr[r + 1] += 1;
    }

    for i in 0..nrow {
        row_ptr[i + 1] += row_ptr[i];
    }

    let mut csr_data = vec![T::default(); nnz];
    let mut csr_data2 = sparse_data.data_2.as_ref().map(|_| vec![U::default(); nnz]);
    let mut csr_col_ind = vec![0; nnz];
    let mut next = row_ptr[..nrow].to_vec();

    for col in 0..(sparse_data.indptr.len() - 1) {
        for idx in sparse_data.indptr[col]..sparse_data.indptr[col + 1] {
            let row = sparse_data.indices[idx];
            let pos = next[row];

            csr_data[pos] = sparse_data.data[idx].clone();
            csr_col_ind[pos] = col;

            // Handle the second layer data
            if let (Some(ref source_data2), Some(ref mut csr_d2)) =
                (&sparse_data.data_2, &mut csr_data2)
            {
                csr_d2[pos] = source_data2[idx].clone();
            }

            next[row] += 1;
        }
    }

    CompressedSparseData {
        data: csr_data,
        indices: csr_col_ind,
        indptr: row_ptr,
        cs_type: CompressedSparseFormat::Csr,
        data_2: csr_data2,
        shape: sparse_data.shape(),
    }
}

/// Transform CSR stored data into CSC stored data
///
/// This version does a full memory copy of the data.
///
/// ### Params
///
/// * `sparse_data` - The CompressedSparseData you want to transform
///
/// ### Returns
///
/// The data in CSC format, i.e., `CscData`
pub fn csr_to_csc<T, U>(sparse_data: &CompressedSparseData<T, U>) -> CompressedSparseData<T, U>
where
    T: Clone + Default + Into<f64> + Sync + Add + PartialEq + Mul,
    U: Clone + Default + Sync + Add + PartialEq + Mul,
{
    let nnz = sparse_data.get_nnz();
    let (_, ncol) = sparse_data.shape();
    let mut col_ptr = vec![0; ncol + 1];

    // Count occurrences per column
    for &c in &sparse_data.indices {
        col_ptr[c + 1] += 1;
    }

    // Cumulative sum to get column pointers
    for i in 0..ncol {
        col_ptr[i + 1] += col_ptr[i];
    }

    let mut csc_data = vec![T::default(); nnz];
    let mut csc_data2 = sparse_data.data_2.as_ref().map(|_| vec![U::default(); nnz]);
    let mut csc_row_ind = vec![0; nnz];
    let mut next = col_ptr[..ncol].to_vec();

    // Iterate through rows and place data in CSC format
    for row in 0..(sparse_data.indptr.len() - 1) {
        for idx in sparse_data.indptr[row]..sparse_data.indptr[row + 1] {
            let col = sparse_data.indices[idx];
            let pos = next[col];

            csc_data[pos] = sparse_data.data[idx].clone();
            csc_row_ind[pos] = row;

            // Handle the second layer data
            if let (Some(ref source_data2), Some(ref mut csc_d2)) =
                (&sparse_data.data_2, &mut csc_data2)
            {
                csc_d2[pos] = source_data2[idx].clone();
            }

            next[col] += 1;
        }
    }

    CompressedSparseData {
        data: csc_data,
        indices: csc_row_ind,
        indptr: col_ptr,
        cs_type: CompressedSparseFormat::Csc,
        data_2: csc_data2,
        shape: sparse_data.shape(),
    }
}

/// Transform COO stored data into CSR
///
/// ### Params
///
/// * `rows` - Row indices
/// * `cols` - Col indices
/// * `vals` - The values to store in the matrix
///
/// ### Returns
///
/// `CompressedSparseData` in CSR format
pub fn coo_to_csr<T>(
    rows: &[usize],
    cols: &[usize],
    vals: &[T],
    shape: (usize, usize),
) -> CompressedSparseData<T>
where
    T: Clone + Default + Into<f64> + Sync + Add<Output = T> + AddAssign + PartialEq + Copy + Mul,
    <T as std::ops::Add>::Output: std::cmp::PartialEq<T>,
{
    let n_rows = shape.0;

    // Sort by (row, col) and merge duplicates
    let mut entries: Vec<(usize, usize, T)> = rows
        .iter()
        .zip(cols.iter())
        .zip(vals.iter())
        .map(|((&r, &c), &v)| (r, c, v))
        .collect();

    entries.sort_unstable_by_key(|&(r, c, _)| (r, c));

    // merge duplicates; can happen during additions
    let mut merged_entries = Vec::new();
    if !entries.is_empty() {
        let mut current = entries[0];

        for &(r, c, v) in &entries[1..] {
            if r == current.0 && c == current.1 {
                current.2 += v;
            } else {
                if current.2 != T::default() {
                    merged_entries.push(current);
                }
                current = (r, c, v);
            }
        }
        if current.2 != T::default() {
            merged_entries.push(current);
        }
    }

    // build CSR from merged entries
    let final_nnz = merged_entries.len();
    let mut data = Vec::with_capacity(final_nnz);
    let mut indices = Vec::with_capacity(final_nnz);
    let mut indptr = vec![0usize; n_rows + 1];

    for &(row, col, val) in &merged_entries {
        data.push(val);
        indices.push(col);
        indptr[row + 1] += 1;
    }

    // Convert counts to cumulative offsets
    for i in 0..n_rows {
        indptr[i + 1] += indptr[i];
    }

    CompressedSparseData::new_csr(&data, &indices, &indptr, None, shape)
}

/// Optimised COO to CSR - assumes input is already sorted by (row, col)
///
/// ### Params
///
/// * `rows` - Row indices (must be sorted by row first, then col)
/// * `cols` - Col indices  
/// * `vals` - Values
/// * `shape` - Matrix dimensions
/// * `is_sorted` - If true, skips sorting step
///
/// ### Returns
///
/// CSR matrix
pub fn coo_to_csr_presorted<T>(
    rows: &[usize],
    cols: &[usize],
    vals: &[T],
    shape: (usize, usize),
) -> CompressedSparseData<T>
where
    T: Clone + Default + Into<f64> + Sync + Add<Output = T> + AddAssign + PartialEq + Copy + Mul,
{
    let n_rows = shape.0;
    let nnz = rows.len();

    let mut data = Vec::with_capacity(nnz);
    let mut indices = Vec::with_capacity(nnz);
    let mut indptr = vec![0usize; n_rows + 1];

    // unsafe to squeeze out performance...
    unsafe {
        data.set_len(nnz);
        indices.set_len(nnz);

        let data_ptr: *mut T = data.as_mut_ptr();
        let indices_ptr: *mut usize = indices.as_mut_ptr();
        let indptr_ptr: *mut usize = indptr.as_mut_ptr();

        for i in 0..nnz {
            *data_ptr.add(i) = *vals.get_unchecked(i);
            *indices_ptr.add(i) = *cols.get_unchecked(i);
            let row = *rows.get_unchecked(i);
            *indptr_ptr.add(row + 1) += 1;
        }

        for i in 0..n_rows {
            *indptr_ptr.add(i + 1) += *indptr_ptr.add(i);
        }
    }

    CompressedSparseData::new_csr(&data, &indices, &indptr, None, shape)
}

/// Add two CSR matrices together
///
/// ### Params
///
/// * `a` - Reference to the first CompressedSparseData (in CSR format!)
/// * `b` - Reference to the second CompressedSparseData (in CSR format!)
///
/// ### Returns
///
/// `CompressedSparseData` with added values between the two.
pub fn sparse_add_csr<T>(
    a: &CompressedSparseData<T>,
    b: &CompressedSparseData<T>,
) -> CompressedSparseData<T>
where
    T: Clone + Default + Into<f64> + Sync + Add<Output = T> + AddAssign + PartialEq + Copy + Mul,
    <T as std::ops::Add>::Output: std::cmp::PartialEq<T>,
{
    assert_eq!(a.shape, b.shape);
    assert!(a.cs_type.is_csr() && b.cs_type.is_csr());

    const EPSILON: f32 = 1e-9;
    let n_rows = a.shape.0;

    let mut rows = Vec::new();
    let mut cols = Vec::new();
    let mut vals = Vec::new();

    for i in 0..n_rows {
        let a_start = a.indptr[i];
        let a_end = a.indptr[i + 1];
        let b_start = b.indptr[i];
        let b_end = b.indptr[i + 1];

        let mut a_idx = a_start;
        let mut b_idx = b_start;

        while a_idx < a_end || b_idx < b_end {
            if a_idx < a_end && (b_idx >= b_end || a.indices[a_idx] < b.indices[b_idx]) {
                rows.push(i);
                cols.push(a.indices[a_idx]);
                vals.push(a.data[a_idx]);
                a_idx += 1;
            } else if b_idx < b_end && (a_idx >= a_end || b.indices[b_idx] < a.indices[a_idx]) {
                rows.push(i);
                cols.push(b.indices[b_idx]);
                vals.push(b.data[b_idx]);
                b_idx += 1;
            } else {
                let val = a.data[a_idx] + b.data[b_idx];
                if val.into().abs() > EPSILON as f64 {
                    rows.push(i);
                    cols.push(a.indices[a_idx]);
                    vals.push(val);
                }
                a_idx += 1;
                b_idx += 1;
            }
        }
    }

    // output is already sorted by (row, col), build CSR directly
    coo_to_csr_presorted(&rows, &cols, &vals, a.shape)
}

/// Scalar multiplication of CSR matrix
///
/// ### Params
///
/// * `a` - Reference to the first CompressedSparseData (in CSR format!)
/// * `scalar` - The scalar value to multiply with
///
/// ### Returns
///
/// `CompressedSparseData` with the data multiplied by the scalar.
pub fn sparse_scalar_multiply_csr<T>(
    a: &CompressedSparseData<T>,
    scalar: T,
) -> CompressedSparseData<T>
where
    T: Clone
        + Default
        + Into<f64>
        + Sync
        + Add<Output = T>
        + AddAssign
        + PartialEq
        + Copy
        + Mul<Output = T>
        + std::marker::Send,
{
    let data: Vec<T> = a.data.par_iter().map(|&v| v * scalar).collect();
    CompressedSparseData::new_csr(&data, &a.indices, &a.indptr, None, a.shape)
}

/// Sparse matrix subtraction
///
/// ### Params
///
/// * `a` - Reference to the first CompressedSparseData (in CSR format!)
/// * `b` - Reference to the second CompressedSparseData (in CSR format!)
///
/// ### Returns
///
/// The subtracted new matrix
pub fn sparse_subtract_csr<T>(
    a: &CompressedSparseData<T>,
    b: &CompressedSparseData<T>,
) -> CompressedSparseData<T>
where
    T: Clone
        + Default
        + Into<f64>
        + Sync
        + Add<Output = T>
        + AddAssign
        + PartialEq
        + Copy
        + Mul<Output = T>
        + std::ops::Sub<Output = T>,
{
    assert_eq!(a.shape, b.shape);
    assert!(a.cs_type.is_csr() && b.cs_type.is_csr());

    const EPSILON: f32 = 1e-9;
    let n_rows = a.shape.0;

    let mut rows = Vec::new();
    let mut cols = Vec::new();
    let mut vals = Vec::new();

    for i in 0..n_rows {
        let a_start = a.indptr[i];
        let a_end = a.indptr[i + 1];
        let b_start = b.indptr[i];
        let b_end = b.indptr[i + 1];

        let mut a_idx = a_start;
        let mut b_idx = b_start;

        while a_idx < a_end || b_idx < b_end {
            if a_idx < a_end && (b_idx >= b_end || a.indices[a_idx] < b.indices[b_idx]) {
                rows.push(i);
                cols.push(a.indices[a_idx]);
                vals.push(a.data[a_idx]);
                a_idx += 1;
            } else if b_idx < b_end && (a_idx >= a_end || b.indices[b_idx] < a.indices[a_idx]) {
                rows.push(i);
                cols.push(b.indices[b_idx]);
                vals.push(T::default() - b.data[b_idx]);
                b_idx += 1;
            } else {
                let val = a.data[a_idx] - b.data[b_idx];
                if val.into().abs() > EPSILON as f64 {
                    rows.push(i);
                    cols.push(a.indices[a_idx]);
                    vals.push(val);
                }
                a_idx += 1;
                b_idx += 1;
            }
        }
    }

    // already sorted
    coo_to_csr_presorted(&rows, &cols, &vals, a.shape)
}

/// Element-wise sparse multiplication
///
/// ### Params
///
/// * `a` - Reference to the first CompressedSparseData (in CSR format!)
/// * `b` - Reference to the second CompressedSparseData (in CSR format!)
///
/// ### Returns
///
/// The multiplied matrix.
pub fn sparse_multiply_elementwise_csr<T>(
    a: &CompressedSparseData<T>,
    b: &CompressedSparseData<T>,
) -> CompressedSparseData<T>
where
    T: Clone
        + Default
        + Into<f64>
        + Sync
        + Add<Output = T>
        + AddAssign
        + PartialEq
        + Copy
        + Mul<Output = T>,
    <T as std::ops::Add>::Output: std::cmp::PartialEq<T>,
{
    assert_eq!(a.shape, b.shape);
    assert!(a.cs_type.is_csr() && b.cs_type.is_csr());
    let n_rows = a.shape.0;
    let mut rows = Vec::new();
    let mut cols = Vec::new();
    let mut vals = Vec::new();
    for i in 0..n_rows {
        let a_start = a.indptr[i];
        let a_end = a.indptr[i + 1];
        let b_start = b.indptr[i];
        let b_end = b.indptr[i + 1];
        let mut a_idx = a_start;
        let mut b_idx = b_start;
        while a_idx < a_end && b_idx < b_end {
            match a.indices[a_idx].cmp(&b.indices[b_idx]) {
                std::cmp::Ordering::Less => {
                    a_idx += 1;
                }
                std::cmp::Ordering::Greater => {
                    b_idx += 1;
                }
                std::cmp::Ordering::Equal => {
                    // Same column - multiply
                    let val = a.data[a_idx] * b.data[b_idx];
                    if val != T::default() {
                        rows.push(i);
                        cols.push(a.indices[a_idx]);
                        vals.push(val);
                    }
                    a_idx += 1;
                    b_idx += 1;
                }
            }
        }
    }
    coo_to_csr(&rows, &cols, &vals, a.shape)
}

/// CSR matrix multiplication
///
/// This function implements the `CSR @ CSR` type math
///
/// ### Params
///
/// * `a` - First matrix in CSR format
/// * `b` - Second matrix in CSR format
///
/// ### Returns
///
/// Result of A @ B in CSR format
#[allow(dead_code)]
pub fn csr_matmul_csr<T>(
    a: &CompressedSparseData<T>,
    b: &CompressedSparseData<T>,
) -> CompressedSparseData<T>
where
    T: Clone
        + Default
        + Into<f64>
        + Sync
        + Add<Output = T>
        + PartialEq
        + Copy
        + Mul<Output = T>
        + AddAssign
        + std::marker::Send,
    <T as std::ops::Add>::Output: std::cmp::PartialEq<T>,
{
    assert!(a.cs_type.is_csr() && b.cs_type.is_csr());
    assert_eq!(a.shape.1, b.shape.0, "Dimension mismatch");

    let nrows = a.shape.0;
    let ncols = b.shape.1;

    // Parallel computation of each row
    let row_results: Vec<Vec<(usize, T)>> = (0..nrows)
        .into_par_iter()
        .map(|i| {
            let mut row_vals = FxHashMap::default();

            let a_row_start = a.indptr[i];
            let a_row_end = a.indptr[i + 1];

            for a_idx in a_row_start..a_row_end {
                let k = a.indices[a_idx];
                let a_val = a.data[a_idx];

                let b_row_start = b.indptr[k];
                let b_row_end = b.indptr[k + 1];

                for b_idx in b_row_start..b_row_end {
                    let j = b.indices[b_idx];
                    let b_val = b.data[b_idx];
                    *row_vals.entry(j).or_insert(T::default()) += a_val * b_val;
                }
            }

            // Sort columns within this row
            let mut sorted_row: Vec<(usize, T)> = row_vals.into_iter().collect();
            sorted_row.sort_unstable_by_key(|(col, _)| *col);
            sorted_row
        })
        .collect();

    // Build CSR directly from row results
    let mut data = Vec::new();
    let mut indices = Vec::new();
    let mut indptr = vec![0; nrows + 1];

    for (i, row) in row_results.iter().enumerate() {
        for &(col, val) in row {
            data.push(val);
            indices.push(col);
        }
        indptr[i + 1] = data.len();
    }

    CompressedSparseData::new_csr(&data, &indices, &indptr, None, (nrows, ncols))
}

/// Normalises the columns of a CSR matrix to a sum of 1 (L1 norm)
///
/// ### Params
///
/// * `csr` - Mutable reference to the CSR matrix (modified in-place)
pub fn normalise_csr_columns_l1<T>(csr: &mut CompressedSparseData<T>)
where
    T: Clone
        + Default
        + Into<f64>
        + Sync
        + Add<Output = T>
        + PartialEq
        + Copy
        + Mul<Output = T>
        + AddAssign
        + DivAssign,
    <T as std::ops::Add>::Output: std::cmp::PartialEq<T>,
{
    assert!(csr.cs_type.is_csr(), "Matrix must be in CSR format");

    let ncols = csr.shape.1;

    let mut col_sums = vec![T::default(); ncols];

    for (idx, &col) in csr.indices.iter().enumerate() {
        col_sums[col] += csr.data[idx]
    }

    for (idx, &col) in csr.indices.iter().enumerate() {
        let sum = col_sums[col];
        if sum.into() > 1e-15 {
            csr.data[idx] /= sum;
        }
    }
}

/// Compute Frobenius norm of sparse matrix
///
/// ### Params
///
/// * `mat` - Sparse matrix in CSR or CSC format
///
/// ### Returns
///
/// Frobenius norm ||A||_F = sqrt(sum(A_ij^2))
pub fn frobenius_norm<T>(mat: &CompressedSparseData<T>) -> f32
where
    T: Clone
        + Default
        + Into<f32>
        + Sync
        + Add<Output = T>
        + PartialEq
        + Copy
        + Mul<Output = T>
        + AddAssign
        + std::iter::Sum,
{
    mat.data
        .iter()
        .map(|&v| {
            let val: f32 = v.into();
            val * val
        })
        .sum::<f32>()
        .sqrt()
}

/// Remove zeros from sparse matrix
///
/// ### Params
///
/// * `mat` - Matrix from which to remove the zeroes
///
/// ### Returns
///
/// The Matrix with 0's removed.
pub fn eliminate_zeros_csr<T>(mat: CompressedSparseData<T>) -> CompressedSparseData<T>
where
    T: Clone
        + Default
        + Into<f64>
        + Sync
        + Add<Output = T>
        + AddAssign
        + PartialEq
        + Copy
        + Mul<Output = T>,
    <T as std::ops::Add>::Output: std::cmp::PartialEq<T>,
{
    let mut rows = Vec::new();
    let mut cols = Vec::new();
    let mut vals = Vec::new();

    let n_rows = mat.shape.0;
    for i in 0..n_rows {
        let start = mat.indptr[i];
        let end = mat.indptr[i + 1];

        for j in start..end {
            if mat.data[j] != T::default() {
                rows.push(i);
                cols.push(mat.indices[j]);
                vals.push(mat.data[j]);
            }
        }
    }

    coo_to_csr(&rows, &cols, &vals, mat.shape)
}

/// Sparse matrix-vector multiplication
///
/// Multiply a sparse CSR matrix with a vector
///
/// ### Params
///
/// * `mat` - The Compressed Sparse matrix in CSR format
/// * `vec` - The Vector to multiply with
///
/// ### Params
///
/// The resulting vector
pub fn csr_matvec<T>(mat: &CompressedSparseData<T>, vec: &[T]) -> Vec<T>
where
    T: Clone
        + Default
        + Into<f64>
        + Sync
        + Add<Output = T>
        + AddAssign
        + PartialEq
        + Copy
        + Mul<Output = T>,
    <T as std::ops::Add>::Output: std::cmp::PartialEq<T>,
{
    let mut result = vec![T::default(); mat.shape.0];
    for i in 0..mat.shape.0 {
        let row_start = mat.indptr[i];
        let row_end = mat.indptr[i + 1];
        let mut sum = T::default();
        for idx in row_start..row_end {
            sum += mat.data[idx] * vec[mat.indices[idx]];
        }
        result[i] = sum;
    }
    result
}

/////////////////////////
// VERY optimised code //
/////////////////////////

/// Sparse accumulator for efficient sparse matrix multiplication
///
/// ### Fields
///
/// * `values` - Vector storing accumulated values for each index
/// * `indices` - Vector of active (non-zero) indices
/// * `flags` - Boolean flags indicating which indices are active
struct SparseAccumulator<T>
where
    T: Copy + Default + AddAssign,
{
    values: Vec<T>,
    indices: Vec<usize>,
    flags: Vec<bool>,
}

impl<T> SparseAccumulator<T>
where
    T: Copy + Default + AddAssign,
{
    /// Create a new sparse accumulator
    ///
    /// ### Params
    ///
    /// * `size` - Maximum number of indices to accumulate
    fn new(size: usize) -> Self {
        Self {
            values: vec![T::default(); size],
            indices: Vec::with_capacity(size / 10),
            flags: vec![false; size],
        }
    }

    /// Add a value to the accumulator at the given index
    ///
    /// ### Params
    ///
    /// * `idx` - Index to accumulate at
    /// * `val` - Value to add
    ///
    /// ### Safety
    ///
    /// `idx` must be less than the size specified during construction
    #[inline]
    unsafe fn add(&mut self, idx: usize, val: T) {
        if !*self.flags.get_unchecked(idx) {
            *self.flags.get_unchecked_mut(idx) = true;
            self.indices.push(idx);
            *self.values.get_unchecked_mut(idx) = val;
        } else {
            *self.values.get_unchecked_mut(idx) += val;
        }
    }

    /// Extract accumulated values as sorted index-value pairs and reset the accumulator
    ///
    /// ### Returns
    ///
    /// Vector of (index, value) pairs sorted by index
    #[inline]
    fn extract_sorted(&mut self) -> Vec<(usize, T)> {
        self.indices.sort_unstable();
        let result: Vec<(usize, T)> = unsafe {
            self.indices
                .iter()
                .map(|&i| (i, *self.values.get_unchecked(i)))
                .collect()
        };
        // Reset for next use
        unsafe {
            for &idx in &self.indices {
                *self.flags.get_unchecked_mut(idx) = false;
                *self.values.get_unchecked_mut(idx) = T::default();
            }
        }
        self.indices.clear();
        result
    }
}

/// Multiply two CSR matrices using sparse accumulators and parallel processing
///
/// ### Params
///
/// * `a` - Left CSR matrix
/// * `b` - Right CSR matrix
///
/// ### Returns
///
/// Product matrix in CSR format
pub fn csr_matmul_csr_optimised<T>(
    a: &CompressedSparseData<T>,
    b: &CompressedSparseData<T>,
) -> CompressedSparseData<T>
where
    T: Clone
        + Default
        + Copy
        + Sync
        + Send
        + Add<Output = T>
        + Mul<Output = T>
        + AddAssign
        + PartialEq
        + Into<f64>,
{
    assert!(a.cs_type.is_csr() && b.cs_type.is_csr());
    assert_eq!(a.shape.1, b.shape.0, "Dimension mismatch");

    let nrows = a.shape.0;
    let ncols = b.shape.1;

    let row_results: Vec<Vec<(usize, T)>> = (0..nrows)
        .into_par_iter()
        .map(|i| {
            let mut acc = SparseAccumulator::new(ncols);

            unsafe {
                let a_indptr = a.indptr.as_ptr();
                let a_indices = a.indices.as_ptr();
                let a_data = a.data.as_ptr();
                let b_indptr = b.indptr.as_ptr();
                let b_indices = b.indices.as_ptr();
                let b_data = b.data.as_ptr();

                let a_start = *a_indptr.add(i);
                let a_end = *a_indptr.add(i + 1);

                for a_idx in a_start..a_end {
                    let k = *a_indices.add(a_idx);
                    let a_val = *a_data.add(a_idx);

                    let b_start = *b_indptr.add(k);
                    let b_end = *b_indptr.add(k + 1);

                    for b_idx in b_start..b_end {
                        let j = *b_indices.add(b_idx);
                        let b_val = *b_data.add(b_idx);
                        acc.add(j, a_val * b_val);
                    }
                }
            }

            acc.extract_sorted()
        })
        .collect();

    // direct CSR construction
    let total_nnz: usize = row_results.iter().map(|r| r.len()).sum();
    let mut data = Vec::with_capacity(total_nnz);
    let mut indices = Vec::with_capacity(total_nnz);
    let mut indptr = Vec::with_capacity(nrows + 1);
    indptr.push(0);

    for row in row_results {
        for (col, val) in row {
            data.push(val);
            indices.push(col);
        }
        indptr.push(data.len());
    }

    CompressedSparseData::new_csr(&data, &indices, &indptr, None, (nrows, ncols))
}

///////////
// Tests //
///////////

#[cfg(test)]
mod tests {
    use super::*;
    use faer::mat;

    #[test]
    fn test_dense_to_sparse_conversion() {
        let dense = mat![[1.0, 0.0, 3.0], [0.0, 2.0, 0.0], [4.0, 0.0, 5.0]];

        let sparse_obj = SparseColumnMatrix::from_dense_matrix(dense.as_ref());

        assert_eq!(sparse_obj.nrow, 3);
        assert_eq!(sparse_obj.ncol, 3);
        assert_eq!(sparse_obj.nnz(), 5);

        let expected_values = vec![1.0, 4.0, 2.0, 3.0, 5.0];

        assert_eq!(sparse_obj.data, expected_values);
    }

    #[test]
    fn test_dense_to_sparse_to_dense_conversion() {
        let dense = mat![[1.0, 0.0, 3.0], [0.0, 2.0, 0.0], [4.0, 0.0, 5.0]];

        let sparse_obj = SparseColumnMatrix::from_dense_matrix(dense.as_ref());

        let redense = sparse_obj.to_dense_matrix();

        assert_eq!(dense, redense);
    }

    #[test]
    fn test_raw_to_dense_conversion() {
        let data = vec![1.0, 4.0, 2.0, 3.0, 5.0];
        let row_indices: Vec<usize> = vec![0, 2, 1, 0, 2];
        let col_ptr: Vec<usize> = vec![0, 2, 3, 5];

        let sparse_obj = SparseColumnMatrix::new(&data, &row_indices, &col_ptr, 3, 3);

        let dense = mat![[1.0, 0.0, 3.0], [0.0, 2.0, 0.0], [4.0, 0.0, 5.0]];

        let redense = sparse_obj.to_dense_matrix();

        assert_eq!(dense, redense);
    }

    #[test]
    fn test_from_upper_triangle_vec_with_diagonal() {
        // Upper triangle with diagonal: [1, 0.8, 0.6, 1, 0.3, 1]
        // Represents correlation matrix:
        // [1.0, 0.8, 0.6]
        // [0.8, 1.0, 0.3]
        // [0.6, 0.3, 1.0]
        let upper_tri = vec![1.0, 0.8, 0.6, 1.0, 0.3, 1.0];
        let sparse = SparseColumnMatrix::from_upper_triangle_sym(&upper_tri, 3, true);

        let dense = sparse.to_dense_matrix();
        let expected = mat![[1.0, 0.8, 0.6], [0.8, 1.0, 0.3], [0.6, 0.3, 1.0]];

        assert_eq!(dense, expected);
    }

    #[test]
    fn test_from_upper_triangle_vec_without_diagonal() {
        // Upper triangle without diagonal: [0.8, 0.6, 0.3]
        // With implied diagonal of 1's, represents:
        // [1.0, 0.8, 0.6]
        // [0.8, 1.0, 0.3]
        // [0.6, 0.3, 1.0]
        let upper_tri = vec![0.8, 0.6, 0.3];
        let sparse = SparseColumnMatrix::from_upper_triangle_sym(&upper_tri, 3, false);

        let dense = sparse.to_dense_matrix();
        let expected = mat![[1.0, 0.8, 0.6], [0.8, 1.0, 0.3], [0.6, 0.3, 1.0]];

        assert_eq!(dense, expected);
    }

    #[test]
    fn test_csc_to_csr_conversion_single_layer() {
        // Test matrix:
        // [1 0 2]
        // [0 3 0]
        // [4 0 5]
        // CSC format
        let data = vec![1, 4, 3, 2, 5];
        let row_ind = vec![0, 2, 1, 0, 2];
        let col_ptr = vec![0, 2, 3, 5];
        let shape = (3, 3);

        let csc_matrix =
            CompressedSparseData::new_csc(&data, &row_ind, &col_ptr, None::<&[i32]>, shape);

        let csr_matrix = csc_to_csr(&csc_matrix);

        // Expected CSR format
        assert_eq!(csr_matrix.data, vec![1, 2, 3, 4, 5]);
        assert_eq!(csr_matrix.indices, vec![0, 2, 1, 0, 2]);
        assert_eq!(csr_matrix.indptr, vec![0, 2, 3, 5]);
        assert_eq!(csr_matrix.data_2, None);
        assert_eq!(csr_matrix.shape, (3, 3));
        assert!(matches!(csr_matrix.cs_type, CompressedSparseFormat::Csr));
    }

    #[test]
    fn test_csc_to_csr_conversion_dual_layer() {
        // Test matrix:
        // [1 0 2] [1.0 0   2.0]
        // [0 3 0] [0   3.0 0  ]
        // [4 0 5] [4.0 0   5.0]
        // CSC format for both layers
        let data = vec![1, 4, 3, 2, 5];
        let data2 = vec![1.0, 4.0, 3.0, 2.0, 5.0];
        let row_ind = vec![0, 2, 1, 0, 2];
        let col_ptr = vec![0, 2, 3, 5];
        let shape = (3, 3);

        let csc_matrix =
            CompressedSparseData::new_csc(&data, &row_ind, &col_ptr, Some(&data2), shape);

        let csr_matrix = csc_to_csr(&csc_matrix);

        // Expected CSR format
        assert_eq!(csr_matrix.data, vec![1, 2, 3, 4, 5]);
        assert_eq!(csr_matrix.indices, vec![0, 2, 1, 0, 2]);
        assert_eq!(csr_matrix.indptr, vec![0, 2, 3, 5]);
        assert_eq!(csr_matrix.data_2, Some(vec![1.0, 2.0, 3.0, 4.0, 5.0]));
        assert_eq!(csr_matrix.shape, (3, 3));
        assert!(matches!(csr_matrix.cs_type, CompressedSparseFormat::Csr));
    }

    #[test]
    fn test_csr_to_csc_conversion_single_layer() {
        // Test matrix:
        // [1 0 2]
        // [0 3 0]
        // [4 0 5]
        // CSR format
        let data = vec![1, 2, 3, 4, 5];
        let col_ind = vec![0, 2, 1, 0, 2];
        let row_ptr = vec![0, 2, 3, 5];
        let shape = (3, 3);

        let csr_matrix =
            CompressedSparseData::new_csr(&data, &col_ind, &row_ptr, None::<&[i32]>, shape);

        let csc_matrix = csr_to_csc(&csr_matrix);

        // Expected CSC format
        assert_eq!(csc_matrix.data, vec![1, 4, 3, 2, 5]);
        assert_eq!(csc_matrix.indices, vec![0, 2, 1, 0, 2]);
        assert_eq!(csc_matrix.indptr, vec![0, 2, 3, 5]);
        assert_eq!(csc_matrix.data_2, None);
        assert_eq!(csc_matrix.shape, (3, 3));
        assert!(matches!(csc_matrix.cs_type, CompressedSparseFormat::Csc));
    }

    #[test]
    fn test_csr_to_csc_conversion_dual_layer() {
        // Test matrix:
        // [1 0 2] [1.0 0   2.0]
        // [0 3 0] [0   3.0 0  ]
        // [4 0 5] [4.0 0   5.0]
        // CSR format for both layers
        let data = vec![1, 2, 3, 4, 5];
        let data2 = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let col_ind = vec![0, 2, 1, 0, 2];
        let row_ptr = vec![0, 2, 3, 5];
        let shape = (3, 3);

        let csr_matrix =
            CompressedSparseData::new_csr(&data, &col_ind, &row_ptr, Some(&data2), shape);

        let csc_matrix = csr_to_csc(&csr_matrix);

        // Expected CSC format
        assert_eq!(csc_matrix.data, vec![1, 4, 3, 2, 5]);
        assert_eq!(csc_matrix.indices, vec![0, 2, 1, 0, 2]);
        assert_eq!(csc_matrix.indptr, vec![0, 2, 3, 5]);
        assert_eq!(csc_matrix.data_2, Some(vec![1.0, 4.0, 3.0, 2.0, 5.0]));
        assert_eq!(csc_matrix.shape, (3, 3));
        assert!(matches!(csc_matrix.cs_type, CompressedSparseFormat::Csc));
    }
}
