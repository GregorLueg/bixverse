//! Some utility functions for the meta cell part. Usually to simplify or cast
//! stuff.

use bixverse_rs::prelude::*;

/// Function to cast types to f32 for meta cells
///
/// ### Params
///
/// * `sparse_mat` - `CompressedSparseData2<f64, f64>` with doubles from R
///
/// ### Returns
///
/// `CompressedSparseData2<f32, f32>`
pub fn cast_compressed_sparse_data(
    sparse_mat: CompressedSparseData2<f64, f64>,
) -> CompressedSparseData2<f32, f32> {
    let data_cast = sparse_mat.data.r_float_convert();
    let data_2_cast = sparse_mat.data_2.map(|x| x.r_float_convert());

    match sparse_mat.cs_type {
        CompressedSparseFormat::Csc => CompressedSparseData2::new_csc(
            &data_cast,
            &sparse_mat.indices,
            &sparse_mat.indptr,
            data_2_cast.as_deref(),
            sparse_mat.shape,
        ),
        CompressedSparseFormat::Csr => CompressedSparseData2::new_csr(
            &data_cast,
            &sparse_mat.indices,
            &sparse_mat.indptr,
            data_2_cast.as_deref(),
            sparse_mat.shape,
        ),
    }
}
