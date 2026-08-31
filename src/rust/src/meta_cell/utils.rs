//! Some utility functions for the meta cell part. Usually to simplify or cast
//! stuff.
//!
//! ### On the second layer
//!
//! The single cell readers want two layers, a raw count and a normalised one,
//! and R only ever hands one matrix over. So the second layer is a cast of the
//! first and the caller picks which of the two the kernel reads by choosing
//! the target type. Deriving it here rather than letting
//! `list_to_sparse_matrix` clone it saves a full nnz-sized `f64` buffer, and
//! taking the `List` rather than an already-parsed matrix is what makes that
//! decision unrepresentable at the call sites.

use bixverse_rs::prelude::*;
use extendr_api::*;

/// Read the meta cell counts into raw `u32` counts and an `f32` second layer
///
/// The layer to read is the raw one, so this is the cast for anything counting
/// rather than scoring. It truncates towards zero, so do not hand it
/// normalised counts.
///
/// ### Params
///
/// * `sparse_data` - The named list from R, needing `data`, `indptr`,
///   `indices`, `nrow`, `ncol` and `cs_type`.
///
/// ### Returns
///
/// The matrix as `CompressedSparseData2<u32, f32>`.
pub fn mc_list_to_sparse_u32(sparse_data: List) -> Result<CompressedSparseData2<u32, f32>> {
    let sparse: CompressedSparseData2<f64, f64> =
        list_to_sparse_matrix(sparse_data, false).to_extendr()?;

    // Two `collect`s over the same buffer rather than one push loop: each is
    // vectorised, the interleaved version is not.
    let raw: Vec<u32> = sparse.data.iter().map(|&x| x as u32).collect();
    let norm: Vec<f32> = sparse.data.iter().map(|&x| x as f32).collect();

    Ok(CompressedSparseData2::from_parts(
        raw,
        sparse.indices,
        sparse.indptr,
        Some(norm),
        sparse.cs_type,
        sparse.shape,
    ))
}

/// Read the meta cell counts into `f32` for both layers
///
/// The layer to read is the second one, so this is the cast for anything
/// scoring normalised counts. Keeping both layers in `f32` is the point:
/// [mc_list_to_sparse_u32] would truncate them.
///
/// ### Params
///
/// * `sparse_data` - The named list from R, needing `data`, `indptr`,
///   `indices`, `nrow`, `ncol` and `cs_type`.
///
/// ### Returns
///
/// The matrix as `CompressedSparseData2<f32, f32>`.
pub fn mc_list_to_sparse_f32(sparse_data: List) -> Result<CompressedSparseData2<f32, f32>> {
    let sparse: CompressedSparseData2<f64, f64> =
        list_to_sparse_matrix(sparse_data, false).to_extendr()?;

    let raw: Vec<f32> = sparse.data.iter().map(|&x| x as f32).collect();
    let norm = raw.clone();

    Ok(CompressedSparseData2::from_parts(
        raw,
        sparse.indices,
        sparse.indptr,
        Some(norm),
        sparse.cs_type,
        sparse.shape,
    ))
}
