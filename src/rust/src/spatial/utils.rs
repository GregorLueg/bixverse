//! Shared conversions for the spatial Rust <> R interface.
//!
//! Three things cross this boundary in more than one place: spot coordinates,
//! adjacency lists and the local/global index distinction. All of it lives
//! here so the four binding modules do not each grow their own version.
//!
//! ### Index spaces
//!
//! | thing | base | space |
//! |---|---|---|
//! | `spots_to_keep` | 0 | global, positions in the count store |
//! | `neighbours` / `weights` | 0 | local, `0..n_spots-1` against `spots_to_keep` |
//! | `labels` | 0 | local |
//! | coordinates | n/a | full-res pixels, `f32`, `(x, y)` |
//!
//! `MoransI::new` asserts `neighbours.len() == spots_to_keep.len()`, which is
//! what makes the neighbour lists local rather than global spot ids.
//! `build_spatial_graph` already returns local 0-based indices, so the two
//! compose directly.

use bixverse_rs::spatial::sp_graph::SpatialAdjacency;
use extendr_api::*;

/////////////////
// Coordinates //
/////////////////

/// Turn an `n x 2` R matrix into the coordinate pairs the crate expects.
///
/// Column 1 is `x`, column 2 is `y`, both in full-res pixels. R matrices are
/// column-major, so the two columns are contiguous blocks.
///
/// ### Params
///
/// * `coords` - Numeric matrix with two columns, one row per spot.
///
/// ### Returns
///
/// The `(x, y)` pairs as `f32`, in row order.
pub fn coords_from_r_matrix(coords: &RMatrix<f64>) -> Result<Vec<(f32, f32)>> {
    if coords.ncols() != 2 {
        return Err(Error::Other(format!(
            "The coordinate matrix needs exactly two columns (x, y), got {}.",
            coords.ncols()
        )));
    }

    let n = coords.nrows();
    let data = coords.data();

    Ok((0..n)
        .map(|i| (data[i] as f32, data[i + n] as f32))
        .collect())
}

///////////////
// Adjacency //
///////////////

/// Turn an R list of integer vectors into local 0-based neighbour lists.
///
/// ### Params
///
/// * `r_list` - One integer vector per spot, holding 0-based local indices.
///
/// ### Returns
///
/// The neighbour lists.
pub fn neighbours_from_r_list(r_list: &List) -> Result<Vec<Vec<usize>>> {
    (0..r_list.len())
        .map(|i| {
            let elt = r_list.elt(i)?;
            let ints = elt.as_integer_slice().ok_or_else(|| {
                Error::Other(format!("Neighbour element {i} is not an integer vector."))
            })?;
            Ok(ints.iter().map(|&x| x as usize).collect())
        })
        .collect()
}

/// Turn an R list of numeric vectors into edge weights.
///
/// ### Params
///
/// * `r_list` - One numeric vector per spot, aligned with the neighbour lists.
///
/// ### Returns
///
/// The weights as `f32`.
pub fn weights_from_r_list(r_list: &List) -> Result<Vec<Vec<f32>>> {
    (0..r_list.len())
        .map(|i| {
            let elt = r_list.elt(i)?;
            let reals = elt.as_real_slice().ok_or_else(|| {
                Error::Other(format!("Weight element {i} is not a numeric vector."))
            })?;
            Ok(reals.iter().map(|&x| x as f32).collect())
        })
        .collect()
}

/// Read an adjacency pair off two R lists, checking they line up.
///
/// ### Params
///
/// * `indices` - Neighbour lists, 0-based local indices.
/// * `weights` - Edge weights, aligned element for element.
///
/// ### Returns
///
/// The `(neighbours, weights)` pair.
pub fn adjacency_from_r_lists(indices: &List, weights: &List) -> Result<SpatialAdjacency> {
    let neighbours = neighbours_from_r_list(indices)?;
    let weights = weights_from_r_list(weights)?;

    if neighbours.len() != weights.len() {
        return Err(Error::Other(format!(
            "`indices` ({}) and `weights` ({}) must have the same length.",
            neighbours.len(),
            weights.len()
        )));
    }
    for (i, (n, w)) in neighbours.iter().zip(weights.iter()).enumerate() {
        if n.len() != w.len() {
            return Err(Error::Other(format!(
                "Element {i}: {} neighbours but {} weights.",
                n.len(),
                w.len()
            )));
        }
    }

    Ok((neighbours, weights))
}

/////////////
// Returns //
/////////////

/// Pack neighbour lists into an R list of integer vectors.
///
/// ### Params
///
/// * `neighbours` - Local 0-based neighbour lists.
///
/// ### Returns
///
/// The R list.
pub fn neighbours_to_r_list(neighbours: &[Vec<usize>]) -> List {
    neighbours
        .iter()
        .map(|nb| nb.iter().map(|&x| x as i32).collect::<Vec<i32>>())
        .collect::<List>()
}

/// Pack edge weights into an R list of numeric vectors.
///
/// ### Params
///
/// * `weights` - Edge weights, aligned with the neighbour lists.
///
/// ### Returns
///
/// The R list.
pub fn weights_to_r_list(weights: &[Vec<f32>]) -> List {
    weights
        .iter()
        .map(|w| w.iter().map(|&x| x as f64).collect::<Vec<f64>>())
        .collect::<List>()
}

/// Pack a gene-major flat matrix into an R matrix.
///
/// The crate stores per-(gene, kernel) results row-major with the gene as the
/// slow axis. R wants column-major, so the transposition happens here rather
/// than leaving a flat vector and a shape convention for the caller to get
/// wrong.
///
/// ### Params
///
/// * `values` - Flat `n_rows * n_cols`, row-major.
/// * `n_rows` - Number of rows.
/// * `n_cols` - Number of columns.
///
/// ### Returns
///
/// An `n_rows x n_cols` R matrix.
pub fn row_major_to_r_matrix(values: &[f64], n_rows: usize, n_cols: usize) -> RMatrix<f64> {
    RMatrix::new_matrix(n_rows, n_cols, |r, c| values[r * n_cols + c])
}
