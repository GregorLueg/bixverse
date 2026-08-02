//! Rust <> R interface for the spatial extras of an h5ad file.
//!
//! Only `obsm/spatial` and the `uns/spatial` group. The counts, `obs` and `var`
//! come out of `load_h5ad()` on the R side and never touch this.

use bixverse_rs::prelude::*;
use bixverse_rs::spatial::sp_data::{
    SpatialH5adData, SpatialH5adParams, SpatialOrientation, read_spatial_h5ad,
};
use extendr_api::*;

////////////////////
// extendr Module //
////////////////////

extendr_module! {
    mod r_sp_data;
    fn rs_sp_read_h5ad_spatial;
}

//////////////
// Reader //
//////////////

/// Read the spatial extras out of an h5ad file
///
/// @description
/// `r lifecycle::badge("experimental")`
/// Pulls `obsm/spatial` and whatever survived in `uns/spatial`, and resolves
/// which column of the coordinate array is `x`.
///
/// **The column order is not fixed by the AnnData spec.** It gets worked out
/// from the file rather than assumed, strongest evidence first:
///
/// \itemize{
///   \item `"image_tissue"` - one order puts the spots on the tissue in the
///   shipped image and the other does not. The only test that measures the
///   thing that has a consequence.
///   \item `"image_frame"` - one order keeps every spot inside the frame the
///   image implies, the other pushes spots off the edge.
///   \item `"obs_pixel_columns"` - `obs` carries `pxl_row_in_fullres` and
///   `pxl_col_in_fullres` and one matched a column exactly. Weaker than it
///   looks: `scanpy.read_visium` swaps those two names relative to Space
///   Ranger, so the labels mean the opposite in a file that went through it.
///   \item `"assumed"` - nothing in the file settled it, so `orientation` was
///   taken as given.
/// }
///
/// Getting it wrong costs nothing statistically. A swap of `x` and `y` is a
/// reflection, and every distance survives it, so the graph, Moran's I and
/// SPARK-X all come out identical. It costs everything for image features,
/// which would then be cut from the transpose of each spot.
///
/// @param h5_path String. Path to the `.h5ad` file.
/// @param library_id String or `NULL`. Which `uns/spatial` library to read.
/// `NULL` takes the only one and errors when there is a choice to make.
/// @param orientation String. `"xy"` or `"yx"`. The column order to fall back
/// on when nothing in the file settles it. `"xy"` is what
/// `scanpy.read_visium` produces.
///
/// @return A list with the following elements
/// \itemize{
///   \item coords - Numeric matrix, spots x 2. Columns `x` then `y`, in
///   full-resolution pixels, already in the contract order.
///   \item orientation - String. The column order that was used.
///   \item evidence - String. What settled it, see above.
///   \item obsm_keys - Character. Every key under `obsm`.
///   \item library_id - String. The `uns/spatial` library read, or `NA`.
///   \item library_ids - Character. Every library the file holds.
///   \item scale_factor_names, scale_factor_values - Character and numeric.
///   The scale factors, passed through untouched.
///   \item image_keys - Character. Images present under the library.
///   \item image_heights, image_widths - Integer. Their pixel dimensions,
///   aligned with `image_keys`.
///   \item metadata_keys - Character. Keys under the library's `metadata`.
///   \item has_array_indices - Boolean. Whether `obs` carries `array_row` and
///   `array_col`, which the hex and square graph layouts need.
///   \item has_pixel_columns - Boolean. Whether `obs` carries the
///   `pxl_*_in_fullres` pair.
/// }
///
/// @export
///
/// @keywords internal
#[extendr]
fn rs_sp_read_h5ad_spatial(
    h5_path: String,
    library_id: Nullable<String>,
    orientation: String,
) -> extendr_api::Result<List> {
    let assumed = SpatialOrientation::parse(&orientation).ok_or_else(|| {
        Error::Other(format!(
            "`orientation` must be 'xy' or 'yx', got '{orientation}'."
        ))
    })?;

    let params = SpatialH5adParams::new(
        match library_id {
            Nullable::NotNull(id) => Some(id),
            Nullable::Null => None,
        },
        assumed,
    );

    let data: SpatialH5adData = read_spatial_h5ad(&h5_path, Some(params)).to_extendr()?;

    let n = data.n_spots();
    let coords = RMatrix::new_matrix(n, 2, |r, c| {
        if c == 0 {
            data.coordinates[r].0
        } else {
            data.coordinates[r].1
        }
    });

    Ok(list!(
        coords = coords,
        orientation = data.orientation.as_str(),
        evidence = data.evidence.as_str(),
        obsm_keys = data.obsm_keys,
        library_id = match data.library_id {
            Some(id) => Robj::from(id),
            None => Robj::from(Rstr::na()),
        },
        library_ids = data.library_ids,
        scale_factor_names = data
            .scale_factors
            .iter()
            .map(|(k, _)| k.clone())
            .collect::<Vec<String>>(),
        scale_factor_values = data
            .scale_factors
            .iter()
            .map(|(_, v)| *v)
            .collect::<Vec<f64>>(),
        image_keys = data
            .images
            .iter()
            .map(|e| e.key.clone())
            .collect::<Vec<String>>(),
        image_heights = data
            .images
            .iter()
            .map(|e| e.height as i32)
            .collect::<Vec<i32>>(),
        image_widths = data
            .images
            .iter()
            .map(|e| e.width as i32)
            .collect::<Vec<i32>>(),
        metadata_keys = data.metadata_keys,
        has_array_indices = data.has_array_indices,
        has_pixel_columns = data.has_pixel_columns
    ))
}
