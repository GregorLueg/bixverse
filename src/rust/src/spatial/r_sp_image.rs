//! Rust <> R interface for histology image features.
//!
//! Needs the `spatial-image` feature on `bixverse-rs`, which in turn needs the
//! OpenSlide C library present at build and run time.

use bixverse_rs::prelude::*;
use bixverse_rs::spatial::sp_image::source::{image_metadata, open_image_source, ImageMetadata};
use bixverse_rs::spatial::sp_image::{spatial_image_features, SpatialImageParams, SpatialImageRes};
use extendr_api::*;
use std::time::Instant;

use crate::spatial::utils::{coords_from_r_matrix, row_major_to_r_matrix};

////////////////////
// extendr Module //
////////////////////

extendr_module! {
    mod r_sp_image;
    fn rs_sp_image_features;
    fn rs_sp_image_metadata;
}

//////////////
// Features //
//////////////

/// Per-spot histology image features
///
/// @description
/// `r lifecycle::badge("experimental")`
/// Cuts one tile per spot out of the histology image and reduces it to a fixed
/// feature vector: first-order statistics per colour channel, mean optical
/// density, and Haralick texture features off a grey-level co-occurrence
/// matrix. Reading and reduction are fused, so peak memory is one tile per
/// thread rather than one per spot.
///
/// Spots whose tile falls entirely off the image are dropped rather than
/// failing the run, which is why `spot_idx` can be shorter than the input.
///
/// The image backend is picked by probing the file, not by its extension.
/// OpenSlide gets first refusal; anything it declines is decoded whole.
///
/// @param image_path String. Path to the histology image.
/// @param image_scalef Float. Scale from full-resolution coordinates into the
/// image's own pixel space, e.g. `tissue_hires_scalef` out of
/// `scalefactors_json.json`. Use `1.0` if the file already is full-res. On a
/// slide this also selects the pyramid level.
/// @param coords Numeric matrix. Two columns, `x` then `y`, in
/// **full-resolution** pixels. One row per spot. The scale factor is applied
/// inside, do not pre-scale.
/// @param spot_diameter_fullres Float. Spot diameter in full-res pixels, i.e.
/// `spot_diameter_fullres` out of `scalefactors_json.json`.
/// @param image_params List. With the following elements:
/// \itemize{
///   \item tile_scale - Float. Multiplier on the spot diameter when cutting a
///   tile. `1.0` takes the spot, larger pulls in context.
///   \item glcm_levels - Integer. Grey levels the GLCM quantises to, 2 to 255.
///   \item glcm_offsets_dy, glcm_offsets_dx - Integer. Neighbour offsets, given
///   as two aligned vectors. Both or neither.
///   \item stain_haem, stain_eosin - Numeric, length 3 each. Stain basis for
///   colour deconvolution. Both or neither; defaults to Ruifrok-Johnston H&E.
/// }
/// @param verbose Integer. `0L` - quiet; `1L` - normal verbosity; `2L` -
/// detailed verbosity.
///
/// @return A list with the following elements
/// \itemize{
///   \item spot_idx - Integer. 0-based positional indices into `coords` for
///   the spots that produced features.
///   \item feature_names - Character. One name per column of `values`.
///   \item values - Numeric matrix, spots x features.
/// }
///
/// @references
/// Ruifrok & Johnston, Anal Quant Cytol Histol, 2001
///
/// Haralick, Shanmugam & Dinstein, IEEE Trans Syst Man Cybern, 1973
///
/// @export
///
/// @keywords internal
#[extendr]
fn rs_sp_image_features(
    image_path: String,
    image_scalef: f64,
    coords: RMatrix<f64>,
    spot_diameter_fullres: f64,
    image_params: List,
    verbose: usize,
) -> extendr_api::Result<List> {
    let verbosity = parse_verbosity_level(verbose);
    let coordinates = coords_from_r_matrix(&coords)?;
    let params = SpatialImageParams::from_r_list(image_params)?;

    // One scale factor covers both backends: it is the fallback scale for a
    // plain image and the target working scale for a slide pyramid.
    let source =
        open_image_source(&image_path, image_scalef as f32, image_scalef as f32).to_extendr()?;

    let start = Instant::now();

    let res: SpatialImageRes = spatial_image_features(
        &source,
        &coordinates,
        spot_diameter_fullres as f32,
        Some(params),
    )
    .to_extendr()?;

    let n_spots = res.n_spots();
    let n_features = res.n_features();

    if verbosity.normal_verbosity() {
        println!(
            "Image features: {} of {} spots, {} features in {:.2?}",
            n_spots,
            coordinates.len(),
            n_features,
            start.elapsed()
        );
    }

    let values: Vec<f64> = res.values.iter().map(|&x| x as f64).collect();

    Ok(list!(
        spot_idx = res
            .spot_indices
            .iter()
            .map(|&x| x as i32)
            .collect::<Vec<i32>>(),
        feature_names = res.feature_names,
        values = row_major_to_r_matrix(&values, n_spots, n_features)
    ))
}

//////////////
// Metadata //
//////////////

/// Describe a histology image without decoding it
///
/// @description
/// `r lifecycle::badge("experimental")`
/// Reads the header, or the pyramid table on a whole-slide image. Cheap enough
/// to call before deciding which resolution to work at. A plain PNG or JPEG
/// reports one level, a downsample of `1.0` and no vendor, so the shape is the
/// same either way.
///
/// @param image_path String. Path to the image.
///
/// @return A list with the following elements
/// \itemize{
///   \item width - Integer. Width of level 0, in pixels.
///   \item height - Integer. Height of level 0, in pixels.
///   \item n_levels - Integer. Number of pyramid levels. `1L` for a plain
///   image.
///   \item level_dims - Integer matrix, levels x 2. Width then height.
///   \item downsamples - Numeric. Downsample per level, relative to level 0.
///   \item vendor - String. OpenSlide vendor, `NA` for a plain image.
/// }
///
/// @export
///
/// @keywords internal
#[extendr]
fn rs_sp_image_metadata(image_path: String) -> extendr_api::Result<List> {
    let meta: ImageMetadata = image_metadata(&image_path).to_extendr()?;

    let n_levels = meta.n_levels;
    let level_dims = RMatrix::new_matrix(n_levels, 2, |r, c| {
        if c == 0 {
            meta.level_dims[r].0 as i32
        } else {
            meta.level_dims[r].1 as i32
        }
    });

    Ok(list!(
        width = meta.width as i32,
        height = meta.height as i32,
        n_levels = n_levels as i32,
        level_dims = level_dims,
        downsamples = meta.downsamples,
        vendor = match meta.vendor {
            Some(v) => Robj::from(v),
            None => Robj::from(Rstr::na()),
        }
    ))
}
