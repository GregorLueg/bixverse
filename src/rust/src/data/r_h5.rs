use extendr_api::prelude::*;

use std::mem::discriminant;

use bixverse_r::prelude::*;
use bixverse_r::single_cell::sc_data::h5_10x_io::{write_tenx_h5, TenxFeatures, TenxVersion};
use bixverse_r::single_cell::sc_data::h5_scan::{detect_raw_count_slot, get_h5ad_dimensions};
use bixverse_r::single_cell::sc_data::h5ad_io::*;

/////////////
// extendR //
/////////////

extendr_module! {
    mod r_h5;
    fn rs_h5ad_data;
    fn rs_h5ad_dimensions;
    fn rs_h5ad_raw_slot;
    fn rs_write_tenx_h5;
}

///////////////
// Fucntions //
///////////////

/// Load in h5ad data via Rust
///
/// @description
/// `r lifecycle::badge("experimental")`
/// Loads in h5ad data within Rust and automatically converts the data into
/// CSR with cells x genes.
///
/// @param f_path File path. The path to the h5ad file.
/// @param cs_type String. One of `c("csr", "csc")`. How the data is stored in
/// the file. Other values raise an error.
/// @param nrows Integer. Number of rows in the file.
/// @param ncols Integer. Number of columns in the file.
/// @param cell_quality List. Specifying the cell quality. Please refer
/// to [bixverse::params_sc_min_quality()].
/// @param slot String. In which slot the raw data can be found. One of
/// `c("X", "raw", "layers.counts")`. Unknown strings default to `"X"`.
/// @param verbose Boolean. Controls verbosity of the function
///
/// @returns A list with the CSR data (cells x genes) of the cells and genes
/// passing `cell_quality`:
/// \itemize{
///   \item data - The counts of the sparse matrix.
///   \item indices - The 0-based gene indices of the sparse matrix.
///   \item indptr - The index pointers of the sparse matrix.
///   \item no_genes - No of genes in the sparse matrix (i.e., ncol).
///   \item no_cells - No of cells in the sparse matrix (i.e., nrow).
/// }
///
/// @export
#[extendr]
fn rs_h5ad_data(
    f_path: String,
    cs_type: String,
    nrows: usize,
    ncols: usize,
    cell_quality: List,
    slot: String,
    verbose: bool,
) -> extendr_api::Result<List> {
    let cell_quality = MinCellQuality::from_r_list(cell_quality)?;

    let raw_slot = parse_raw_slot(&slot).unwrap_or_else(|| {
        println!(
            "The provided string ({:?}) could not be matched. Defaulting to X",
            slot
        );
        RawDataSlot::default()
    });

    let file_format = parse_compressed_sparse_format(&cs_type)
        .ok_or_else(|| ScDataErrors::UnknownSparseFormat(cs_type.to_string()))
        .to_extendr()?;

    let file_quality = match file_format {
        CompressedSparseFormat::Csr => {
            parse_h5_csr_quality(&f_path, (nrows, ncols), &raw_slot, &cell_quality, verbose)
                .to_extendr()?
        }
        CompressedSparseFormat::Csc => {
            parse_h5_csc_quality(&f_path, (nrows, ncols), &cell_quality, &raw_slot, verbose)
                .to_extendr()?
        }
    };

    let file_data: CompressedSparseData2<u32> = match file_format {
        CompressedSparseFormat::Csr => {
            read_h5ad_x_data_csr(&f_path, &file_quality, &raw_slot, verbose).to_extendr()?
        }
        CompressedSparseFormat::Csc => {
            let data =
                read_h5ad_x_data_csc(&f_path, &file_quality, &raw_slot, verbose).to_extendr()?;
            data.transpose_and_convert()
        }
    };

    Ok(list!(
        data = file_data
            .data
            .iter()
            .map(|x| *x as i32)
            .collect::<Vec<i32>>(),
        indices = file_data
            .indices
            .iter()
            .map(|x| *x as i32)
            .collect::<Vec<i32>>(),
        indptr = file_data
            .indptr
            .iter()
            .map(|x| *x as i32)
            .collect::<Vec<i32>>(),
        no_genes = file_quality.genes_to_keep.len(),
        no_cells = file_quality.cells_to_keep.len()
    ))
}

/// Dimensions and storage layout of an h5ad file
///
/// @description
/// `r lifecycle::badge("experimental")`
/// Reads the obs and var lengths and how `X` is stored, without touching the
/// counts. Sparse `X` takes its layout from the `encoding-type` attribute,
/// dense `X` from comparing its shape against the obs and var lengths.
///
/// @param f_path String. Path to the h5ad file.
///
/// @returns A list with
/// \itemize{
///   \item n_obs - Integer. Number of cells (rows of obs).
///   \item n_var - Integer. Number of genes (rows of var).
///   \item format - String. One of `c("CSR", "CSC", "DENSE_ROW",
///   "DENSE_COL")`.
/// }
///
/// @export
///
/// @keywords internal
#[extendr]
fn rs_h5ad_dimensions(f_path: &str) -> extendr_api::Result<List> {
    let dims = get_h5ad_dimensions(f_path).to_extendr()?;

    Ok(list!(
        n_obs = dims.n_obs,
        n_var = dims.n_var,
        format = dims.format.as_str()
    ))
}

/// Detect which slot of an h5ad file holds raw integer counts
///
/// @description
/// `r lifecycle::badge("experimental")`
/// Samples each candidate slot in turn and returns the first whose non-zero
/// values are whole numbers at a rate of at least `threshold`. A missing slot,
/// or one sampled to all zeros, is skipped.
///
/// @param f_path String. Path to the h5ad file.
/// @param candidates Character vector. Slots to test, in priority order. Any of
/// `c("layers.counts", "raw.X", "X")`.
/// @param n_sample Integer. Values to sample per slot.
/// @param threshold Numeric. Minimum fraction of whole numbers among the
/// non-zero values.
///
/// @returns The first qualifying entry of `candidates`, or `NULL`.
///
/// @export
///
/// @keywords internal
#[extendr]
fn rs_h5ad_raw_slot(
    f_path: &str,
    candidates: Vec<String>,
    n_sample: usize,
    threshold: f64,
) -> extendr_api::Result<Robj> {
    let slots = candidates
        .iter()
        .map(|c| {
            parse_raw_slot(c)
                .ok_or_else(|| Error::Other(format!("Unknown raw count slot '{c}'.")))
        })
        .collect::<extendr_api::Result<Vec<RawDataSlot>>>()?;

    let found = detect_raw_count_slot(f_path, &slots, n_sample, threshold).to_extendr()?;

    // RawDataSlot has no PartialEq; the variant identifies the candidate
    Ok(found
        .and_then(|hit| {
            slots
                .iter()
                .position(|s| discriminant(s) == discriminant(&hit))
        })
        .map_or_else(|| NULL.into(), |i| candidates[i].clone().into()))
}

/// Write a cells x features count matrix as a CellRanger style h5 file
///
/// @description
/// `r lifecycle::badge("experimental")`
/// The cells x features CSR buffers are the features x cells CSC buffers 10x
/// stores, so they are written as they are. Counts, indices, pointers and
/// shape go out as `int32`, strings as variable-length UTF-8. An existing file
/// is truncated.
///
/// @param f_path String. Path of the file to create.
/// @param version String. One of `c("v3", "v2")`.
/// @param indptr Integer vector. Row pointers of the cells x features CSR
/// matrix (0-indexed).
/// @param indices Integer vector. Feature indices (0-indexed).
/// @param data Integer vector. The counts. Must be non-negative.
/// @param n_cells Integer. Number of cells.
/// @param n_features Integer. Number of features.
/// @param barcodes Character vector. One barcode per cell.
/// @param feature_ids Character vector. One id per feature.
/// @param feature_names Character vector. One name per feature.
/// @param feature_types Optional character vector. One type per feature, v3
/// only. `NULL` writes `"Gene Expression"` for every feature.
///
/// @returns Invisible `NULL`.
///
/// @export
///
/// @keywords internal
#[extendr]
#[allow(clippy::too_many_arguments)]
fn rs_write_tenx_h5(
    f_path: &str,
    version: &str,
    indptr: Vec<i32>,
    indices: Vec<i32>,
    data: Vec<i32>,
    n_cells: usize,
    n_features: usize,
    barcodes: Vec<String>,
    feature_ids: Vec<String>,
    feature_names: Vec<String>,
    feature_types: Nullable<Vec<String>>,
) -> extendr_api::Result<()> {
    let version = match version {
        "v3" => TenxVersion::V3,
        "v2" => TenxVersion::V2,
        other => {
            return Err(Error::Other(format!(
                "Unknown 10x version '{other}'. Use 'v3' or 'v2'."
            )))
        }
    };

    if data.iter().any(|&x| x < 0) {
        return Err(Error::Other("Counts must be non-negative.".into()));
    }

    let data: Vec<u32> = data.iter().map(|&x| x as u32).collect();
    let indices: Vec<u32> = indices.iter().map(|&x| x as u32).collect();
    let indptr: Vec<u32> = indptr.iter().map(|&x| x as u32).collect();
    let counts: CompressedSparseData2<u32> =
        CompressedSparseData2::new_csr(&data, &indices, &indptr, None, (n_cells, n_features));

    let feature_types = match feature_types {
        Nullable::NotNull(x) => Some(x),
        Nullable::Null => None,
    };
    let features = TenxFeatures {
        ids: &feature_ids,
        names: &feature_names,
        feature_types: feature_types.as_deref(),
    };

    write_tenx_h5(f_path, version, &counts, &barcodes, &features).to_extendr()
}
