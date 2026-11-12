//! Parquet export of the obs / var metadata for DuckDB to ingest.

use bixverse_r::prelude::*;
use bixverse_r::single_cell::sc_data::h5_10x_io::parse_tenx_version;
use bixverse_r::single_cell::sc_data::obs_var_io::{
    h5ad_frame_to_parquet, parse_frame_kind, tenx_frame_to_parquet, FrameKind, ParquetColumnType,
    ParquetFrameInfo,
};
use extendr_api::prelude::*;

/////////////
// extendR //
/////////////

extendr_module! {
    mod r_sc_obs_var;
    fn rs_h5ad_frame_to_parquet;
    fn rs_tenx_frame_to_parquet;
}

/////////////
// Helpers //
/////////////

/// Parse the frame argument
///
/// ### Params
///
/// * `frame` - `"obs"` or `"var"`
///
/// ### Returns
///
/// The [`FrameKind`], or an error naming the bad value.
fn frame_from_r(frame: &str) -> extendr_api::Result<FrameKind> {
    parse_frame_kind(frame).ok_or_else(|| {
        Error::Other(format!(
            "frame must be one of 'obs' or 'var', got '{frame}'."
        ))
    })
}

/// Turn the optional 0-based R row selection into Rust indices
///
/// ### Params
///
/// * `keep` - Optional 0-based integer vector
///
/// ### Returns
///
/// The indices, or an error for a negative or missing entry.
fn keep_from_r(keep: Nullable<Vec<i32>>) -> extendr_api::Result<Option<Vec<usize>>> {
    match keep {
        Nullable::Null => Ok(None),
        Nullable::NotNull(k) => k
            .into_iter()
            .map(|x| {
                if x < 0 {
                    Err(Error::Other(format!(
                        "keep holds {x}; it must be 0-based and non-missing."
                    )))
                } else {
                    Ok(x as usize)
                }
            })
            .collect::<extendr_api::Result<Vec<usize>>>()
            .map(Some),
    }
}

/// Describe a written table for R
///
/// ### Params
///
/// * `info` - The [`ParquetFrameInfo`] of the written file
///
/// ### Returns
///
/// The R list documented on `rs_h5ad_frame_to_parquet()`.
fn info_to_r(info: ParquetFrameInfo) -> List {
    let names: Vec<String> = info.columns.iter().map(|c| c.name.clone()).collect();
    let types: Vec<&str> = info
        .columns
        .iter()
        .map(|c| match c.col_type {
            ParquetColumnType::Integer => "integer",
            ParquetColumnType::Double => "double",
            ParquetColumnType::Boolean => "boolean",
            ParquetColumnType::String => "string",
            ParquetColumnType::Categorical => "categorical",
        })
        .collect();
    let categories: Vec<Robj> = info
        .columns
        .into_iter()
        .map(|c| match c.categories {
            Some(cats) => cats.into(),
            None => r!(NULL),
        })
        .collect();

    list!(
        names = names,
        types = types,
        categories = List::from_values(categories),
        skipped = info.skipped,
        n_rows = info.n_rows as f64
    )
}

///////////////
// Functions //
///////////////

/// Write the obs or var dataframe of an h5ad file to Parquet
///
/// @description
/// `r lifecycle::badge("experimental")`
/// Writes the table the obs / var DuckDB tables are built from: `cell_idx` /
/// `gene_idx` (1-based after filtering), `cell_id` / `gene_id`, then the
/// remaining columns snake_cased. Categoricals are written as strings; their
/// levels come back for the caller to cast to an ENUM.
///
/// @param h5_path String. Path to the h5ad file.
/// @param frame String. One of `c("obs", "var")`.
/// @param keep Optional integer vector. 0-indexed(!) rows to keep, in output
/// order. `NULL` keeps all.
/// @param id_col Optional string. Column to use as identifier when the group
/// has no index; the first column otherwise.
/// @param out_path String. Parquet file to write.
///
/// @returns A list with
/// \itemize{
///   \item names - Character vector. Column names in file order.
///   \item types - Character vector. One of `"integer"`, `"double"`,
///   `"boolean"`, `"string"` or `"categorical"` per column.
///   \item categories - List. The levels in order for a categorical column,
///   `NULL` otherwise.
///   \item skipped - Character vector. Source columns left out for an
///   unsupported encoding.
///   \item n_rows - Numeric. Number of rows written.
/// }
///
/// @export
///
/// @keywords internal
#[extendr]
fn rs_h5ad_frame_to_parquet(
    h5_path: &str,
    frame: &str,
    keep: Nullable<Vec<i32>>,
    id_col: Nullable<String>,
    out_path: &str,
) -> extendr_api::Result<List> {
    let kind = frame_from_r(frame)?;
    let keep = keep_from_r(keep)?;
    let id_col = match id_col {
        Nullable::NotNull(s) => Some(s),
        Nullable::Null => None,
    };

    let info = h5ad_frame_to_parquet(h5_path, kind, keep.as_deref(), id_col.as_deref(), out_path)
        .to_extendr()?;

    Ok(info_to_r(info))
}

/// Write the barcodes or features of a 10x h5 file to Parquet
///
/// @description
/// `r lifecycle::badge("experimental")`
/// Obs is `cell_idx`, `cell_id`. Var is `gene_idx`, `gene_id`, `gene_name`,
/// plus `feature_type` for v3. Features are read in file order across all
/// modalities, so `keep` indexes into that order.
///
/// @param h5_path String. Path to the 10x h5 file.
/// @param version String. One of `c("v2", "v3")`.
/// @param frame String. One of `c("obs", "var")`.
/// @param keep Optional integer vector. 0-indexed(!) rows to keep, in output
/// order. `NULL` keeps all.
/// @param out_path String. Parquet file to write.
///
/// @returns See `rs_h5ad_frame_to_parquet()`.
///
/// @export
///
/// @keywords internal
#[extendr]
fn rs_tenx_frame_to_parquet(
    h5_path: &str,
    version: &str,
    frame: &str,
    keep: Nullable<Vec<i32>>,
    out_path: &str,
) -> extendr_api::Result<List> {
    let kind = frame_from_r(frame)?;
    let keep = keep_from_r(keep)?;
    let version = parse_tenx_version(version).ok_or_else(|| {
        Error::Other(format!(
            "version must be one of 'v2' or 'v3', got '{version}'."
        ))
    })?;

    let info =
        tenx_frame_to_parquet(h5_path, version, kind, keep.as_deref(), out_path).to_extendr()?;

    Ok(info_to_r(info))
}
