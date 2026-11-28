//! Header-level helpers of the v5 single cell store: feature and norm kind,
//! caller metadata, checksums and the cell-ID hash.

use bixverse_r::prelude::*;
use bixverse_r::single_cell::sc_data::store::v5::{
    check_store_pair_paths, peek_store_info, set_cell_ids, set_store_metadata, verify_store,
    NormKind,
};
use extendr_api::prelude::*;
use extendr_api::Result;

////////////////////
// extendr Module //
////////////////////

extendr_module! {
    mod r_sc_store;
    fn rs_store_feature_kind;
    fn rs_store_norm_kind;
    fn rs_store_metadata;
    fn rs_set_store_metadata;
    fn rs_verify_store;
    fn rs_check_store_pair;
    fn rs_set_cell_ids;
    fn rs_cell_id_hash;
}

///////////
// Kinds //
///////////

/// Read what the features of a store are
///
/// @description
/// `r lifecycle::badge("experimental")`
/// Reads the header only.
///
/// @param path String. A cell- or gene-based binary (v4 or v5).
///
/// @returns One of `"gene"`, `"peak"`, `"tile"`, `"snp"`, `"protein"` or
/// `"gene_activity"`. A v4 file is `"gene"`.
///
/// @export
///
/// @keywords internal
#[extendr]
fn rs_store_feature_kind(path: &str) -> Result<String> {
    let info = peek_store_info(path).to_extendr()?;
    Ok(info.feature_kind.name().to_string())
}

/// Read how the norm of a store derives
///
/// @description
/// `r lifecycle::badge("experimental")`
/// Reads the header only.
///
/// @param path String. A cell- or gene-based binary (v4 or v5).
///
/// @returns One of `"none"`, `"log_cpm"` or `"tf_idf"`. A v4 file is
/// `"log_cpm"`.
///
/// @export
///
/// @keywords internal
#[extendr]
fn rs_store_norm_kind(path: &str) -> Result<String> {
    let info = peek_store_info(path).to_extendr()?;
    Ok(match info.norm_kind {
        NormKind::None => "none",
        NormKind::LogCpm => "log_cpm",
        NormKind::TfIdf => "tf_idf",
    }
    .to_string())
}

//////////////
// Metadata //
//////////////

/// Caller metadata of a store
///
/// @description
/// `r lifecycle::badge("experimental")`
///
/// @param path String. A cell- or gene-based binary.
///
/// @returns A named character vector of the `key = value` pairs, empty for a
/// v4 file.
///
/// @export
///
/// @keywords internal
#[extendr]
fn rs_store_metadata(path: &str) -> Result<Robj> {
    let info = peek_store_info(path).to_extendr()?;
    let (keys, values): (Vec<String>, Vec<String>) = info.metadata.into_iter().unzip();
    let mut out: Robj = values.into_robj();
    out.set_names(keys)?;
    Ok(out)
}

/// Set caller metadata on both files of a v5 store
///
/// @description
/// `r lifecycle::badge("experimental")`
/// Given keys replace, all others are kept.
///
/// @param f_path_cells String. Cell-based binary.
/// @param f_path_genes String. Gene-based binary.
/// @param keys Character vector. The keys.
/// @param values Character vector. The values, same length as `keys`.
///
/// @returns Invisible `NULL`. Errors for a multi-line key or value, or a v4
/// store.
///
/// @export
///
/// @keywords internal
#[extendr]
fn rs_set_store_metadata(
    f_path_cells: &str,
    f_path_genes: &str,
    keys: Vec<String>,
    values: Vec<String>,
) -> Result<()> {
    if keys.len() != values.len() {
        return Err(Error::Other("keys and values differ in length".into()));
    }
    let pairs: Vec<(&str, &str)> = keys
        .iter()
        .zip(&values)
        .map(|(k, v)| (k.as_str(), v.as_str()))
        .collect();
    set_store_metadata(&[f_path_cells, f_path_genes], &pairs).to_extendr()
}

///////////////
// Integrity //
///////////////

/// Check every chunk of a store file against its checksums
///
/// @description
/// `r lifecycle::badge("experimental")`
///
/// @param path String. A cell- or gene-based binary.
///
/// @returns `TRUE` when every chunk matched, `FALSE` for a file without
/// checksums (v4). Errors on a corrupt chunk.
///
/// @export
///
/// @keywords internal
#[extendr]
fn rs_verify_store(path: &str) -> Result<bool> {
    verify_store(path).to_extendr()
}

/// Error unless a cell file and a gene file belong to the same store
///
/// @description
/// `r lifecycle::badge("experimental")`
/// Reads the headers only.
///
/// @param f_path_cells String. Cell-based binary.
/// @param f_path_genes String. Gene-based binary.
///
/// @returns Invisible `NULL`, or an error naming the field that differs.
///
/// @export
///
/// @keywords internal
#[extendr]
fn rs_check_store_pair(f_path_cells: &str, f_path_genes: &str) -> Result<()> {
    check_store_pair_paths(f_path_cells, f_path_genes).to_extendr()
}

//////////////
// Cell IDs //
//////////////

/// Record the cell IDs of a v5 store in both headers
///
/// @description
/// `r lifecycle::badge("experimental")`
/// Only their hash is stored. Stores that share a cell axis (a gene activity
/// store and its source) carry the same one.
///
/// @param f_path_cells String. Cell-based binary.
/// @param f_path_genes String. Gene-based binary.
/// @param ids Character vector. Cell IDs in store order.
///
/// @returns Invisible `NULL`. Errors for a store that is not v5 or a wrong
/// number of IDs.
///
/// @export
///
/// @keywords internal
#[extendr]
fn rs_set_cell_ids(f_path_cells: &str, f_path_genes: &str, ids: Vec<String>) -> Result<()> {
    set_cell_ids(&[f_path_cells, f_path_genes], &ids).to_extendr()
}

/// The cell-ID hash in a binary's header
///
/// @description
/// `r lifecycle::badge("experimental")`
///
/// @param path String. A cell- or gene-based binary.
///
/// @returns The 16-byte hash as hex, `NULL` when the store records no cell
/// IDs.
///
/// @export
///
/// @keywords internal
#[extendr]
fn rs_cell_id_hash(path: &str) -> Result<Option<String>> {
    let info = peek_store_info(path).to_extendr()?;
    Ok(info
        .cell_id_hash
        .map(|h| h.iter().map(|b| format!("{b:02x}")).collect()))
}
