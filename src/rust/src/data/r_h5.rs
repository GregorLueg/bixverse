use extendr_api::prelude::*;

use std::mem::discriminant;

use bixverse_r::prelude::*;
use bixverse_r::single_cell::sc_data::h5_10x_io::{
    detect_tenx_version, get_tenx_dimensions, write_tenx_h5, TenxFeatures, TenxVersion,
};
use bixverse_r::single_cell::sc_data::h5_scan::{detect_raw_count_slot, get_h5ad_dimensions};
use bixverse_r::single_cell::sc_data::h5ad_frame_io::{
    list_h5_group, read_h5_attrs, read_h5_column, read_h5_matrix, read_h5ad_frame,
    write_h5ad_plain, ColumnData, PlainX,
};
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
    fn rs_h5ad_frame;
    fn rs_h5_read;
    fn rs_h5_attrs;
    fn rs_h5_ls;
    fn rs_tenx_h5_metadata;
    fn rs_write_h5ad_plain;
}

/////////////
// Helpers //
/////////////

/// Turn a decoded column into an R vector; categoricals become factors
///
/// ### Params
///
/// * `col` - The column
///
/// ### Returns
///
/// The R vector.
fn column_to_r(col: ColumnData) -> Robj {
    match col {
        ColumnData::Integer(v) => Integers::from_values(v.into_iter().map(Rint::from)).into(),
        ColumnData::Double(v) => Doubles::from_values(v.into_iter().map(Rfloat::from)).into(),
        ColumnData::Boolean(v) => Logicals::from_values(v.into_iter().map(|x| match x {
            Some(b) => Rbool::from(b),
            None => Rbool::na(),
        }))
        .into(),
        ColumnData::String(v) => Strings::from_values(v.into_iter().map(|x| match x {
            Some(s) => Rstr::from(s),
            None => Rstr::na(),
        }))
        .into(),
        ColumnData::Categorical { codes, categories } => {
            let n = categories.len() as i32;
            let codes = Integers::from_values(codes.into_iter().map(|c| {
                if (0..n).contains(&c) {
                    Rint::from(c + 1)
                } else {
                    Rint::na()
                }
            }));
            let mut x: Robj = codes.into();
            x.set_attrib("levels", categories).unwrap();
            x.set_class(&["factor"]).unwrap();
            x
        }
    }
}

/// Named R list from `(name, column)` pairs
///
/// ### Params
///
/// * `cols` - The columns
///
/// ### Returns
///
/// The named list.
fn columns_to_r(cols: Vec<(String, ColumnData)>) -> List {
    let (names, vals): (Vec<String>, Vec<Robj>) =
        cols.into_iter().map(|(n, c)| (n, column_to_r(c))).unzip();
    let mut l = List::from_values(vals);
    l.set_names(names).unwrap();
    l
}

/// Turn a named R list of vectors into columns
///
/// ### Params
///
/// * `cols` - Named list of integer, double, logical, character or factor
///   vectors
///
/// ### Returns
///
/// The columns in list order.
fn r_to_columns(cols: List) -> extendr_api::Result<Vec<(String, ColumnData)>> {
    cols.iter()
        .map(|(name, x)| {
            let col = match x.rtype() {
                Rtype::Integers if x.inherits("factor") => ColumnData::Categorical {
                    codes: x
                        .as_integer_slice()
                        .unwrap()
                        .iter()
                        .map(|&c| if c == i32::MIN { -1 } else { c - 1 })
                        .collect(),
                    categories: x
                        .get_attrib("levels")
                        .and_then(|l| l.as_string_vector())
                        .unwrap_or_default(),
                },
                Rtype::Integers => ColumnData::Integer(
                    x.as_integer_slice()
                        .unwrap()
                        .iter()
                        .map(|&c| (c != i32::MIN).then_some(c))
                        .collect(),
                ),
                Rtype::Doubles => ColumnData::Double(
                    x.as_real_iter()
                        .unwrap()
                        .map(|v| (!v.is_na()).then_some(*v))
                        .collect(),
                ),
                Rtype::Logicals => ColumnData::Boolean(
                    x.as_logical_iter()
                        .unwrap()
                        .map(|b| (!b.is_na()).then(|| b.is_true()))
                        .collect(),
                ),
                Rtype::Strings => ColumnData::String(
                    x.as_str_iter()
                        .unwrap()
                        .map(|s| (!s.is_na()).then(|| s.to_string()))
                        .collect(),
                ),
                other => {
                    return Err(Error::Other(format!(
                        "Column '{name}' has unsupported type {other:?}."
                    )))
                }
            };
            Ok((name.to_string(), col))
        })
        .collect()
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
            parse_raw_slot(c).ok_or_else(|| Error::Other(format!("Unknown raw count slot '{c}'.")))
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

/// Read an h5ad dataframe group into R
///
/// @description
/// Decodes plain, categorical and nullable columns, old-style
/// `__categories` and the `_index` attribute or dataset.
///
/// @param f_path String. Path to the h5ad file.
/// @param group String. The dataframe group, e.g. `"obs"`.
///
/// @returns A list with `cols` (named list of vectors, categoricals as
/// factors, in HDF5 name order with column groups first), `idx` (the index,
/// or `NULL`) and `skipped` (columns with an unsupported encoding).
///
/// @export
///
/// @keywords internal
#[extendr]
fn rs_h5ad_frame(f_path: &str, group: &str) -> extendr_api::Result<List> {
    let frame = read_h5ad_frame(f_path, group).to_extendr()?;
    let idx: Robj = frame.idx.map(column_to_r).unwrap_or(().into());
    Ok(list!(
        cols = columns_to_r(frame.cols),
        idx = idx,
        skipped = frame.skipped
    ))
}

/// Read one HDF5 dataset into R
///
/// @description
/// 1-D datasets and encoded columns (categorical, nullable) come back as
/// vectors. A 2-D numeric dataset comes back as a double matrix with the
/// dimensions reversed (the on-disk row-major `rows x cols` read as R's
/// column-major `cols x rows`).
///
/// @param f_path String. Path to the HDF5 file.
/// @param path String. Path of the dataset inside the file.
///
/// @returns The vector or matrix.
///
/// @export
///
/// @keywords internal
#[extendr]
fn rs_h5_read(f_path: &str, path: &str) -> extendr_api::Result<Robj> {
    if let Some(col) = read_h5_column(f_path, path).to_extendr()? {
        return Ok(column_to_r(col));
    }
    let (data, rows, cols) = read_h5_matrix(f_path, path).to_extendr()?;
    let mut x: Robj = Doubles::from_values(data).into();
    x.set_attrib("dim", [cols as i32, rows as i32])?;
    Ok(x)
}

/// Read the attributes of an HDF5 group or dataset into R
///
/// @param f_path String. Path to the HDF5 file.
/// @param path String. Path of the group or dataset; `"/"` for the root.
///
/// @returns A named list of vectors. Unsupported dtypes are left out.
///
/// @export
///
/// @keywords internal
#[extendr]
fn rs_h5_attrs(f_path: &str, path: &str) -> extendr_api::Result<List> {
    Ok(columns_to_r(read_h5_attrs(f_path, path).to_extendr()?))
}

/// Version and dimensions of a 10x CellRanger h5 file
///
/// @param f_path String. Path to the 10x h5 file.
///
/// @returns A list with `version` (`"v2"` or `"v3"`), `n_genes` and
/// `n_cells`.
///
/// @export
///
/// @keywords internal
#[extendr]
fn rs_tenx_h5_metadata(f_path: &str) -> extendr_api::Result<List> {
    let version = detect_tenx_version(f_path).to_extendr()?;
    let (n_cells, n_genes) = get_tenx_dimensions(f_path, Some(version)).to_extendr()?;
    let version = match version {
        TenxVersion::V3 => "v3",
        TenxVersion::V2 => "v2",
    };
    Ok(list!(
        version = version,
        n_genes = n_genes as i32,
        n_cells = n_cells as i32
    ))
}

/// Write a minimal h5ad file
///
/// @description
/// Writes `X` plus `obs` / `var` as plain datasets without anndata encoding
/// attributes. The first column of `obs` and `var` becomes `_index`.
///
/// @param f_path String. Output file, overwritten.
/// @param x_type String. One of `c("csr", "csc", "dense")`.
/// @param data Numeric. Non-zero values, or the dense values cells x genes in
/// row-major order.
/// @param indices Integer. 0-based minor indices; ignored for dense.
/// @param indptr Integer. 0-based major offsets; ignored for dense.
/// @param n_obs Integer. Number of cells.
/// @param n_var Integer. Number of genes.
/// @param obs Named list of cell columns, identifier first.
/// @param var Named list of gene columns, identifier first.
///
/// @export
///
/// @keywords internal
#[extendr]
#[allow(clippy::too_many_arguments)]
fn rs_write_h5ad_plain(
    f_path: &str,
    x_type: &str,
    data: &[f64],
    indices: &[i32],
    indptr: &[i32],
    n_obs: usize,
    n_var: usize,
    obs: List,
    var: List,
) -> extendr_api::Result<()> {
    let x = match x_type {
        "csr" => PlainX::Csr {
            data,
            indices,
            indptr,
        },
        "csc" => PlainX::Csc {
            data,
            indices,
            indptr,
        },
        "dense" => PlainX::Dense(data),
        other => {
            return Err(Error::Other(format!(
                "x_type must be one of 'csr', 'csc' or 'dense', got '{other}'."
            )))
        }
    };
    write_h5ad_plain(
        f_path,
        x,
        n_obs,
        n_var,
        &r_to_columns(obs)?,
        &r_to_columns(var)?,
    )
    .to_extendr()
}

/// List the members of an HDF5 group
///
/// @param f_path String. Path to the HDF5 file.
/// @param group String. Path of the group; `"/"` for the root.
///
/// @returns A list with `name`, `otype` (`"group"` or `"dataset"`) and `rank`
/// (number of dimensions, 0 for groups and scalars).
///
/// @export
///
/// @keywords internal
#[extendr]
fn rs_h5_ls(f_path: &str, group: &str) -> extendr_api::Result<List> {
    let members = list_h5_group(f_path, group).to_extendr()?;
    let otype: Vec<&str> = members
        .iter()
        .map(|m| if m.is_group { "group" } else { "dataset" })
        .collect();
    let rank: Vec<i32> = members.iter().map(|m| m.shape.len() as i32).collect();
    let name: Vec<String> = members.into_iter().map(|m| m.name).collect();
    Ok(list!(name = name, otype = otype, rank = rank))
}
