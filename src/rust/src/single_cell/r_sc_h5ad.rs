use extendr_api::*;
use std::collections::HashMap;
use std::path::Path;

use bixverse_rs::prelude::*;
use scx_core::api::write::{H5AdBuilder, H5AdOptions};
use scx_core::ir::{
    Column, ColumnData, DenseMatrix, Embeddings, ObsTable, UnsTable, VarTable, Varm,
};
use sprs::CsMatI;

use crate::single_cell::utils::norm_counts_to_f32;

/////////////
// extendR //
/////////////

extendr_module! {
    mod r_sc_h5ad;
    fn rs_save_h5ad;
}

/////////////
// Helpers //
/////////////

/// Error unless a length matches the axis it belongs to
///
/// ### Params
///
/// * `what` - Name of the input, for the error message.
/// * `got` - Its length.
/// * `expected` - The length of the axis.
///
/// ### Returns
///
/// `Ok(())` if the lengths match.
fn check_len(what: &str, got: usize, expected: usize) -> Result<()> {
    if got == expected {
        Ok(())
    } else {
        Err(Error::Other(format!(
            "{what} has length {got}, expected {expected}."
        )))
    }
}

/// Unwrap the data of an R column, missing values included
///
/// ### Params
///
/// * `name` - Column name, for the error message.
/// * `data` - The result of one of the `as_*` accessors.
///
/// ### Returns
///
/// The data. Unlike the `TryFrom` conversions, which reject any `NA`, this
/// keeps them: `NA` codes, `NaN` and `"NA"` strings are all valid h5ad values.
fn na_ok<T>(name: &str, data: Option<T>) -> Result<T> {
    data.ok_or_else(|| Error::Other(format!("Column '{name}' could not be read.")))
}

/// Translate an R factor into an scx categorical
///
/// ### Params
///
/// * `name` - Column name, for error messages.
/// * `x` - The factor.
///
/// ### Returns
///
/// The categorical column data.
fn factor_to_categorical(name: &str, x: &Robj) -> Result<ColumnData> {
    let levels = x
        .get_attrib("levels")
        .and_then(|l| l.as_string_vector())
        .ok_or_else(|| Error::Other(format!("Factor '{name}' has no levels.")))?;
    let codes = na_ok(name, x.as_integer_slice())?
        .iter()
        // `NA_integer_` becomes u32::MAX, which scx narrows to -1, the anndata
        // code for a missing value, in every code width
        .map(|&c| {
            if c == i32::MIN {
                u32::MAX
            } else {
                (c - 1) as u32
            }
        })
        .collect();
    Ok(ColumnData::Categorical { codes, levels })
}

/// Translate a named R list of columns into scx columns
///
/// ### Params
///
/// * `cols` - Named list of factor, double, integer, logical or character
///   vectors. The R side decides the encoding; this only maps the types.
/// * `n` - Length every column must have.
///
/// ### Returns
///
/// The columns in list order.
fn r_list_to_columns(cols: List, n: usize) -> Result<Vec<Column>> {
    cols.iter()
        .map(|(name, x)| {
            let data = match x.rtype() {
                Rtype::Integers if x.inherits("factor") => factor_to_categorical(name, &x)?,
                Rtype::Integers => ColumnData::Int(na_ok(name, x.as_integer_slice())?.to_vec()),
                Rtype::Doubles => ColumnData::Float(na_ok(name, x.as_real_slice())?.to_vec()),
                Rtype::Logicals => ColumnData::Bool(
                    na_ok(name, x.as_logical_slice())?
                        .iter()
                        .map(|b| b.is_true())
                        .collect(),
                ),
                Rtype::Strings => ColumnData::String(na_ok(name, x.as_string_vector())?),
                _ => {
                    return Err(Error::Other(format!(
                        "Column '{name}' has a type with no h5ad encoding."
                    )))
                }
            };
            check_len(&format!("Column '{name}'"), data.len(), n)?;
            Ok(Column {
                name: name.to_string(),
                data,
            })
        })
        .collect()
}

/// Translate a named R list of double matrices into row-major scx matrices
///
/// ### Params
///
/// * `mats` - Named list of numeric matrices.
/// * `n` - Number of rows every matrix must have.
///
/// ### Returns
///
/// A map from name to the row-major matrix.
fn r_list_to_dense(mats: List, n: usize) -> Result<HashMap<String, DenseMatrix>> {
    mats.iter()
        .map(|(name, x)| {
            let mat = RMatrix::<f64>::try_from(x)?;
            let (nrow, ncol) = (mat.nrows(), mat.ncols());
            check_len(&format!("Matrix '{name}'"), nrow, n)?;
            let col_major = mat.data();
            let data = (0..nrow)
                .flat_map(|i| (0..ncol).map(move |j| col_major[i + j * nrow]))
                .collect();
            Ok((
                name.to_string(),
                DenseMatrix {
                    shape: (nrow, ncol),
                    data,
                },
            ))
        })
        .collect()
}

/// Translate an R CSR list into a square sprs matrix
///
/// ### Params
///
/// * `name` - Matrix name, for error messages.
/// * `x` - List with `indptr`, `indices` (0-indexed) and `data`.
/// * `n` - Number of rows and columns.
///
/// ### Returns
///
/// The validated CSR matrix.
fn r_list_to_csr(name: &str, x: Robj, n: usize) -> Result<CsMatI<f32, u32>> {
    let csr = List::try_from(x)?;
    let get_u32 = |key: &str| -> Result<Vec<u32>> {
        Ok(<&[i32]>::try_from(&csr.dollar(key)?)?
            .iter()
            .map(|&i| i as u32)
            .collect())
    };
    let data = <&[f64]>::try_from(&csr.dollar("data")?)?
        .iter()
        .map(|&v| v as f32)
        .collect();
    CsMatI::try_new((n, n), get_u32("indptr")?, get_u32("indices")?, data)
        .map_err(|(_, _, _, e)| Error::Other(format!("obsp '{name}': {e}")))
}

//////////
// h5ad //
//////////

/// Write a single cell experiment to h5ad
///
/// @description
/// Streams the counts out of the cell-based binary file, cell batch by cell
/// batch, into a spec-compliant h5ad file via
/// [scx-core](https://github.com/btraven00/scx), together with the obs and
/// var tables, dense embeddings and sparse cell x cell graphs supplied from
/// R. Counts are stored as `float32`. All inputs are checked before the file
/// is created.
///
/// @param f_path_cells String. Path to the `counts_cells.bin` file.
/// @param h5_path String. Path of the h5ad file to create.
/// @param cell_indices Integer vector. The cells to write (0-indexed!), in
/// the order they shall appear in the file.
/// @param norm Boolean. Write the normalised instead of the raw counts.
/// @param obs_index Character vector. One name per cell.
/// @param obs Named list. The obs columns, see the R wrapper for the types.
/// @param var_index Character vector. One name per gene.
/// @param var Named list. The var columns.
/// @param obsm Named list of numeric matrices with one row per cell.
/// @param varm Named list of numeric matrices with one row per gene.
/// @param obsp Named list of CSR matrices, each a list with `indptr`,
/// `indices` (0-indexed) and `data`.
/// @param uns_json String. JSON object written to `uns`.
/// @param chunk_size Integer. Number of cells per streaming batch.
///
/// @returns Invisible `NULL`.
///
/// @export
///
/// @keywords internal
#[extendr]
#[allow(clippy::too_many_arguments)]
fn rs_save_h5ad(
    f_path_cells: &str,
    h5_path: &str,
    cell_indices: &[i32],
    norm: bool,
    obs_index: Vec<String>,
    obs: List,
    var_index: Vec<String>,
    var: List,
    obsm: List,
    varm: List,
    obsp: List,
    uns_json: &str,
    chunk_size: usize,
) -> Result<()> {
    let cell_indices = cell_indices.r_int_convert();
    let n_obs = cell_indices.len();
    let n_vars = var_index.len();
    let chunk_size = chunk_size.max(1);

    // convert and check everything before the file exists, so that bad input
    // fails here and not after the counts have been streamed
    let reader = ParallelSparseReader::new(f_path_cells).to_extendr()?;
    check_len("obs_index", obs_index.len(), n_obs)?;
    let obs = ObsTable {
        index: obs_index,
        columns: r_list_to_columns(obs, n_obs)?,
    };
    let var = VarTable {
        index: var_index,
        columns: r_list_to_columns(var, n_vars)?,
    };
    let obsm = Embeddings {
        map: r_list_to_dense(obsm, n_obs)?,
    };
    let varm = Varm {
        map: r_list_to_dense(varm, n_vars)?,
    };
    let obsp = obsp
        .iter()
        .map(|(name, x)| Ok((name, r_list_to_csr(name, x, n_obs)?)))
        .collect::<Result<Vec<_>>>()?;
    let uns = UnsTable {
        raw: serde_json::from_str(uns_json).to_extendr()?,
    };

    let opts = H5AdOptions {
        // ponytail: fixed gzip level, the one the rhdf5 writers used; expose
        // it if anyone needs faster or smaller files
        compression: Some(6),
        chunk_size: Some(chunk_size),
    };
    let mut builder = H5AdBuilder::new(Path::new(h5_path), n_obs, n_vars, &opts).to_extendr()?;

    for (i, batch) in cell_indices.chunks(chunk_size).enumerate() {
        let cells = reader.read_cells_parallel(batch).to_extendr()?;
        let nnz = cells.iter().map(|c| c.indices.len()).sum();
        let mut indptr: Vec<u64> = Vec::with_capacity(cells.len() + 1);
        let mut indices: Vec<u32> = Vec::with_capacity(nnz);
        let mut data: Vec<f32> = Vec::with_capacity(nnz);
        indptr.push(0);
        for cell in &cells {
            indices.extend_from_slice(&cell.indices);
            indptr.push(indices.len() as u64);
        }
        if norm {
            data.extend(cells.iter().flat_map(|c| norm_counts_to_f32(&c.data_norm)));
        } else {
            data.extend(
                cells
                    .iter()
                    .flat_map(|c| c.data_raw.iter().map(|x| x as f32)),
            );
        }
        // the vectors move into scx, no copy
        builder
            .push_x_csr_chunk(i * chunk_size, indptr, indices, data)
            .to_extendr()?;
    }

    builder.obs(obs).to_extendr()?;
    builder.var(var).to_extendr()?;
    builder.add_obsm(obsm).to_extendr()?;
    builder.add_varm(varm).to_extendr()?;
    for (name, mat) in &obsp {
        builder.add_obsp_csr(name, mat.view()).to_extendr()?;
    }
    builder.add_uns(uns).to_extendr()?;
    builder.finalize().to_extendr()?;

    Ok(())
}
