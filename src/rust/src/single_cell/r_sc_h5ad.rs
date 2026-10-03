use extendr_api::*;
use std::collections::HashMap;
use std::path::Path;

use bixverse_rs::prelude::*;
use scx_core::api::write::{H5AdBuilder, H5AdOptions};
use scx_core::ir::{
    Column, ColumnData, DenseMatrix, Embeddings, ObsTable, UnsTable, VarTable, Varm,
};
use sprs::CsMatI;

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

/// Translate a named R list of columns into scx columns
///
/// ### Params
///
/// * `cols` - Named list of factor, double, integer, logical or character
///   vectors. The R side decides the encoding; this only maps the types.
///
/// ### Returns
///
/// The columns in list order.
fn r_list_to_columns(cols: List) -> Result<Vec<Column>> {
    cols.iter()
        .map(|(name, x)| {
            let data = if x.inherits("factor") {
                let codes = x
                    .as_integer_slice()
                    .ok_or_else(|| Error::Other(format!("Factor '{name}' has no codes.")))?
                    .iter()
                    // `NA_integer_` becomes u32::MAX, which scx narrows to -1,
                    // the anndata code for a missing value, in every code width
                    .map(|&c| {
                        if c == i32::MIN {
                            u32::MAX
                        } else {
                            (c - 1) as u32
                        }
                    })
                    .collect();
                let levels = x
                    .get_attrib("levels")
                    .and_then(|l| l.as_string_vector())
                    .unwrap_or_default();
                ColumnData::Categorical { codes, levels }
            } else if let Some(v) = x.as_real_slice() {
                ColumnData::Float(v.to_vec())
            } else if let Some(v) = x.as_integer_slice() {
                ColumnData::Int(v.to_vec())
            } else if let Some(v) = x.as_logical_slice() {
                ColumnData::Bool(v.iter().map(|b| b.is_true()).collect())
            } else if let Some(v) = x.as_string_vector() {
                ColumnData::String(v)
            } else {
                return Err(Error::Other(format!(
                    "Column '{name}' has a type with no h5ad encoding."
                )));
            };
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
///
/// ### Returns
///
/// A map from name to the row-major matrix.
fn r_list_to_dense(mats: List) -> Result<HashMap<String, DenseMatrix>> {
    mats.iter()
        .map(|(name, x)| {
            let mat = RMatrix::<f64>::try_from(x)?;
            let (nrow, ncol) = (mat.nrows(), mat.ncols());
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
/// R. Counts are stored as `float32`.
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

    let opts = H5AdOptions {
        // ponytail: fixed gzip level, the one the rhdf5 writers used; expose
        // it if anyone needs faster or smaller files
        compression: Some(6),
        chunk_size: Some(chunk_size),
    };
    let mut builder = H5AdBuilder::new(Path::new(h5_path), n_obs, n_vars, &opts).to_extendr()?;

    let reader = ParallelSparseReader::new(f_path_cells).to_extendr()?;
    for (i, batch) in cell_indices.chunks(chunk_size).enumerate() {
        let cells = reader.read_cells_parallel(batch).to_extendr()?;
        let mut indptr: Vec<u64> = Vec::with_capacity(batch.len() + 1);
        let mut indices: Vec<u32> = Vec::new();
        let mut data: Vec<f32> = Vec::new();
        indptr.push(0);
        for cell in &cells {
            indices.extend_from_slice(&cell.indices);
            if norm {
                data.extend(cell.data_norm.iter().map(|&x| {
                    let v: half::f16 = x.into();
                    v.to_f32()
                }));
            } else {
                data.extend(cell.data_raw.iter().map(|x| x as f32));
            }
            indptr.push(indices.len() as u64);
        }
        builder
            .push_x_csr_chunk(i * chunk_size, &indptr, &indices, &data)
            .to_extendr()?;
    }

    builder
        .obs(ObsTable {
            index: obs_index,
            columns: r_list_to_columns(obs)?,
        })
        .to_extendr()?;
    builder
        .var(VarTable {
            index: var_index,
            columns: r_list_to_columns(var)?,
        })
        .to_extendr()?;
    builder
        .add_obsm(Embeddings {
            map: r_list_to_dense(obsm)?,
        })
        .to_extendr()?;
    builder
        .add_varm(Varm {
            map: r_list_to_dense(varm)?,
        })
        .to_extendr()?;

    for (name, x) in obsp.iter() {
        let csr = List::try_from(x)?;
        let get_u32 = |key: &str| -> Result<Vec<u32>> {
            let v = csr.dollar(key)?;
            let v = v
                .as_integer_slice()
                .ok_or_else(|| Error::Other(format!("obsp '{name}': {key} is not integer.")))?;
            Ok(v.iter().map(|&i| i as u32).collect())
        };
        let data = csr.dollar("data")?;
        let data: Vec<f32> = data
            .as_real_slice()
            .ok_or_else(|| Error::Other(format!("obsp '{name}': data is not double.")))?
            .iter()
            .map(|&v| v as f32)
            .collect();
        let mat = CsMatI::<f32, u32>::try_new(
            (n_obs, n_obs),
            get_u32("indptr")?,
            get_u32("indices")?,
            data,
        )
        .map_err(|(_, _, _, e)| Error::Other(format!("obsp '{name}': {e}")))?;
        builder.add_obsp_csr(name, mat.view()).to_extendr()?;
    }

    let uns: serde_json::Value = serde_json::from_str(uns_json).to_extendr()?;
    builder.add_uns(UnsTable { raw: uns }).to_extendr()?;

    builder.finalize().to_extendr()?;

    Ok(())
}
