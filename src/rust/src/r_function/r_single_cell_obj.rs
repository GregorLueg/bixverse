use extendr_api::prelude::*;

use crate::{helpers::structs_sparse::*, utils::general::flatten_vector};
use rayon::prelude::*;

// Extendr unfortunately cannot do Roxygen2 manipulation of R6 type
// classes. This will have to be done manually in R... Documentation
// still here to make it easier.

///////////
// Enums //
///////////

/// Enum for the count type to return
#[derive(Clone, Debug)]
pub enum AssayType {
    /// Return the `Raw` counts
    Raw,
    /// Return the log `Norm` counts
    Norm,
}

/// Enum for mixed types
#[derive(Clone, Debug)]
enum AssayData {
    /// The `Raw` data as `i32`.
    Raw(Vec<i32>),
    /// The `Norm` data as `f32`.
    Norm(Vec<f32>),
}

impl AssayData {
    /// Get the length of the vector
    ///
    /// ### Returns
    ///
    /// The length of the vector
    fn len(&self) -> usize {
        match self {
            AssayData::Raw(data) => data.len(),
            AssayData::Norm(data) => data.len(),
        }
    }

    /// Flatten the data into an R vector
    ///
    /// ### Params
    ///
    /// * `data` - A vector of AssayData that shall be flattened and transformed
    ///   into an R object
    ///
    /// ### Returns
    ///
    /// The flattened vector as an `Robj`
    fn flatten_into_r_vector(data: Vec<AssayData>) -> Robj {
        if data.is_empty() {
            return Robj::from(Vec::<f64>::new());
        }

        match &data[0] {
            AssayData::Raw(_) => {
                let flattened: Vec<i32> = data
                    .into_iter()
                    .flat_map(|d| match d {
                        AssayData::Raw(vec) => vec,
                        AssayData::Norm(_) => unreachable!(),
                    })
                    .collect();
                Robj::from(flattened)
            }
            AssayData::Norm(_) => {
                let flattened: Vec<f64> = data
                    .into_iter()
                    .flat_map(|d| match d {
                        AssayData::Norm(vec) => {
                            vec.into_iter().map(|x| x as f64).collect::<Vec<_>>()
                        }
                        AssayData::Raw(_) => unreachable!(),
                    })
                    .collect();
                Robj::from(flattened)
            }
        }
    }
}

/////////////
// Helpers //
/////////////

/// Parsing the count types
///
/// ### Params
///
/// * `s` - String defining the count type to return
///
/// ### Returns
///
/// The `AssayType`.
pub fn parse_count_type(s: &str) -> Option<AssayType> {
    match s.to_lowercase().as_str() {
        "raw" => Some(AssayType::Raw),
        "norm" => Some(AssayType::Norm),
        _ => None,
    }
}

////////////////
// Structures //
////////////////

/// A class for handling single cell count data
///
/// ### Params
///
/// * `f_path_cells` - Path to the .bin file for the cells.
/// * `f_path_genes` - Path to the .bin file for the genes.
#[extendr]
struct SingeCellCountData {
    pub f_path_cells: String,
    pub f_path_genes: String,
}

#[extendr]
impl SingeCellCountData {
    /// Create new instance of the class
    ///
    /// ### Params
    ///
    /// * `f_path_cells` - Path to the .bin file for the cells.
    /// * `f_path_genes` - Path to the .bin file for the genes.
    pub fn new(f_path_cells: String, f_path_genes: String) -> Self {
        Self {
            f_path_cells,
            f_path_genes,
        }
    }

    /// Write data from R CSC to disk
    ///
    /// Helper function to write CSC matrices from R to disk.
    ///
    /// ### Params
    ///
    /// * `no_cells` - Number of cells, i.e., columns in the data (`csc_matrix@Dim[2]`).
    /// * `no_genes` - Number of genes, i.e., rows in the data (`csc_matrix@Dim[1]`).
    /// * `data` - Slice of the data (`csc_matrix@x`).
    /// * `col_ptr` - The column pointers of the data (`csc_matrix@p`).
    /// * `row_idx` - The row indices of the data (`csc_matrix@i`).
    /// * `target_size` - To which target size to normalise to. 1e6 ->
    ///   CPM normalisation
    pub fn r_csc_mat_to_file(
        &self,
        no_cells: usize,
        no_genes: usize,
        data: &[i32],
        col_ptr: &[i32],
        row_idx: &[i32],
        target_size: f64,
    ) {
        let ptr_range = 1..col_ptr.len();

        let mut writer =
            CellGeneSparseWriter::new(&self.f_path_cells, true, no_cells, no_genes).unwrap();

        for i in ptr_range {
            // get the index position
            let end_i = col_ptr[i] as usize;
            let start_i = col_ptr[i - 1] as usize;
            // create chunk from raw data
            let chunk_i = CscCellChunk::from_r_data(
                &data[start_i..end_i],
                &row_idx[start_i..end_i],
                i - 1,
                target_size as f32,
            );

            writer.write_cell_chunk(chunk_i).unwrap();
        }

        writer.finalise().unwrap();
    }

    /// Returns the full data as a List
    ///
    /// ### Params
    ///
    /// * `assay` - String. Return the raw counts or log-normalised counts. One
    /// of `"raw"` or `"norm"`.
    ///
    /// ### Returns
    ///
    /// An R list with all of the info that was stored in the .bin file
    pub fn file_to_r_csc_mat(&self, assay: &str) -> List {
        let mut reader = StreamingSparseReader::new(&self.f_path_cells).unwrap();

        let assay_type = parse_count_type(assay).unwrap();

        let no_cells = reader.get_header().total_cells;
        let no_genes = reader.get_header().total_genes;

        let mut data: Vec<AssayData> = Vec::new();
        let mut row_idx: Vec<Vec<u16>> = Vec::new();
        let mut col_ptr: Vec<usize> = Vec::new();

        let all_cells: Vec<_> = reader.iter_cells().map(|r| r.unwrap()).collect();

        let mut current_pol_ptr = 0_usize;
        col_ptr.push(current_pol_ptr);

        for cell in all_cells {
            let data_i = match assay_type {
                AssayType::Raw => AssayData::Raw(cell.data_raw.iter().map(|&x| x as i32).collect()),
                AssayType::Norm => AssayData::Norm(
                    cell.data_norm
                        .iter()
                        .map(|x| {
                            let f16_val: half::f16 = (*x).into();
                            f16_val.to_f32()
                        })
                        .collect(),
                ),
            };
            let len_data_i = data_i.len();
            current_pol_ptr += len_data_i;
            // Add data
            data.push(data_i);
            row_idx.push(cell.row_indices);
            col_ptr.push(current_pol_ptr);
        }

        let data = AssayData::flatten_into_r_vector(data);
        let row_idx = flatten_vector(row_idx);

        let row_idx_final = row_idx.par_iter().map(|x| *x as i32).collect::<Vec<i32>>();

        list!(
            col_ptr = col_ptr,
            row_idx = row_idx_final,
            data = data,
            no_cells = no_cells,
            no_genes = no_genes
        )
    }

    /// Returns the data for the specified indices
    ///
    /// ### Params
    ///
    /// * `indices` - Slice of the index positions for the cells
    /// * `assay` - String. Return the raw counts or log-normalised counts. One
    /// of `"raw"` or `"norm"`
    ///
    /// ### Returns
    ///
    /// An R list with the selected cell data that was stored in the .bin file
    pub fn get_cells_by_indices(&self, indices: &[i32], assay: &str) -> List {
        // TODO better error handling here
        let mut reader = StreamingSparseReader::new(&self.f_path_cells).unwrap();

        let assay_type = parse_count_type(assay).unwrap();

        let no_genes = reader.get_header().total_genes;

        let indices = indices
            .iter()
            .map(|x| (*x - 1_i32) as usize)
            .collect::<Vec<usize>>();
        let cells = reader.read_cells_by_indices(&indices).unwrap();

        let mut data: Vec<AssayData> = Vec::new();
        let mut row_idx: Vec<Vec<u16>> = Vec::new();
        let mut col_ptr: Vec<usize> = Vec::new();

        let mut current_pol_ptr = 0_usize;
        col_ptr.push(current_pol_ptr);

        for cell in &cells {
            let data_i = match assay_type {
                AssayType::Raw => AssayData::Raw(cell.data_raw.iter().map(|&x| x as i32).collect()),
                AssayType::Norm => AssayData::Norm(
                    cell.data_norm
                        .iter()
                        .map(|x| {
                            let f16_val: half::f16 = (*x).into();
                            f16_val.to_f32()
                        })
                        .collect(),
                ),
            };
            let len_data_i = data_i.len();
            current_pol_ptr += len_data_i;
            // Add data
            data.push(data_i);
            row_idx.push(cell.row_indices.clone());
            col_ptr.push(current_pol_ptr);
        }

        let data = AssayData::flatten_into_r_vector(data);
        let row_idx = flatten_vector(row_idx);

        let row_idx_final = row_idx.par_iter().map(|x| *x as i32).collect::<Vec<i32>>();

        list!(
            col_ptr = col_ptr,
            row_idx = row_idx_final,
            data = data,
            no_cells = cells.len(),
            no_genes = no_genes
        )
    }
}

extendr_module! {
    mod r_single_cell_obj;
    impl SingeCellCountData;
}
