use extendr_api::prelude::*;
use rayon::prelude::*;

use crate::core::data::sparse_io::*;
use crate::core::data::sparse_io_h5::*;
use crate::core::data::sparse_io_mtx::*;
use crate::core::data::sparse_structures::*;
use crate::single_cell::processing::*;
use crate::utils::general::flatten_vector;
use crate::utils::traits::F16;

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
    pub n_cells: usize,
    pub n_genes: usize,
    pub gene_mask: Vec<u16>,
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
            n_cells: usize::default(),
            n_genes: usize::default(),
            gene_mask: Vec::new(),
        }
    }

    /////////////
    // Helpers //
    /////////////

    /// Get the shape
    pub fn get_shape(&mut self) -> Vec<usize> {
        vec![self.n_cells, self.n_genes]
    }

    ///////////
    // Cells //
    ///////////

    /// Write data from R CSR to disk
    ///
    /// Helper function to write CSR matrices from R to disk.
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
    /// * `min_genes` - Number of minimum genes needed to be considered included
    ///   in analysis.
    ///
    /// ### Returns
    ///
    /// The mask of cells that have the number of genes expected.
    #[allow(clippy::too_many_arguments)]
    pub fn r_csr_mat_to_file(
        &mut self,
        no_cells: usize,
        no_genes: usize,
        data: &[i32],
        row_ptr: &[i32],
        col_idx: &[i32],
        target_size: f64,
        min_genes: usize,
    ) -> List {
        self.n_cells = no_cells;
        self.n_genes = no_genes;
        let row_ptr_usize = row_ptr.iter().map(|x| *x as usize).collect::<Vec<usize>>();

        let (nnz, cell_mask) = filter_by_nnz(&row_ptr_usize, min_genes);

        let mut writer =
            CellGeneSparseWriter::new(&self.f_path_cells, true, no_cells, no_genes).unwrap();

        for i in 0..self.n_cells {
            // get the index position
            let end_i = row_ptr_usize[i + 1];
            let start_i = row_ptr_usize[i];
            let to_keep_i = cell_mask[i];

            // create chunk from raw data
            let chunk_i = CsrCellChunk::from_data(
                &data[start_i..end_i],
                &col_idx[start_i..end_i],
                i,
                target_size as f32,
                to_keep_i,
            );

            writer.write_cell_chunk(chunk_i).unwrap();
        }

        writer.finalise().unwrap();

        list!(nnz = nnz, cell_mask = cell_mask)
    }

    /// Save h5 to file
    ///
    /// ### Params
    ///
    /// * `cs_type` - How was the h5 data saved. CSC or CSR.
    /// * `h5_path` - Path to the h5 file.
    /// * `no_cells` - Number of cells in the h5 file.
    /// * `no_genes` - Number of genes in the h5 file.
    /// * `qc_params` - List with the quality control parameters.
    pub fn h5_to_file(
        &mut self,
        cs_type: String,
        h5_path: String,
        no_cells: usize,
        no_genes: usize,
        qc_params: List,
    ) -> List {
        // assign the values internally
        self.n_cells = no_cells;
        self.n_genes = no_genes;

        let qc_params = MinCellQuality::from_r_list(qc_params);

        let cell_qc: CellQuality = write_h5_counts(
            &h5_path,
            &self.f_path_cells,
            &cs_type,
            no_cells,
            no_genes,
            qc_params,
        );

        list!(
            cell_mask = cell_qc.to_keep,
            lib_size = cell_qc.lib_size.unwrap_or_default(),
            nnz = cell_qc.no_genes.unwrap_or_default()
        )
    }

    pub fn mtx_to_file(&mut self, mtx_path: String, qc_params: List) -> extendr_api::Result<List> {
        let qc_params = MinCellQuality::from_r_list(qc_params);

        let mtx_reader = MtxReader::new(&mtx_path, qc_params)
            .map_err(|e| extendr_api::Error::from(Box::new(e) as Box<dyn std::error::Error>))?;

        let mtx_res = mtx_reader
            .process_mtx_and_write_bin(&self.f_path_cells)
            .map_err(|e| extendr_api::Error::from(Box::new(e) as Box<dyn std::error::Error>))?;

        self.n_cells = mtx_res.no_cells;
        self.n_genes = mtx_res.no_genes;

        Ok(list!(cells_kept = mtx_res.to_keep))
    }

    /// Returns the full matrix
    ///
    /// ### Params
    ///
    /// * `assay` - String. Return the raw counts or log-normalised counts. One
    ///   of `"raw"` or `"norm"`.
    /// * `cell_based` - Boolean. Shall the data be returned in CSR or CSC.
    /// * `verbose` - Boolean. Verbosity of the function.
    ///
    /// ### Returns
    ///
    /// An R list with all of the info that was stored in the .bin file
    pub fn return_full_mat(&self, assay: &str, cell_based: bool, verbose: bool) -> List {
        let mut data: Vec<AssayData> = Vec::new();
        let mut indices: Vec<Vec<i32>> = Vec::new();
        let mut indptr: Vec<usize> = Vec::new();
        let assay_type = parse_count_type(assay).unwrap();

        if cell_based {
            let reader = ParallelSparseReader::new(&self.f_path_cells).unwrap();
            let cell_chunks = reader.get_all_cells();

            if verbose {
                println!("All cells loaded in successfully.")
            }

            let mut current_ptr = 0_usize;
            indptr.push(current_ptr);

            for cell in cell_chunks {
                let data_i = match assay_type {
                    AssayType::Raw => {
                        AssayData::Raw(cell.data_raw.iter().map(|&x| x as i32).collect())
                    }
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
                current_ptr += len_data_i;
                // Add data
                data.push(data_i);
                indices.push(cell.col_indices.iter().map(|x| *x as i32).collect());
                indptr.push(current_ptr);
            }
        } else {
            let reader = ParallelSparseReader::new(&self.f_path_genes).unwrap();
            let gene_chunks = reader.get_all_genes();

            if verbose {
                println!("All cells loaded in successfully.")
            }

            let mut current_ptr = 0_usize;
            indptr.push(current_ptr);

            for gene in gene_chunks {
                let data_i = match assay_type {
                    AssayType::Raw => {
                        AssayData::Raw(gene.data_raw.iter().map(|&x| x as i32).collect())
                    }
                    AssayType::Norm => AssayData::Norm(
                        gene.data_norm
                            .iter()
                            .map(|x| {
                                let f16_val: half::f16 = (*x).into();
                                f16_val.to_f32()
                            })
                            .collect(),
                    ),
                };
                let len_data_i = data_i.len();
                current_ptr += len_data_i;
                data.push(data_i);
                indices.push(gene.row_indices.iter().map(|x| *x as i32).collect());
                indptr.push(current_ptr);
            }
        };

        let data: Robj = AssayData::flatten_into_r_vector(data);
        let indices = flatten_vector(indices);

        list!(
            indptr = indptr,
            indices = indices,
            data = data,
            no_cells = self.n_cells,
            no_genes = self.n_genes
        )
    }

    /// To write
    pub fn get_cells_by_indices(&self, indices: &[i32], assay: &str) -> List {
        let reader = ParallelSparseReader::new(&self.f_path_cells).unwrap();
        let assay_type = parse_count_type(assay).unwrap();

        let indices: Vec<usize> = indices.iter().map(|x| (*x - 1) as usize).collect();

        // Parallel read all cells at once
        let cells = reader.read_cells_parallel(&indices);

        // Parallel processing of results
        let results: Vec<(AssayData, Vec<u16>)> = cells
            .par_iter()
            .map(|cell| {
                let data_i = match assay_type {
                    AssayType::Raw => {
                        AssayData::Raw(cell.data_raw.iter().map(|&x| x as i32).collect())
                    }
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

                // let (data_i, indices_i) = if !self.gene_mask.is_empty() {

                // } else {

                // }

                (data_i, cell.col_indices.clone())
            })
            .collect();

        // Sequential assembly (required for row_ptr)
        let mut data = Vec::new();
        let mut col_idx = Vec::new();
        let mut row_ptr = vec![0];

        for (data_i, indices) in results {
            let len = data_i.len();
            row_ptr.push(row_ptr.last().unwrap() + len);
            data.push(data_i);
            col_idx.push(indices);
        }

        let data = AssayData::flatten_into_r_vector(data);
        let col_idx = flatten_vector(col_idx)
            .par_iter()
            .map(|x| *x as i32)
            .collect::<Vec<i32>>();

        list!(
            indptr = row_ptr,
            indices = col_idx,
            data = data,
            no_cells = cells.len(),
            no_genes = self.n_genes
        )
    }

    ///////////
    // Genes //
    ///////////

    /// Transforms already written cell data also into the gene data
    ///
    /// This function will read the .bin file at `self.f_path_cells` and
    /// transform the data into the gene-based file format. This happens for
    /// the full data set in memory and might cause memory pressure, pending
    /// the size of the data.
    ///
    /// ### Params
    ///
    /// ### Returns
    ///
    ///
    pub fn generate_gene_based_data(&mut self, qc_params: List) -> List {
        let reader = ParallelSparseReader::new(&self.f_path_cells).unwrap();

        let qc_params = MinCellQuality::from_r_list(qc_params);

        let no_cells = reader.get_header().total_cells;
        let no_genes = reader.get_header().total_genes;

        // Extract all data
        let all_cells: Vec<_> = reader.get_all_cells();

        let mut data: Vec<Vec<u16>> = Vec::new();
        let mut data_2: Vec<Vec<F16>> = Vec::new();
        let mut col_idx: Vec<Vec<u16>> = Vec::new();
        let mut row_ptr: Vec<usize> = Vec::new();

        let mut current_row_ptr = 0_usize;
        row_ptr.push(current_row_ptr);

        for cell in all_cells {
            let data_i = cell.data_raw;
            let data_norm_i = cell.data_norm;

            let len_data_i = data_i.len();
            current_row_ptr += len_data_i;
            // Add data
            data.push(data_i);
            data_2.push(data_norm_i);
            col_idx.push(cell.col_indices);
            row_ptr.push(current_row_ptr);
        }

        let data = flatten_vector(data);
        let data_2 = flatten_vector(data_2);
        let col_idx = flatten_vector(col_idx);
        let col_idx = col_idx.iter().map(|x| *x as usize).collect::<Vec<usize>>();

        let sparse_data = CompressedSparseData::new_csr(
            &data,
            &col_idx,
            &row_ptr,
            Some(&data_2),
            (no_cells, no_genes),
        );

        let sparse_data = sparse_data.transform();

        let (nnz, gene_mask) = filter_by_nnz(&sparse_data.indptr, qc_params.min_cells);

        // Write the data in CSC format to disk
        let mut writer =
            CellGeneSparseWriter::new(&self.f_path_genes, false, no_cells, no_genes).unwrap();

        let data_2 = sparse_data.get_data2_unsafe();

        let mut number_genes_written = 0_usize;

        #[allow(clippy::needless_range_loop)]
        for i in 0..no_genes {
            if gene_mask[i] {
                // get the index position
                let start_i = sparse_data.indptr[i];
                let end_i = sparse_data.indptr[i + 1];
                let to_keep_i = gene_mask[i];
                // create the chunk and write to disk
                let chunk_i = CscGeneChunk::from_conversion(
                    &sparse_data.data[start_i..end_i],
                    &data_2[start_i..end_i],
                    &sparse_data.indices[start_i..end_i],
                    i,
                    to_keep_i,
                );

                writer.write_gene_chunk(chunk_i).unwrap();
                number_genes_written += 1;
            }
        }

        writer.update_header_no_genes(number_genes_written);
        writer.finalise().unwrap();

        let gene_indices: Vec<u16> = gene_mask
            .iter()
            .enumerate()
            .filter_map(|(i, &val)| if val { Some(i as u16) } else { None })
            .collect();

        self.gene_mask = gene_indices;
        self.n_genes = number_genes_written;

        list!(nnz = nnz, gene_mask = gene_mask)
    }

    pub fn get_genes_by_indices(&self, indices: &[i32], assay: &str) -> List {
        let reader = ParallelSparseReader::new(&self.f_path_genes).unwrap();

        let assay_type = parse_count_type(assay).unwrap();

        let no_cells = reader.get_header().total_cells;

        let indices = indices
            .iter()
            .map(|x| (*x - 1) as usize)
            .collect::<Vec<usize>>();

        let genes = reader.read_gene_parallel(&indices);

        let mut data: Vec<AssayData> = Vec::new();
        let mut col_idx: Vec<Vec<u32>> = Vec::new();
        let mut row_ptr: Vec<usize> = Vec::new();

        let mut current_col_ptr = 0_usize;
        row_ptr.push(current_col_ptr);

        for gene in &genes {
            let data_i = match assay_type {
                AssayType::Raw => AssayData::Raw(gene.data_raw.iter().map(|&x| x as i32).collect()),
                AssayType::Norm => AssayData::Norm(
                    gene.data_norm
                        .iter()
                        .map(|x| {
                            let f16_val: half::f16 = (*x).into();
                            f16_val.to_f32()
                        })
                        .collect(),
                ),
            };
            let len_data_i = data_i.len();
            current_col_ptr += len_data_i;
            data.push(data_i);
            col_idx.push(gene.row_indices.clone());
            row_ptr.push(current_col_ptr);
        }

        let data = AssayData::flatten_into_r_vector(data);
        let col_idx = flatten_vector(col_idx)
            .par_iter()
            .map(|x| *x as i32)
            .collect::<Vec<i32>>();

        list!(
            indptr = row_ptr,
            indices = col_idx,
            data = data,
            no_cells = no_cells,
            no_genes = indices.len()
        )
    }
}

extendr_module! {
    mod r_count_obj;
    impl SingeCellCountData;
}
