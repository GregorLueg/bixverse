use extendr_api::prelude::*;
use rayon::prelude::*;
use rustc_hash::{FxHashMap, FxHashSet};
use std::sync::Arc;
use std::time::Instant;
use thousands::Separable;

use crate::core::data::sparse_io::*;
use crate::core::data::sparse_io_h5::*;
use crate::core::data::sparse_io_mtx::*;
use crate::core::data::sparse_io_r_obj::*;
use crate::core::data::sparse_structures::*;
use crate::single_cell::processing::*;
use crate::utils::general::flatten_vector;
use crate::utils::r_rust_interface::list_to_sparse_matrix;
use crate::utils::traits::F16;

// Extendr unfortunately cannot do Roxygen2 manipulation of R6 type
// classes. This will have to be done manually in R... Documentation
// still here to make it easier.

//////////////////
// Type aliases //
//////////////////

/// Type to store the GeneData during streaming
type GeneData = Vec<FxHashMap<u16, Vec<(u32, u16, F16)>>>;

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

/// Helper function to retrieve and optionally filter cell data
///
/// ### Params
///
/// * `indices` - The original indices (in this case column).
/// * `data_raw`- The raw counts for the cell.
/// * `data_norm` - The normalised counts for the cell.
/// * `assay_type` - Which assay to return
///
/// ### Returns
///
/// Tuple of the indices and respective data
fn get_cell_data(
    indices: &[u16],
    data_raw: &[u16],
    data_norm: &[F16],
    assay_type: &AssayType,
) -> (Vec<i32>, AssayData) {
    let all_indices: Vec<i32> = indices.iter().map(|&x| x as i32).collect();

    let data = match assay_type {
        AssayType::Raw => AssayData::Raw(data_raw.iter().map(|&x| x as i32).collect()),
        AssayType::Norm => {
            let norm_data: Vec<f32> = data_norm
                .iter()
                .map(|&x| {
                    let f16_val: half::f16 = x.into();
                    f16_val.to_f32()
                })
                .collect();
            AssayData::Norm(norm_data)
        }
    };

    (all_indices, data)
}

/// Helper function to retrieve the gene data
///
/// ### Params
///
/// * `indices` - The original indices (in this case column).
/// * `data_raw`- The raw counts for the cell.
/// * `data_norm` - The normalised counts for the cell.
/// * `assay_type` - Which assay to return
///
/// ### Returns
///
/// Tuple of the indices and respective data
fn get_gene_data(
    indices: &[u32],
    data_raw: &[u16],
    data_norm: &[F16],
    assay_type: &AssayType,
) -> (Vec<i32>, AssayData) {
    let all_indices: Vec<i32> = indices.iter().map(|&x| x as i32).collect();

    let data = match assay_type {
        AssayType::Raw => AssayData::Raw(data_raw.iter().map(|&x| x as i32).collect()),
        AssayType::Norm => {
            let norm_data: Vec<f32> = data_norm
                .iter()
                .map(|&x| {
                    let f16_val: half::f16 = x.into();
                    f16_val.to_f32()
                })
                .collect();
            AssayData::Norm(norm_data)
        }
    };

    (all_indices, data)
}

/// Parsing the count types
///
/// ### Params
///
/// * `s` - String defining the count type to return
///
/// ### Returns
///
/// The `AssayType`.
fn parse_count_type(s: &str) -> Option<AssayType> {
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
/// * `n_cells` - No of cells represented in the data.
/// * `n_genes` - No of genes represented in the data.
/// * `cell_mask` - HashSet for which cells to keep.
#[extendr]
struct SingeCellCountData {
    pub f_path_cells: String,
    pub f_path_genes: String,
    pub n_cells: usize,
    pub n_genes: usize,
    pub cell_mask: FxHashSet<u32>,
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
            cell_mask: FxHashSet::default(),
        }
    }

    /////////////
    // Helpers //
    /////////////

    /// Get the shape
    ///
    /// ### Returns
    ///
    /// Vector with rows x cells
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
    /// * `r_data` - A list that can be transformed into CompressedSparseData.
    ///   Needs to have following elements: `"indptr"`, `"indices"`, `"data"`,
    ///   `"nrow"`, `"ncol"` and `"format"`.
    /// * `qc_params` - List with the quality control parameters.
    /// * `verbose` - Controls verbosity of the function.
    ///
    /// ### Returns
    ///
    /// A list with QC parameters.
    pub fn r_data_to_file(&mut self, r_data: List, qc_params: List, verbose: bool) -> List {
        let qc_params = MinCellQuality::from_r_list(qc_params);

        let compressed_data: CompressedSparseData<u32> = list_to_sparse_matrix(r_data);

        let (no_cells, no_genes, cell_qc): (usize, usize, CellQuality) =
            write_r_counts(&self.f_path_cells, compressed_data, qc_params, verbose);

        self.n_cells = no_cells;
        self.n_genes = no_genes;

        list!(
            cell_indices = cell_qc.cell_indices,
            gene_indices = cell_qc.gene_indices,
            lib_size = cell_qc.lib_size,
            nnz = cell_qc.no_genes
        )
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
    /// * `verbose` - Controls verbosity of the function.
    ///
    /// ### Returns
    ///
    /// A list with qc parameters.
    pub fn h5_to_file(
        &mut self,
        cs_type: String,
        h5_path: String,
        no_cells: usize,
        no_genes: usize,
        qc_params: List,
        verbose: bool,
    ) -> List {
        let qc_params = MinCellQuality::from_r_list(qc_params);

        let (no_cells, no_genes, cell_qc): (usize, usize, CellQuality) = write_h5_counts(
            &h5_path,
            &self.f_path_cells,
            &cs_type,
            no_cells,
            no_genes,
            qc_params,
            verbose,
        );

        self.n_cells = no_cells;
        self.n_genes = no_genes;

        list!(
            cell_indices = cell_qc.cell_indices,
            gene_indices = cell_qc.gene_indices,
            lib_size = cell_qc.lib_size,
            nnz = cell_qc.no_genes
        )
    }

    /// Save h5 to file
    ///
    /// Slower version that is less memory heavy and will make usage of
    /// streaming where possible.
    ///
    /// ### Params
    ///
    /// * `cs_type` - How was the h5 data saved. CSC or CSR.
    /// * `h5_path` - Path to the h5 file.
    /// * `no_cells` - Number of cells in the h5 file.
    /// * `no_genes` - Number of genes in the h5 file.
    /// * `qc_params` - List with the quality control parameters.
    /// * `verbose` - Controls verbosity of the function.
    ///
    /// ### Returns
    ///
    /// A list with qc parameters.
    pub fn h5_to_file_streaming(
        &mut self,
        cs_type: String,
        h5_path: String,
        no_cells: usize,
        no_genes: usize,
        qc_params: List,
        verbose: bool,
    ) -> List {
        let qc_params = MinCellQuality::from_r_list(qc_params);

        let (no_cells, no_genes, cell_qc): (usize, usize, CellQuality) = stream_h5_counts(
            &h5_path,
            &self.f_path_cells,
            &cs_type,
            no_cells,
            no_genes,
            qc_params,
            verbose,
        );

        self.n_cells = no_cells;
        self.n_genes = no_genes;

        list!(
            cell_indices = cell_qc.cell_indices,
            gene_indices = cell_qc.gene_indices,
            lib_size = cell_qc.lib_size,
            nnz = cell_qc.no_genes
        )
    }

    /// Save mtx to file
    ///
    /// ### Params
    ///
    /// * `mtx_path` - Path to the mtx file.
    /// * `qc_params` - List with the quality control parameters.
    /// * `cells_as_rows` - Do the cells represent rows (= TRUE) or columns.
    /// * `verbose` - Controls verbosity of the function.
    ///
    /// ### Returns
    ///
    /// A list with qc parameters.
    pub fn mtx_to_file(
        &mut self,
        mtx_path: String,
        qc_params: List,
        cells_as_rows: bool,
        verbose: bool,
    ) -> extendr_api::Result<List> {
        let qc_params = MinCellQuality::from_r_list(qc_params);

        let mut mtx_reader = MtxReader::new(&mtx_path, qc_params, cells_as_rows)
            .map_err(|e| extendr_api::Error::from(Box::new(e) as Box<dyn std::error::Error>))?;

        let mtx_quality_data = mtx_reader
            .parse_mtx_quality(verbose)
            .map_err(|e| extendr_api::Error::from(Box::new(e) as Box<dyn std::error::Error>))?;

        let mtx_res: MtxFinalData = mtx_reader
            .process_mtx_and_write_bin(&self.f_path_cells, &mtx_quality_data, verbose)
            .map_err(|e| extendr_api::Error::from(Box::new(e) as Box<dyn std::error::Error>))?;

        self.n_cells = mtx_res.no_cells;
        self.n_genes = mtx_res.no_genes;

        Ok(list!(
            cell_indices = mtx_res.cell_qc.cell_indices,
            gene_indices = mtx_res.cell_qc.gene_indices,
            lib_size = mtx_res.cell_qc.lib_size,
            nnz = mtx_res.cell_qc.no_genes
        ))
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
                let (indices_i, data_i) =
                    get_cell_data(&cell.indices, &cell.data_raw, &cell.data_norm, &assay_type);

                let len_data_i = data_i.len();
                current_ptr += len_data_i;
                // Add data
                data.push(data_i);
                indices.push(indices_i);
                indptr.push(current_ptr);
            }
        } else {
            let reader = ParallelSparseReader::new(&self.f_path_genes).unwrap();
            let gene_chunks = reader.get_all_genes();

            if verbose {
                println!("All genes loaded in successfully.")
            }

            let mut current_ptr = 0_usize;
            indptr.push(current_ptr);

            for gene in gene_chunks {
                let (indices_i, data_i) =
                    get_gene_data(&gene.indices, &gene.data_raw, &gene.data_norm, &assay_type);
                let len_data_i = data_i.len();
                current_ptr += len_data_i;
                data.push(data_i);
                indices.push(indices_i);
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

    /// Return cells by index positions
    ///
    /// Leverages the CSR-stored data for fast cell retrieval
    ///
    /// ### Params
    ///
    /// * `indices` - The cell indices which to return
    /// * `assay` - Shall the raw or norm counts be returned
    ///
    /// ### Returns
    ///
    /// A list that can be parsed into a CSR matrix in R
    pub fn get_cells_by_indices(&self, indices: &[i32], assay: &str) -> List {
        let reader = ParallelSparseReader::new(&self.f_path_cells).unwrap();
        let assay_type = parse_count_type(assay).unwrap();

        let indices: Vec<usize> = indices.iter().map(|x| (*x - 1) as usize).collect();

        // Parallel read all cells at once
        let cells = reader.read_cells_parallel(&indices);

        // Parallel processing of results
        let results: Vec<(Vec<i32>, AssayData)> = cells
            .par_iter()
            .map(|cell| get_cell_data(&cell.indices, &cell.data_raw, &cell.data_norm, &assay_type))
            .collect();

        // Sequential assembly (required for row_ptr)
        let mut data = Vec::new();
        let mut col_idx = Vec::new();
        let mut row_ptr = vec![0];

        for (indices, data_i) in results {
            let len = data_i.len();
            row_ptr.push(row_ptr.last().unwrap() + len);
            data.push(data_i);
            col_idx.push(indices);
        }

        let data = AssayData::flatten_into_r_vector(data);
        let col_idx = flatten_vector(col_idx);

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
    /// * `qc_params` - List with the quality control parameters.
    /// * `verbose` - Controls verbosity of the function.
    pub fn generate_gene_based_data(&mut self, verbose: bool) {
        let reader = ParallelSparseReader::new(&self.f_path_cells).unwrap();

        let no_cells = reader.get_header().total_cells;
        let no_genes = reader.get_header().total_genes;

        let start_conversion = Instant::now();

        if verbose {
            println!("Loading in the cell data and saving it into a gene-friendly file")
        }

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
            col_idx.push(cell.indices);
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

        // get the data
        let sparse_data = sparse_data.transform();
        let data_2 = sparse_data.get_data2_unsafe();

        // write the data in CSC format to disk
        let mut writer =
            CellGeneSparseWriter::new(&self.f_path_genes, false, no_cells, no_genes).unwrap();

        #[allow(clippy::needless_range_loop)]
        for i in 0..no_genes {
            // get the index position
            let start_i = sparse_data.indptr[i];
            let end_i = sparse_data.indptr[i + 1];
            // create the chunk and write to disk
            let chunk_i = CscGeneChunk::from_conversion(
                &sparse_data.data[start_i..end_i],
                &data_2[start_i..end_i],
                &sparse_data.indices[start_i..end_i],
                i,
                true,
            );

            writer.write_gene_chunk(chunk_i).unwrap();
        }

        let end_conversion = start_conversion.elapsed();

        if verbose {
            println!(
                "Convertion data into gene-friendly format done: {:.2?}",
                end_conversion
            );
        }

        writer.finalise().unwrap();
    }

    /// Generate gene-based data with streaming to reduce memory pressure
    ///
    /// This approach builds CSC format directly without creating intermediate
    /// CSR structures. Ideal for very large data sets.
    ///
    /// ### Params
    ///
    /// * `batch_size` - Size of the batch to process in one go. The larger, the
    ///   more memory pressure will occur.
    /// * `verbose` - Controls verbosity of the function.
    pub fn generate_gene_based_data_streaming(&mut self, batch_size: usize, verbose: bool) {
        let reader = ParallelSparseReader::new(&self.f_path_cells).unwrap();
        let header = reader.get_header();
        let no_cells = header.total_cells;
        let no_genes = header.total_genes;

        if verbose {
            println!("Streaming cell data to gene-friendly format with reduced memory usage.");
        }

        let start_conversion = Instant::now();

        // Use HashMap to accumulate gene data - only keep current batch in memory
        let mut gene_data_map: FxHashMap<u16, Vec<(u32, u16, F16)>> = FxHashMap::default();

        // Process cells in batches to control memory usage
        let total_batches = no_cells.div_ceil(batch_size);

        let mut writer =
            CellGeneSparseWriter::new(&self.f_path_genes, false, no_cells, no_genes).unwrap();

        for batch_idx in 0..total_batches {
            let start_cell = batch_idx * batch_size;
            let end_cell = std::cmp::min(start_cell + batch_size, no_cells);
            let cell_indices: Vec<usize> = (start_cell..end_cell).collect();

            if verbose {
                let progress = (batch_idx + 1) as f32 / total_batches as f32 * 100.0;
                if batch_idx % (total_batches / 10).max(1) == 0 || batch_idx == total_batches - 1 {
                    println!(
                        " Progress: {:.1}% ( {} / {} batches)",
                        progress,
                        batch_idx + 1,
                        total_batches
                    );
                }
            }

            // Load current batch of cells
            let cell_batch = reader.read_cells_parallel(&cell_indices);

            // Distribute cell data to genes
            for cell in cell_batch {
                let cell_id = cell.original_index as u32;

                for (idx, &gene_id) in cell.indices.iter().enumerate() {
                    let raw_count = cell.data_raw[idx];
                    let norm_count = cell.data_norm[idx];

                    gene_data_map
                        .entry(gene_id)
                        .or_default()
                        .push((cell_id, raw_count, norm_count));
                }
            }
        }

        if verbose {
            println!("All cells processed, writing gene chunks to disk");
        }

        // Write gene data to disk
        for gene_id in 0..no_genes {
            let gene_id_u16 = gene_id as u16;

            if let Some(mut gene_entries) = gene_data_map.remove(&gene_id_u16) {
                // Sort by cell_id to maintain order
                gene_entries.sort_by_key(|&(cell_id, _, _)| cell_id);

                // Extract data for this gene
                let data_raw: Vec<u16> = gene_entries.iter().map(|(_, raw, _)| *raw).collect();
                let data_norm: Vec<F16> = gene_entries.iter().map(|(_, _, norm)| *norm).collect();
                let row_indices: Vec<usize> = gene_entries
                    .iter()
                    .map(|(cell_id, _, _)| *cell_id as usize)
                    .collect();

                let chunk = CscGeneChunk::from_conversion(
                    &data_raw,
                    &data_norm,
                    &row_indices,
                    gene_id,
                    true,
                );

                writer.write_gene_chunk(chunk).unwrap();
            } else {
                // Gene has no expression - write empty chunk
                let chunk = CscGeneChunk::from_conversion(&[], &[], &[], gene_id, true);
                writer.write_gene_chunk(chunk).unwrap();
            }
        }

        let end_conversion = start_conversion.elapsed();

        if verbose {
            println!(
                "Convertion data into gene-friendly format done: {:.2?}",
                end_conversion
            );
        }

        writer.finalise().unwrap();
    }

    /// Generate gene-based data with memory-bounded accumulation
    ///
    /// This processes genes in phases to limit memory usage. Each phase:
    ///
    /// 1. Reads all cells (unavoidable for CSC conversion)
    /// 2. Only accumulates data for genes in current phase
    /// 3. Writes those genes to disk
    /// 4. Clears memory and moves to next phase
    ///
    /// ### Params
    ///
    /// * `max_genes_in_memory` - Maximum genes to accumulate at once
    ///   (e.g., 2000)
    /// * `cell_batch_size` - How many cells to process at once
    ///   (e.g., 100000)
    /// * `verbose` - Controls verbosity
    pub fn generate_gene_based_data_memory_bounded(
        &mut self,
        max_genes_in_memory: usize,
        cell_batch_size: usize,
        verbose: bool,
    ) {
        let reader = Arc::new(ParallelSparseReader::new(&self.f_path_cells).unwrap());
        let header = reader.get_header();
        let no_cells = header.total_cells;
        let no_genes = header.total_genes;

        if verbose {
            println!("Converting CSR to CSC with memory-bounded approach:");
            println!("  Total cells: {}", no_cells.separate_with_underscores());
            println!("  Total genes: {}", no_genes.separate_with_underscores());
            println!(
                "  Processing {} genes per phase",
                max_genes_in_memory.separate_with_underscores()
            );
        }

        let start_conversion = Instant::now();

        let mut writer =
            CellGeneSparseWriter::new(&self.f_path_genes, false, no_cells, no_genes).unwrap();

        let num_phases = no_genes.div_ceil(max_genes_in_memory);

        for phase in 0..num_phases {
            let gene_phase_start = phase * max_genes_in_memory;
            let gene_phase_end = ((phase + 1) * max_genes_in_memory).min(no_genes);

            if verbose {
                println!(
                    "Phase {}/{}: Processing genes {}-{}",
                    phase + 1,
                    num_phases,
                    gene_phase_start,
                    gene_phase_end
                );
            }

            let num_cell_batches = no_cells.div_ceil(cell_batch_size);
            let cell_batches: Vec<(usize, usize)> = (0..num_cell_batches)
                .map(|i| {
                    let start = i * cell_batch_size;
                    let end = ((i + 1) * cell_batch_size).min(no_cells);
                    (start, end)
                })
                .collect();

            // Parallel accumulation
            let gene_data_parts: GeneData = cell_batches
                .par_iter()
                .map(|&(cell_start, cell_end)| {
                    let mut local_gene_data: FxHashMap<u16, Vec<(u32, u16, F16)>> =
                        FxHashMap::default();

                    // Create a separate reader for this thread
                    if let Ok(local_reader) = ParallelSparseReader::new(&self.f_path_cells) {
                        let cells = local_reader.read_cells_range(cell_start, cell_end);

                        for cell in cells {
                            let cell_id = cell.original_index as u32;

                            for (idx, &gene_id) in cell.indices.iter().enumerate() {
                                if (gene_id as usize) >= gene_phase_start
                                    && (gene_id as usize) < gene_phase_end
                                {
                                    let raw_count = cell.data_raw[idx];
                                    let norm_count = cell.data_norm[idx];

                                    local_gene_data
                                        .entry(gene_id)
                                        .or_insert_with(|| Vec::with_capacity(100))
                                        .push((cell_id, raw_count, norm_count));
                                }
                            }
                        }
                    }

                    local_gene_data
                })
                .collect();

            if verbose {
                println!("  Merging parallel results...");
            }

            // Merge results
            let mut gene_data: FxHashMap<u16, Vec<(u32, u16, F16)>> = FxHashMap::default();
            for local_data in gene_data_parts {
                for (gene_id, mut entries) in local_data {
                    gene_data.entry(gene_id).or_default().append(&mut entries);
                }
            }

            if verbose {
                println!("  Writing {} genes to disk...", gene_data.len());
            }

            // Write genes in sorted order
            for gene_id in gene_phase_start..gene_phase_end {
                if let Some(mut entries) = gene_data.remove(&(gene_id as u16)) {
                    entries.sort_unstable_by_key(|&(cell_id, _, _)| cell_id);

                    let data_raw: Vec<u16> = entries.iter().map(|(_, raw, _)| *raw).collect();
                    let data_norm: Vec<F16> = entries.iter().map(|(_, _, norm)| *norm).collect();
                    let row_indices: Vec<usize> = entries
                        .iter()
                        .map(|(cell_id, _, _)| *cell_id as usize)
                        .collect();

                    let chunk = CscGeneChunk::from_conversion(
                        &data_raw,
                        &data_norm,
                        &row_indices,
                        gene_id,
                        true,
                    );
                    writer.write_gene_chunk(chunk).unwrap();
                } else {
                    let empty_chunk = CscGeneChunk::from_conversion(&[], &[], &[], gene_id, true);
                    writer.write_gene_chunk(empty_chunk).unwrap();
                }
            }

            if verbose {
                println!("Phase {}/{} complete", phase + 1, num_phases);
            }
        }

        writer.finalise().unwrap();

        let end_conversion = start_conversion.elapsed();

        if verbose {
            println!(
                "Conversion into gene-friendly format done: {:.2?}",
                end_conversion
            );
        }
    }

    /// Return genes by index positions
    ///
    /// Leverages the CSC-stored data for fast gene retrieval
    ///
    /// ### Params
    ///
    /// * `indices` - The gene indices which to return
    /// * `assay` - Shall the raw or norm counts be returned
    ///
    /// ### Returns
    ///
    /// A list that can be parsed into a CSC matrix in R
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
        let mut col_idx: Vec<Vec<i32>> = Vec::new();
        let mut row_ptr: Vec<usize> = Vec::new();

        let mut current_col_ptr = 0_usize;
        row_ptr.push(current_col_ptr);

        for gene in &genes {
            let (indices_i, data_i) =
                get_gene_data(&gene.indices, &gene.data_raw, &gene.data_norm, &assay_type);
            let len_data_i = data_i.len();
            current_col_ptr += len_data_i;
            data.push(data_i);
            col_idx.push(indices_i);
            row_ptr.push(current_col_ptr);
        }

        let data = AssayData::flatten_into_r_vector(data);
        let col_idx = flatten_vector(col_idx);

        list!(
            indptr = row_ptr,
            indices = col_idx,
            data = data,
            no_cells = no_cells,
            no_genes = indices.len()
        )
    }

    /// Add cell indices
    ///
    /// ### Param
    ///
    /// * `cell_idx` - The cell indices to keep
    pub fn add_cells_to_keep(&mut self, cell_idx: Vec<i32>) {
        let cell_idx: Vec<u32> = cell_idx.into_iter().map(|x| x as u32).collect();
        let cell_idx_set: FxHashSet<u32> = cell_idx.into_iter().collect();
        self.cell_mask = cell_idx_set;
    }

    /// Set cell numbers and genes
    ///
    /// ### Params
    ///
    /// * `cell_no` - No of cells
    /// * `gene_no` - No of genes
    pub fn set_from_file(&mut self) {
        let reader = ParallelSparseReader::new(&self.f_path_cells).unwrap();
        let header = reader.get_header();
        self.n_cells = header.total_cells;
        self.n_genes = header.total_genes;
    }
}

extendr_module! {
    mod r_count_obj;
    impl SingeCellCountData;
}
