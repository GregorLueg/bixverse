use extendr_api::prelude::*;
use rayon::prelude::*;
use rustc_hash::FxHashMap;
use std::sync::Arc;
use std::time::Instant;
use thousands::Separable;

use bixverse_rs::prelude::*;
use bixverse_rs::single_cell::sc_data::{
    bin_merge_io::*, data_io::*, h5ad_io::*, h5ad_multifile_io::*, mtx_io::*, mtx_multifile_io::*,
    r_obj_io::*,
};

// Extendr unfortunately cannot do Roxygen2 manipulation of R6 type
// classes. This will have to be done manually in R... Documentation
// still here to make it easier.

/////////////
// extendR //
/////////////

extendr_module! {
    mod r_count_obj;
    impl SingleCellCountData;
}

//////////////////
// Type aliases //
//////////////////

/// Type to store the GeneData during streaming
type GeneData = Vec<FxHashMap<u32, Vec<(u32, u32, F16)>>>;

/// Type to store the Results structure
type GeneDataRes = std::result::Result<FxHashMap<u32, Vec<(u32, u32, F16)>>, BixverseErrors>;

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
    fn len(&self) -> usize {
        match self {
            AssayData::Raw(data) => data.len(),
            AssayData::Norm(data) => data.len(),
        }
    }

    /// Flatten the data into an R vector
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

/// Retrieve cell data for a given assay type
fn get_cell_data(
    indices: &[u32],
    data_raw: &RawCounts,
    data_norm: &[F16],
    assay_type: &AssayType,
) -> (Vec<i32>, AssayData) {
    let all_indices: Vec<i32> = indices.iter().map(|&x| x as i32).collect();

    let data = match assay_type {
        AssayType::Raw => AssayData::Raw(data_raw.iter().map(|x| x as i32).collect()),
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

/// Retrieve gene data for a given assay type
fn get_gene_data(
    indices: &[u32],
    data_raw: &RawCounts,
    data_norm: &[F16],
    assay_type: &AssayType,
) -> (Vec<i32>, AssayData) {
    let all_indices: Vec<i32> = indices.iter().map(|&x| x as i32).collect();

    let data = match assay_type {
        AssayType::Raw => AssayData::Raw(data_raw.iter().map(|x| x as i32).collect()),
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

/// Parse a count type string into an `AssayType`
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

/// Single cell count data handler
///
/// @description
/// A class for handling single cell count data stored on disk in two
/// complementary binary representations: a CSR-like layout (`f_path_cells`)
/// for fast cell-wise access and a CSC-like layout (`f_path_genes`) for fast
/// gene-wise access. Both raw counts and log-normalised counts are stored
/// side by side. Provides methods for ingesting data from R, `h5ad`, and
/// `mtx` sources (including multi-file workflows), converting between
/// layouts, retrieving slices of the matrix, and merging existing binary
/// objects.
///
/// @usage NULL
/// @format NULL
///
/// @param f_path_cells (`character`)\cr
/// Path to the `.bin` file for the cell-based (CSR-like) representation.
/// @param f_path_genes (`character`)\cr
/// Path to the `.bin` file for the gene-based (CSC-like) representation.
/// @param n_cells (`integer`)\cr
/// Number of cells represented in the data.
/// @param n_genes (`integer`)\cr
/// Number of genes represented in the data.
///
/// @return A new instance of the `SingleCellCountData` class.
///
/// @export
#[extendr]
struct SingleCellCountData {
    pub f_path_cells: String,
    pub f_path_genes: String,
    pub n_cells: usize,
    pub n_genes: usize,
}

#[extendr]
impl SingleCellCountData {
    /// Create a new instance of the class
    ///
    /// @param f_path_cells (`character`)\cr
    /// Path to the `.bin` file for the cell-based representation.
    /// @param f_path_genes (`character`)\cr
    /// Path to the `.bin` file for the gene-based representation.
    ///
    /// @return A new `SingleCellCountData` instance with `n_cells` and
    /// `n_genes` initialised to zero.
    pub fn new(f_path_cells: String, f_path_genes: String) -> Self {
        Self {
            f_path_cells,
            f_path_genes,
            n_cells: usize::default(),
            n_genes: usize::default(),
        }
    }

    /////////////
    // Helpers //
    /////////////

    /// Get the shape of the matrix
    ///
    /// @return An integer vector `c(n_cells, n_genes)`.
    pub fn get_shape(&mut self) -> Vec<usize> {
        vec![self.n_cells, self.n_genes]
    }

    /// Populate `n_cells` and `n_genes` from the cells binary file
    ///
    /// @description
    /// Reads the header of the file at `f_path_cells` and updates the
    /// `n_cells` and `n_genes` fields accordingly. Useful when reconnecting
    /// to an existing object on disk.
    ///
    /// @return Invisible `NULL`.
    pub fn set_from_file(&mut self) -> Result<(), extendr_api::Error> {
        let reader = ParallelSparseReader::new(&self.f_path_cells).to_extendr()?;
        let header = reader.get_header();
        self.n_cells = header.total_cells;
        self.n_genes = header.total_genes;

        Ok(())
    }

    /////////////////////
    // Cells ingestion //
    /////////////////////

    ////////////
    // From R //
    ////////////

    /// Write a CSR matrix from R to the cells binary file
    ///
    /// @description
    /// Ingest a sparse matrix passed in from R, apply per-cell QC, and write
    /// the result to `f_path_cells`.
    ///
    /// @param r_data (`list`)\cr
    /// A list convertible into `CompressedSparseData2`. Must contain the
    /// elements `"indptr"`, `"indices"`, `"data"`, `"nrow"`, `"ncol"` and
    /// `"format"`.
    /// @param qc_params (`list`)\cr
    /// Quality control parameters parseable into `MinCellQuality`.
    /// @param verbose (`logical`)\cr
    /// Controls verbosity of the function.
    ///
    /// @return A list with `cell_indices`, `gene_indices`, `lib_size` and
    /// `nnz`.
    pub fn r_data_to_file(
        &mut self,
        r_data: List,
        qc_params: List,
        verbose: bool,
    ) -> Result<List, extendr_api::Error> {
        let qc_params = MinCellQuality::from_r_list(qc_params)?;

        let compressed_data: CompressedSparseData2<u32> =
            list_to_sparse_matrix(r_data, false).to_extendr()?;

        let (no_cells, no_genes, cell_qc): (usize, usize, CellQuality) =
            write_r_counts(&self.f_path_cells, compressed_data, qc_params, verbose);

        self.n_cells = no_cells;
        self.n_genes = no_genes;

        Ok(list!(
            cell_indices = cell_qc.cell_indices,
            gene_indices = cell_qc.gene_indices,
            lib_size = cell_qc.lib_size,
            nnz = cell_qc.nnz
        ))
    }

    ///////////////
    // From h5ad //
    ///////////////

    /// Write an h5ad file to the cells binary file
    ///
    /// @param cs_type (`character`)\cr
    /// Storage layout of the h5ad data. One of `"CSC"` or `"CSR"`.
    /// @param h5_path (`character`)\cr
    /// Path to the h5ad file.
    /// @param no_cells (`integer`)\cr
    /// Number of cells in the h5ad file.
    /// @param no_genes (`integer`)\cr
    /// Number of genes in the h5ad file.
    /// @param qc_params (`list`)\cr
    /// Quality control parameters parseable into `MinCellQuality`.
    /// @param verbose (`logical`)\cr
    /// Controls verbosity of the function.
    ///
    /// @return A list with `cell_indices`, `gene_indices`, `lib_size` and
    /// `nnz`.
    pub fn h5ad_to_file(
        &mut self,
        cs_type: String,
        h5_path: String,
        no_cells: usize,
        no_genes: usize,
        qc_params: List,
        verbose: bool,
    ) -> extendr_api::Result<List> {
        let qc_params = MinCellQuality::from_r_list(qc_params)?;

        let (no_cells, no_genes, cell_qc) = write_h5_counts(
            &h5_path,
            &self.f_path_cells,
            &cs_type,
            no_cells,
            no_genes,
            qc_params,
            verbose,
        )
        .to_extendr()?;

        self.n_cells = no_cells;
        self.n_genes = no_genes;

        Ok(list!(
            cell_indices = cell_qc.cell_indices,
            gene_indices = cell_qc.gene_indices,
            lib_size = cell_qc.lib_size,
            nnz = cell_qc.nnz
        ))
    }

    /// Write an h5ad file with normalised counts to the cells binary file
    ///
    /// @description
    /// For data sets where only normalised counts are available in `X`.
    /// Reads library sizes from a specified `obs` column to reconstruct raw
    /// counts before writing.
    ///
    /// @param cs_type (`character`)\cr
    /// Storage layout of the h5 data. One of `"CSC"` or `"CSR"`.
    /// @param h5_path (`character`)\cr
    /// Path to the h5 file.
    /// @param no_cells (`integer`)\cr
    /// Number of cells in the h5 file.
    /// @param no_genes (`integer`)\cr
    /// Number of genes in the h5 file.
    /// @param obs_lib_size_col (`character`)\cr
    /// Name of the `obs` column containing total counts per cell
    /// (e.g. `"nCount_RNA"`).
    /// @param target_size (`numeric`)\cr
    /// Target size used in the original normalisation (e.g. `1e4`).
    /// @param qc_params (`list`)\cr
    /// Quality control parameters parseable into `MinCellQuality`.
    /// @param verbose (`logical`)\cr
    /// Controls verbosity of the function.
    ///
    /// @return A list with `cell_indices`, `gene_indices`, `lib_size` and
    /// `nnz`.
    #[allow(clippy::too_many_arguments)]
    pub fn norm_h5ad_to_file(
        &mut self,
        cs_type: String,
        h5_path: String,
        no_cells: usize,
        no_genes: usize,
        obs_lib_size_col: String,
        target_size: f64,
        qc_params: List,
        verbose: bool,
    ) -> extendr_api::Result<List> {
        let qc_params = MinCellQuality::from_r_list(qc_params)?;

        let (no_cells, no_genes, cell_qc) = write_h5_normalised_counts(
            &h5_path,
            &self.f_path_cells,
            &cs_type,
            no_cells,
            no_genes,
            &obs_lib_size_col,
            target_size as f32,
            qc_params,
            verbose,
        )
        .to_extendr()?;

        self.n_cells = no_cells;
        self.n_genes = no_genes;

        Ok(list!(
            cell_indices = cell_qc.cell_indices,
            gene_indices = cell_qc.gene_indices,
            lib_size = cell_qc.lib_size,
            nnz = cell_qc.nnz
        ))
    }

    /// Write an h5ad file to disk using streaming
    ///
    /// @description
    /// Slower but lighter on memory than `h5ad_to_file`; streams the input
    /// where possible.
    ///
    /// @param cs_type (`character`)\cr
    /// Storage layout of the h5 data. One of `"CSC"` or `"CSR"`.
    /// @param h5_path (`character`)\cr
    /// Path to the h5 file.
    /// @param no_cells (`integer`)\cr
    /// Number of cells in the h5 file.
    /// @param no_genes (`integer`)\cr
    /// Number of genes in the h5 file.
    /// @param qc_params (`list`)\cr
    /// Quality control parameters parseable into `MinCellQuality`.
    /// @param verbose (`logical`)\cr
    /// Controls verbosity of the function.
    ///
    /// @return A list with `cell_indices`, `gene_indices`, `lib_size` and
    /// `nnz`.
    pub fn h5ad_to_file_streaming(
        &mut self,
        cs_type: String,
        h5_path: String,
        no_cells: usize,
        no_genes: usize,
        qc_params: List,
        verbose: bool,
    ) -> extendr_api::Result<List> {
        let qc_params = MinCellQuality::from_r_list(qc_params)?;

        let (no_cells, no_genes, cell_qc) = stream_h5_counts(
            &h5_path,
            &self.f_path_cells,
            &cs_type,
            no_cells,
            no_genes,
            qc_params,
            verbose,
        )
        .to_extendr()?;

        self.n_cells = no_cells;
        self.n_genes = no_genes;

        Ok(list!(
            cell_indices = cell_qc.cell_indices,
            gene_indices = cell_qc.gene_indices,
            lib_size = cell_qc.lib_size,
            nnz = cell_qc.nnz
        ))
    }

    /// Load multiple h5ad files into a single binary
    ///
    /// @param file_tasks (`list`)\cr
    /// A list of lists, each produced by the R prescan function. Each inner
    /// list must contain `exp_id`, `h5_path`, `cs_type`, `no_cells`,
    /// `no_genes` and `gene_local_to_universe`.
    /// @param universe_size (`integer`)\cr
    /// Total number of genes in the universe.
    /// @param qc_params (`list`)\cr
    /// Quality control parameters (`min_unique_genes`, `min_lib_size`,
    /// `min_cells`, `target_size`).
    /// @param verbose (`logical`)\cr
    /// Controls verbosity.
    ///
    /// @return A list with `global_gene_indices`, `total_cells`,
    /// `total_genes` and `per_file` (a list of lists with `exp_id`,
    /// `cell_indices`, `lib_size`, `nnz`).
    pub fn multi_h5ad_to_file(
        &mut self,
        file_tasks: List,
        universe_size: i32,
        qc_params: List,
        verbose: bool,
    ) -> Result<List, extendr_api::Error> {
        let qc = MinCellQuality::from_r_list(qc_params)?;

        let tasks: Vec<H5adFileTask> = file_tasks
            .into_iter()
            .map(|(_, robj)| {
                let inner_list = List::try_from(robj).expect("Each file_task must be a list");
                H5adFileTask::from_r_list(inner_list)
            })
            .collect::<Result<Vec<_>, _>>()?;

        let result = multi_h5ad_to_file(
            &tasks,
            &self.f_path_cells,
            universe_size as usize,
            &qc,
            verbose,
        )
        .to_extendr()?;

        self.n_cells = result.total_cells;
        self.n_genes = result.total_genes;

        let per_file: List = result
            .per_file
            .into_iter()
            .map(|f| {
                list!(
                    exp_id = f.exp_id,
                    cell_indices = f.cells_to_keep,
                    lib_size = f.lib_size,
                    nnz = f.nnz
                )
            })
            .collect::<List>();

        Ok(list!(
            global_gene_indices = result.global_gene_indices,
            total_cells = result.total_cells,
            total_genes = result.total_genes,
            per_file = per_file
        ))
    }

    //////////////
    // From mtx //
    //////////////

    /// Write an mtx file to the cells binary file
    ///
    /// @param mtx_path (`character`)\cr
    /// Path to the mtx file.
    /// @param qc_params (`list`)\cr
    /// Quality control parameters parseable into `MinCellQuality`.
    /// @param cells_as_rows (`logical`)\cr
    /// `TRUE` if cells are rows in the mtx file, `FALSE` if cells are
    /// columns.
    /// @param verbose (`logical`)\cr
    /// Controls verbosity of the function.
    ///
    /// @return A list with `cell_indices`, `gene_indices`, `lib_size` and
    /// `nnz`.
    pub fn mtx_to_file(
        &mut self,
        mtx_path: String,
        qc_params: List,
        cells_as_rows: bool,
        verbose: bool,
    ) -> extendr_api::Result<List> {
        let qc_params = MinCellQuality::from_r_list(qc_params)?;

        let mut mtx_reader = MtxReader::new(&mtx_path, qc_params, cells_as_rows).to_extendr()?;

        let mtx_quality_data = mtx_reader.parse_mtx_quality(verbose).to_extendr()?;

        let mtx_res: MtxFinalData = mtx_reader
            .process_mtx_and_write_bin(&self.f_path_cells, &mtx_quality_data, verbose)
            .to_extendr()?;

        self.n_cells = mtx_res.no_cells;
        self.n_genes = mtx_res.no_genes;

        Ok(list!(
            cell_indices = mtx_res.cell_qc.cell_indices,
            gene_indices = mtx_res.cell_qc.gene_indices,
            lib_size = mtx_res.cell_qc.lib_size,
            nnz = mtx_res.cell_qc.nnz
        ))
    }

    /// Write an mtx file to the cells binary file using streaming
    ///
    /// @param mtx_path (`character`)\cr
    /// Path to the mtx file.
    /// @param qc_params (`list`)\cr
    /// Quality control parameters parseable into `MinCellQuality`.
    /// @param cells_as_rows (`logical`)\cr
    /// `TRUE` if cells are rows in the mtx file, `FALSE` if cells are
    /// columns.
    /// @param verbose (`logical`)\cr
    /// Controls verbosity of the function.
    ///
    /// @return A list with `cell_indices`, `gene_indices`, `lib_size` and
    /// `nnz`.
    pub fn mtx_to_file_streaming(
        &mut self,
        mtx_path: String,
        qc_params: List,
        cells_as_rows: bool,
        verbose: bool,
    ) -> extendr_api::Result<List> {
        let qc_params = MinCellQuality::from_r_list(qc_params)?;

        let mut mtx_reader = MtxReader::new(&mtx_path, qc_params, cells_as_rows).to_extendr()?;

        let mtx_quality_data = mtx_reader.parse_mtx_quality(verbose).to_extendr()?;

        let mtx_res: MtxFinalData = mtx_reader
            .process_mtx_and_write_bin_streaming(&self.f_path_cells, &mtx_quality_data, verbose)
            .to_extendr()?;

        self.n_cells = mtx_res.no_cells;
        self.n_genes = mtx_res.no_genes;

        Ok(list!(
            cell_indices = mtx_res.cell_qc.cell_indices,
            gene_indices = mtx_res.cell_qc.gene_indices,
            lib_size = mtx_res.cell_qc.lib_size,
            nnz = mtx_res.cell_qc.nnz
        ))
    }

    /// Load multiple mtx files into a single binary
    ///
    /// @param file_tasks (`list`)\cr
    /// A list of lists, each containing `exp_id`, `mtx_path`,
    /// `cells_as_rows` and `gene_local_to_universe` (integer vector, `NA`
    /// for unmapped genes).
    /// @param universe_size (`integer`)\cr
    /// Number of genes in the intersection universe.
    /// @param qc_params (`list`)\cr
    /// Quality control parameters parseable into `MinCellQuality`.
    /// @param verbose (`logical`)\cr
    /// Controls verbosity.
    ///
    /// @return A list with `global_gene_indices`, `total_cells`,
    /// `total_genes` and `per_file` (a list of lists with `exp_id`,
    /// `cell_indices`, `lib_size`, `nnz`).
    pub fn multi_mtx_to_file(
        &mut self,
        file_tasks: List,
        universe_size: i32,
        qc_params: List,
        verbose: bool,
    ) -> Result<List, extendr_api::Error> {
        let qc = MinCellQuality::from_r_list(qc_params)?;

        let tasks: Vec<MtxFileTask> = file_tasks
            .into_iter()
            .map(|(_, robj)| {
                let inner = List::try_from(robj).expect("Each file_task must be a list");
                MtxFileTask::from_r_list(inner)
            })
            .collect::<Result<Vec<_>, _>>()?;

        let result = multi_mtx_to_file(
            &tasks,
            &self.f_path_cells,
            universe_size as usize,
            &qc,
            verbose,
        )
        .to_extendr()?;

        self.n_cells = result.total_cells;
        self.n_genes = result.total_genes;

        let per_file: List = result
            .per_file
            .into_iter()
            .map(|f| {
                list!(
                    exp_id = f.exp_id,
                    cell_indices = f.cells_to_keep,
                    lib_size = f.lib_size,
                    nnz = f.nnz
                )
            })
            .collect::<List>();

        Ok(list!(
            global_gene_indices = result.global_gene_indices,
            total_cells = result.total_cells,
            total_genes = result.total_genes,
            per_file = per_file
        ))
    }

    //////////////////////////////
    // Return cell-based counts //
    //////////////////////////////

    /// Return the full count matrix
    ///
    /// @param assay (`character`)\cr
    /// One of `"raw"` or `"norm"`. Selects whether raw counts or
    /// log-normalised counts are returned.
    /// @param cell_based (`logical`)\cr
    /// If `TRUE`, the data is returned in CSR layout (cells as rows). If
    /// `FALSE`, the data is returned in CSC layout (genes as columns).
    /// @param verbose (`logical`)\cr
    /// Controls verbosity of the function.
    ///
    /// @return A list with `indptr`, `indices`, `data`, `no_cells` and
    /// `no_genes`, parseable into a sparse matrix in R.
    pub fn return_full_mat(
        &self,
        assay: &str,
        cell_based: bool,
        verbose: bool,
    ) -> Result<List, extendr_api::Error> {
        let mut data: Vec<AssayData> = Vec::new();
        let mut indices: Vec<Vec<i32>> = Vec::new();
        let mut indptr: Vec<usize> = Vec::new();
        let assay_type = parse_count_type(assay).unwrap();

        if cell_based {
            let reader = ParallelSparseReader::new(&self.f_path_cells).to_extendr()?;
            let cell_chunks = reader.get_all_cells().to_extendr()?;

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
                data.push(data_i);
                indices.push(indices_i);
                indptr.push(current_ptr);
            }
        } else {
            let reader = ParallelSparseReader::new(&self.f_path_genes).to_extendr()?;
            let gene_chunks = reader.get_all_genes().to_extendr()?;

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

        Ok(list!(
            indptr = indptr,
            indices = indices,
            data = data,
            no_cells = self.n_cells,
            no_genes = self.n_genes
        ))
    }

    /// Return cells by index positions
    ///
    /// @description
    /// Leverages the CSR-stored data for fast cell retrieval.
    ///
    /// @param indices (`integer`)\cr
    /// The cell indices to return (1-indexed).
    /// @param assay (`character`)\cr
    /// One of `"raw"` or `"norm"`.
    ///
    /// @return A list with `indptr`, `indices`, `data`, `no_cells` and
    /// `no_genes`, parseable into a CSR matrix in R.
    pub fn get_cells_by_indices(
        &self,
        indices: &[i32],
        assay: &str,
    ) -> Result<List, extendr_api::Error> {
        let reader = ParallelSparseReader::new(&self.f_path_cells).to_extendr()?;
        let assay_type = parse_count_type(assay).unwrap();

        let indices: Vec<usize> = indices.iter().map(|x| (*x - 1) as usize).collect();

        let cells = reader.read_cells_parallel(&indices).to_extendr()?;

        let results: Vec<(Vec<i32>, AssayData)> = cells
            .par_iter()
            .map(|cell| get_cell_data(&cell.indices, &cell.data_raw, &cell.data_norm, &assay_type))
            .collect();

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

        Ok(list!(
            indptr = row_ptr,
            indices = col_idx,
            data = data,
            no_cells = cells.len(),
            no_genes = self.n_genes
        ))
    }

    ///////////
    // Genes //
    ///////////

    //////////////////////////////
    // Transform the CSR to CSC //
    //////////////////////////////

    /// Generate gene-based data from the cells binary file
    ///
    /// @description
    /// Reads the `.bin` file at `f_path_cells` and writes a gene-friendly
    /// (CSC) representation to `f_path_genes`. The conversion happens fully
    /// in memory and may cause memory pressure on large data sets; see
    /// `generate_gene_based_data_streaming` or
    /// `generate_gene_based_data_memory_bounded` for lighter alternatives.
    ///
    /// @param verbose (`logical`)\cr
    /// Controls verbosity of the function.
    ///
    /// @return Invisible `NULL`.
    pub fn generate_gene_based_data(&mut self, verbose: bool) -> Result<(), extendr_api::Error> {
        let reader = ParallelSparseReader::new(&self.f_path_cells).to_extendr()?;

        let no_cells = reader.get_header().total_cells;
        let no_genes = reader.get_header().total_genes;

        let start_conversion = Instant::now();

        if verbose {
            println!("Loading in the cell data and saving it into a gene-friendly file")
        }

        let all_cells: Vec<_> = reader.get_all_cells().to_extendr()?;

        let mut data: Vec<Vec<u32>> = Vec::new();
        let mut data_2: Vec<Vec<F16>> = Vec::new();
        let mut col_idx: Vec<Vec<u32>> = Vec::new();
        let mut row_ptr: Vec<usize> = Vec::new();

        let mut current_row_ptr = 0_usize;
        row_ptr.push(current_row_ptr);

        for cell in all_cells {
            let data_i: Vec<u32> = cell.data_raw.iter().collect();
            let data_norm_i = cell.data_norm;

            let len_data_i = data_i.len();
            current_row_ptr += len_data_i;
            data.push(data_i);
            data_2.push(data_norm_i);
            col_idx.push(cell.indices);
            row_ptr.push(current_row_ptr);
        }

        let data = flatten_vector(data);
        let data_2 = flatten_vector(data_2);
        let col_idx = flatten_vector(col_idx);
        let col_idx = col_idx.iter().map(|x| *x as usize).collect::<Vec<usize>>();

        let sparse_data = CompressedSparseData2::new_csr(
            &data,
            &col_idx,
            &row_ptr,
            Some(&data_2),
            (no_cells, no_genes),
        );

        let sparse_data = sparse_data.transform();
        let data_2 = sparse_data.get_data2_unsafe();

        let mut writer = CellGeneSparseWriter::new(&self.f_path_genes, false, no_cells, no_genes)
            .to_extendr()?;

        for i in 0..no_genes {
            let start_i = sparse_data.indptr[i];
            let end_i = sparse_data.indptr[i + 1];

            let raw_slice = &sparse_data.data[start_i..end_i];
            let chunk_i = CscGeneChunk::from_conversion(
                RawCounts::from_u32_auto(raw_slice),
                &data_2[start_i..end_i],
                &sparse_data.indices[start_i..end_i],
                i,
                true,
            );

            writer.write_gene_chunk(chunk_i).to_extendr()?;
        }

        let end_conversion = start_conversion.elapsed();

        if verbose {
            println!(
                "Convertion data into gene-friendly format done: {:.2?}",
                end_conversion
            );
        }

        writer.finalise().to_extendr()?;

        Ok(())
    }

    /// Generate gene-based data with streaming
    ///
    /// @description
    /// Builds the CSC representation directly without creating intermediate
    /// CSR structures. Suitable for very large data sets where the all-in-
    /// memory path is too costly.
    ///
    /// @param batch_size (`integer`)\cr
    /// Number of cells processed per batch. Larger values increase memory
    /// pressure but reduce overhead.
    /// @param verbose (`logical`)\cr
    /// Controls verbosity of the function.
    ///
    /// @return Invisible `NULL`.
    pub fn generate_gene_based_data_streaming(
        &mut self,
        batch_size: usize,
        verbose: bool,
    ) -> Result<(), extendr_api::Error> {
        let reader = ParallelSparseReader::new(&self.f_path_cells).to_extendr()?;
        let header = reader.get_header();
        let no_cells = header.total_cells;
        let no_genes = header.total_genes;

        if verbose {
            println!("Streaming cell data to gene-friendly format with reduced memory usage.");
        }

        let start_conversion = Instant::now();

        let mut gene_data_map: FxHashMap<u32, Vec<(u32, u32, F16)>> = FxHashMap::default();

        let total_batches = no_cells.div_ceil(batch_size);

        let mut writer = CellGeneSparseWriter::new(&self.f_path_genes, false, no_cells, no_genes)
            .to_extendr()?;

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

            let cell_batch = reader.read_cells_parallel(&cell_indices).to_extendr()?;

            for cell in cell_batch {
                let cell_id = cell.original_index as u32;

                for (idx, &gene_id) in cell.indices.iter().enumerate() {
                    let raw_count = cell.data_raw.get(idx);
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

        for gene_id in 0..no_genes {
            let gene_id_u32 = gene_id as u32;

            if let Some(mut gene_entries) = gene_data_map.remove(&gene_id_u32) {
                gene_entries.sort_by_key(|&(cell_id, _, _)| cell_id);

                let data_raw: Vec<u32> = gene_entries.iter().map(|(_, raw, _)| *raw).collect();
                let data_norm: Vec<F16> = gene_entries.iter().map(|(_, _, norm)| *norm).collect();
                let row_indices: Vec<usize> = gene_entries
                    .iter()
                    .map(|(cell_id, _, _)| *cell_id as usize)
                    .collect();

                let chunk = CscGeneChunk::from_conversion(
                    RawCounts::from_u32_auto(&data_raw),
                    &data_norm,
                    &row_indices,
                    gene_id,
                    true,
                );

                writer.write_gene_chunk(chunk).to_extendr()?;
            } else {
                let chunk = CscGeneChunk::from_conversion(
                    RawCounts::U16(Vec::new()),
                    &[],
                    &[],
                    gene_id,
                    true,
                );
                writer.write_gene_chunk(chunk).to_extendr()?;
            }
        }

        let end_conversion = start_conversion.elapsed();

        if verbose {
            println!(
                "Conversion of count data into gene-friendly format done: {:.2?}",
                end_conversion
            );
        }

        writer.finalise().to_extendr()?;

        Ok(())
    }

    /// Generate gene-based data with memory-bounded accumulation
    ///
    /// @description
    /// Processes genes in phases to cap memory usage. Each phase:
    /// \enumerate{
    ///   \item reads all cells (unavoidable for CSC conversion);
    ///   \item only accumulates data for genes in the current phase;
    ///   \item writes those genes to disk;
    ///   \item clears memory and moves to the next phase.
    /// }
    ///
    /// @param max_genes_in_memory (`integer`)\cr
    /// Maximum number of genes to accumulate at once (e.g. `2000`).
    /// @param cell_batch_size (`integer`)\cr
    /// Number of cells to process at once (e.g. `100000`).
    /// @param verbose (`logical`)\cr
    /// Controls verbosity.
    ///
    /// @return Invisible `NULL`.
    pub fn generate_gene_based_data_memory_bounded(
        &mut self,
        max_genes_in_memory: usize,
        cell_batch_size: usize,
        verbose: bool,
    ) -> Result<(), extendr_api::Error> {
        let reader = Arc::new(ParallelSparseReader::new(&self.f_path_cells).to_extendr()?);
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

        let mut writer = CellGeneSparseWriter::new(&self.f_path_genes, false, no_cells, no_genes)
            .to_extendr()?;

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

            let gene_data_parts: GeneData = cell_batches
                .par_iter()
                .map(|&(cell_start, cell_end)| -> GeneDataRes {
                    let mut local_gene_data: FxHashMap<u32, Vec<(u32, u32, F16)>> =
                        FxHashMap::default();

                    let local_reader = ParallelSparseReader::new(&self.f_path_cells)?;
                    let cells = local_reader.read_cells_range(cell_start, cell_end)?;

                    for cell in cells {
                        let cell_id = cell.original_index as u32;
                        for (idx, &gene_id) in cell.indices.iter().enumerate() {
                            if (gene_id as usize) >= gene_phase_start
                                && (gene_id as usize) < gene_phase_end
                            {
                                let raw_count = cell.data_raw.get(idx);
                                let norm_count = cell.data_norm[idx];
                                local_gene_data
                                    .entry(gene_id)
                                    .or_insert_with(|| Vec::with_capacity(100))
                                    .push((cell_id, raw_count, norm_count));
                            }
                        }
                    }

                    Ok(local_gene_data)
                })
                .collect::<std::result::Result<Vec<_>, BixverseErrors>>()
                .to_extendr()?;

            if verbose {
                println!("  Merging parallel results...");
            }

            let mut gene_data: FxHashMap<u32, Vec<(u32, u32, F16)>> = FxHashMap::default();
            for local_data in gene_data_parts {
                for (gene_id, mut entries) in local_data {
                    gene_data.entry(gene_id).or_default().append(&mut entries);
                }
            }

            if verbose {
                println!("  Writing {} genes to disk...", gene_data.len());
            }

            for gene_id in gene_phase_start..gene_phase_end {
                if let Some(mut entries) = gene_data.remove(&(gene_id as u32)) {
                    entries.sort_unstable_by_key(|&(cell_id, _, _)| cell_id);

                    let data_raw: Vec<u32> = entries.iter().map(|(_, raw, _)| *raw).collect();
                    let data_norm: Vec<F16> = entries.iter().map(|(_, _, norm)| *norm).collect();
                    let row_indices: Vec<usize> = entries
                        .iter()
                        .map(|(cell_id, _, _)| *cell_id as usize)
                        .collect();

                    let chunk = CscGeneChunk::from_conversion(
                        RawCounts::from_u32_auto(&data_raw),
                        &data_norm,
                        &row_indices,
                        gene_id,
                        true,
                    );
                    writer.write_gene_chunk(chunk).unwrap();
                } else {
                    let empty_chunk = CscGeneChunk::from_conversion(
                        RawCounts::U16(Vec::new()),
                        &[],
                        &[],
                        gene_id,
                        true,
                    );
                    writer.write_gene_chunk(empty_chunk).to_extendr()?;
                }
            }

            if verbose {
                println!("Phase {}/{} complete", phase + 1, num_phases);
            }
        }

        writer.finalise().to_extendr()?;

        let end_conversion = start_conversion.elapsed();

        if verbose {
            println!(
                "Conversion into gene-friendly format done: {:.2?}",
                end_conversion
            );
        }

        Ok(())
    }

    //////////////////////////////
    // Return gene-based counts //
    //////////////////////////////

    /// Return genes by index positions
    ///
    /// @description
    /// Leverages the CSC-stored data for fast gene retrieval.
    ///
    /// @param indices (`integer`)\cr
    /// The gene indices to return (1-indexed).
    /// @param assay (`character`)\cr
    /// One of `"raw"` or `"norm"`.
    ///
    /// @return A list with `indptr`, `indices`, `data`, `no_cells` and
    /// `no_genes`, parseable into a CSC matrix in R.
    pub fn get_genes_by_indices(
        &self,
        indices: &[i32],
        assay: &str,
    ) -> Result<List, extendr_api::Error> {
        let reader = ParallelSparseReader::new(&self.f_path_genes).to_extendr()?;

        let assay_type = parse_count_type(assay).unwrap();

        let no_cells = reader.get_header().total_cells;

        let indices = indices
            .iter()
            .map(|x| (*x - 1) as usize)
            .collect::<Vec<usize>>();

        let genes = reader.read_gene_parallel(&indices).to_extendr()?;

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

        Ok(list!(
            indptr = row_ptr,
            indices = col_idx,
            data = data,
            no_cells = no_cells,
            no_genes = indices.len()
        ))
    }

    /// Get the number of cells expressing each gene
    ///
    /// @param gene_indices (`integer` or `NULL`)\cr
    /// Optional 1-indexed gene indices. If `NULL`, results are returned for
    /// all genes.
    ///
    /// @return An integer vector of NNZ counts for the requested genes.
    pub fn get_nnz_genes(
        &mut self,
        gene_indices: Option<&[i32]>,
    ) -> Result<Vec<i32>, extendr_api::Error> {
        let reader = ParallelSparseReader::new(&self.f_path_genes).to_extendr()?;

        let nnz = match gene_indices {
            Some(indices) => {
                let gene_indices = indices
                    .iter()
                    .map(|&x| (x - 1) as usize)
                    .collect::<Vec<usize>>();
                reader.read_gene_nnz(&gene_indices).to_extendr()?
            }
            None => reader.get_all_gene_nnz().to_extendr()?,
        };

        Ok(nnz.r_int_convert())
    }

    ///////////////////////
    // Combining objects //
    ///////////////////////

    /// Merge multiple existing bin files into the cells binary file
    ///
    /// @param merge_tasks (`list`)\cr
    /// A list of lists. Each inner list must contain `exp_id`,
    /// `bin_cells_path`, `cells_to_keep` (0-indexed integer vector) and
    /// `gene_local_to_universe` (integer vector, `-1` for genes absent from
    /// the universe).
    /// @param universe_size (`integer`)\cr
    /// Number of genes in the intersection universe.
    /// @param renormalise (`logical`)\cr
    /// If `TRUE`, recompute `data_norm` against `target_size` using each
    /// cell's surviving raw counts. If `FALSE`, pass `data_norm` through
    /// untouched; the caller must guarantee all inputs were normalised
    /// against the same `target_size`.
    /// @param target_size (`numeric`)\cr
    /// Target library size for renormalisation. Ignored when
    /// `renormalise = FALSE`.
    /// @param verbose (`logical`)\cr
    /// Controls verbosity.
    ///
    /// @return A list with `total_cells`, `total_genes` and `per_file` (a
    /// list of lists with `exp_id`, `lib_size`, `nnz`).
    pub fn merge_sc_files(
        &mut self,
        merge_tasks: List,
        universe_size: i32,
        renormalise: bool,
        target_size: f64,
        verbose: bool,
    ) -> Result<List, extendr_api::Error> {
        let tasks: Vec<BinMergeTask> = merge_tasks
            .into_iter()
            .map(|(_, robj)| {
                let inner = List::try_from(robj).expect("Each merge_task must be a list");
                BinMergeTask::from_r_list(inner)
            })
            .collect::<Result<Vec<_>, _>>()?;

        let result = merge_sc_bin_files(
            &tasks,
            &self.f_path_cells,
            universe_size as usize,
            renormalise,
            target_size as f32,
            verbose,
        )
        .to_extendr()?;

        self.n_cells = result.total_cells;
        self.n_genes = result.total_genes;

        let per_file: List = result
            .per_file
            .into_iter()
            .map(|f| list!(exp_id = f.exp_id, lib_size = f.lib_size, nnz = f.nnz))
            .collect::<List>();

        Ok(list!(
            total_cells = result.total_cells,
            total_genes = result.total_genes,
            per_file = per_file
        ))
    }
}
