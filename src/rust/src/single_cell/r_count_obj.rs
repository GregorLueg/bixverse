use extendr_api::prelude::*;
use rayon::prelude::*;
use std::time::Instant;

use bixverse_rs::prelude::*;
use bixverse_rs::single_cell::sc_data::{
    archive_io::*, bin_merge_io::*, gene_file_io::*, h5_10x_io::*, h5_10x_multifile_io::*,
    h5ad_io::*, h5ad_multifile_io::*, mtx_io::*, mtx_multifile_io::*, r_obj_io::*,
};
use bixverse_rs::single_cell::sc_processing::cellsweep::{
    run_cellsweep, CellSweepParams, CellSweepSample,
};
use bixverse_rs::single_cell::sc_r_wrappers::cellsweep_sample_from_r_list;

/////////////
// extendR //
/////////////

extendr_module! {
    mod r_count_obj;
    impl SingleCellCountData;
}

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
    ///
    /// ### Params
    ///
    /// * `data` - Chunks of one variant; the first chunk decides the type.
    ///
    /// ### Returns
    ///
    /// An integer vector for `Raw`, a double vector for `Norm`, or an empty
    /// double vector if `data` is empty.
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

/// Retrieve one cell's data for a given assay type
///
/// ### Params
///
/// * `indices` - Gene indices (0-indexed) of the cell's non-zero entries
/// * `data_raw` - The raw counts, aligned to `indices`.
/// * `data_norm` - The normalised counts, aligned to `indices`.
/// * `assay_type` - Which assay type to return, see [AssayType]
///
/// ### Returns
///
/// A tuple of the gene indices as `i32` and the assay data.
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

/// Retrieve one gene's data for a given assay type
///
/// ### Params
///
/// * `indices` - Cell indices (0-indexed) of the gene's non-zero entries
/// * `data_raw` - The raw counts, aligned to `indices`.
/// * `data_norm` - The normalised counts, aligned to `indices`.
/// * `assay_type` - Which assay type to return, see [AssayType]
///
/// ### Returns
///
/// A tuple of the cell indices as `i32` and the assay data.
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
///
/// ### Params
///
/// * `s` - String to parse. One of `"raw"` or `"norm"`, case-insensitive.
///
/// ### Returns
///
/// The optional [AssayType]; `None` for anything else.
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
/// `r lifecycle::badge("experimental")`
/// A class for handling single cell count data stored on disk in two
/// complementary binary representations: a CSR-like layout (`f_path_cells`)
/// for fast cell-wise access and a CSC-like layout (`f_path_genes`) for fast
/// gene-wise access. Both raw counts and log-normalised counts are stored
/// side by side. Provides methods for ingesting data from R, `h5ad`, `mtx`
/// and 10x CellRanger h5 sources (including multi-file workflows),
/// converting between layouts, retrieving slices of the matrix, merging
/// existing binary objects and writing CellSweep-denoised counts.
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
/// @returns A new instance of the `SingleCellCountData` class.
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
    /// @returns A new `SingleCellCountData` instance with `n_cells` and
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
    /// @returns An integer vector `c(n_cells, n_genes)`.
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
    /// @returns Invisible `NULL`.
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
    /// A named list convertible into `CompressedSparseData2`. Must contain
    /// `"indptr"`, `"indices"` (0-indexed), `"data"`, `"nrow"`, `"ncol"` and
    /// `"cs_type"` (`"csr"` or `"csc"`).
    /// @param qc_params (`list`)\cr
    /// Quality control parameters parseable into `MinCellQuality`.
    /// @param verbose (`logical`)\cr
    /// Controls verbosity of the function.
    ///
    /// @returns A list with `cell_indices` and `gene_indices` (0-indexed,
    /// surviving QC), `lib_size` and `nnz`.
    pub fn r_data_to_file(
        &mut self,
        r_data: List,
        qc_params: List,
        verbose: bool,
    ) -> Result<List, extendr_api::Error> {
        let start = Instant::now();

        let qc_params = MinCellQuality::from_r_list(qc_params)?;

        if verbose {
            println!("Transforming R data into compressed sparse data.")
        }

        let compressed_data: CompressedSparseData2<u32> =
            list_to_sparse_matrix(r_data, false).to_extendr()?;

        if verbose {
            println!(" Done in {:.2?}", start.elapsed())
        }

        if verbose {
            println!("Preparing generation of binary files.")
        }

        let (no_cells, no_genes, cell_qc): (usize, usize, CellQuality) =
            write_r_counts(&self.f_path_cells, compressed_data, qc_params, verbose).to_extendr()?;

        if verbose {
            println!(" Done in {:.2?}", start.elapsed())
        }

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
    /// @param slot (`character`)\cr
    /// Where to find the raw counts. One of `"X"`, `"raw.X"` (or `"raw"`, for
    /// CellXGene data) or `"layers.counts"`. Unmatched values fall back to
    /// `"X"`.
    /// @param verbose (`logical`)\cr
    /// Controls verbosity of the function.
    ///
    /// @returns A list with `cell_indices` and `gene_indices` (0-indexed,
    /// surviving QC), `lib_size` and `nnz`.
    #[allow(clippy::too_many_arguments)]
    pub fn h5ad_to_file(
        &mut self,
        cs_type: String,
        h5_path: String,
        no_cells: usize,
        no_genes: usize,
        qc_params: List,
        slot: String,
        verbose: bool,
    ) -> extendr_api::Result<List> {
        let qc_params = MinCellQuality::from_r_list(qc_params)?;

        let raw_slot = parse_raw_slot(&slot).unwrap_or_else(|| {
            println!(
                "The provided string ({:?}) could not be matched. Defaulting to X",
                slot
            );
            RawDataSlot::default()
        });

        let (no_cells, no_genes, cell_qc) = write_h5_counts(
            &h5_path,
            &self.f_path_cells,
            &cs_type,
            no_cells,
            no_genes,
            qc_params,
            raw_slot,
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
    /// @returns A list with `cell_indices` and `gene_indices` (0-indexed,
    /// surviving QC), `lib_size` and `nnz`.
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
    /// @param slot (`character`)\cr
    /// Where to find the raw counts. One of `"X"`, `"raw.X"` (or `"raw"`, for
    /// CellXGene data) or `"layers.counts"`. Unmatched values fall back to
    /// `"X"`.
    /// @param verbose (`logical`)\cr
    /// Controls verbosity of the function.
    ///
    /// @returns A list with `cell_indices` and `gene_indices` (0-indexed,
    /// surviving QC), `lib_size` and `nnz`.
    #[allow(clippy::too_many_arguments)]
    pub fn h5ad_to_file_streaming(
        &mut self,
        cs_type: String,
        h5_path: String,
        no_cells: usize,
        no_genes: usize,
        qc_params: List,
        slot: String,
        verbose: bool,
    ) -> extendr_api::Result<List> {
        let qc_params = MinCellQuality::from_r_list(qc_params)?;

        let raw_slot = parse_raw_slot(&slot).unwrap_or_else(|| {
            println!(
                "The provided string ({:?}) could not be matched. Defaulting to X",
                slot
            );
            RawDataSlot::default()
        });

        let (no_cells, no_genes, cell_qc) = stream_h5_counts(
            &h5_path,
            &self.f_path_cells,
            &cs_type,
            no_cells,
            no_genes,
            qc_params,
            raw_slot,
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
    /// `no_genes` and `gene_local_to_universe` (0-indexed integer vector,
    /// `NA` for unmapped genes), plus an optional `raw_slot`.
    /// @param universe_size (`integer`)\cr
    /// Total number of genes in the universe.
    /// @param qc_params (`list`)\cr
    /// Quality control parameters (`min_unique_genes`, `min_lib_size`,
    /// `min_cells`, `target_size`).
    /// @param verbose (`logical`)\cr
    /// Controls verbosity.
    ///
    /// @returns A list with `global_gene_indices`, `total_cells`,
    /// `total_genes` and `per_file` (a list of lists with `exp_id`,
    /// `cell_indices` (file-local, 0-indexed), `lib_size`, `nnz`).
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
                let inner_list = List::try_from(robj)?;
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
    /// @returns A list with `cell_indices` and `gene_indices` (0-indexed,
    /// surviving QC), `lib_size` and `nnz`.
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
    /// @returns A list with `cell_indices` and `gene_indices` (0-indexed,
    /// surviving QC), `lib_size` and `nnz`.
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
    /// @returns A list with `global_gene_indices`, `total_cells`,
    /// `total_genes` and `per_file` (a list of lists with `exp_id`,
    /// `cell_indices` (file-local, 0-indexed), `lib_size`, `nnz`).
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
                let inner = List::try_from(robj)?;
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

    /////////////////////////
    // From h5 10x outputs //
    /////////////////////////

    /// Write a 10x CellRanger h5 file to the cells binary file using streaming
    ///
    /// @description
    /// Ingests the gene-expression modality from a CellRanger v2/v3 h5 file.
    /// Other modalities (e.g. Antibody Capture) are filtered out via
    /// `feature_type`.
    ///
    /// @param h5_path (`character`)\cr
    /// Path to the 10x h5 file.
    /// @param version (`character`)\cr
    /// One of `"auto"`, `"v2"` or `"v3"`. `"auto"`, and any unmatched value,
    /// detects the layout from the file.
    /// @param no_cells (`integer`)\cr
    /// Number of cells (columns) in the file.
    /// @param no_genes (`integer`)\cr
    /// Number of features (rows), including all modalities.
    /// @param qc_params (`list`)\cr
    /// Quality control parameters parseable into `MinCellQuality`.
    /// @param feature_type (`character` or `NULL`)\cr
    /// Target modality for v3. Defaults to `"Gene Expression"` when `NULL`.
    /// @param verbose (`logical`)\cr
    /// Controls verbosity of the function.
    ///
    /// @returns A list with `cell_indices` and `gene_indices` (0-indexed,
    /// surviving QC), `lib_size` and `nnz`.
    #[allow(clippy::too_many_arguments)]
    pub fn tenx_h5_to_file_streaming(
        &mut self,
        h5_path: String,
        version: String,
        no_cells: usize,
        no_genes: usize,
        qc_params: List,
        feature_type: Option<String>,
        verbose: bool,
    ) -> extendr_api::Result<List> {
        let qc_params = MinCellQuality::from_r_list(qc_params)?;

        let version = match version.to_lowercase().as_str() {
            "auto" => None,
            other => match parse_tenx_version(other) {
                Some(v) => Some(v),
                None => {
                    println!(
                            "The provided version ({:?}) could not be matched. Falling back to auto-detection.",
                            version
                        );
                    None
                }
            },
        };

        let (no_cells, no_genes, cell_qc) = stream_h5_tenx_counts(
            &h5_path,
            &self.f_path_cells,
            version,
            no_cells,
            no_genes,
            qc_params,
            feature_type.as_deref(),
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

    /// Load multiple 10x CellRanger h5 files into a single binary
    ///
    /// @param file_tasks (`list`)\cr
    /// A list of lists, each produced by the R prescan function. Each inner
    /// list must contain `exp_id`, `h5_path`, `version` (`"v2"` or `"v3"`),
    /// `no_cells`, `no_genes`, `gene_local_to_universe` (integer vector, `NA`
    /// for unmapped / non-gene features) and `feature_type` (optional string,
    /// defaults to `"Gene Expression"`).
    /// @param universe_size (`integer`)\cr
    /// Total number of genes in the universe.
    /// @param qc_params (`list`)\cr
    /// Quality control parameters (`min_unique_genes`, `min_lib_size`,
    /// `min_cells`, `target_size`).
    /// @param verbose (`logical`)\cr
    /// Controls verbosity.
    ///
    /// @returns A list with `global_gene_indices`, `total_cells`, `total_genes`
    /// and `per_file` (a list of lists with `exp_id`, `cell_indices`
    /// (file-local, 0-indexed), `lib_size`, `nnz`).
    pub fn multi_tenx_h5_to_file(
        &mut self,
        file_tasks: List,
        universe_size: i32,
        qc_params: List,
        verbose: bool,
    ) -> Result<List, extendr_api::Error> {
        let qc = MinCellQuality::from_r_list(qc_params)?;

        let tasks: Vec<TenxFileTask> = file_tasks
            .into_iter()
            .map(|(_, robj)| {
                let inner = List::try_from(robj).expect("Each file_task must be a list");
                TenxFileTask::from_r_list(inner)
            })
            .collect::<Result<Vec<_>, _>>()?;

        let result = multi_10x_h5_to_file(
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
    /// @returns A list with `indptr`, `indices` (0-indexed), `data`,
    /// `no_cells` and `no_genes`, parseable into a sparse matrix in R. `data`
    /// is integer for `"raw"` and double for `"norm"`.
    pub fn return_full_mat(
        &self,
        assay: &str,
        cell_based: bool,
        verbose: bool,
    ) -> Result<List, extendr_api::Error> {
        let mut data: Vec<AssayData> = Vec::new();
        let mut indices: Vec<Vec<i32>> = Vec::new();
        let mut indptr: Vec<usize> = Vec::new();
        let assay_type = parse_count_type(assay).ok_or_else(|| {
            extendr_api::Error::Other(format!("Invalid assay '{assay}'. Use 'raw' or 'norm'."))
        })?;

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
    /// @returns A list with `indptr`, `indices` (0-indexed gene positions),
    /// `data`, `no_cells` (number of returned cells) and `no_genes`,
    /// parseable into a CSR matrix in R.
    pub fn get_cells_by_indices(
        &self,
        indices: &[i32],
        assay: &str,
    ) -> Result<List, extendr_api::Error> {
        let reader = ParallelSparseReader::new(&self.f_path_cells).to_extendr()?;
        let assay_type = parse_count_type(assay).ok_or_else(|| {
            extendr_api::Error::Other(format!("Invalid assay '{assay}'. Use 'raw' or 'norm'."))
        })?;

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
    /// (CSC) representation to `f_path_genes`. A parallel counting-sort
    /// transpose: one pass counts the non-zeros per gene, then genes are
    /// converted in phases that fit `max_mem_gb`. Each phase re-reads the
    /// cell file, so a tighter budget trades time for memory.
    ///
    /// @param max_mem_gb (`numeric` or `NULL`)\cr
    /// Memory for the conversion buffers in GB, at 10 bytes per non-zero.
    /// `NULL` converts in a single phase.
    /// @param verbose (`logical`)\cr
    /// Controls verbosity of the function.
    ///
    /// @returns Invisible `NULL`.
    pub fn generate_gene_based_data(
        &mut self,
        max_mem_gb: Option<f64>,
        verbose: bool,
    ) -> Result<(), extendr_api::Error> {
        let max_nnz = max_mem_gb.map(|gb| ((gb * 1e9) as usize / GENE_FILE_BYTES_PER_NNZ).max(1));

        write_gene_file(&self.f_path_cells, &self.f_path_genes, max_nnz, verbose).to_extendr()
    }

    /////////////
    // Archive //
    /////////////

    /// Archive the cell-based binary for cold storage
    ///
    /// @description
    /// Writes a zstd-compressed archive of `f_path_cells`. The gene-based
    /// file is not archived; `restore_archive()` rebuilds it. Normalised
    /// values are only stored for cells where they cannot be recomputed from
    /// the raw counts.
    ///
    /// @param f_path_archive (`character`)\cr
    /// Path of the archive to write.
    /// @param level (`integer`)\cr
    /// zstd compression level, 1 to 22.
    /// @param verbose (`logical`)\cr
    /// Controls verbosity of the function.
    ///
    /// @returns A list with `n_cells`, `nnz`, `n_norm_stored` and
    /// `archive_bytes`.
    pub fn archive(
        &self,
        f_path_archive: &str,
        level: i32,
        verbose: bool,
    ) -> Result<List, extendr_api::Error> {
        let stats =
            archive_cell_file(&self.f_path_cells, f_path_archive, level, verbose).to_extendr()?;

        Ok(list!(
            n_cells = stats.n_cells,
            nnz = stats.nnz as f64,
            n_norm_stored = stats.n_norm_stored,
            archive_bytes = stats.archive_bytes as f64
        ))
    }

    /// Restore both binaries from an archive
    ///
    /// @description
    /// Rebuilds `f_path_cells` from the archive, then generates
    /// `f_path_genes` from it.
    ///
    /// @param f_path_archive (`character`)\cr
    /// Path to the archive.
    /// @param max_mem_gb (`numeric` or `NULL`)\cr
    /// Memory for the gene file conversion buffers in GB. `NULL` converts in
    /// a single phase.
    /// @param verbose (`logical`)\cr
    /// Controls verbosity of the function.
    ///
    /// @returns Invisible `NULL`.
    pub fn restore_archive(
        &mut self,
        f_path_archive: &str,
        max_mem_gb: Option<f64>,
        verbose: bool,
    ) -> Result<(), extendr_api::Error> {
        let max_nnz = max_mem_gb.map(|gb| ((gb * 1e9) as usize / GENE_FILE_BYTES_PER_NNZ).max(1));

        restore_archive(
            f_path_archive,
            &self.f_path_cells,
            &self.f_path_genes,
            max_nnz,
            verbose,
        )
        .to_extendr()
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
    /// @returns A list with `indptr`, `indices` (0-indexed cell positions),
    /// `data`, `no_cells` and `no_genes` (number of returned genes),
    /// parseable into a CSC matrix in R.
    pub fn get_genes_by_indices(
        &self,
        indices: &[i32],
        assay: &str,
    ) -> Result<List, extendr_api::Error> {
        let reader = ParallelSparseReader::new(&self.f_path_genes).to_extendr()?;

        let assay_type = parse_count_type(assay).ok_or_else(|| {
            extendr_api::Error::Other(format!("Invalid assay '{assay}'. Use 'raw' or 'norm'."))
        })?;

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
    /// @returns An integer vector of NNZ counts for the requested genes.
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
    /// `gene_local_to_universe` (0-indexed integer vector, `NA` for genes
    /// absent from the universe).
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
    /// @returns A list with `total_cells`, `total_genes` and `per_file` (a
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

    ///////////////
    // CellSweep //
    ///////////////

    /// Run CellSweep and write the denoised barcodes into the cells binary
    ///
    /// One independent EM fit per sample, since the ambient profile is a
    /// property of a single emulsion. Only the real barcodes are written: the
    /// empty droplets exist to train the ambient profile, and barcodes that
    /// are neither empty nor annotated are not part of the model.
    ///
    /// Writes the cell-based file only; regenerate the gene-based companion
    /// afterwards, as for a merge.
    ///
    /// @param f_path_source (`character`)\cr
    /// Path to the raw `counts_cells.bin`, which must still contain the empty
    /// droplets.
    /// @param samples (`list`)\cr
    /// A list of lists. Each inner list must contain `sample_id`,
    /// `real_cells` and `empty_cells` (0-indexed integer vectors of store
    /// indices), `celltype_idx` (0-indexed integer vector, one entry per
    /// `real_cells` entry) and `n_celltypes`.
    /// @param cellsweep_params (`list`)\cr
    /// The CellSweep model parameters. Missing entries fall back to the
    /// reference implementation's defaults.
    /// @param target_size (`numeric`)\cr
    /// Library size the normalised layer is scaled to.
    /// @param verbose (`integer`)\cr
    /// `0` silent, `1` per-sample progress, `2` per-EM-iteration.
    ///
    /// @returns A list with `cell_order` (0-indexed source indices in output
    /// order), `lib_size`, `nnz` and `fits` (one list per sample with
    /// `sample_id`, `alpha`, `z_hat` (1-indexed), `beta`, `ambient`,
    /// `celltype_profiles`, `n_celltypes`, `log_likelihood`, `n_iter` and
    /// `converged`).
    pub fn cellsweep(
        &mut self,
        f_path_source: String,
        samples: List,
        cellsweep_params: List,
        target_size: f64,
        verbose: usize,
    ) -> Result<List, extendr_api::Error> {
        let sample_specs: Vec<CellSweepSample> = samples
            .into_iter()
            .map(|(_, robj)| {
                let inner = List::try_from(robj).expect("Each sample must be a list");
                cellsweep_sample_from_r_list(inner)
            })
            .collect::<Result<Vec<_>, _>>()?;

        let params = CellSweepParams::from_r_list(cellsweep_params)?;

        let reader = ParallelSparseReader::new(&f_path_source).to_extendr()?;
        let result = run_cellsweep(
            &reader,
            &sample_specs,
            params,
            &self.f_path_cells,
            target_size as f32,
            verbose,
        )
        .to_extendr()?;

        self.n_cells = result.cell_order.len();
        self.n_genes = reader.get_header().total_genes;

        let fits: List = result
            .fits
            .iter()
            .zip(&sample_specs)
            .map(|(fit, sample)| fit.to_r_list(&sample.sample_id))
            .collect::<List>();

        Ok(list!(
            cell_order = result.cell_order.r_int_convert(),
            lib_size = result.library_size.r_int_convert(),
            nnz = result.nnz.r_int_convert(),
            fits = fits
        ))
    }
}
