use bixverse_rs::prelude::*;
use bixverse_rs::single_cell::mc_analysis::metacell_density::*;
use bixverse_rs::single_cell::mc_generation::cell_aggregation_utils::*;
use bixverse_rs::single_cell::mc_generation::hdwgcna_meta_cells::*;
use bixverse_rs::single_cell::mc_generation::seacells::*;
use bixverse_rs::single_cell::mc_generation::super_cells::*;
use bixverse_rs::single_cell::sc_r_wrappers::{assignments_to_r_list, metacells_to_r_list};
use extendr_api::*;
use faer::Mat;
use rustc_hash::FxHashMap;
use std::time::Instant;

use crate::single_cell::utils::{knn_data_to_rust, knn_indices_processing};

///////////////////
// Index mapping //
///////////////////

/// Resolved mapping between working space and the original count file.
///
/// Two index spaces meet in the meta cell entry points and conflating them
/// silently aggregates the wrong cells. *Working space* is the row space of the
/// embedding or kNN graph handed over from R, holding one row per cell that
/// passed QC. *Original space* is the row space of `counts_cells.bin`, which is
/// what [`aggregate_meta_cells`] indexes into and what R reports back as
/// `original_cell_idx`. The two coincide only when every cell in the file
/// passed QC, so memberships are always translated through this mapping before
/// they are aggregated or returned.
pub struct CellMapping {
    /// Working rows to retain, in output order. `None` means every working row
    /// is used and the embedding needs no subsetting.
    pub rows_to_use: Option<Vec<usize>>,
    /// Maps a working row onto its row in the original count file.
    pub subset_to_orig: Vec<usize>,
    /// Number of cells in the original count file, i.e. the size of the space
    /// `subset_to_orig` points into.
    pub n_total: usize,
}

impl CellMapping {
    /// Whether the incoming embedding or kNN graph has to be subset before use.
    ///
    /// ### Returns
    ///
    /// `true` when a narrowed set of cells was requested, which also implies the
    /// kNN graph has to be rebuilt on those cells.
    pub fn needs_subsetting(&self) -> bool {
        self.rows_to_use.is_some()
    }
}

/// Work out how working-space rows map back onto the original count file.
///
/// Resolves the identity case (no QC filtering), the QC-filtered case, and the
/// narrowed case where only a subset of the QC-passing cells is requested.
/// Requested cells that are absent from `cells_to_keep` are dropped silently,
/// since the caller may legitimately ask for cells that failed QC.
///
/// ### Params
///
/// * `cells_to_keep` - Original row indices the embedding was built from, in
///   embedding row order, 0-indexed. `None` assumes the embedding covers the
///   count file one-to-one.
/// * `cells_to_use` - Optional further narrowing, given in original space and
///   0-indexed. `None` keeps every working row.
/// * `n_working_rows` - Number of embedding or kNN rows supplied by R.
/// * `n_total` - Number of cells in the count file.
///
/// ### Returns
///
/// A [`CellMapping`], or an error when `cells_to_keep` disagrees with
/// `n_working_rows` or when none of the requested cells survive the mapping.
fn resolve_cell_mapping(
    cells_to_keep: Option<&[i32]>,
    cells_to_use: Option<&[i32]>,
    n_working_rows: usize,
    n_total: usize,
) -> extendr_api::Result<CellMapping> {
    // Working row i corresponds to cells_to_keep[i]. Without it we can only
    // assume the embedding covers the file one-to-one.
    let keep: Vec<usize> = match cells_to_keep {
        Some(keep) => {
            if keep.len() != n_working_rows {
                return Err(format!(
                    "'cells_to_keep' has {} entries but the embedding/kNN has {} rows",
                    keep.len(),
                    n_working_rows
                )
                .into());
            }
            keep.iter().map(|&x| x as usize).collect()
        }
        None => (0..n_working_rows).collect(),
    };

    match cells_to_use {
        Some(use_cells) => {
            let orig_to_working: FxHashMap<usize, usize> = keep
                .iter()
                .enumerate()
                .map(|(working_row, &orig_idx)| (orig_idx, working_row))
                .collect();

            let mut rows_to_use = Vec::with_capacity(use_cells.len());
            let mut subset_to_orig = Vec::with_capacity(use_cells.len());

            for &orig_idx in use_cells {
                let orig_idx = orig_idx as usize;
                if let Some(&working_row) = orig_to_working.get(&orig_idx) {
                    rows_to_use.push(working_row);
                    subset_to_orig.push(orig_idx);
                }
            }

            if subset_to_orig.is_empty() {
                return Err(
                    "None of the requested 'cells_to_use' are part of 'cells_to_keep'".into(),
                );
            }

            Ok(CellMapping {
                rows_to_use: Some(rows_to_use),
                subset_to_orig,
                n_total,
            })
        }
        // Still a remap whenever cells_to_keep is not the identity.
        None => Ok(CellMapping {
            rows_to_use: None,
            subset_to_orig: keep,
            n_total,
        }),
    }
}

/// Gather selected rows of an R embedding matrix into a faer matrix.
///
/// R matrices are column-major, so the source element for row `i`, column `j`
/// sits at `rows_to_use[i] + j * nrow`. Downcast to `f32` matches the precision
/// the meta cell algorithms work in; the embedding is a PCA projection, so the
/// dynamic range is small enough for this to be safe.
///
/// ### Params
///
/// * `embd` - Column-major embedding matrix as handed over from R, with one row
///   per working-space cell.
/// * `rows_to_use` - Working row indices to gather, in output order.
///
/// ### Returns
///
/// A `rows_to_use.len()` by `embd.ncols()` matrix holding the selected rows.
fn subset_embedding(embd: &RMatrix<f64>, rows_to_use: &[usize]) -> Mat<f32> {
    let ncol = embd.ncols();
    let nrow = embd.nrows();
    let data = embd.data();

    Mat::from_fn(rows_to_use.len(), ncol, |i, j| {
        data[rows_to_use[i] + j * nrow] as f32
    })
}

/// Read the number of cells recorded in a `counts_cells.bin` header.
///
/// Opens the file purely to size the original index space. The aggregation step
/// opens its own reader later, so this handle is dropped immediately.
///
/// ### Params
///
/// * `f_path` - Path to the `counts_cells.bin` file.
///
/// ### Returns
///
/// The total number of cells stored in the file, or an error when the file
/// cannot be opened.
fn count_file_n_cells(f_path: &str) -> extendr_api::Result<usize> {
    let reader = ParallelSparseReader::new(f_path)
        .map_err(|e| Error::Other(format!("Could not open '{}': {:?}", f_path, e)))?;
    Ok(reader.get_header().total_cells)
}

/////////////
// extendR //
/////////////

extendr_module! {
    mod r_sc_metacells;
    // meta cells
    fn rs_get_metacells_bootstrapped;
    fn rs_get_seacells;
    fn rs_supercell;
    // meta cell metrics
    fn rs_metacell_density;
    fn rs_metacell_compactness;
    fn rs_metacell_separation;
    // pseudo-bulking
    fn rs_pseudobulk_cells_dense;
    fn rs_pseudobulk_cells_sparse;
}

////////////////
// Meta cells //
////////////////

//////////////////
// Bootstrapped //
//////////////////

/// Generate meta cells (hdWGCNA method)
///
/// @description
/// `r lifecycle::badge("experimental")`
/// This function implements the approach from Morabito, et al. to generate meta
/// cells. You can provide an already pre-computed kNN matrix or an embedding to
/// regenerate the kNN matrix with specified parameters in the meta_cell_params.
/// If `knn_mat` is provided, this one will be used. You need to at least
/// provide `knn_mat` or `embd`!
///
/// @param f_path String. Path to the `counts_cells.bin` file.
/// @param knn_mat Optional integer matrix. The kNN matrix you wish to use
/// for the generation of the meta cells. This function expects 0-indices!
/// @param embd Optional numerical matrix. The embedding matrix (for example
/// PCA embedding) you wish to use for the generation of the kNN graph that
/// is used subsequently for aggregation of the meta cells.
/// @param cells_to_keep Optional indices of the cells to keep, i.e., the
/// cells used for the generation of the embedding.
/// @param cells_to_use Optional indices of cells to use for meta cell
/// generation. Useful if you wish to generate meta cells in specific cell
/// types. If this is provided, the kNN graph will be regenerated.
/// @param meta_cell_params A list containing the meta cell parameters.
/// @param target_size Numeric. Target library size for re-normalisation of
/// the meta cells. Typically `1e4`.
/// @param verbose Integer. `0L` - quiet; `1L` - normal verbosity; `2L` -
/// detailed verbosity.
/// @param verbose Boolean. Controls verbosity of the function.
///
/// @returns A list with the following elements:
/// \itemize{
///  \item assignments - A list containing assignment information with elements:
///    assignments (vector), metacells (list), unassigned (vector), n_metacells,
///    n_cells, n_unassigned
///  \item aggregated - A list with indptr, indices, raw_counts, norm_counts,
///    nrow, ncol in sparse format.
/// }
///
/// @export
///
/// @keywords internal
#[extendr]
#[allow(clippy::too_many_arguments)]
fn rs_get_metacells_bootstrapped(
    f_path: String,
    knn_mat: Option<RMatrix<i32>>,
    embd: Option<RMatrix<f64>>,
    cells_to_keep: Option<Vec<i32>>,
    cells_to_use: Option<Vec<i32>>,
    meta_cell_params: List,
    target_size: f64,
    seed: usize,
    verbose: usize,
) -> extendr_api::Result<List> {
    let meta_cell_params = BootstrappedMetaCellParams::from_r_list(meta_cell_params)?;

    let verbosity = parse_verbosity_level(verbose);

    if cells_to_use.is_some() && (cells_to_keep.is_none() || embd.is_none()) {
        return Err(
            "When using 'cells_to_use', both 'cells_to_keep' and 'embd' must be provided".into(),
        );
    }

    let n_working_rows = match (&knn_mat, &embd) {
        (Some(mat), _) => mat.nrows(),
        (_, Some(em)) => em.nrows(),
        _ => return Err("Must provide either 'knn_mat' or 'embd' parameter".into()),
    };

    let mapping = resolve_cell_mapping(
        cells_to_keep.as_deref(),
        cells_to_use.as_deref(),
        n_working_rows,
        count_file_n_cells(&f_path)?,
    )?;

    let knn_graph = match &mapping.rows_to_use {
        // Narrowed to a set of cells, so the kNN has to be rebuilt on them.
        Some(rows_to_use) => {
            let embd = embd.as_ref().ok_or_else(|| {
                Error::Other("When using 'cells_to_use', 'embd' must be provided".into())
            })?;

            if verbosity.normal_verbosity() {
                println!(
                    "Subsetting to {} cells (from {} QC-passing cells) and regenerating kNN graph",
                    rows_to_use.len(),
                    n_working_rows
                );
            }

            let embd_subset = subset_embedding(embd, rows_to_use);

            let (knn, _) = generate_knn_with_dist(
                embd_subset.as_ref(),
                &meta_cell_params.knn_params,
                false,
                false,
                seed,
                verbosity.detailed_verbosity(),
            )
            .to_extendr()?;

            knn
        }
        None => match (knn_mat, embd) {
            (Some(knn_mat), _) => {
                if verbosity.normal_verbosity() {
                    println!("Using provided kNN matrix");
                }
                knn_indices_processing(knn_mat)
            }
            (None, Some(embd)) => {
                if verbosity.normal_verbosity() {
                    println!("Calculating the kNN matrix from the provided data.");
                }

                let embd = r_matrix_to_faer_fp32(&embd);

                let (knn, _) = generate_knn_with_dist(
                    embd.as_ref(),
                    &meta_cell_params.knn_params,
                    false,
                    false,
                    seed,
                    verbosity.detailed_verbosity(),
                )
                .to_extendr()?;

                knn
            }
            (None, None) => {
                return Err("Must provide either 'knn_mat' or 'embd' parameter".into());
            }
        },
    };

    let subset_to_orig = mapping.subset_to_orig;
    let n_total_cells = mapping.n_total;

    let nn_map = build_nn_map(&knn_graph);

    if verbosity.normal_verbosity() {
        println!("Identifying meta cells.");
    }

    let centres = identify_meta_cells(
        &nn_map,
        meta_cell_params.max_shared,
        meta_cell_params.target_no_metacells,
        meta_cell_params.max_iter,
        seed,
        verbose,
    );

    let metacells_subset: Vec<Vec<usize>> = assign_bootstrapped_meta_cells(&centres, &nn_map);

    // Always remap: the aggregation reader and R both work in original space.
    let metacells_original: Vec<Vec<usize>> = remap_metacells_to_original(
        &metacells_subset
            .iter()
            .map(|v| v.as_slice())
            .collect::<Vec<_>>(),
        &subset_to_orig,
    );

    let assignment_list: List = metacells_to_r_list(&metacells_original, n_total_cells);

    if verbosity.normal_verbosity() {
        println!("Aggregating meta cells.");
    }

    let reader = ParallelSparseReader::new(&f_path).unwrap();
    let n_genes = reader.get_header().total_genes;

    let metacells_refs: Vec<&[usize]> = metacells_original.iter().map(|v| v.as_slice()).collect();

    let aggregated: CompressedSparseData2<u32, f32> =
        aggregate_meta_cells(&reader, &metacells_refs, target_size as f32, n_genes).to_extendr()?;

    Ok(list!(
        assignments = assignment_list,
        aggregated = list!(
            indptr = aggregated.indptr.r_int_convert(),
            indices = aggregated.indices.r_int_convert(),
            raw_counts = aggregated.data.r_int_convert(),
            norm_counts = aggregated.data_2.unwrap().r_float_convert(),
            nrow = aggregated.shape.0,
            ncol = aggregated.shape.1
        )
    ))
}

//////////////
// SEACells //
//////////////

/// Generate SEACells
///
/// @description
/// `r lifecycle::badge("experimental")`
/// This function implements the SEACells algorithm for generating meta cells
/// from Persad et al. An embedding matrix must be provided which is used to
/// construct the kNN graph and kernel matrix for the SEACells algorithm. This
/// version is highly memory and speed-optimised and will truncate small values
/// during matrix operations which can affect convergence.
///
/// @param f_path String. Path to the `counts_cells.bin` file.
/// @param embd Numerical matrix. The embedding matrix (for example PCA embedding)
/// used for the generation of the kNN graph and kernel matrix.
/// @param cells_to_keep Optional indices of the cells to keep, i.e., the
/// cells used for the generation of the embedding.
/// @param cells_to_use Optional indices of cells to use for meta cell
/// generation. Useful if you wish to generate meta cells in specific cell
/// types.
/// @param knn_data Optional list. This contains pre-computed kNN data
/// (including distances). The user has to ensure consistency!
/// @param seacells_params A list containing the SEACells parameters.
/// @param target_size Numeric. Target library size for re-normalisation of
/// the meta cells. Typically `1e4`.
/// @param seed Integer. For reproducibility purposes.
/// @param verbose Integer. `0L` - quiet; `1L` - normal verbosity; `2L` -
/// detailed verbosity.
///
/// @returns A list with the following elements:
/// \itemize{
///  \item assignments - A list containing assignment information with elements:
///    assignments (vector), metacells (list), unassigned (vector), n_metacells,
///    n_cells, n_unassigned
///  \item aggregated - A list with indptr, indices, raw_counts, norm_counts,
///    nrow, ncol in sparse format.
///  \item rss - Vector of RSS values from each iteration.
///  \item archetypes - Vector of cell indices selected as archetypes.
/// }
///
/// @export
///
/// @references Persad, et al., Nat. Biotechnol., 2023.
///
/// @keywords internal
#[extendr]
#[allow(clippy::too_many_arguments)]
fn rs_get_seacells(
    f_path: String,
    embd: RMatrix<f64>,
    cells_to_keep: Option<Vec<i32>>,
    cells_to_use: Option<Vec<i32>>,
    knn_data: Nullable<List>,
    seacells_params: List,
    target_size: f64,
    seed: usize,
    verbose: usize,
) -> extendr_api::Result<List> {
    let start_seacell = Instant::now();

    let seacells_params = SEACellsParams::from_r_list(seacells_params)?;
    let verbosity = parse_verbosity_level(verbose);
    let knn_provided = knn_data != extendr_api::Nullable::Null;

    if cells_to_use.is_some() && knn_provided {
        println!(
            "[WARNING!] 'knn_data' is ignored when 'cells_to_use' is set; the kNN graph will be regenerated on the subset"
        );
    }

    let mapping = resolve_cell_mapping(
        cells_to_keep.as_deref(),
        cells_to_use.as_deref(),
        embd.nrows(),
        count_file_n_cells(&f_path)?,
    )?;

    let embd_mat = match &mapping.rows_to_use {
        Some(rows_to_use) => {
            if verbosity.normal_verbosity() {
                println!(
                    "Subsetting to {} cells (from {} QC-passing cells)",
                    rows_to_use.len(),
                    embd.nrows()
                );
            }
            subset_embedding(&embd, rows_to_use)
        }
        None => r_matrix_to_faer_fp32(&embd),
    };

    let is_subset = mapping.needs_subsetting();
    let subset_to_orig = mapping.subset_to_orig;
    let n_total_cells = mapping.n_total;

    let (knn_indices, knn_dist, dist_squared) = if knn_provided && !is_subset {
        if verbosity.normal_verbosity() {
            println!("Using provided kNN graph.")
        }

        let knn_data = knn_data
            .into_robj()
            .as_list()
            .ok_or_else(|| Error::Other("'knn_data' is not a list".into()))?;
        let (knn_indices, knn_dist, _, distance) = knn_data_to_rust(knn_data)?;

        if knn_indices.len() != embd_mat.nrows() {
            return Err(format!(
                "kNN indices have {} rows but embedding has {}",
                knn_indices.len(),
                embd_mat.nrows()
            )
            .into());
        }

        let dist_squared = distance == "euclidean";
        (knn_indices, knn_dist, dist_squared)
    } else {
        let start_knn = Instant::now();

        if verbosity.normal_verbosity() {
            println!("Regenerating kNN graph.")
        }

        let (knn_indices, knn_dist) = generate_knn_with_dist(
            embd_mat.as_ref(),
            &seacells_params.knn_params,
            true,
            false,
            seed,
            verbosity.detailed_verbosity(),
        )
        .to_extendr()?;
        let knn_dist = knn_dist.unwrap();
        let dist_squared = seacells_params.knn_params.ann_dist == "euclidean";

        if verbosity.normal_verbosity() {
            println!(
                "kNN generation done in : {:.2?} with {}",
                start_knn.elapsed(),
                seacells_params.knn_params.knn_method
            );
        }

        (knn_indices, knn_dist, dist_squared)
    };

    let mut seacell = SEACells::new(embd_mat.nrows(), &seacells_params);

    seacell.construct_kernel_mat(embd_mat.as_ref(), &knn_indices, &knn_dist, verbose);

    if let Some(n_landmarks) = seacells_params.n_landmarks {
        seacell
            .initialise_archetypes_landmark(
                embd_mat.as_ref(),
                &knn_indices,
                &knn_dist,
                dist_squared,
                n_landmarks,
                verbose,
                seed as u64,
            )
            .to_extendr()?;
    } else {
        seacell
            .initialise_archetypes(&knn_indices, &knn_dist, verbose, dist_squared, seed as u64)
            .to_extendr()?;
    }

    seacell.fit(seed, verbose).to_extendr()?;

    let assignments_raw = seacell.get_hard_assignments().to_extendr()?;
    let archetypes_raw = seacell.get_archetypes().to_extendr()?;
    let k = seacells_params.n_sea_cells;

    let rss = seacell.get_rss_history();

    let groups_raw = assignments_to_metacells(&assignments_raw, k);

    let mut id_remap: Vec<Option<usize>> = vec![None; k];
    let mut groups_kept: Vec<Vec<usize>> = Vec::new();
    let mut archetypes_kept: Vec<usize> = Vec::new();

    for (old_id, group) in groups_raw.into_iter().enumerate() {
        if !group.is_empty() {
            id_remap[old_id] = Some(groups_kept.len());
            archetypes_kept.push(archetypes_raw[old_id]);
            groups_kept.push(group);
        }
    }

    if verbosity.normal_verbosity() && groups_kept.len() < k {
        println!(
            "Dropped {} empty archetype(s); keeping {} of {} requested",
            k - groups_kept.len(),
            groups_kept.len(),
            k
        );
    }

    let assignments_subset_opt: Vec<Option<usize>> =
        assignments_raw.iter().map(|&old| id_remap[old]).collect();

    // Always remap: the aggregation reader and R both work in original space.
    let assignments_full: Vec<Option<usize>> =
        remap_assignments_to_original(&assignments_subset_opt, &subset_to_orig, n_total_cells);

    let archetypes_original: Vec<usize> = archetypes_kept
        .iter()
        .map(|&idx| subset_to_orig[idx])
        .collect();

    let assignment_list = assignments_to_r_list(&assignments_full, n_total_cells);

    if verbosity.normal_verbosity() {
        println!("Aggregating meta cells.");
    }

    let metacells_original: Vec<Vec<usize>> = remap_metacells_to_original(
        &groups_kept.iter().map(|v| v.as_slice()).collect::<Vec<_>>(),
        &subset_to_orig,
    );

    let metacells_refs: Vec<&[usize]> = metacells_original.iter().map(|v| v.as_slice()).collect();

    let reader = ParallelSparseReader::new(&f_path).unwrap();
    let n_genes = reader.get_header().total_genes;

    let aggregated: CompressedSparseData2<u32, f32> =
        aggregate_meta_cells(&reader, &metacells_refs, target_size as f32, n_genes).to_extendr()?;

    if verbosity.normal_verbosity() {
        println!("SEACells found in : {:.2?}", start_seacell.elapsed());
    }

    Ok(list!(
        assignments = assignment_list,
        aggregated = list!(
            indptr = aggregated.indptr.r_int_convert(),
            indices = aggregated.indices.r_int_convert(),
            raw_counts = aggregated.data.r_int_convert(),
            norm_counts = aggregated.data_2.unwrap().r_float_convert(),
            nrow = aggregated.shape.0,
            ncol = aggregated.shape.1
        ),
        rss = rss.to_vec(),
        archetypes = archetypes_original.r_int_convert()
    ))
}

////////////////
// Supercells //
////////////////

/// Generate SuperCells.
///
/// @description
/// `r lifecycle::badge("experimental")`
/// This function implements the approach from Bilous, et al.
/// to generate meta cells or called here SuperCells. You can provide
/// pre-computed kNN data (indices + distances) via `knn_data`, or an
/// embedding via `embd` from which the kNN graph will be generated. You
/// need to at least provide `knn_data` or `embd`. When `cells_to_use` is
/// supplied, the kNN graph is always regenerated on the subset and any
/// `knn_data` is ignored. Distances are required when the SuperCell
/// parameters request the kernel-weighted graph.
///
/// @param f_path String. Path to the `counts_cells.bin` file.
/// @param embd Optional numerical matrix. The embedding matrix (for example
/// PCA embedding) used for the generation of the kNN graph. Required when
/// `knn_data` is not provided, and required when using `cells_to_use`.
/// @param cells_to_keep Optional indices of the cells to keep, i.e., the
/// cells used for the generation of the embedding.
/// @param cells_to_use Optional indices of cells to use for meta cell
/// generation. Useful if you wish to generate meta cells in specific cell
/// types. If this is provided, `embd` and `cells_to_keep` are required and
/// the kNN graph will be regenerated on the subset.
/// @param knn_data Optional list. This contains pre-computed kNN data
/// (including distances). The user has to ensure consistency! Ignored when
/// `cells_to_use` is set.
/// @param supercell_params A list containing the SuperCell parameters.
/// @param target_size Numeric. Target library size for re-normalisation of
/// the meta cells. Typically `1e4`.
/// @param seed Integer. For reproducibility purposes.
/// @param verbose Integer. `0L` - quiet; `1L` - normal verbosity; `2L` -
/// detailed verbosity.
///
/// @returns A list with the following elements:
/// \itemize{
///  \item assignments - A list containing assignment information with elements:
///    assignments (vector), metacells (list), unassigned (vector), n_metacells,
///    n_cells, n_unassigned
///  \item aggregated - A list with indptr, indices, raw_counts, norm_counts,
///    nrow, ncol in sparse format.
/// }
///
/// @export
///
/// @keywords internal
#[extendr]
#[allow(clippy::too_many_arguments)]
fn rs_supercell(
    f_path: String,
    embd: Option<RMatrix<f64>>,
    cells_to_keep: Option<Vec<i32>>,
    cells_to_use: Option<Vec<i32>>,
    knn_data: Nullable<List>,
    supercell_params: List,
    target_size: f64,
    seed: usize,
    verbose: usize,
) -> extendr_api::Result<List> {
    let supercell_params = SuperCellParams::from_r_list(supercell_params)?;
    let verbosity = parse_verbosity_level(verbose);

    if cells_to_use.is_some() && (cells_to_keep.is_none() || embd.is_none()) {
        return Err(
            "When using 'cells_to_use', both 'cells_to_keep' and 'embd' must be provided".into(),
        );
    }

    let knn_provided = knn_data != extendr_api::Nullable::Null;
    let is_subset = cells_to_use.is_some();

    if is_subset && knn_provided {
        println!(
            "[WARNING!] 'knn_data' is ignored when 'cells_to_use' is set; the kNN graph will be regenerated on the subset"
        );
    }

    if !knn_provided && embd.is_none() {
        return Err("Must provide either 'knn_data' or 'embd'".into());
    }

    // kNN rows define the working space when no embedding was handed over.
    let n_working_rows = match (&knn_data, &embd) {
        (_, Some(em)) => em.nrows(),
        _ => {
            let knn_list = knn_data
                .clone()
                .into_robj()
                .as_list()
                .ok_or_else(|| Error::Other("'knn_data' is not a list".into()))?;
            knn_data_to_rust(knn_list)?.0.len()
        }
    };

    let mapping = resolve_cell_mapping(
        cells_to_keep.as_deref(),
        cells_to_use.as_deref(),
        n_working_rows,
        count_file_n_cells(&f_path)?,
    )?;

    let (knn_indices, knn_dist, dist_squared) = match &mapping.rows_to_use {
        Some(rows_to_use) => {
            let embd = embd.as_ref().ok_or_else(|| {
                Error::Other("When using 'cells_to_use', 'embd' must be provided".into())
            })?;

            if verbosity.normal_verbosity() {
                println!(
                    "Subsetting to {} cells and regenerating kNN graph",
                    rows_to_use.len()
                );
            }

            let embd_subset = subset_embedding(embd, rows_to_use);

            let start_knn = Instant::now();
            let (knn_idx, knn_d) = generate_knn_with_dist(
                embd_subset.as_ref(),
                &supercell_params.knn_params,
                true,
                false,
                seed,
                verbosity.detailed_verbosity(),
            )
            .to_extendr()?;
            let knn_d = knn_d.unwrap();
            let dist_sq = supercell_params.knn_params.ann_dist == "euclidean";

            if verbosity.normal_verbosity() {
                println!(
                    "kNN generation done in: {:.2?} with {}",
                    start_knn.elapsed(),
                    supercell_params.knn_params.knn_method
                );
            }

            (knn_idx, knn_d, dist_sq)
        }
        None => {
            if knn_provided {
                let knn_list = knn_data
                    .into_robj()
                    .as_list()
                    .ok_or_else(|| Error::Other("'knn_data' is not a list".into()))?;
                let (knn_idx, knn_d, _, distance) = knn_data_to_rust(knn_list)?;

                if let Some(ref em) = embd {
                    if knn_idx.len() != em.nrows() {
                        return Err(format!(
                            "kNN indices have {} rows but embedding has {}",
                            knn_idx.len(),
                            em.nrows()
                        )
                        .into());
                    }
                }

                (knn_idx, knn_d, distance == "euclidean")
            } else {
                let embd = embd.unwrap();
                let embd_mat = r_matrix_to_faer_fp32(&embd);

                if verbosity.normal_verbosity() {
                    println!("Calculating kNN matrix from provided data");
                }

                let start_knn = Instant::now();
                let (knn_idx, knn_d) = generate_knn_with_dist(
                    embd_mat.as_ref(),
                    &supercell_params.knn_params,
                    true,
                    false,
                    seed,
                    verbosity.detailed_verbosity(),
                )
                .to_extendr()?;
                let knn_d = knn_d.unwrap();
                let dist_sq = supercell_params.knn_params.ann_dist == "euclidean";

                if verbosity.normal_verbosity() {
                    println!(
                        "kNN generation done in: {:.2?} with {}",
                        start_knn.elapsed(),
                        supercell_params.knn_params.knn_method
                    );
                }

                (knn_idx, knn_d, dist_sq)
            }
        }
    };

    let subset_to_orig = mapping.subset_to_orig;
    let n_total_cells = mapping.n_total;

    let n_meta_cells =
        (knn_indices.len() as f64 / supercell_params.graining_factor).ceil() as usize;

    if verbosity.normal_verbosity() {
        println!("Running SuperCell with Walktrap");
    }

    let membership_subset = supercell(
        &knn_indices,
        &knn_dist,
        &supercell_params,
        dist_squared,
        n_meta_cells,
        verbose,
    );

    let assignments_subset_opt: Vec<Option<usize>> =
        membership_subset.iter().map(|&x| Some(x)).collect();

    // Always remap: the aggregation reader and R both work in original space.
    let assignments_full =
        remap_assignments_to_original(&assignments_subset_opt, &subset_to_orig, n_total_cells);

    let assignment_list = assignments_to_r_list(&assignments_full, n_total_cells);

    if verbosity.normal_verbosity() {
        println!("Aggregating metacells");
    }

    let meta_cell_indices = assignments_to_metacells(&membership_subset, n_meta_cells);

    let metacells_original: Vec<Vec<usize>> = remap_metacells_to_original(
        &meta_cell_indices
            .iter()
            .map(|v| v.as_slice())
            .collect::<Vec<_>>(),
        &subset_to_orig,
    );

    let metacells_refs: Vec<&[usize]> = metacells_original.iter().map(|v| v.as_slice()).collect();

    let reader = ParallelSparseReader::new(&f_path).unwrap();
    let n_genes = reader.get_header().total_genes;

    let aggregated: CompressedSparseData2<u32, f32> =
        aggregate_meta_cells(&reader, &metacells_refs, target_size as f32, n_genes).to_extendr()?;

    Ok(list!(
        assignments = assignment_list,
        aggregated = list!(
            indptr = aggregated.indptr.r_int_convert(),
            indices = aggregated.indices.r_int_convert(),
            raw_counts = aggregated.data.r_int_convert(),
            norm_counts = aggregated.data_2.unwrap().r_float_convert(),
            nrow = aggregated.shape.0,
            ncol = aggregated.shape.1
        )
    ))
}

////////////////
// MetaCells2 //
////////////////

//////////////////////
// MetaCell density //
//////////////////////

/// Calculates diffusion maps for density calculations for meta cells
///
/// @description
/// `r lifecycle::badge("experimental")`
/// Generates diffusion maps and identifies in which density region a given
/// cell sits (defined as distance to k-nearest neighbours quite).
///
/// @param knn_data Named list. Needs to have the relevant data from the kNN
/// graph.
/// @param n_dcs Integer. The number of diffusion coordinates to return.
/// Typically `10`.
/// @param k_density Integer. The k-nearest neighbour to use for the density
/// estimation. Typically `150`.
/// @param knn_params List. The kNN parameters defined by
/// [params_sc_neighbours()].
/// @param verbose Integer. `0L` - quiet; `1L` - normal verbosity; `2L` -
/// detailed verbosity.
/// @param seed Integer. For reproducibility.
///
/// @return A list with the following items
/// \itemize{
///   \item dcs - Density coordinates
///   \item density_distances - Density distances at `k_density` neighbours.
///   \item regions - Region of the manifold where this given cell is.
/// }
///
/// @export
///
/// @references Persad, et al., Nat. Biotechnol., 2023.
///
/// @keywords internal
#[extendr]
fn rs_metacell_density(
    knn_data: List,
    n_dcs: usize,
    k_density: usize,
    knn_params: List,
    verbose: usize,
    seed: usize,
) -> Result<List> {
    let (knn_indices, knn_distances, original_k, distance) = knn_data_to_rust(knn_data)?;
    let knn_params = KnnParams::from_r_list(knn_params)?;
    let squared_dist = distance == "euclidean";

    let density_res: DiffusionDensity = compute_diffusion_density(
        &knn_indices,
        &knn_distances,
        squared_dist,
        original_k,
        n_dcs,
        k_density,
        &knn_params,
        seed as u64,
        verbose,
    )
    .to_extendr()?;

    let density_to_string = |x: &DensityRegion| -> String {
        match x {
            DensityRegion::High => "high".to_string(),
            DensityRegion::Mid => "mid".to_string(),
            DensityRegion::Low => "low".to_string(),
        }
    };

    let dcs = faer_to_r_matrix(density_res.dcs.as_ref());
    let density_distances = density_res.density_distances.r_float_convert();
    let regions: Vec<String> = density_res.regions.iter().map(density_to_string).collect();

    Ok(list![
        dcs = dcs,
        density_distances = density_distances,
        regions = regions
    ])
}

/// Calculates the compactness of the MetaCells based on diffusion map
/// coordinates
///
/// @description
/// `r lifecycle::badge("experimental")`
/// Calculates the meta cell compactness based on the diffusion map coordinates.
///
/// @param dc Numerical matrix. The diffusion map coordinates.
/// @param meta_cells List. The cell indices of the meta cells.
///
/// @returns The compactness results
///
/// @export
///
/// @keywords internal
#[extendr]
fn rs_metacell_compactness(dc: RMatrix<f64>, meta_cells: List) -> Result<Vec<f64>> {
    let dc = r_matrix_to_faer_fp32(&dc);

    let mut meta_cell_indices: Vec<Vec<usize>> = Vec::with_capacity(meta_cells.len());

    for i in 0..meta_cells.len() {
        let indices = meta_cells.elt(i)?.as_integer_vector().ok_or_else(|| {
            Error::Other("Could not convert the meta cell indices to Rust usize".into())
        })?;

        meta_cell_indices.push(indices.r_int_convert_shift());
    }

    let compactness = compute_compactness(dc.as_ref(), &meta_cell_indices);

    Ok(compactness.r_float_convert())
}

/// Calculates the separation of the centroids of the MetaCells based on
/// diffusion map coordinates.
///
/// @description
/// `r lifecycle::badge("experimental")`
/// Calculates the separation of the single cells of a given meta cell based on
/// the diffusion map.
///
/// @param dc Numerical matrix. The diffusion map coordinates.
/// @param meta_cells List. The cell indices of the meta cells.
///
/// @returns The separation results
///
/// @export
///
/// @keywords internal
#[extendr]
fn rs_metacell_separation(dc: RMatrix<f64>, meta_cells: List) -> Result<Vec<f64>> {
    let dc = r_matrix_to_faer_fp32(&dc);

    let mut meta_cell_indices: Vec<Vec<usize>> = Vec::with_capacity(meta_cells.len());

    for i in 0..meta_cells.len() {
        let indices = meta_cells.elt(i)?.as_integer_vector().ok_or_else(|| {
            Error::Other("Could not convert the meta cell indices to Rust usize".into())
        })?;

        meta_cell_indices.push(indices.r_int_convert_shift());
    }

    let compactness = compute_separation(dc.as_ref(), &meta_cell_indices);

    Ok(compactness.r_float_convert())
}

////////////////////
// Pseudo bulking //
////////////////////

/// Pseudo-bulk a set of cells (dense)
///
/// @description
/// `r lifecycle::badge("experimental")`
/// This function will return a dense matrix of
/// `length(cell_indices_ls) x number of genes`. The function has the option
/// to return the sum of the sum of the raw counts or the average of the
/// normalised counts.
///
/// @param f_path String. Path to the `counts_cells.bin` file.
/// @param cell_indices_ls List. Must contains 0-indexed positions of the
/// cells to aggregate per element.
/// @param assay String. One of `c("raw", "norm")`. Which counts to normalise.
/// @param verbose Integer. `0L` - quiet; `1L` - normal verbosity; `2L` -
/// detailed verbosity.
///
/// @returns A dense matrix with the pseudo-bulked data.
///
/// @export
///
/// @keywords internal
#[extendr]
fn rs_pseudobulk_cells_dense(
    f_path: String,
    cell_indices_ls: List,
    assay: String,
    verbose: usize,
) -> extendr_api::Result<RArray<f64, 2>> {
    let bulk_type = parse_pseudo_bulk(&assay).unwrap_or_default();

    let mut cell_indices: Vec<Vec<usize>> = Vec::with_capacity(cell_indices_ls.len());

    for i in 0..cell_indices_ls.len() {
        let element = cell_indices_ls.elt(i)?;
        let vec_i = element.as_integer_vector().unwrap().r_int_convert();
        cell_indices.push(vec_i);
    }

    let data =
        get_pseudo_bulked_counts_dense(&f_path, &cell_indices, bulk_type, verbose).to_extendr()?;

    Ok(faer_to_r_matrix(data.as_ref()))
}

/// Pseudo-bulk a set of cells (sparse)
///
/// @description
/// `r lifecycle::badge("experimental")`
/// This function will return a sparse matrix of
/// `length(cell_indices_ls) x number of genes` (in list form in CSR).
/// The function has the option to return the sum of the sum of the raw counts
/// or the average of the normalised counts.
///
/// @param f_path String. Path to the `counts_cells.bin` file.
/// @param cell_indices_ls List. Must contains 0-indexed positions of the
/// cells to aggregate per element.
/// @param assay String. One of `c("raw", "norm")`. Which counts to normalise.
/// @param verbose Integer. `0L` - quiet; `1L` - normal verbosity; `2L` -
/// detailed verbosity.
///
/// @returns A list with the following elements (easy to convert into CSR in R)
/// \itemize{
///   \item indptr - The index pointers (representing cells)
///   \item indices - The indices (representing genes)
///   \item data - The pseudo-bulked data
///   \item nrow - Number of rows (i.e., pseudo-bulked samples)
///   \item ncol - Number of columns
/// }
///
/// @export
///
/// @keywords internal
#[extendr]
fn rs_pseudobulk_cells_sparse(
    f_path: String,
    cell_indices_ls: List,
    assay: String,
    verbose: usize,
) -> extendr_api::Result<List> {
    let bulk_type = parse_pseudo_bulk(&assay).unwrap_or_default();

    let mut cell_indices: Vec<Vec<usize>> = Vec::with_capacity(cell_indices_ls.len());

    for i in 0..cell_indices_ls.len() {
        let element = cell_indices_ls.elt(i).unwrap();
        let vec_i = element.as_integer_vector().unwrap().r_int_convert();
        cell_indices.push(vec_i);
    }

    let data: CompressedSparseData2<f64> =
        get_pseudo_bulked_counts_sparse(&f_path, &cell_indices, bulk_type, verbose).to_extendr()?;

    Ok(list!(
        indptr = data.indptr,
        indices = data.indices,
        data = data.data,
        nrow = data.shape.0,
        ncol = data.shape.1
    ))
}
