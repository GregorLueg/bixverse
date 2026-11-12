use bixverse_r::prelude::*;
use bixverse_r::sc_r_wrappers::{assignments_to_r_list, metacells_to_r_list};
use bixverse_r::single_cell::mc_analysis::metacell_density::*;
use bixverse_r::single_cell::mc_generation::cell_aggregation_utils::*;
use bixverse_r::single_cell::mc_generation::generate::{
    bootstrapped_metacells_on_store, metacells2_on_store, seacells_on_store,
    supercells_on_store, MetaCellsRes,
};
use bixverse_r::single_cell::mc_generation::hdwgcna_meta_cells::*;
use bixverse_r::single_cell::mc_generation::metacells2::{DacParams, MetacellsParams};
use bixverse_r::single_cell::mc_generation::seacells::*;
use bixverse_r::single_cell::mc_generation::super_cells::*;
use extendr_api::*;

use crate::single_cell::utils::{
    knn_data_to_rust, knn_indices_processing, optional_knn, r_list_to_index_vecs,
};

/////////////
// Helpers //
/////////////

/// The aggregated meta cell counts as the R list the wrappers expect.
///
/// ### Params
///
/// * `aggregated` - Meta cells x genes, CSR, raw counts in `data` and the
///   normalised counts in `data_2`.
///
/// ### Returns
///
/// A list with `indptr`, `indices`, `raw_counts`, `norm_counts`, `nrow` and
/// `ncol`.
fn aggregated_to_r_list(aggregated: CompressedSparseData2<u32, f32>) -> Result<List> {
    let norm = aggregated
        .data_2
        .ok_or_else(|| Error::Other("Aggregated meta cells carry no normalised layer".into()))?;

    Ok(list!(
        indptr = aggregated.indptr.r_int_convert(),
        indices = aggregated.indices.r_int_convert(),
        raw_counts = aggregated.data.r_int_convert(),
        norm_counts = norm.r_float_convert(),
        nrow = aggregated.shape.0,
        ncol = aggregated.shape.1
    ))
}

/// Per-cell meta cell id for non-overlapping meta cells.
///
/// ### Params
///
/// * `res` - The meta cells, store indices per meta cell.
///
/// ### Returns
///
/// One entry per store cell, `None` for cells in no meta cell.
fn exclusive_assignments(res: &MetaCellsRes) -> Vec<Option<usize>> {
    let mut out = vec![None; res.n_total];
    for (mc, cells) in res.metacells.iter().enumerate() {
        for &cell in cells {
            out[cell] = Some(mc);
        }
    }
    out
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
    fn rs_metacells2;
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
/// @param cells_to_keep Optional integer vector. Original cell indices
/// (0-indexed!) of the rows of `knn_mat` / `embd`, in row order. If `NULL`,
/// the rows are assumed to map one-to-one onto the count file.
/// @param cells_to_use Optional integer vector. Original cell indices
/// (0-indexed!) to restrict the meta cell generation to, e.g. specific cell
/// types. Needs `cells_to_keep` and `embd`; the kNN graph is then regenerated
/// on the subset. Cells not in `cells_to_keep` are dropped silently.
/// @param meta_cell_params A list containing the meta cell parameters.
/// @param target_size Numeric. Target library size for re-normalisation of
/// the meta cells. Typically `1e4`.
/// @param seed Integer. For reproducibility.
/// @param verbose Integer. `0L` - quiet; `1L` - normal verbosity; `2L` -
/// detailed verbosity.
///
/// @returns A list with the following elements:
/// \itemize{
///  \item assignments - A list with `assignments` (list with one integer
///    vector of 1-indexed meta cell ids per cell, as meta cells can overlap),
///    `metacells` (list of 1-indexed original cell indices per meta cell),
///    `unassigned` (1-indexed cells in no meta cell), `n_metacells`,
///    `n_cells` and `n_unassigned`.
///  \item aggregated - A CSR list (meta cells x genes) with indptr, indices,
///    raw_counts, norm_counts, nrow and ncol.
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
    let knn = knn_mat.map(knn_indices_processing);
    let embd = embd.map(|e| r_matrix_to_faer_fp32(&e));
    let cells_to_keep = cells_to_keep.map(|v| v.r_int_convert());
    let cells_to_use = cells_to_use.map(|v| v.r_int_convert());

    let res = bootstrapped_metacells_on_store(
        &f_path,
        knn.as_deref(),
        embd.as_ref().map(|e| e.as_ref()),
        cells_to_keep.as_deref(),
        cells_to_use.as_deref(),
        &meta_cell_params,
        target_size as f32,
        seed,
        verbose,
    )
    .to_extendr()?;

    Ok(list!(
        assignments = metacells_to_r_list(&res.metacells, res.n_total),
        aggregated = aggregated_to_r_list(res.aggregated)?
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
/// @param embd Numerical matrix. The embedding matrix (for example PCA
/// embedding) used for the generation of the kNN graph and kernel matrix.
/// @param cells_to_keep Optional integer vector. Original cell indices
/// (0-indexed!) of the rows of `embd`, in row order. If `NULL`, the rows are
/// assumed to map one-to-one onto the count file.
/// @param cells_to_use Optional integer vector. Original cell indices
/// (0-indexed!) to restrict the meta cell generation to, e.g. specific cell
/// types. The kNN graph is then regenerated on the subset. Cells not in
/// `cells_to_keep` are dropped silently.
/// @param knn_data Optional list. This contains pre-computed kNN data
/// (including distances). The user has to ensure consistency! Ignored when
/// `cells_to_use` is set.
/// @param seacells_params A list containing the SEACells parameters.
/// @param target_size Numeric. Target library size for re-normalisation of
/// the meta cells. Typically `1e4`.
/// @param seed Integer. For reproducibility purposes.
/// @param verbose Integer. `0L` - quiet; `1L` - normal verbosity; `2L` -
/// detailed verbosity.
///
/// @returns A list with the following elements:
/// \itemize{
///  \item assignments - A list with `assignments` (integer vector with the
///    1-indexed meta cell id per original cell, `-1` if unassigned),
///    `metacells` (list of 1-indexed original cell indices per meta cell),
///    `unassigned` (1-indexed), `n_metacells`, `n_cells` and `n_unassigned`.
///    Empty archetypes are dropped and the ids renumbered.
///  \item aggregated - A CSR list (meta cells x genes) with indptr, indices,
///    raw_counts, norm_counts, nrow and ncol.
///  \item rss - Numerical vector of RSS values from each iteration.
///  \item archetypes - Integer vector with the original cell indices
///    (0-indexed!) of the archetypes of the retained meta cells.
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
    let seacells_params = SEACellsParams::from_r_list(seacells_params)?;
    let embd = r_matrix_to_faer_fp32(&embd);
    let knn = optional_knn(knn_data)?;
    let cells_to_keep = cells_to_keep.map(|v| v.r_int_convert());
    let cells_to_use = cells_to_use.map(|v| v.r_int_convert());

    let res = seacells_on_store(
        &f_path,
        embd.as_ref(),
        knn.as_ref().map(|(i, d)| (i.as_slice(), d.as_slice())),
        cells_to_keep.as_deref(),
        cells_to_use.as_deref(),
        &seacells_params,
        target_size as f32,
        seed,
        verbose,
    )
    .to_extendr()?;

    Ok(list!(
        assignments = assignments_to_r_list(&exclusive_assignments(&res.base), res.base.n_total),
        aggregated = aggregated_to_r_list(res.base.aggregated)?,
        rss = res.rss,
        archetypes = res.archetypes.r_int_convert()
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
/// @param cells_to_keep Optional integer vector. Original cell indices
/// (0-indexed!) of the rows of `embd` / the kNN data, in row order. If
/// `NULL`, the rows are assumed to map one-to-one onto the count file.
/// @param cells_to_use Optional integer vector. Original cell indices
/// (0-indexed!) to restrict the meta cell generation to, e.g. specific cell
/// types. If this is provided, `embd` and `cells_to_keep` are required and
/// the kNN graph will be regenerated on the subset. Cells not in
/// `cells_to_keep` are dropped silently.
/// @param knn_data Optional list. This contains pre-computed kNN data
/// (including distances). The user has to ensure consistency! Ignored when
/// `cells_to_use` is set.
/// @param supercell_params A list containing the SuperCell parameters. The
/// number of meta cells is `ceiling(n_cells / graining_factor)`.
/// @param target_size Numeric. Target library size for re-normalisation of
/// the meta cells. Typically `1e4`.
/// @param seed Integer. For reproducibility purposes.
/// @param verbose Integer. `0L` - quiet; `1L` - normal verbosity; `2L` -
/// detailed verbosity.
///
/// @returns A list with the following elements:
/// \itemize{
///  \item assignments - A list with `assignments` (integer vector with the
///    1-indexed meta cell id per original cell, `-1` if unassigned),
///    `metacells` (list of 1-indexed original cell indices per meta cell),
///    `unassigned` (1-indexed), `n_metacells`, `n_cells` and `n_unassigned`.
///  \item aggregated - A CSR list (meta cells x genes) with indptr, indices,
///    raw_counts, norm_counts, nrow and ncol.
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
    let embd = embd.map(|e| r_matrix_to_faer_fp32(&e));
    let knn = optional_knn(knn_data)?;
    let cells_to_keep = cells_to_keep.map(|v| v.r_int_convert());
    let cells_to_use = cells_to_use.map(|v| v.r_int_convert());

    let res = supercells_on_store(
        &f_path,
        knn.as_ref().map(|(i, d)| (i.as_slice(), d.as_slice())),
        embd.as_ref().map(|e| e.as_ref()),
        cells_to_keep.as_deref(),
        cells_to_use.as_deref(),
        &supercell_params,
        target_size as f32,
        seed,
        verbose,
    )
    .to_extendr()?;

    Ok(list!(
        assignments = assignments_to_r_list(&exclusive_assignments(&res), res.n_total),
        aggregated = aggregated_to_r_list(res.aggregated)?
    ))
}

////////////////
// MetaCells2 //
////////////////

/// Generate meta cells via MetaCells2 (divide and conquer)
///
/// @description
/// `r lifecycle::badge("experimental")`
/// Rust port of the metacells divide and conquer pipeline. Works on the raw
/// counts directly, so no embedding or kNN graph is needed. Cells are split
/// into piles, each pile gets its own feature selection, similarity graph and
/// partition, and outliers are re-run in further rounds. Cells sharing a rare
/// gene module are piled together first unless switched off in `dac_params`.
///
/// @param f_path String. Path to the `counts_cells.bin` file.
/// @param cells_to_keep Integer vector. 0-indexed positions in the count file
/// of the cells to include.
/// @param mc2_params List. The per-pile parameters, see
/// [params_sc_metacells2()].
/// @param dac_params List. The divide and conquer parameters, see
/// [params_sc_metacells2_dac()].
/// @param target_size Numeric. Library size to normalise the aggregated counts
/// to.
/// @param seed Integer. For reproducibility.
/// @param verbose Boolean. Controls verbosity.
///
/// @returns A list with the following elements
/// \itemize{
///   \item assignments - The meta cell assignments in the original cell space.
///   \item aggregated - The aggregated raw and normalised counts in CSR form.
///   \item dissolved - Logical per included cell. Was the cell's candidate
///   dissolved in its last pile.
///   \item level - Integer per included cell. Outlier round that placed it.
///   \item rare_module_of_cell - Integer per included cell. 0-indexed rare
///   gene module, `-1` for none.
///   \item rare_module_of_gene - Integer per gene. 0-indexed rare gene module,
///   `-1` for none.
///   \item selected_genes - Integer. 1-indexed genes selected in any pile plus
///   the rare module genes.
/// }
///
/// @export
///
/// @keywords internal
///
/// @references Ben-Kiki, et al., Genome Biol., 2022.
#[extendr]
fn rs_metacells2(
    f_path: String,
    cells_to_keep: Vec<i32>,
    mc2_params: List,
    dac_params: List,
    target_size: f64,
    seed: usize,
    verbose: bool,
) -> extendr_api::Result<List> {
    let mc2_params = MetacellsParams::from_r_list(mc2_params)?;
    let dac_params = DacParams::from_r_list(dac_params)?;

    let res = metacells2_on_store(
        &f_path,
        &cells_to_keep.r_int_convert(),
        &mc2_params,
        &dac_params,
        target_size as f32,
        seed,
        verbose,
    )
    .to_extendr()?;

    let selected_genes: Vec<i32> = res.selected_genes.iter().map(|&g| g as i32 + 1).collect();

    Ok(list!(
        assignments = assignments_to_r_list(&exclusive_assignments(&res.base), res.base.n_total),
        aggregated = aggregated_to_r_list(res.base.aggregated)?,
        dissolved = res.dissolved,
        level = res.level,
        rare_module_of_cell = res.rare_module_of_cell,
        rare_module_of_gene = res.rare_module_of_gene,
        selected_genes = selected_genes
    ))
}

//////////////////////
// MetaCell density //
//////////////////////

/// Calculates diffusion maps for density calculations for meta cells
///
/// @description
/// `r lifecycle::badge("experimental")`
/// Builds multiscale diffusion components from the kNN graph and uses the
/// distance to the `k_density`-th neighbour in that space as a density proxy.
/// The lower quartile of these distances is tagged high density, the upper
/// quartile low density, the rest mid.
///
/// @param knn_data Named list. The kNN data with `indices` (0-indexed!),
/// `dist`, `k` and `dist_metric`.
/// @param n_dcs Integer. The number of diffusion coordinates to return.
/// Typically `10`.
/// @param k_density Integer. The k-nearest neighbour to use for the density
/// estimation. Typically `150`.
/// @param knn_params List. The kNN parameters defined by
/// [params_sc_neighbours()], used for the search in diffusion space.
/// @param verbose Integer. `0L` - quiet; `1L` - normal verbosity; `2L` -
/// detailed verbosity.
/// @param seed Integer. For reproducibility.
///
/// @returns A list with the following items
/// \itemize{
///   \item dcs - Numerical matrix of cells x `n_dcs` with the multiscale
///   diffusion components.
///   \item density_distances - Numerical vector. Distance to the
///   `k_density`-th neighbour in diffusion space per cell.
///   \item regions - Character vector. `"high"`, `"mid"` or `"low"` density
///   per cell.
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
    let (knn_indices, knn_distances, _, _) = knn_data_to_rust(knn_data)?;
    let knn_params = KnnParams::from_r_list(knn_params)?;

    let density_res: DiffusionDensity = compute_diffusion_density(
        &knn_indices,
        &knn_distances,
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
/// Calculates the meta cell compactness, i.e. the average variance across the
/// diffusion components over the cells of each meta cell. Lower is better.
///
/// @param dc Numerical matrix. The diffusion map coordinates, cells x
/// components.
/// @param meta_cells List. Per meta cell, an integer vector with the row
/// indices (1-indexed!) into `dc`.
///
/// @returns Numerical vector with one compactness value per meta cell. Empty
/// meta cells yield `NaN`.
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
/// Calculates the separation, i.e. the Euclidean distance from each meta cell
/// centroid in diffusion space to the nearest other meta cell centroid.
/// Higher is better.
///
/// @param dc Numerical matrix. The diffusion map coordinates, cells x
/// components.
/// @param meta_cells List. Per meta cell, an integer vector with the row
/// indices (1-indexed!) into `dc`.
///
/// @returns Numerical vector with one separation value per meta cell. Empty
/// meta cells yield `NaN`; a lone non-empty meta cell yields `Inf`.
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
/// to return the sum of the raw counts or the average of the normalised
/// counts.
///
/// @param f_path String. Path to the `counts_cells.bin` file.
/// @param cell_indices_ls List. Each element contains the 0-indexed positions
/// of the cells to aggregate.
/// @param assay String. One of `c("raw", "norm")`. `"raw"` sums the raw
/// counts, `"norm"` averages the normalised counts.
/// @param verbose Integer. `0L` - quiet; `1L` - normal verbosity; `2L` -
/// detailed verbosity.
///
/// @returns A dense numerical matrix of pseudo-bulked samples x genes.
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
    let bulk_type = parse_pseudo_bulk(&assay)
        .ok_or_else(|| Error::Other(format!("Unknown assay '{assay}'. Use 'raw' or 'norm'.")))?;
    let cell_indices = r_list_to_index_vecs(cell_indices_ls, "cell_indices_ls")?;

    let reader = ParallelSparseReader::new(&f_path).to_extendr()?;

    let data =
        get_pseudo_bulked_counts_dense(&reader, &cell_indices, bulk_type, verbose).to_extendr()?;

    Ok(faer_to_r_matrix(data.as_ref()))
}

/// Pseudo-bulk a set of cells (sparse)
///
/// @description
/// `r lifecycle::badge("experimental")`
/// This function will return a sparse matrix of
/// `length(cell_indices_ls) x number of genes` (in list form in CSR).
/// The function has the option to return the sum of the raw counts or the
/// average of the normalised counts.
///
/// @param f_path String. Path to the `counts_cells.bin` file.
/// @param cell_indices_ls List. Each element contains the 0-indexed positions
/// of the cells to aggregate.
/// @param assay String. One of `c("raw", "norm")`. `"raw"` sums the raw
/// counts, `"norm"` averages the normalised counts.
/// @param verbose Integer. `0L` - quiet; `1L` - normal verbosity; `2L` -
/// detailed verbosity.
///
/// @returns A list with the following elements (easy to convert into CSR in R)
/// \itemize{
///   \item indptr - The index pointers (representing pseudo-bulked samples)
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
    let bulk_type = parse_pseudo_bulk(&assay)
        .ok_or_else(|| Error::Other(format!("Unknown assay '{assay}'. Use 'raw' or 'norm'.")))?;
    let cell_indices = r_list_to_index_vecs(cell_indices_ls, "cell_indices_ls")?;

    let reader = ParallelSparseReader::new(&f_path).to_extendr()?;

    let data: CompressedSparseData2<f64> =
        get_pseudo_bulked_counts_sparse(&reader, &cell_indices, bulk_type, verbose).to_extendr()?;

    Ok(list!(
        indptr = data.indptr,
        indices = data.indices,
        data = data.data,
        nrow = data.shape.0,
        ncol = data.shape.1
    ))
}
