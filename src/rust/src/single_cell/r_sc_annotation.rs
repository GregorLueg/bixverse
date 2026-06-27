//! Functions related to cell type label transfer or cell type annotations

use bixverse_rs::prelude::*;
use bixverse_rs::single_cell::sc_annotation::sc_type::*;
use bixverse_rs::single_cell::sc_annotation::symphony::*;
use bixverse_rs::single_cell::sc_batch_correction::harmony::HarmonyParams;
use bixverse_rs::single_cell::sc_batch_correction::harmony_v2::HarmonyParamsV2;
use bixverse_rs::single_cell::sc_processing::pca::SingleCellPcaParams;
use extendr_api::*;

use crate::single_cell::utils::process_cell_markers;

////////////////////
// extendr Module //
////////////////////

extendr_module! {
    // module
    mod r_sc_annotation;
    // sctype
    fn rs_sc_type;
    fn rs_sc_type_cluster_assignment;
    // symphony
    fn rs_build_symphony_ref;
    fn rs_symphony_map_query;
    fn rs_transfer_labels_symphony;
}

////////////
// ScType //
////////////

/// Run the ScType scoring approach
///
/// @description
/// `r lifecycle::badge("experimental")`
/// This Rust function implements the cell type scoring approach from Ianevski
/// et al. (2022).
///
/// @param f_path String. Path to the `counts_genes.bin` file.
/// @param cell_indices Integer vector. 0-indexed(!) positions of cells to
/// include in the analysis
/// @param cell_markers A list with the cell marker gene indices.
/// @param sensitivity Boolean. Shall a sensitivity correction be applied that
/// downweights common cell type markers.
/// @param weight_floor Optional numeric. If `sensitivity = TRUE`, what is
/// the weight floor. If not provided, defaults to `0.1`.
/// @param verbose Integer. `0L` - quiet; `1L` - normal verbosity; `2L` -
/// detailed verbosity.
///
/// @returns A list with
/// \itemize{
///   \item cell_types - String vector. The cell types
///   \item scores - Row-major scores (cells x cell_types).
///   \item n_cells - Number of cells
///   \item n_cell_types - Number of cell types
/// }
///
/// @export
///
/// @keywords internal
#[extendr]
fn rs_sc_type(
    f_path: &str,
    cell_indices: Vec<i32>,
    cell_markers: List,
    sensitivity: bool,
    weight_floor: Option<f64>,
    verbose: usize,
) -> Result<List> {
    let cell_markers: Vec<CellTypeMarkers> = process_cell_markers(cell_markers)?;
    let cell_indices = cell_indices.r_int_convert();

    let res: SctypeRes = run_sctype(
        f_path,
        &cell_indices,
        &cell_markers,
        sensitivity,
        None,
        weight_floor.map(|x| x as f32),
        verbose,
    )
    .to_extendr()?;

    Ok(list!(
        cell_types = res.cell_types,
        scores = res.scores.r_float_convert(),
        n_cells = res.n_cells,
        n_cell_types = res.n_cell_types
    ))
}

/// Score the individual clusters based on ScType
///
/// @description
/// `r lifecycle::badge("experimental")`
/// This Rust function implements the cell type scoring approach from Ianevski
/// et al. (2022).
///
/// @param sc_type_res List. The ScType results.
/// @param cluster_labels Integer. Cluster assignment. Needs to be of length
/// of scored cells.
///
/// @returns A list with
/// \itemize{
///  \item cluster_id - The cluster id/integer
///  \item cell_type - String; the predicted cell type
///  \item score - The final score for the clsuter.
///  \item n_cells - The number of cells in the cluster.
/// }
///
/// @export
///
/// @keywords internal
#[extendr]
fn rs_sc_type_cluster_assignment(sc_type_res: List, cluster_labels: Vec<i32>) -> Result<List> {
    let cluster_labels = cluster_labels.r_int_convert();

    let sc_res = SctypeRes::from_r_list(sc_type_res)?;

    let cluster_assignments: Vec<ScTypeClusterAssignment> =
        assign_clusters(&sc_res, &cluster_labels).to_extendr()?;

    let mut cluster_id: Vec<usize> = Vec::with_capacity(cluster_assignments.len());
    let mut cell_type: Vec<String> = Vec::with_capacity(cluster_assignments.len());
    let mut scores: Vec<f64> = Vec::with_capacity(cluster_assignments.len());
    let mut n_cells: Vec<usize> = Vec::with_capacity(cluster_assignments.len());

    for res in cluster_assignments {
        cluster_id.push(res.cluster);
        cell_type.push(res.cell_type);
        scores.push(res.score as f64);
        n_cells.push(res.n_cells);
    }

    Ok(list!(
        cluster_id = cluster_id.r_int_convert(),
        cell_type = cell_type,
        scores = scores,
        n_cells = n_cells.r_int_convert()
    ))
}

//////////////
// Symphony //
//////////////

/// Build a Symphony reference (Rust)
///
/// @description
/// `r lifecycle::badge("experimental")`
/// Builds the Symphony reference in Rust, see Kang et al.
///
/// @param f_path_gene String. Path to the gene-based binary file.
/// @param f_path_cell String. Path to the cell-based binary file.
/// @param cell_indices Integer vector. 0-based cell indices.
/// @param hvg_indices Integer vector. 0-based HVG indices.
/// @param batch_labels List of 0-indexed integer vectors (one per batch
/// variable).
/// @param pca_params List. Output of `params_sc_pca()`.
/// @param no_pcs Integer.
/// @param harmony_params List. Output of `params_sc_harmony()` or
/// `params_sc_harmony_v2()`.
/// @param harmony_version String. "v1" or "v2".
/// @param seed Integer.
/// @param verbose Integer. `0L` - quiet; `1L` - normal verbosity; `2L` -
/// detailed verbosity.
///
/// @return A list with gene_means, gene_sds, loadings, z_orig, z_corr, r,
/// centroids, nr, c.
///
/// @references
/// Kang et al., Nat Comm, 2021.
///
/// @export
///
/// @keywords internal
#[extendr]
#[allow(clippy::too_many_arguments)]
fn rs_build_symphony_ref(
    f_path_gene: String,
    f_path_cell: &str,
    cell_indices: Vec<i32>,
    hvg_indices: Vec<i32>,
    batch_labels: List,
    pca_params: List,
    no_pcs: usize,
    harmony_params: List,
    harmony_version: String,
    seed: usize,
    verbose: usize,
) -> extendr_api::Result<List> {
    let verbosity = parse_verbosity_level(verbose);
    let cell_indices: Vec<usize> = cell_indices.r_int_convert();
    let hvg_indices: Vec<usize> = hvg_indices.r_int_convert();

    let mut batch_indices: Vec<Vec<usize>> = Vec::new();
    for i in 0..batch_labels.len() {
        let batch_indices_i = batch_labels.elt(i)?;
        let batch_indices_i = batch_indices_i.as_integer_vector().unwrap();
        batch_indices.push(batch_indices_i.r_int_convert());
    }

    let pca_params = SingleCellPcaParams::from_r_list(pca_params)?;

    let harmony_backend = match harmony_version.as_str() {
        "v1" => HarmonyBackend::V1(HarmonyParams::from_r_list(harmony_params)?),
        "v2" => HarmonyBackend::V2(HarmonyParamsV2::from_r_list(harmony_params)?),
        _ => {
            return Err(extendr_api::Error::Other(
                "harmony_version must be 'v1' or 'v2'".into(),
            ))
        }
    };

    let offsets = if pca_params.clr {
        if verbosity.normal_verbosity() {
            println!("PFlogPF-transformation requested. Loading offsets from disk.")
        }

        let reader = ParallelSparseReader::new(f_path_cell).to_extendr()?;

        let offsets = reader.get_clr_offsets(&cell_indices, None).to_extendr()?;

        Some(offsets)
    } else {
        None
    };

    let reference = build_symphony_reference(
        &f_path_gene,
        &cell_indices,
        &hvg_indices,
        &batch_indices,
        &pca_params,
        no_pcs,
        harmony_backend,
        offsets.as_deref(),
        seed,
        verbose,
    )
    .to_extendr()?;

    let gene_means: Vec<f64> = reference.gene_means.iter().map(|&x| x as f64).collect();
    let gene_sds: Vec<f64> = reference.gene_sds.iter().map(|&x| x as f64).collect();
    let nr: Vec<f64> = reference.nr.iter().map(|&x| x as f64).collect();

    Ok(list!(
        gene_means = gene_means,
        gene_sds = gene_sds,
        loadings = faer_to_r_matrix(reference.loadings.as_ref()),
        z_orig = faer_to_r_matrix(reference.z_orig.as_ref()),
        z_corr = faer_to_r_matrix(reference.z_corr.as_ref()),
        r = faer_to_r_matrix(reference.r.as_ref()),
        centroids = faer_to_r_matrix(reference.centroids.as_ref()),
        nr = nr,
        c = faer_to_r_matrix(reference.c.as_ref())
    ))
}

/// Map a query onto a Symphony reference (Rust)
///
/// @description
/// `r lifecycle::badge("experimental")`
/// Maps a given other SingleCells data set to the Symphony reference, see
/// Kang et al.
///
/// @param f_path_query String. Path to the query gene-based binary file.
/// @param cell_indices_query Integer vector. 0-based query cell indices.
/// @param gene_means,gene_sds Numerical vectors. Reference per-HVG stats.
/// @param loadings Reference PCA loadings (n_hvgs x d).
/// @param centroids Reference centroids (K x d).
/// @param nr Reference cluster sizes (length K).
/// @param c_cache Reference compression term R*Z_corr (K x d).
/// @param ref_to_query_gene_map Integer vector. For each reference HVG slot,
/// the 0-based query gene index, or `NA_integer_` if absent.
/// @param batch_labels_query List of 0-indexed integer vectors (empty = no
/// batch correction).
/// @param params_symphony Named list. Contains the parameters for the referemce
/// generation.
/// @param verbose Integer. `0L` - quiet; `1L` - normal verbosity; `2L` -
/// detailed verbosity.
///
/// @return A list with z_pca, z_corr, r.
///
/// @references
/// Kang et al., Nat Comm, 2021.
///
/// @export
///
/// @keywords internal
#[extendr]
#[allow(clippy::too_many_arguments)]
fn rs_symphony_map_query(
    f_path_query: String,
    cell_indices_query: Vec<i32>,
    gene_means: Vec<f64>,
    gene_sds: Vec<f64>,
    loadings: RMatrix<f64>,
    centroids: RMatrix<f64>,
    nr: Vec<f64>,
    c_cache: RMatrix<f64>,
    ref_to_query_gene_map: Vec<i32>,
    batch_labels_query: List,
    params_symphony: List,
    verbose: usize,
) -> extendr_api::Result<List> {
    let cell_indices_query: Vec<usize> = cell_indices_query.r_int_convert();

    let na_int = i32::MIN;
    let gene_map: Vec<Option<usize>> = ref_to_query_gene_map
        .iter()
        .map(|&x| {
            if x == na_int || x < 0 {
                None
            } else {
                Some(x as usize)
            }
        })
        .collect();

    let mut batch_indices: Vec<Vec<usize>> = Vec::new();
    for i in 0..batch_labels_query.len() {
        let batch_indices_i = batch_labels_query.elt(i)?;
        let batch_indices_i = batch_indices_i.as_integer_vector().unwrap();
        batch_indices.push(batch_indices_i.r_int_convert());
    }

    let gene_means_f32: Vec<f32> = gene_means.iter().map(|&x| x as f32).collect();
    let gene_sds_f32: Vec<f32> = gene_sds.iter().map(|&x| x as f32).collect();
    let nr_f32: Vec<f32> = nr.iter().map(|&x| x as f32).collect();
    let loadings = r_matrix_to_faer_fp32(&loadings);
    let centroids = r_matrix_to_faer_fp32(&centroids);
    let c_cache = r_matrix_to_faer_fp32(&c_cache);

    let params = SymphonyMapParams::from_r_list(params_symphony)?;

    let res = symphony_map_query_parts(
        &f_path_query,
        &cell_indices_query,
        &gene_means_f32,
        &gene_sds_f32,
        loadings.as_ref(),
        centroids.as_ref(),
        &nr_f32,
        c_cache.as_ref(),
        &gene_map,
        &batch_indices,
        &params,
        verbose,
    )
    .to_extendr()?;

    Ok(list!(
        z_pca = faer_to_r_matrix(res.z_pca.as_ref()),
        z_corr = faer_to_r_matrix(res.z_corr.as_ref()),
        r = faer_to_r_matrix(res.r.as_ref())
    ))
}

/// Transfer labels from a Symphony reference to a query via kNN vote
///
/// @description
/// `r lifecycle::badge("experimental")`
///
/// @param reference_z_corr Reference Harmony-corrected embedding (N_ref x d).
/// @param query_z_corr Query Symphony-corrected embedding (N_q x d).
/// @param reference_labels 0-based integer-encoded reference labels.
/// @param n_labels Number of distinct labels.
/// @param knn_params List. Output of `params_sc_knn()`.
/// @param seed Integer.
/// @param verbose Integer. 0/1/2.
///
/// @return A list with `predicted` (0-based integer per query cell) and
/// `confidence` (vote share of the winning label).
///
/// @export
///
/// @keywords internal
#[extendr]
fn rs_transfer_labels_symphony(
    reference_z_corr: RMatrix<f64>,
    query_z_corr: RMatrix<f64>,
    reference_labels: Vec<i32>,
    n_labels: usize,
    knn_params: List,
    seed: usize,
    verbose: usize,
) -> extendr_api::Result<List> {
    let ref_mat = r_matrix_to_faer_fp32(&reference_z_corr);
    let qry_mat = r_matrix_to_faer_fp32(&query_z_corr);
    let ref_labels: Vec<usize> = reference_labels.r_int_convert();
    let knn_params = KnnParams::from_r_list(knn_params)?;
    let verbosity = parse_verbosity_level(verbose);

    let (knn_idx, _knn_dist) = generate_knn_cross_with_dist(
        ref_mat.as_ref(),
        qry_mat.as_ref(),
        &knn_params,
        seed,
        verbosity.normal_verbosity(),
    )
    .to_extendr()?;

    let (preds, conf) = knn_majority_vote(&knn_idx, &ref_labels, n_labels);

    let preds_i32: Vec<i32> = preds.iter().map(|&x| x as i32).collect();
    let conf_f64: Vec<f64> = conf.iter().map(|&x| x as f64).collect();

    Ok(list!(predicted = preds_i32, confidence = conf_f64))
}
