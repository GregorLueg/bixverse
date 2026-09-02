//! Helpers for the single cell part, especially around transforming structures
//! into R lists.

use bixverse_rs::prelude::*;
use bixverse_rs::single_cell::sc_analysis::dialogue::{DialogueResult, ProgrammeSignature};
use bixverse_rs::single_cell::sc_analysis::fast_clusters::{
    FastLouvainGridResult, FastLouvainResults,
};
use bixverse_rs::single_cell::sc_analysis::nebula::NebulaScRes;
use bixverse_rs::single_cell::sc_annotation::sc_type::CellTypeMarkers;
use bixverse_rs::single_cell::sc_processing::hvg::HvgDispersionRes;
use either::Either;
use std::collections::HashMap;

use extendr_api::*;

///////////
// Types //
///////////

/// Results data from the neighbours class
///
/// ### Fields
///
/// * `0` - The kNN indices
/// * `1` - The kNN distances
/// * `2` - Number of neighbours
/// * `3` - Distance metric
pub type NeighboursData = Result<(Vec<Vec<usize>>, Vec<Vec<f32>>, usize, String)>;

////////////////////
// Dispersion res //
////////////////////

/// Helper function to flatten the dispersion results to a List
///
/// ### Params
///
/// * `results` - A vector of `HvgDispersionRes`.
///
/// ### Returns
///
/// The results list
pub fn flatten_dispersion_batches(results: Vec<HvgDispersionRes>) -> List {
    let n_genes = results.first().map_or(0, |res| res.mean.len());
    let total_len = n_genes * results.len();
    let mut mean_flat = Vec::with_capacity(total_len);
    let mut disp_flat = Vec::with_capacity(total_len);
    let mut disp_scaled_flat = Vec::with_capacity(total_len);
    let mut bin_flat = Vec::with_capacity(total_len);
    let mut batch_idx = Vec::with_capacity(total_len);
    let mut gene_idx = Vec::with_capacity(total_len);

    for (batch, res) in results.into_iter().enumerate() {
        mean_flat.extend(res.mean);
        disp_flat.extend(res.dispersion);
        disp_scaled_flat.extend(res.dispersion_scaled);
        bin_flat.extend(res.bin);
        batch_idx.extend(vec![batch as i32; n_genes]);
        gene_idx.extend(0..n_genes as i32);
    }

    list!(
        mean = mean_flat,
        dispersion = disp_flat,
        dispersion_scaled = disp_scaled_flat,
        bin = bin_flat,
        batch = batch_idx,
        gene_idx = gene_idx
    )
}

/////////
// kNN //
/////////

/// Process R KNN indices to the Rust variant
///
/// ### Params
///
/// * `knn_mat` - Samples x indices of the k-nearest neighbours (0-indexed!)
///
/// ### Returns
///
/// A `Vec<Vec<usize>>`
pub fn knn_indices_processing(knn_mat: RMatrix<i32>) -> Vec<Vec<usize>> {
    let nrow = knn_mat.nrows();
    let ncol = knn_mat.ncols();
    let data = knn_mat.data();

    let mut out: Vec<Vec<usize>> = (0..nrow).map(|_| Vec::with_capacity(ncol)).collect();
    for i in 0..ncol {
        let col_offset = i * nrow;
        for j in 0..nrow {
            out[j].push(data[col_offset + j] as usize);
        }
    }
    out
}

/// Process R KNN distances to the Rust variant
///
/// ### Params
///
/// * `knn_dist` - Samples x indices of the k-nearest neighbours (0-indexed!)
///
/// ### Returns
///
/// A `Vec<Vec<f32>>`
pub fn knn_distances_processing(knn_dist: RMatrix<f64>) -> Vec<Vec<f32>> {
    let nrow = knn_dist.nrows();
    let ncol = knn_dist.ncols();
    let data = knn_dist.data();

    let mut out: Vec<Vec<f32>> = (0..nrow).map(|_| Vec::with_capacity(ncol)).collect();
    for i in 0..ncol {
        let col_offset = i * nrow;
        for j in 0..nrow {
            out[j].push(data[col_offset + j] as f32);
        }
    }
    out
}

/// Transform R kNN data to Rust data
///
/// ### Params
///
/// * `knn_data` - R list with the kNN data
///
/// ###
///
/// The [NeighboursData] or an error.
pub fn knn_data_to_rust(knn_data: List) -> NeighboursData {
    let data: HashMap<&str, Robj> = knn_data.try_into()?;

    let knn_indices: RArray<i32, 2> = data
        .get("indices")
        .ok_or_else(|| Error::Other("missing 'indices'".into()))?
        .as_matrix()
        .ok_or_else(|| Error::Other("'indices' is not a matrix".into()))?;

    let knn_dist: RArray<f64, 2> = data
        .get("dist")
        .ok_or_else(|| Error::Other("missing 'dist'".into()))?
        .as_matrix()
        .ok_or_else(|| Error::Other("'dist' is not a matrix".into()))?;

    let dist_metric = data
        .get("dist_metric")
        .ok_or_else(|| Error::Other("missing 'dist_metric'".into()))?
        .as_str()
        .ok_or_else(|| Error::Other("'dist_metric' is not a string".into()))?
        .to_string();

    let k = data
        .get("k")
        .ok_or_else(|| Error::Other("missing 'k'".into()))?
        .as_integer()
        .ok_or_else(|| Error::Other("'k' is not an integer".into()))? as usize;

    let knn_indices = knn_indices_processing(knn_indices);
    let knn_dist = knn_distances_processing(knn_dist);

    Ok((knn_indices, knn_dist, k, dist_metric))
}

/////////////////
// FastCluster //
/////////////////

/// Type for the fast clustering single results
///
/// ### Fields
///
/// * `0` - The Louvain membership per resolution
/// * `1` - Optional k-means cluster membership
/// * `2` - Optional k-means centroids
pub type FastClusterSingle = Result<(Vec<Vec<usize>>, Option<Vec<i32>>, Option<RMatrix<f64>>)>;

/// Type for the fast clustering single results
///
/// ### Fields
///
/// * `0` - The [FastLouvainGridResult] per resolution
/// * `1` - Optional k-means cluster membership
/// * `2` - Optional k-means centroids
pub type FastClusterGrid = Result<(
    Vec<FastLouvainGridResult>,
    Option<Vec<i32>>,
    Option<RMatrix<f64>>,
)>;

/// Helper function to extract (single) results from the FastClustering
///
/// ### Params
///
/// * `res` - The [FastLouvainResults]
/// * `return_km` - Shall the km results be returned
///
/// ### Returns
///
/// [FastClusterSingle]
pub fn fast_cluster_unwrap_single(res: FastLouvainResults, return_km: bool) -> FastClusterSingle {
    let memberships = match res.get_assignments() {
        Either::Left(v) => v,
        Either::Right(_) => panic!("expected Single variant"),
    };

    let (k_means_membership, centroids) = if return_km {
        let k_means_clusters = res.get_k_mean_clusters().to_extendr()?;
        let centroids = res.get_centroids().to_extendr()?;

        (
            Some(k_means_clusters.r_int_convert()),
            Some(faer_to_r_matrix(centroids.as_ref())),
        )
    } else {
        (None, None)
    };

    Ok((memberships, k_means_membership, centroids))
}

/// Helper function to extract (single) results from the FastClustering
///
/// ### Params
///
/// * `res` - The [FastLouvainResults]
/// * `return_km` - Shall the km results be returned
///
/// ### Returns
///
/// [FastClusterSingle]
pub fn fast_cluster_unwrap_multiple(res: FastLouvainResults, return_km: bool) -> FastClusterGrid {
    let memberships = match res.get_assignments() {
        Either::Left(_) => panic!("expected Grid variant"),
        Either::Right(v) => v,
    };

    let (k_means_membership, centroids) = if return_km {
        let k_means_clusters = res.get_k_mean_clusters().to_extendr()?;
        let centroids = res.get_centroids().to_extendr()?;

        (
            Some(k_means_clusters.r_int_convert()),
            Some(faer_to_r_matrix(centroids.as_ref())),
        )
    } else {
        (None, None)
    };

    Ok((memberships, k_means_membership, centroids))
}

/// Process the fast cluster Louvain results
///
/// ### Params
///
/// * `results` - Vector of [FastLouvainGridResult]
///
/// ### Retuns
///
/// A list with the results
pub fn process_fc_louvain_results(results: Vec<FastLouvainGridResult>) -> Result<List> {
    let mut mean_ari: Vec<f32> = Vec::with_capacity(results.len());
    let mut median_ari: Vec<f32> = Vec::with_capacity(results.len());
    let mut mean_conductance: Vec<f32> = Vec::with_capacity(results.len());
    let mut median_conductance: Vec<f32> = Vec::with_capacity(results.len());
    let mut mean_n_comms: Vec<f32> = Vec::with_capacity(results.len());

    let mut res = List::new(results.len());

    for (index, louvain_res) in results.iter().enumerate() {
        let membership = louvain_res.best_labels.clone().r_int_convert();
        res.set_elt(index, Robj::from(membership))?;

        mean_ari.push(louvain_res.mean_ari);
        median_ari.push(louvain_res.median_ari);
        mean_conductance.push(louvain_res.mean_conductance);
        median_conductance.push(louvain_res.median_conductance);
        mean_n_comms.push(louvain_res.mean_n_communities);
    }

    let stats = list![
        mean_ari = mean_ari.r_float_convert(),
        median_ari = median_ari.r_float_convert(),
        mean_conductance = mean_conductance.r_float_convert(),
        median_conductance = median_conductance.r_float_convert(),
        mean_n_comms = mean_n_comms.r_float_convert(),
    ];

    Ok(list![memberships = res, stats = stats])
}

////////////
// ScType //
////////////

/// Process a list of cell markers
///
/// ### Params
///
/// * `r_list` - The R list to parse
///
/// ### Returns
///
/// A vector of [CellTypeMarkers].
pub fn process_cell_markers(r_list: List) -> Result<Vec<CellTypeMarkers>> {
    let mut res: Vec<CellTypeMarkers> = Vec::with_capacity(r_list.len());

    let iterator = 0..r_list.len();

    for i in iterator {
        let element = r_list
            .elt(i)?
            .as_list()
            .ok_or_else(|| Error::Other(format!("Cell marker element {} is not a list", i)))?;
        let markers = CellTypeMarkers::from_r_list(element)?;

        res.push(markers);
    }

    Ok(res)
}

///////////////////////////////////////
// Hotspot working size calculations //
///////////////////////////////////////

/// Derive a panel size from a working-memory budget for the streaming pair
/// path.
///
/// Working set is roughly `20 * n_cells * panel_size` bytes (two resident
/// panels plus transient load buffers); the `n_genes^2` output matrices sit on
/// top and are not controlled here. Clamps to `[512, n_genes]`: below ~512 the
/// GEMM blocks are too small to be efficient, and `panel_size == n_genes` is
/// the single-panel (no-reload) case.
///
/// ### Parmas
///
/// * `working_mem_gb` - Working memory to allocate here
/// * `n_cells` - Number of cells
/// * `n_genes` - Number of genes
///
/// ### Returns
///
/// The panel size
pub fn panel_size_from_mem(working_mem_gb: f64, n_cells: usize, n_genes: usize) -> usize {
    const BYTES_PER_GENE_PER_CELL: f64 = 20.0;
    let budget = working_mem_gb * (1u64 << 30) as f64;
    let raw = (budget / (BYTES_PER_GENE_PER_CELL * n_cells.max(1) as f64)).floor() as usize;
    let lo = 512.min(n_genes.max(1));
    raw.clamp(lo, n_genes.max(1))
}

//////////////
// NicheNet //
//////////////

/// NicheNet network data
///
/// ### Fields
///
/// * `0` - The from indices
/// * `1` - The to indices
/// * `2` - The edge weights
pub type NicheNetNetwork = Result<(Vec<u32>, Vec<u32>, Vec<f64>)>;

/// Transform network data into a [NicheNetNetwork]
///
/// ### Params
///
/// * `network_data` - R list. Needs to have the `"from"`, `"to"` and
///   `"weights"` of the network.
///
/// ### Returns
///
/// The [NicheNetNetwork]
pub fn prep_nichenet_network(network_data: List) -> NicheNetNetwork {
    let network_data: HashMap<&str, Robj> = network_data.try_into()?;

    let from = network_data
        .get("from")
        .ok_or_else(|| Error::Other("missing 'from'".into()))?
        .as_integer_vector()
        .ok_or_else(|| Error::Other("'from' is not a integer vector".into()))?;
    let from = from.iter().map(|x| *x as u32).collect();

    let to = network_data
        .get("to")
        .ok_or_else(|| Error::Other("missing 'to'".into()))?
        .as_integer_vector()
        .ok_or_else(|| Error::Other("'to' is not a integer vector".into()))?;
    let to = to.iter().map(|x| *x as u32).collect();

    let weight = network_data
        .get("weight")
        .ok_or_else(|| Error::Other("missing 'weight'".into()))?
        .as_real_vector()
        .ok_or_else(|| Error::Other("'weight' is not a real vector".into()))?;

    Ok((from, to, weight))
}

//////////////
// DIALOGUE //
//////////////

/// The per-cell-type inputs DIALOGUE needs, owned so the borrows stay alive.
///
/// ### Fields
///
/// * `0` - Global cell indices per cell type
/// * `1` - Feature matrices per cell type, still in R's column-major storage
pub type DialogueInputs = Result<(Vec<Vec<usize>>, Vec<RMatrix<f64>>)>;

/// Pulls the per-cell-type cell indices and feature matrices off two R lists.
///
/// The matrices are returned owned rather than as `MatRef`, because
/// `r_matrix_to_faer` borrows and the caller has to keep the storage alive for
/// the length of the run.
///
/// ### Params
///
/// * `cell_type_indices` - List of integer vectors, 0-indexed global cell
///   positions per cell type
/// * `features` - List of numeric matrices, one per cell type
///
/// ### Returns
///
/// The [DialogueInputs], or an error naming the offending element. Whether the
/// two line up in length and shape is left to `dialogue_run`, which checks it
/// anyway and reports it in the same terms for every caller.
pub fn dialogue_inputs_to_rust(cell_type_indices: List, features: List) -> DialogueInputs {
    let mut cells: Vec<Vec<usize>> = Vec::with_capacity(cell_type_indices.len());
    for i in 0..cell_type_indices.len() {
        let elem = cell_type_indices.elt(i)?;
        let idx = elem.as_integer_vector().ok_or_else(|| {
            Error::Other(format!("cell type {} indices are not an integer vector", i))
        })?;
        if let Some(bad) = idx.iter().find(|&&v| v < 0) {
            return Err(Error::Other(format!(
                "cell type {i} holds a negative cell index ({bad}); these are 0-indexed"
            )));
        }
        cells.push(idx.iter().map(|&v| v as usize).collect());
    }

    let mut mats: Vec<RMatrix<f64>> = Vec::with_capacity(features.len());
    for i in 0..features.len() {
        let elem = features.elt(i)?;
        let mat: RMatrix<f64> = elem
            .as_matrix()
            .ok_or_else(|| Error::Other(format!("features {} is not a numeric matrix", i)))?;
        mats.push(mat);
    }

    Ok((cells, mats))
}

/// Flattens a nested `[cell_type][programme]` signature list into an R list.
///
/// ### Params
///
/// * `signatures` - The per-cell-type, per-programme gene lists
///
/// ### Returns
///
/// A list of lists, innermost holding `up` and `down` as 0-indexed gene
/// positions.
fn signatures_to_r_list(signatures: &[Vec<ProgrammeSignature>]) -> List {
    List::from_values(signatures.iter().map(|per_type| {
        List::from_values(per_type.iter().map(|sig| {
            list!(
                up = sig.up.clone().r_int_convert(),
                down = sig.down.clone().r_int_convert()
            )
        }))
    }))
}

/// Flattens the DIALOGUE result into an R list.
///
/// Everything index-like stays 0-indexed and the R wrappers add the one, the
/// same as every other single cell binding. Shapes: `emp_p` and `pair_cor` are
/// `k x n_pairs` with the columns in `combn(n_cell_types, 2)` order,
/// `refit_fidelity` is `n_cell_types x k`.
///
/// ### Params
///
/// * `res` - The finished [DialogueResult]
///
/// ### Returns
///
/// The results list.
pub fn dialogue_res_to_r_list(res: &DialogueResult) -> List {
    let verdicts = list!(
        cell_type = res
            .verdicts
            .iter()
            .map(|v| v.cell_type as i32)
            .collect::<Vec<i32>>(),
        programme = res
            .verdicts
            .iter()
            .map(|v| v.programme as i32)
            .collect::<Vec<i32>>(),
        gene = res
            .verdicts
            .iter()
            .map(|v| v.gene as i32)
            .collect::<Vec<i32>>(),
        up = res.verdicts.iter().map(|v| v.up).collect::<Vec<bool>>(),
        n_supporting = res
            .verdicts
            .iter()
            .map(|v| v.n_supporting as i32)
            .collect::<Vec<i32>>(),
        support_fraction = res
            .verdicts
            .iter()
            .map(|v| v.support_fraction)
            .collect::<Vec<f64>>(),
        p_up = res.verdicts.iter().map(|v| v.p_up).collect::<Vec<f64>>(),
        p_down = res.verdicts.iter().map(|v| v.p_down).collect::<Vec<f64>>(),
        coefficient = res
            .verdicts
            .iter()
            .map(|v| v.coefficient)
            .collect::<Vec<f64>>()
    );

    list!(
        shared_samples = res.shared_samples.clone().r_int_convert(),
        kept_features =
            List::from_values(res.kept_features.iter().map(|f| f.clone().r_int_convert())),
        mcp_cell_types =
            List::from_values(res.mcp_cell_types.iter().map(|t| t.clone().r_int_convert())),
        ws = List::from_values(res.ws.iter().map(|m| faer_to_r_matrix(m.as_ref()))),
        scores = List::from_values(res.scores.iter().map(|m| faer_to_r_matrix(m.as_ref()))),
        cca_scores = List::from_values(res.cca_scores.iter().map(|m| faer_to_r_matrix(m.as_ref()))),
        emp_p = faer_to_r_matrix(res.emp_p.as_ref()),
        pair_cor = faer_to_r_matrix(res.pair_cor.as_ref()),
        refit_fidelity = faer_to_r_matrix(res.refit_fidelity.as_ref()),
        verdicts = verdicts,
        permissive = signatures_to_r_list(&res.permissive),
        strict = signatures_to_r_list(&res.strict)
    )
}

////////////
// NEBULA //
////////////

/// Flatten the NEBULA fits into an R list
///
/// Shared by the single cell and meta cell entry points, which run the same
/// kernel and only differ in where the counts come from.
///
/// `coefficients` and `se` arrive row-major over `n_kept * n_coef` and go back
/// as matrices of genes x coefficients, so R never has to know the stride.
///
/// ### Params
///
/// * `res` - The `NebulaScRes` from the fit.
///
/// ### Returns
///
/// The results as an R list.
pub fn nebula_res_to_r_list(res: NebulaScRes) -> List {
    let n_kept = res.gene_idx.len();
    let n_coef = res.n_coef;

    let coefficients = RMatrix::new_matrix(n_kept, n_coef, |r, c| res.coefficients[r * n_coef + c]);
    let se = RMatrix::new_matrix(n_kept, n_coef, |r, c| res.se[r * n_coef + c]);

    // `cell_overdispersion_shrunk` is absent when the shrinkage was turned off.
    // `Nullable` keeps that an R `NULL` rather than a zero-length vector, which
    // would be indistinguishable from "every gene shrank to nothing".
    let shrunk = match res.cell_overdispersion_shrunk {
        Some(v) => Nullable::NotNull(v),
        None => Nullable::Null,
    };

    list!(
        gene_idx = res.gene_idx.r_int_convert(),
        coefficients = coefficients,
        se = se,
        subject_overdispersion = res.subject_overdispersion,
        cell_overdispersion = res.cell_overdispersion,
        cell_overdispersion_shrunk = shrunk,
        convergence = res.convergence,
        sigma_at_bound = res.sigma_at_bound,
        log_fc = res.log_fc,
        effect_se = res.effect_se,
        z = res.z,
        p_values = res.p_val,
        fdr = res.fdr
    )
}
