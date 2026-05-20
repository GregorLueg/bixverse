//! Helpers for the single cell part, especially around transforming structures
//! into R lists.

use bixverse_rs::prelude::*;
use bixverse_rs::single_cell::sc_analysis::fast_clusters::{
    FastLouvainGridResult, FastLouvainResults,
};
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
    let n_genes = results[0].mean.len();
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
/// * `knn_mat` - Samples x indices of the k-nearest neighbours (1-indexed!)
///
/// ### Returns
///
/// A `Vec<Vec<usize>>`
pub fn knn_indices_processing(knn_mat: RMatrix<i32>) -> Vec<Vec<usize>> {
    let ncol = knn_mat.ncols();
    let nrow = knn_mat.nrows();
    let data = knn_mat.data();

    (0..nrow)
        .map(|j| (0..ncol).map(|i| data[j + i * nrow] as usize).collect())
        .collect()
}

/// Process R KNN distances to the Rust variant
///
/// ### Params
///
/// * `knn_dist` - Samples x indices of the k-nearest neighbours (1-indexed!)
///
/// ### Returns
///
/// A `Vec<Vec<f32>>`
pub fn knn_distances_processing(knn_dist: RMatrix<f64>) -> Vec<Vec<f32>> {
    let ncol = knn_dist.ncols();
    let nrow = knn_dist.nrows();
    let data: Vec<f32> = knn_dist.data().iter().map(|x| *x as f32).collect();

    (0..nrow)
        .map(|j| (0..ncol).map(|i| data[j + i * nrow]).collect())
        .collect()
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
            .ok_or_else(|| Error::Other("missing 'k'".into()))?;
        let markers = CellTypeMarkers::from_r_list(element)?;

        res.push(markers);
    }

    Ok(res)
}
