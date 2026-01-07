//! Module containing functions specifically designed for VERY large
//! data sets.

use ann_search_rs::*;
use extendr_api::List;
use faer::MatRef;
use std::time::Instant;

////////////////
// Neighbours //
////////////////

////////////
// Params //
////////////

/// KnnParams
///
/// ### Fields
///
/// **General**
///
/// * `quantisation_method` - Which of the quantisation methods to use. A choice
///   of `"SQ8"`, `"PQ"` or `"OPQ"`
/// * `ann_dist` - Approximate nearest neighbour distance measure. One of
///   `"euclidean"` or `"cosine"`.
/// * `k` - Number of neighbours to search.
/// * `n_centroids` - Number of centroids to use.
/// * `n_probes` - Number of centroids to probe.
/// * `n_query` - Number of queries to push into the GPU in one go.
#[derive(Clone, Debug)]
pub struct KnnParamsBigData {
    pub ann_dist: String,
    pub k: usize,
    pub n_centroids: Option<usize>,
    pub n_probes: Option<usize>,
    pub n_query: Option<usize>,
}

impl KnnParamsBigData {
    /// Generate KnnParamsBigData from an R list
    ///
    /// Should values not be found within the List, the parameters will default
    /// to sensible defaults based on heuristics.
    ///
    /// ### Params
    ///
    /// * `r_list` - The list with the kNN parameters.
    ///
    /// ### Returns
    ///
    /// The `KnnParamsBigData` with all parameters set.
    pub fn from_r_list(r_list: List) -> Self {
        let params_list = r_list.into_hashmap();

        let ann_dist = std::string::String::from(
            params_list
                .get("ann_dist")
                .and_then(|v| v.as_str())
                .unwrap_or("cosine"),
        );

        let k = params_list
            .get("k")
            .and_then(|v| v.as_integer())
            .unwrap_or(15) as usize;

        // ivf
        let n_centroids = params_list
            .get("n_centroids")
            .and_then(|v| v.as_integer())
            .map(|v| v as usize);

        let n_probes = params_list
            .get("n_probes")
            .and_then(|v| v.as_integer())
            .map(|v| v as usize);

        // gpu
        let n_query = params_list
            .get("n_query")
            .and_then(|v| v.as_integer())
            .map(|v| v as usize);

        Self {
            ann_dist,
            k,
            n_centroids,
            n_probes,
            n_query,
        }
    }
}

//////////
// Main //
//////////

/// Wrapper function to use quantised indices
///
/// ### Params
///
/// * `mat` - Matrix in which rows represent the samples and columns the
///   respective embeddings for that sample
/// * `knn_params` - The `KnnParamsBigData` with all of the parameters to use
///   for the quantised indices
/// * `seed` - Seed for reproducibility
/// * `verbose` - Controls verbosity.
///
/// ### Returns
///
/// The nearest neighbour indices
#[allow(clippy::too_many_arguments)]
pub fn generate_knn_quantised(
    mat: MatRef<f32>,
    knn_params: &KnnParamsBigData,
    seed: usize,
    verbose: bool,
) -> Vec<Vec<usize>> {
    let device: cubecl::wgpu::WgpuDevice = Default::default();

    let start_index = Instant::now();
    let ivf_gpu_idx = build_ivf_index_gpu::<f32, cubecl::wgpu::WgpuRuntime>(
        mat,
        knn_params.n_centroids,
        None,
        &knn_params.ann_dist,
        seed,
        verbose,
        device.clone(),
    );
    let end_index = start_index.elapsed();
    if verbose {
        println!("Generated IVF-GPU index: {:.2?}", end_index);
    }

    let start_search = Instant::now();
    let (indices, _) = query_ivf_index_gpu_self(
        &ivf_gpu_idx,
        knn_params.k + 1,
        knn_params.n_probes,
        knn_params.n_query,
        false,
        true,
    );
    let end_search = start_search.elapsed();

    if verbose {
        println!(
            "Identified approximate nearest neighbours via IVF-GPU: {:.2?}.",
            end_search
        );
    }

    let res: Vec<Vec<usize>> = indices
        .into_iter()
        .map(|mut v| {
            v.remove(0);
            v
        })
        .collect();

    res
}
