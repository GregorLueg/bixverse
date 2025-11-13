use extendr_api::{Conversions, List};
use faer::{Mat, RowRef};
use indexmap::IndexSet;
use rayon::prelude::*;
use std::time::Instant;

use std::collections::HashSet;

use crate::assert_same_len;
use crate::core::base::linear_algebra::linear_regression;
use crate::core::base::stats::{calc_fdr, inv_logit, logit, z_scores_to_pval};
use crate::core::base::utils::rank_vector;
use crate::core::data::sparse_io::*;

////////////
// VISION //
////////////

/// Structure to store the indices of the SignatureGenes
///
/// ### Fields
///
/// * `positive` - The gene indices of the positive genes
/// * `negative` - The gene indices of the negative genes
#[derive(Clone, Debug)]
pub struct SignatureGenes {
    pub positive: Vec<usize>,
    pub negative: Vec<usize>,
}

impl SignatureGenes {
    /// Generate a SignatureGenes from an R list
    ///
    /// ### Params
    ///
    /// * `r_list` - An R list that is expected to have `"pos"` and `"neg"` with
    ///   0-index positions of the gene for this gene set.
    pub fn from_r_list(r_list: List) -> Self {
        let r_list = r_list.into_hashmap();

        let positive: Vec<usize> = r_list
            .get("pos")
            .and_then(|v| v.as_integer_vector())
            .unwrap_or_default()
            .iter()
            .map(|x| *x as usize)
            .collect();

        let negative: Vec<usize> = r_list
            .get("neg")
            .and_then(|v| v.as_integer_vector())
            .unwrap_or_default()
            .iter()
            .map(|x| *x as usize)
            .collect();

        Self { positive, negative }
    }
}

/// Helper function to transform an R gene set list to `Vec<SignatureGenes>`
///
/// ### Params
///
/// * `gs_list` - Initial R list with the gene sets for VISION
///
/// ### Returns
///
/// The vector of SignatureGenes
pub fn r_list_to_sig_genes(gs_list: List) -> extendr_api::Result<Vec<SignatureGenes>> {
    let mut gene_signatures: Vec<SignatureGenes> = Vec::with_capacity(gs_list.len());

    for i in 0..gs_list.len() {
        let r_obj = gs_list.elt(i)?;
        let gs_list_i = r_obj.as_list().ok_or_else(|| {
            extendr_api::Error::from(
                "The lists in the gs_list could not be converted. Please check!",
            )
        })?;
        gene_signatures.push(SignatureGenes::from_r_list(gs_list_i));
    }

    Ok(gene_signatures)
}

/////////////
// Helpers //
/////////////

/// Calculate VISION signature scores for a single cell
///
/// ### Params
///
/// * `cell` - The CsrCellChunk
/// * `signatures` - Slice of `SignatureGenes` to calculate the scores for
/// * `total_genes` - Total number of represented genes
///
/// ### Returns
///
/// A Vec<f32> with a score for each of the `SignatureGenes` that were supplied
fn calculate_vision_scores_for_cell(
    cell: &CsrCellChunk,
    signatures: &[SignatureGenes],
    total_genes: usize,
) -> Vec<f32> {
    // helper
    let get_expr = |gene_idx: usize| -> f32 {
        match cell.indices.binary_search(&(gene_idx as u16)) {
            Ok(pos) => cell.data_norm[pos].to_f32(),
            Err(_) => 0.0,
        }
    };

    // general cell statistics
    let sum: f32 = cell.data_norm.iter().map(|x| x.to_f32()).sum();
    let sum_sq: f32 = cell
        .data_norm
        .iter()
        .map(|x| {
            let v = x.to_f32();
            v * v
        })
        .sum();

    let mu_j = sum / total_genes as f32;
    let sigma_sq_j = (sum_sq / total_genes as f32) - (mu_j * mu_j);

    // score the signatures
    signatures
        .iter()
        .map(|sig| {
            let sum_pos: f32 = sig.positive.iter().map(|&idx| get_expr(idx)).sum();
            let sum_neg: f32 = sig.negative.iter().map(|&idx| get_expr(idx)).sum();

            let n = sig.positive.len() as f32;
            let m = sig.negative.len() as f32;
            let total = n + m;

            let s_j = (sum_pos - sum_neg) / total;
            let e_rj = ((n - m) / total) * mu_j;
            let var_rj = sigma_sq_j / total;

            if var_rj > 0.0 {
                (s_j - e_rj) / var_rj.sqrt()
            } else {
                0.0
            }
        })
        .collect()
}

/// Calculate Geary's C for a single signature (pathway)
///
/// This implements a modified approach akin to VISION.
///
/// ### Params
///
/// * `scores` - Ranked signature scores for all cells
/// * `knn_indices` - KNN indices matrix (cells x k)
/// * `knn_weights` - KNN weights matrix (cells x k)
///
/// ### Returns
///
/// Geary's C statistic
fn geary_c(scores: &[f64], knn_indices: &[Vec<usize>], knn_weights: &[Vec<f32>]) -> f64 {
    let n = scores.len();

    let mean: f64 = scores.iter().sum::<f64>() / n as f64;
    let variance: f64 = scores.iter().map(|x| (x - mean).powi(2)).sum::<f64>();

    if variance == 0.0 {
        return 0.0;
    }

    let mut numerator = 0.0;
    let mut total_weight = 0.0;

    // unsafe unchecked access in hot loop
    for (i, (indices, weights)) in knn_indices.iter().zip(knn_weights.iter()).enumerate() {
        let xi = unsafe { *scores.get_unchecked(i) };
        for (&j, &w) in indices.iter().zip(weights.iter()) {
            let xj = unsafe { *scores.get_unchecked(j) };
            numerator += w as f64 * (xi - xj).powi(2);
            total_weight += w as f64;
        }
    }

    let norm = 2.0 * total_weight * variance / (n as f64 - 1.0);

    numerator / norm
}

/// Calculate KNN weights using exponential kernel
///
/// ### Params
///
/// * `knn_indices` - KNN indices (cells x k)
/// * `knn_distances` - KNN squared distances (cells x k) - use Euclidean here!
///
/// ### Returns
///
/// KNN weights matrix (cells x k)
fn calc_knn_weights(knn_indices: &[Vec<usize>], knn_distances: &[Vec<f32>]) -> Vec<Vec<f32>> {
    knn_indices
        .par_iter() // Parallel iterator
        .zip(knn_distances.par_iter())
        .map(|(indices, distances)| {
            if distances.is_empty() {
                return vec![];
            }

            let sigma_sq = distances.last().copied().unwrap_or(1.0);

            if sigma_sq == 0.0 {
                return vec![1.0; indices.len()];
            }

            distances
                .iter()
                .map(|&d_sq| (-d_sq / sigma_sq).exp())
                .collect()
        })
        .collect()
}

//////////
// Main //
//////////

/// Calculate VISION signature scores across a set of cells
///
/// ### Params
///
/// * `f_path` -  File path to the cell-based binary file.
/// * `signatures` - Slice of `SignatureGenes` to calculate the scores for
/// * `cells_to_keep` - Vector of indices with the cells to keep.
/// * `verbose` - Controls verbosity of the function
///
/// ### Returns
///
/// A Vec<Vec<f32>> with cells x scores per gene set (pair)
pub fn calculate_vision(
    f_path: &str,
    gene_signs: &[SignatureGenes],
    cells_to_keep: &[usize],
    verbose: bool,
) -> Vec<Vec<f32>> {
    let start_read = Instant::now();
    let reader = ParallelSparseReader::new(f_path).unwrap();
    let no_genes = reader.get_header().total_genes;
    let cell_chunks: Vec<CsrCellChunk> = reader.read_cells_parallel(cells_to_keep);
    let end_read = start_read.elapsed();

    if verbose {
        println!("Loaded in data: {:.2?}", end_read);
    }

    let start_signatures = Instant::now();
    let signature_scores: Vec<Vec<f32>> = cell_chunks
        .par_iter()
        .map(|chunk| calculate_vision_scores_for_cell(chunk, gene_signs, no_genes))
        .collect();
    let end_signatures = start_signatures.elapsed();

    if verbose {
        println!("Calculated VISION scores: {:.2?}", end_signatures);
    }

    signature_scores // cells x signatures
}

/// Calculate VISION signature scores across a set of cells (streaming)
///
/// The streaming version of the function.
///
/// ### Params
///
/// * `f_path` -  File path to the cell-based binary file.
/// * `signatures` - Slice of `SignatureGenes` to calculate the scores for
/// * `cells_to_keep` - Vector of indices with the cells to keep.
/// * `verbose` - Controls verbosity of the function
///
/// ### Returns
///
/// A Vec<Vec<f32>> with cells x scores per gene set (pair)
pub fn calculate_vision_streaming(
    f_path: &str,
    gene_signs: &[SignatureGenes],
    cells_to_keep: &[usize],
    verbose: bool,
) -> Vec<Vec<f32>> {
    const CHUNK_SIZE: usize = 50000;

    let total_chunks = cells_to_keep.len().div_ceil(CHUNK_SIZE);
    let reader = ParallelSparseReader::new(f_path).unwrap();
    let no_genes = reader.get_header().total_genes;

    let mut all_results: Vec<Vec<f32>> = Vec::with_capacity(cells_to_keep.len());

    for (chunk_idx, cell_indices_chunk) in cells_to_keep.chunks(CHUNK_SIZE).enumerate() {
        let start_chunk = Instant::now();

        let cell_chunks = reader.read_cells_parallel(cell_indices_chunk);

        let chunk_scores: Vec<Vec<f32>> = cell_chunks
            .par_iter()
            .map(|chunk| calculate_vision_scores_for_cell(chunk, gene_signs, no_genes))
            .collect();

        all_results.extend(chunk_scores);

        if verbose {
            let elapsed = start_chunk.elapsed();
            let pct_complete = ((chunk_idx + 1) as f32 / total_chunks as f32) * 100.0;
            println!(
                "Processing chunk {} out of {} (took {:.2?}, completed {:.1}%)",
                chunk_idx + 1,
                total_chunks,
                elapsed,
                pct_complete
            );
        }
    }

    all_results
}

/// Calculate VISION local autocorrelation scores
///
/// ### Params
///
/// * `pathway_scores` - Vector representing the actual vision scores:
///   `cells x pathways
/// * `random_scores_by_cluster` - Vector representing the random scores by
///   cluster: `clusters -> (cells x sigs)`
/// * `cluster_membership` - Vector representing to which cluster a given
///   gene set belongs.
/// * `knn_indices` - KNN indices from embedding (cells x k)
/// * `knn_distances` - KNN squared distances (cells x k)
/// * `verbose` - Print progress
///
/// ### Returns
///
/// Tuple of (consistency_scores, p_values) for each pathway
pub fn calc_autocorr_with_clusters(
    pathway_scores: &[Vec<f32>],
    random_scores_by_cluster: &[Vec<Vec<f32>>],
    cluster_membership: &[usize],
    knn_indices: Vec<Vec<usize>>,
    knn_distances: Vec<Vec<f32>>,
    verbose: bool,
) -> (Vec<f64>, Vec<f64>) {
    let start = Instant::now();

    let knn_weights = calc_knn_weights(&knn_indices, &knn_distances);

    if verbose {
        println!("Computed KNN weights: {:.2?}", start.elapsed());
    }

    let n_pathways = pathway_scores[0].len();

    // calculate Geary's C for actual pathways
    let pathway_consistency: Vec<f64> = (0..n_pathways)
        .into_par_iter()
        .map(|pathway_idx| {
            let scores: Vec<f32> = pathway_scores
                .iter()
                .map(|cell| cell[pathway_idx])
                .collect();

            let ranks = rank_vector(&scores);
            let c = geary_c(&ranks, &knn_indices, &knn_weights);
            1.0 - c
        })
        .collect();

    if verbose {
        println!("Calculated pathway consistency: {:.2?}", start.elapsed());
    }

    let cluster_bg_consistency: Vec<Vec<f64>> = random_scores_by_cluster
        .iter()
        .map(|cluster_scores| {
            let n_random = cluster_scores[0].len();

            // Single level of parallelism
            (0..n_random)
                .into_par_iter() // Only this one is parallel
                .map(|sig_idx| {
                    let scores: Vec<f32> =
                        cluster_scores.iter().map(|cell| cell[sig_idx]).collect();

                    let ranks = rank_vector(&scores);
                    let c = geary_c(&ranks, &knn_indices, &knn_weights);
                    1.0 - c
                })
                .collect()
        })
        .collect();

    // calculate p vals
    let p_vals: Vec<f64> = pathway_consistency
        .iter()
        .enumerate()
        .map(|(i, &fg_c)| {
            let cluster_idx = cluster_membership[i];
            let bg_dist = &cluster_bg_consistency[cluster_idx];

            let n = bg_dist.len();
            let num_greater_equal = bg_dist.iter().filter(|&&bg_c| bg_c >= fg_c).count();

            (num_greater_equal + 1) as f64 / (n + 1) as f64
        })
        .collect();

    if verbose {
        println!("Calculated p-values: {:.2?}", start.elapsed());
    }

    (pathway_consistency, p_vals)
}

/////////////
// Hotspot //
/////////////

////////////
// Params //
////////////

/// HotSpot parameters
///
/// ### Fields
///
/// **General**
///
/// * `knn_method` - Which of the kNN methods to use. One of `"annoy"`, `"hnsw"`
///   or `"nndescent"`.
/// * `ann_dist` - Approximate nearest neighbour distance measure. One of
///   `"euclidean"` or `"cosine"`.
/// * `k` - Number of neighbours to search
///
/// **Annoy**
///
/// * `n_tree` - Number of trees for the generation of the index
/// * `search_budget` - Search budget during querying
///
/// **NN Descent**
///
/// * `max_iter` - Maximum iterations for the algorithm
/// * `rho` - Sampling rate for the algorithm
/// * `delta` - Early termination criterium
///
/// **HotSpot**
///
/// * `model` - The model to use for modelling the GEX. Choice of
///   `"danb"`, `"bernoulli"` or `"normal"`.
/// * `normalise` - Shall the data be normalised.
pub struct HotSpotParams {
    // general params
    pub knn_method: String,
    pub ann_dist: String,
    pub k: usize,
    // annoy params
    pub n_tree: usize,
    pub search_budget: usize,
    // nn descent params
    pub max_iter: usize,
    pub rho: f32,
    pub delta: f32,
    // hotspot parameters
    pub model: String,
    pub normalise: bool,
}

impl HotSpotParams {
    /// Generate HotSpotParams from an R list
    ///
    /// Should values not be found within the List, the parameters will default
    /// to sensible defaults based on heuristics.
    ///
    /// ### Params
    ///
    /// * `r_list` - The list with the Boost parameters.
    ///
    /// ### Returns
    ///
    /// The `HotSpotParams` with all parameters set.
    pub fn from_r_list(r_list: List) -> Self {
        let params_list = r_list.into_hashmap();

        // general
        let knn_method = std::string::String::from(
            params_list
                .get("knn_algorithm")
                .and_then(|v| v.as_str())
                .unwrap_or("annoy"),
        );

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

        // annoy
        let n_tree = params_list
            .get("n_tree")
            .and_then(|v| v.as_integer())
            .unwrap_or(100) as usize;

        let search_budget = params_list
            .get("search_budget")
            .and_then(|v| v.as_integer())
            .unwrap_or(100) as usize;

        // nn descent
        let max_iter = params_list
            .get("nn_max_iter")
            .and_then(|v| v.as_integer())
            .unwrap_or(25) as usize;

        let rho = params_list
            .get("rho")
            .and_then(|v| v.as_real())
            .unwrap_or(1.0) as f32;

        let delta = params_list
            .get("delta")
            .and_then(|v| v.as_real())
            .unwrap_or(0.001) as f32;

        // hotspot
        let model = std::string::String::from(
            params_list
                .get("model")
                .and_then(|v| v.as_str())
                .unwrap_or("normal"),
        );

        let normalise = params_list
            .get("normalise")
            .and_then(|v| v.as_bool())
            .unwrap_or(true);

        Self {
            knn_method,
            ann_dist,
            k,
            n_tree,
            search_budget,
            max_iter,
            rho,
            delta,
            model,
            normalise,
        }
    }
}

/////////////
// Helpers //
/////////////

#[derive(Debug, Clone)]
pub enum GexModel {
    /// Use depth-adjusted negative binomial model
    DephAdjustNegBinom,
    /// Uses Bernoulli distribution to model prediction probability
    Bernoulli,
    /// Use depth-adjusted normal model
    Normal,
}

/// Parse the model to use gene expression
///
/// ### Params
///
/// * `s` - Type of model to use the model
///
/// ### Returns
///
/// Option of the GexModel to use (some not yet implemented)
pub fn parse_gex_model(s: &str) -> Option<GexModel> {
    match s.to_lowercase().as_str() {
        "danb" => Some(GexModel::DephAdjustNegBinom),
        "bernoulli" => Some(GexModel::Bernoulli),
        "normal" => Some(GexModel::Normal),
        _ => None,
    }
}

/// Structure for the gene results
///
/// ### Fields
///
/// * `gene_idx` - Gene index of the analysed gene
/// * `c` - Geary's C statistic for this gene
/// * `z` - Z-score for this gene
/// * `pval` - P-value based on the Z-score
/// * `fdr` - False discovery corrected pvals
#[derive(Debug, Clone)]
pub struct HotSpotGeneRes {
    pub gene_idx: Vec<usize>,
    pub c: Vec<f64>,
    pub z: Vec<f64>,
    pub pval: Vec<f64>,
    pub fdr: Vec<f64>,
}

/// Structure for pair-wise correlations
///
/// ### Fields
///
/// * `cor` - Symmetric matrix with cor coefficients (N_genes x N_genes)
/// * `z_scores` - Symmetric matrix with Z scores (N_genex x N_genes)
#[derive(Debug, Clone)]
pub struct HotSpotPairRes {
    pub cor: Mat<f32>,
    pub z_scores: Mat<f32>,
}

/// Compute momentum weights
///
/// Calculates the expected value (EG) and expected squared value (EG2) of the
/// local covariance statistic under the null hypothesis of no spatial
/// autocorrelation.
///
/// ### Params
///
/// * `mu` - Mean expression values for each cell
/// * `x2` - Second moment (variance + mean²) for each cell
/// * `neighbours` - Neighbour indices for each cell
/// * `weights` - Edge weights for each neighbour connection
///
/// ### Returns
///
/// Tuple of (EG, EG2) where:
///
/// - EG: Expected value of the spatial covariance statistic
/// - EG2: Expected value of the squared spatial covariance statistic
fn compute_moments_weights(
    mu: &[f32],
    x2: &[f32],
    neighbours: &[Vec<usize>],
    weights: &[Vec<f32>],
) -> (f32, f32) {
    let n = neighbours.len();
    let mu_sq: Vec<f32> = mu.iter().map(|&m| m * m).collect();

    let mut eg = 0_f32;
    let mut t1 = vec![0_f32; n];
    let mut t2 = vec![0_f32; n];

    for i in 0..n {
        let mu_i = mu[i];

        for (k, &j) in neighbours[i].iter().enumerate() {
            let wij = weights[i][k];
            let mu_j = mu[j];

            eg += wij * mu_i * mu_j;

            t1[i] += wij * mu_j;
            let wij_sq = wij * wij;
            t2[i] += wij_sq * mu_j * mu_j;

            // Add these back:
            t1[j] += wij * mu_i;
            t2[j] += wij_sq * mu_i * mu_i;
        }
    }

    let mut eg2 = 0_f32;

    for i in 0..n {
        eg2 += (x2[i] - mu_sq[i]) * (t1[i] * t1[i] - t2[i]);
    }

    for i in 0..n {
        let x2_i = x2[i];
        let mu_sq_i = mu_sq[i];

        for (k, &j) in neighbours[i].iter().enumerate() {
            let wij = weights[i][k];
            eg2 += wij * wij * (x2_i * x2[j] - mu_sq_i * mu_sq[j]);
        }
    }

    eg2 += eg * eg;

    (eg, eg2)
}

/// Remove redundancy in bidirectional edge weights
///
/// Consolidates weights from bidirectional edges by accumulating both
/// directions into the lower-indexed node's edge and zeroing the
/// higher-indexed node's reciprocal edge.
///
/// ### Params
///
/// * `neighbours` - Neighbour indices for each node
/// * `weights` - Edge weights for each neighbour connection
///
/// ### Returns
///
/// Modified weights with redundant edges zeroed
fn make_weights_non_redundant(neighbours: &[Vec<usize>], weights: &[Vec<f32>]) -> Vec<Vec<f32>> {
    let mut w_no_redundant = weights.to_vec();

    for i in 0..neighbours.len() {
        for k in 0..neighbours[i].len() {
            let j = neighbours[i][k];

            if j < i {
                continue;
            }

            // check if j has i as a neighbour
            for k2 in 0..neighbours[j].len() {
                if neighbours[j][k2] == i {
                    let w_ji = w_no_redundant[j][k2];
                    w_no_redundant[j][k2] = 0.0;
                    w_no_redundant[i][k] += w_ji;
                    break;
                }
            }
        }
    }

    w_no_redundant
}

/// Compute node degree from edge weights
///
/// Calculates the degree (sum of incident edge weights) for each node.
/// Each edge contributes to the degree of both its endpoints.
///
/// ### Params
///
/// * `neighbours` - Neighbour indices for each node
/// * `weights` - Edge weights for each neighbour connection
///
/// ### Returns
///
/// Vector of degree values for each node
fn compute_node_degree(neighbours: &[Vec<usize>], weights: &[Vec<f32>]) -> Vec<f32> {
    let mut d = vec![0.0_f32; neighbours.len()];

    for i in 0..neighbours.len() {
        for k in 0..neighbours[i].len() {
            let j = neighbours[i][k];
            let w_ij = weights[i][k];

            d[i] += w_ij;
            d[j] += w_ij;
        }
    }

    d
}

/// Compute local covariance using edge weights
///
/// Calculates the weighted local covariance statistic for spatial
/// autocorrelation. This is the numerator of Geary's C statistic.
///
/// ### Params
///
/// * `vals` - Gene expression values for each cell
/// * `neighbours` - Neighbour indices for each cell
/// * `weights` - Edge weights for each neighbour connection
///
/// ### Returns
///
/// The local covariance statistic
fn local_cov_weights(vals: &[f32], neighbours: &[Vec<usize>], weights: &[Vec<f32>]) -> f32 {
    let mut out = 0.0;

    for i in 0..vals.len() {
        let xi = vals[i];
        if xi == 0.0 {
            continue;
        }

        for (k, &j) in neighbours[i].iter().enumerate() {
            let xj = vals[j];
            let wij = weights[i][k];

            if xj != 0.0 && wij != 0.0 {
                out += xi * xj * wij;
            }
        }
    }

    out
}

/// Compute maximum possible local covariance
///
/// Calculates the theoretical maximum value of the local covariance statistic
/// given the node degrees and expression values. Used to normalise Geary's C.
///
/// ### Params
///
/// * `node_degrees` - Sum of edge weights for each node
/// * `vals` - Gene expression values for each cell
///
/// ### Returns
///
/// Maximum possible local covariance
fn compute_local_cov_max(node_degrees: &[f32], vals: &[f32]) -> f32 {
    let mut tot = 0.0;
    for i in 0..node_degrees.len() {
        tot += node_degrees[i] * vals[i] * vals[i];
    }
    tot / 2.0
}

/// Center (Z-score) the values
///
/// Transforms values to have zero means and unit variance of one using the
/// provided stats.
///
/// ### Params
///
/// * `vals` - Mutable reference to the values to scale
/// * `mu` - The mean values
/// * `var` - The variance of the values
fn center_values(vals: &mut [f32], mu: &[f32], var: &[f32]) {
    assert_same_len!(vals, mu, var);

    for i in 0..vals.len() {
        vals[i] = (vals[i] - mu[i]) / var[i].sqrt();
    }
}

//////////////////
// Corr helpers //
//////////////////

/// Compute local covariance for gene pairs
///
/// Test statistic for local pairwise autocorrelation. Calculates the weighted
/// covariance between two genes across neighbouring cells.
///
/// ### Params
///
/// * `x` - RowRef for first gene.
/// * `y` - RowRef for second gene.
/// * `neighbours` - Neighbour indices for each cell
/// * `weights` - Edge weights for each neighbour connection
///
/// ### Returns
///
/// Local covariance statistic
fn local_cov_pair(
    x: RowRef<f32>,
    y: RowRef<f32>,
    neighbours: &[Vec<usize>],
    weights: &[Vec<f32>],
) -> f32 {
    let mut out = 0.0;

    for i in 0..x.ncols() {
        let xi = x[i];
        let yi = y[i];
        if xi == 0.0 && yi == 0.0 {
            continue;
        }
        for k in 0..neighbours[i].len() {
            let j = neighbours[i][k];
            let w_ij = weights[i][k];

            let xj = x[j];
            let yj = y[j];

            out += w_ij * (xi * yj + yi * xj) / 2.0;
        }
    }

    out
}

/// Compute local covariance for gene pairs
///
/// Test statistic for local pairwise autocorrelation. Calculates the weighted
/// covariance between two genes across neighbouring cells.
///
/// ### Params
///
/// * `x` - Slice for first gene.
/// * `y` - Slice for second gene.
/// * `neighbours` - Neighbour indices for each cell
/// * `weights` - Edge weights for each neighbour connection
///
/// ### Returns
///
/// Local covariance statistic
fn local_cov_pair_vec(
    x: &[f32],
    y: &[f32],
    neighbours: &[Vec<usize>],
    weights: &[Vec<f32>],
) -> f32 {
    neighbours
        .iter()
        .zip(weights.iter())
        .enumerate()
        .map(|(i, (neighs, ws))| {
            let xi = x[i];
            let yi = y[i];
            neighs
                .iter()
                .zip(ws.iter())
                .map(|(&j, &w)| {
                    let xj = x[j];
                    let yj = y[j];
                    w * (xi * yj + yi * xj) / 2.0
                })
                .sum::<f32>()
        })
        .sum()
}

/// Compute conditional EG2 for correlation
///
/// Calculates the expected value of G_square for the conditional correlation
/// statistic, assuming standardised variables.
///
/// ### Params
///
/// * `x` - Standardised expression values for a gene
/// * `neighbors` - Neighbour indices for each cell
/// * `weights` - Edge weights for each neighbour connection
///
/// ### Returns
///
/// Expected value of G²
fn conditional_eg2(x: &[f32], neighbours: &[Vec<usize>], weights: &[Vec<f32>]) -> f32 {
    let n = neighbours.len();

    let mut t1x = vec![0_f32; n];

    for i in 0..n {
        for k in 0..neighbours[i].len() {
            let j = neighbours[i][k];
            let wij = weights[i][k];

            if wij == 0.0 {
                continue;
            }

            t1x[i] += wij * x[j];
            t1x[j] += wij * x[i];
        }
    }

    t1x.iter().map(|&t| t * t).sum()
}

/// Compute maximum possible pairwise local covariance
///
/// Calculates the theoretical maximum for pairwise correlation normalisation.
///
/// ### Params
///
/// * `node_degrees` - Sum of edge weights for each node
/// * `counts` - Centred gene expression matrix (genes × cells)
///
/// ### Returns
///
/// Matrix of maximum covariances (genes × genes)
fn compute_local_cov_pairs_max(node_degrees: &[f32], counts: &Mat<f32>) -> Mat<f32> {
    let n_genes = counts.nrows();

    let gene_maxs: Vec<f32> = (0..n_genes)
        .into_par_iter()
        .map(|i| {
            let row = counts.row(i);
            let row_vec = row.iter().copied().collect::<Vec<f32>>();
            compute_local_cov_max(node_degrees, &row_vec)
        })
        .collect();

    let values: Vec<f32> = (0..n_genes * n_genes)
        .into_par_iter()
        .map(|idx| {
            let i = idx / n_genes;
            let j = idx % n_genes;
            (gene_maxs[i] + gene_maxs[j]) / 2.0
        })
        .collect();

    let mut result = Mat::zeros(n_genes, n_genes);
    for (idx, &val) in values.iter().enumerate() {
        result[(idx / n_genes, idx % n_genes)] = val;
    }

    result
}

/// Centre gene counts for correlation computation
///
/// Standardises gene expression using the specified model, transforming to
/// zero mean and unit variance.
///
/// ### Params
///
/// * `gene` - Reference to gene expression data
/// * `umi_counts` - Total UMI counts per cell
/// * `n_cells` - Number of cells
/// * `model` - Statistical model to use
///
/// ### Returns
///
/// Vector of centred expression values
fn create_centered_counts_gene(
    gene: &CscGeneChunk,
    umi_counts: &[f32],
    n_cells: usize,
    model: &GexModel,
) -> Vec<f32> {
    let (mu, var, _) = match model {
        GexModel::DephAdjustNegBinom => danb_model(gene, umi_counts, n_cells),
        GexModel::Bernoulli => bernoulli_model(gene, umi_counts, n_cells),
        GexModel::Normal => normal_model(gene, umi_counts, n_cells),
    };

    let mut vals = vec![0_f32; n_cells];
    for (&idx, &val) in gene.indices.iter().zip(&gene.data_raw) {
        vals[idx as usize] = val as f32;
    }

    center_values(&mut vals, &mu, &var);

    vals
}

////////////////
// DANB model //
////////////////

/// Depth-adjusted negative binomial (DANB) model
///
/// Fits a negative binomial distribution to gene expression data, adjusting
/// for sequencing depth differences between cells.
///
/// ### Params
///
/// * `gene` - Reference to the CscGeneChunk on which to apply the model.
/// * `umi_counts` - Slice of the UMI counts across these cells (i.e.,
///   sequencing depth).
/// * `n_cells` - Total number of cells
///
/// ### Returns
///
/// Tuple of (mu, var, x2) where:
/// - mu: Mean expression for each cell
/// - var: Variance for each cell
/// - x2: Second moment (var + mu²) for each cell
fn danb_model(
    gene: &CscGeneChunk,
    umi_counts: &[f32],
    n_cells: usize,
) -> (Vec<f32>, Vec<f32>, Vec<f32>) {
    let n = n_cells as f32;
    let total: f32 = umi_counts.iter().sum();

    let tj: f32 = gene.data_raw.iter().map(|&x| x as f32).sum();

    let mu = umi_counts
        .iter()
        .map(|&ti| tj * ti / total)
        .collect::<Vec<f32>>();

    let mut sum_sq = 0_f32;
    for (&idx, &val) in gene.indices.iter().zip(&gene.data_raw) {
        let diff = val as f32 - mu[idx as usize];
        sum_sq += diff * diff;
    }

    for i in 0..n_cells {
        if !gene.indices.contains(&(i as u32)) {
            sum_sq += mu[i] * mu[i];
        }
    }
    let vv = sum_sq / (n - 1.0);

    let tis_sq_sum: f32 = umi_counts.iter().map(|&ti| ti.powi(2)).sum();
    let mut size = ((tj * tj) / total) * (tis_sq_sum / total) / ((n - 1.0) * vv - tj);

    let min_size = 1e-10;
    if size < 0.0 {
        size = 1e9;
    }
    if size < min_size && size >= 0.0 {
        size = min_size;
    }

    let var: Vec<f32> = mu.iter().map(|&m| m * (1.0 + m / size)).collect();
    let x2: Vec<f32> = var.iter().zip(&mu).map(|(&v, &m)| v + m * m).collect();

    (mu, var, x2)
}

/////////////////////
// Bernoulli model //
/////////////////////

/// Bin gene detections by UMI count bins
///
/// Calculates the detection rate within each bin, applying Laplace smoothing
/// to handle edge cases (0% or 100% detection).
///
/// ### Params
///
/// * `detected_gene` - Binary detection indicators (0 or 1) for each cell
/// * `umi_count_bins` - Bin assignment for each cell
/// * `n_bins` - Total number of bins
///
/// ### Returns
///
/// Vector of detection rates per bin (with Laplace smoothing)
fn bin_gene_detection(detected_gene: &[f32], umi_count_bins: &[usize], n_bins: usize) -> Vec<f32> {
    let mut bin_detects = vec![0_f32; n_bins];
    let mut bin_totals = vec![0_f32; n_bins];

    for i in 0..detected_gene.len() {
        let bin_i = umi_count_bins[i];
        bin_detects[bin_i] += detected_gene[i];
        bin_totals[bin_i] += 1.0;
    }

    // laplace smoothing
    bin_detects
        .iter()
        .zip(&bin_totals)
        .map(|(&d, &t)| (d + 1.0) / (t + 2.0))
        .collect()
}

/// Quantile-based binning with duplicate edge handling
///
/// Generates quantile-based bins from data, dropping duplicate bin edges when
/// they would result in empty bins.
///
/// ### Params
///
/// * `data` - Input data to bin
/// * `n_bins` - Target number of bins
///
/// ### Returns
///
/// Tuple of (bin_assignments, bin_edges) where:
/// - bin_assignments: Vector of bin indices for each data point
/// - bin_edges: Vector of bin edge values (length = n_bins + 1)
fn quantile_cut(data: &[f32], n_bins: usize) -> (Vec<usize>, Vec<f32>) {
    let mut data_sorted = data.to_vec();
    data_sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());

    let n = data_sorted.len();
    let mut edges = vec![data_sorted[0]];

    for i in 1..n_bins {
        let idx = (i * n) / n_bins;
        let value = data_sorted[idx.min(n - 1)];

        if value > *edges.last().unwrap() {
            edges.push(value);
        }
    }

    let max_val = data_sorted[n - 1];
    if max_val > *edges.last().unwrap() {
        edges.push(max_val + 1e-6);
    } else {
        *edges.last_mut().unwrap() += 1e-6;
    }

    let n_actual_bins = edges.len() - 1;

    // binary search is faster here...
    let bin_assignments: Vec<usize> = data
        .iter()
        .map(|&x| {
            edges
                .partition_point(|&edge| edge <= x)
                .saturating_sub(1)
                .min(n_actual_bins - 1)
        })
        .collect();

    (bin_assignments, edges)
}

/// Bernoulli model for gene expression
///
/// Models the probability of detecting gene expression using a Bernoulli
/// distribution. Fits a logistic regression model on binned UMI counts to
/// predict detection probability.
///
/// ### Params
///
/// * `gene` - Reference to the CscGeneChunk containing gene expression data
/// * `umi_counts` - Total UMI counts per cell
/// * `n_cells` - Total number of cells
///
/// ### Returns
///
/// Tuple of (mu, var, x2) where:
/// - mu: Detection probability for each cell
/// - var: Variance (p * (1-p)) for each cell
/// - x2: Second moment (equal to mu for Bernoulli)
fn bernoulli_model(
    gene: &CscGeneChunk,
    umi_counts: &[f32],
    n_cells: usize,
) -> (Vec<f32>, Vec<f32>, Vec<f32>) {
    const N_BIN_TARGET: usize = 30;

    let mut detected_gene = vec![0_f32; n_cells];
    for idx in &gene.indices {
        detected_gene[*idx as usize] = 1.0;
    }

    let log_umi: Vec<f32> = umi_counts
        .iter()
        .map(|&x| if x > 0.0 { x.log10() } else { 0.0 })
        .collect();

    let (umi_count_bins, bin_edges) = quantile_cut(&log_umi, N_BIN_TARGET);
    let n_bins = bin_edges.len() - 1;

    let bin_centers: Vec<f32> = (0..n_bins)
        .map(|i| (bin_edges[i] + bin_edges[i + 1]) / 2.0)
        .collect();

    let bin_detects = bin_gene_detection(&detected_gene, &umi_count_bins, n_bins);

    let lbin_detects: Vec<f32> = bin_detects.iter().map(|&p| logit(p)).collect();
    let coef = linear_regression(&bin_centers, &lbin_detects);

    let mu: Vec<f32> = log_umi
        .iter()
        .map(|&log_u| inv_logit(coef.0 + coef.1 * log_u))
        .collect();

    let var: Vec<f32> = mu.iter().map(|&p| p * (1.0 - p)).collect();
    let x2: Vec<f32> = mu.clone();

    (mu, var, x2)
}

//////////////////
// Normal model //
//////////////////

/// Normal model for gene expression
///
/// Simplest model just using the normalised counts in the data.
///
/// ### Params
///
/// * `gene` - Reference to the CscGeneChunk containing gene expression data
/// * `n_cells` - Total number of cells
///
/// ### Returns
///
/// Tuple of (mu, var, x2) where:
/// - mu: Mean expression for each cell (from linear regression)
/// - var: Residual variance (constant across cells)
/// - x2: Second moment (var + mu²) for each cell
fn normal_model(
    gene: &CscGeneChunk,
    umi_counts: &[f32],
    n_cells: usize,
) -> (Vec<f32>, Vec<f32>, Vec<f32>) {
    let mut gene_raw = vec![0_f32; n_cells];
    for (&idx, &val) in gene.indices.iter().zip(&gene.data_raw) {
        gene_raw[idx as usize] = val as f32;
    }

    // Fit linear regression: expression ~ log(umi_counts)
    let log_umi: Vec<f32> = umi_counts
        .iter()
        .map(|&x| if x > 0.0 { x.ln() } else { 0.0 })
        .collect();

    let (intercept, slope) = linear_regression(&log_umi, &gene_raw);

    // Cell-specific mu from regression
    let mu: Vec<f32> = log_umi.iter().map(|&x| intercept + slope * x).collect();

    // Residual variance (constant across cells)
    let residuals_sq: f32 = gene_raw
        .iter()
        .zip(&mu)
        .map(|(&obs, &pred)| (obs - pred).powi(2))
        .sum();
    let var_val = residuals_sq / (n_cells as f32 - 2.0);

    let var = vec![var_val; n_cells];
    let x2: Vec<f32> = mu.iter().map(|&m| var_val + m * m).collect();

    (mu, var, x2)
}

//////////
// Main //
//////////

/// HotSpot structure
///
/// Main structure for computing spatial autocorrelation and gene <> gene
/// correlations in spatially-resolved transcriptomics data.
///
/// ### Fields
///
/// * `f_path_gene` - File path to the gene-based binary file.
/// * `f_path_cell` - File path to the cell-based binary file.
/// * `neigbours` - Slice if the indices of the cells to include in this
///   analysis.
/// * `weights` - Slice of the distances to the neighbours of a given cell.
/// * `cells_to_keep` - Slice of cells to analyse/keep in this analysis.
/// * `node_degrees` - Pre-computed node-degree for each cell based on the
///   weights.
/// * `umi_counts` - Optional vector with the total UMI counts per cell
/// * `wtot2` -
/// * `n_cells` - Total number of cells analysed in the experiment.
#[derive(Clone, Debug)]
pub struct Hotspot<'a> {
    f_path_gene: String,
    f_path_cell: String,
    neighbours: &'a [Vec<usize>],
    weights: Vec<Vec<f32>>,
    cells_to_keep: &'a [usize],
    node_degrees: Vec<f32>,
    umi_counts: Option<Vec<f32>>,
    wtot2: f32,
    n_cells: usize,
}

impl<'a> Hotspot<'a> {
    /// Initialise a new instance
    ///
    /// ### Params
    ///
    /// * `f_path_gene` - File path to the gene-based binary file.
    /// * `f_path_cell` - File path to the cell-based binary file.
    /// * `cells_to_keep` - Slice if the indices of the cells to include in this
    ///   analysis.
    /// * `neighbours` - Slice of the indices of the neighbours of the given
    ///   cell.
    /// * `weights` - Slice of the distances to the neighbours of a given cell.
    ///
    /// ### Return
    ///
    /// Initialised `HotSpot` class.
    pub fn new(
        f_path_gene: String,
        f_path_cell: String,
        cells_to_keep: &'a [usize],
        neighbours: &'a [Vec<usize>],
        weights: &mut [Vec<f32>],
    ) -> Self {
        let n_cells = neighbours.len();

        let weights = make_weights_non_redundant(neighbours, weights);

        let node_degrees = compute_node_degree(neighbours, &weights);

        let wtot2: f32 = weights.iter().flatten().map(|&w| w * w).sum();

        Self {
            f_path_gene,
            f_path_cell,
            neighbours,
            weights,
            cells_to_keep,
            node_degrees,
            umi_counts: None,
            wtot2,
            n_cells,
        }
    }

    /// Compute spatial autocorrelation for all specified genes
    ///
    /// Calculates Geary's C statistic and Z-scores for spatial autocorrelation
    /// across the specified genes.
    ///
    /// ### Params
    ///
    /// * `gene_indices` - Indices of genes to analyse
    /// * `model` - Statistical model to use ("danb", "bernoulli", or "normal")
    /// * `centered` - Whether to centre the data before computing statistics
    /// * `verbose` - Whether to print progress information
    ///
    /// ### Returns
    ///
    /// `Result<HotSpotGeneRes>` with gene indices, Geary's C, Z-scores, derived
    /// p-values and FDR.
    pub fn compute_all_genes(
        &mut self,
        gene_indices: &[usize],
        model: &str,
        centered: bool,
        verbose: bool,
    ) -> Result<HotSpotGeneRes, String> {
        let gex_model =
            parse_gex_model(model).ok_or_else(|| format!("Invalid model type: {}", model))?;

        self.populate_umi_counts();

        let cell_set: IndexSet<u32> = self.cells_to_keep.iter().map(|&x| x as u32).collect();

        let start_reading = Instant::now();

        let reader = ParallelSparseReader::new(&self.f_path_gene).unwrap();
        let mut gene_chunks: Vec<CscGeneChunk> = reader.read_gene_parallel(gene_indices);

        gene_chunks.par_iter_mut().for_each(|chunk| {
            chunk.filter_selected_cells(&cell_set);
        });

        let end_reading = start_reading.elapsed();

        if verbose {
            println!("Loaded in data: {:.2?}", end_reading);
        }

        let start_calculation = Instant::now();

        let res: Vec<(usize, f32, f32)> = gene_chunks
            .par_iter()
            .map(|chunk| self.compute_single_gene(chunk, &gex_model, centered))
            .collect();

        let mut gene_indices: Vec<usize> = Vec::with_capacity(res.len());
        let mut gaery_c: Vec<f64> = Vec::with_capacity(res.len());
        let mut z_scores: Vec<f64> = Vec::with_capacity(res.len());

        for (idx, c, z) in res {
            gene_indices.push(idx);
            gaery_c.push(c as f64);
            z_scores.push(z as f64);
        }

        let end_calculations = start_calculation.elapsed();

        if verbose {
            println!("Finsished the calculations: {:.2?}", end_calculations);
        }

        let p_vals = z_scores_to_pval(&z_scores, "twosided")?;
        let fdrs = calc_fdr(&p_vals);

        Ok(HotSpotGeneRes {
            gene_idx: gene_indices,
            c: gaery_c,
            z: z_scores,
            pval: p_vals,
            fdr: fdrs,
        })
    }

    /// Compute spatial autocorrelation with streaming (memory-efficient)
    ///
    /// Processes genes in batches to reduce memory usage for large datasets.
    ///
    /// ### Params
    ///
    /// * `gene_indices` - Indices of genes to analyse
    /// * `model` - Statistical model to use ("danb", "bernoulli", or "normal")
    /// * `centered` - Whether to centre the data before computing statistics
    /// * `verbose` - Whether to print progress information
    ///
    /// ### Returns
    ///
    /// `Result<HotSpotGeneRes>` with gene indices, Geary's C, Z-scores, derived
    /// p-values and FDR.
    pub fn compute_all_genes_streaming(
        &mut self,
        gene_indices: &[usize],
        model: &str,
        centered: bool,
        verbose: bool,
    ) -> Result<HotSpotGeneRes, String> {
        const GENE_BATCH_SIZE: usize = 1000;

        let no_genes = gene_indices.len();
        let no_batches = no_genes.div_ceil(GENE_BATCH_SIZE);
        let cell_set: IndexSet<u32> = self.cells_to_keep.iter().map(|&x| x as u32).collect();
        let reader = ParallelSparseReader::new(&self.f_path_gene).unwrap();

        let gex_model =
            parse_gex_model(model).ok_or_else(|| format!("Invalid model type: {}", model))?;

        self.populate_umi_counts();

        let mut results: Vec<(Vec<usize>, Vec<f64>, Vec<f64>)> = Vec::with_capacity(no_batches);

        for batch_idx in 0..no_batches {
            if verbose && batch_idx % 5 == 0 {
                let progress = (batch_idx + 1) as f32 / no_batches as f32 * 100.0;
                println!("  Progress: {:.1}%", progress);
            }

            let start_gene = batch_idx * GENE_BATCH_SIZE;
            let end_gene = ((batch_idx + 1) * GENE_BATCH_SIZE).min(no_genes);

            let gene_indices: Vec<usize> = (start_gene..end_gene).collect();

            let start_loading = Instant::now();

            let mut gene_chunks = reader.read_gene_parallel(&gene_indices);

            gene_chunks.par_iter_mut().for_each(|chunk| {
                chunk.filter_selected_cells(&cell_set);
            });

            let end_loading = start_loading.elapsed();

            if verbose {
                println!("   Loaded batch in: {:.2?}.", end_loading);
            }

            let start_calc = Instant::now();

            let batch_res: Vec<(usize, f32, f32)> = gene_chunks
                .par_iter()
                .map(|chunk| self.compute_single_gene(chunk, &gex_model, centered))
                .collect();

            let mut batch_gene_indices: Vec<usize> = Vec::with_capacity(batch_res.len());
            let mut batch_gaery_c: Vec<f64> = Vec::with_capacity(batch_res.len());
            let mut batch_z_scores: Vec<f64> = Vec::with_capacity(batch_res.len());

            for (idx, c, z) in batch_res {
                batch_gene_indices.push(idx);
                batch_gaery_c.push(c as f64);
                batch_z_scores.push(z as f64);
            }

            let end_calc = start_calc.elapsed();

            if verbose {
                println!("   Finished calculations in: {:.2?}.", end_calc);
            }

            results.push((batch_gene_indices, batch_gaery_c, batch_z_scores));
        }

        let mut gene_indices: Vec<usize> = Vec::new();
        let mut gaery_c: Vec<f64> = Vec::new();
        let mut z_scores: Vec<f64> = Vec::new();

        for (idx, c, z) in results {
            gene_indices.extend(idx);
            gaery_c.extend(c);
            z_scores.extend(z);
        }

        let p_vals = z_scores_to_pval(&z_scores, "twosided")?;
        let fdrs = calc_fdr(&p_vals);

        Ok(HotSpotGeneRes {
            gene_idx: gene_indices,
            c: gaery_c,
            z: z_scores,
            pval: p_vals,
            fdr: fdrs,
        })
    }

    /// Compute a single gene's spatial autocorrelation
    ///
    /// Internal method for calculating Geary's C and Z-score for one gene.
    ///
    /// ### Params
    ///
    /// * `gene_chunk` - Gene expression data
    /// * `gex_model` - Statistical model to apply
    /// * `centered` - Whether to centre the data
    ///
    /// ### Returns
    ///
    /// Tuple of (gene_index, Geary's C, Z-score)
    fn compute_single_gene(
        &self,
        gene_chunk: &CscGeneChunk,
        gex_model: &GexModel,
        centered: bool,
    ) -> (usize, f32, f32) {
        assert!(
            self.umi_counts.is_some(),
            "The internal UMI counts need to be populated"
        );

        let (mu, var, x2) = match gex_model {
            GexModel::DephAdjustNegBinom => {
                danb_model(gene_chunk, self.umi_counts.as_ref().unwrap(), self.n_cells)
            }
            GexModel::Bernoulli => {
                bernoulli_model(gene_chunk, self.umi_counts.as_ref().unwrap(), self.n_cells)
            }
            GexModel::Normal => {
                normal_model(gene_chunk, self.umi_counts.as_ref().unwrap(), self.n_cells)
            }
        };

        let mut vals = vec![0_f32; self.n_cells];
        for (&idx, &val) in gene_chunk.indices.iter().zip(&gene_chunk.data_raw) {
            vals[idx as usize] = val as f32;
        }

        if centered {
            center_values(&mut vals, &mu, &var);
        }

        let g = local_cov_weights(&vals, self.neighbours, &self.weights);

        let (eg, eg2) = if centered {
            (0.0, self.wtot2)
        } else {
            compute_moments_weights(&mu, &x2, self.neighbours, &self.weights)
        };

        let std_g = (eg2 - eg * eg).sqrt();
        let z = (g - eg) / std_g;

        let g_max = compute_local_cov_max(&self.node_degrees, &vals);
        let c = (g - eg) / g_max;

        (gene_chunk.original_index, c, z)
    }

    /// Compute pairwise gene correlations (in-memory version)
    ///
    /// Calculates local spatial correlations between all pairs of specified
    /// genes. Loads all gene data into memory for faster computation.
    /// WARNING: This will create a dense matrix of size n_cells x n_genes
    /// in memory! Should only be used for small data sets or very selected
    /// number of genes!
    ///
    /// ### Params
    ///
    /// * `gene_indices` - Indices of genes to analyse
    /// * `model` - Statistical model to use ("danb", "bernoulli", or "normal")
    /// * `verbose` - Whether to print progress information
    ///
    /// ### Returns
    ///
    /// Result containing HotSpotPairRes with correlation and Z-score matrices
    pub fn compute_gene_cor(
        &mut self,
        gene_indices: &[usize],
        model: &str,
        verbose: bool,
    ) -> Result<HotSpotPairRes, String> {
        let gex_model =
            parse_gex_model(model).ok_or_else(|| format!("Invalid model type: {}", model))?;

        self.populate_umi_counts();

        let cell_set: IndexSet<u32> = self.cells_to_keep.iter().map(|&x| x as u32).collect();

        if verbose {
            println!("Loading {} genes...", gene_indices.len());
        }

        let start_loading = Instant::now();
        let reader = ParallelSparseReader::new(&self.f_path_gene).unwrap();
        let mut gene_chunks: Vec<CscGeneChunk> = reader.read_gene_parallel(gene_indices);

        gene_chunks.par_iter_mut().for_each(|chunk| {
            chunk.filter_selected_cells(&cell_set);
        });

        if verbose {
            println!("Loaded data in {:.2?}", start_loading.elapsed());
            println!("Centering gene expression...");
        }

        let start_center = Instant::now();
        let centered_counts: Vec<Vec<f32>> = gene_chunks
            .par_iter()
            .map(|gene| {
                create_centered_counts_gene(
                    gene,
                    self.umi_counts.as_ref().unwrap(),
                    self.n_cells,
                    &gex_model,
                )
            })
            .collect();

        let n_genes = centered_counts.len();
        let mut counts_mat = Mat::zeros(n_genes, self.n_cells);
        for (i, gene_vec) in centered_counts.iter().enumerate() {
            for (j, &val) in gene_vec.iter().enumerate() {
                counts_mat[(i, j)] = val;
            }
        }

        if verbose {
            println!("Centered in {:.2?}", start_center.elapsed());
            println!("Computing conditional EG2 values...");
        }

        let start_eg2 = Instant::now();
        let eg2s: Vec<f32> = (0..n_genes)
            .into_par_iter()
            .map(|i| {
                let row = counts_mat.row(i);
                let row_vec = row.iter().copied().collect::<Vec<f32>>();
                conditional_eg2(&row_vec, self.neighbours, &self.weights)
            })
            .collect();

        if verbose {
            println!("Computed EG2 in {:.2?}", start_eg2.elapsed());
            println!("Computing pairwise correlations...");
        }

        let start_pairs = Instant::now();
        let n_pairs = (n_genes * (n_genes - 1)) / 2;

        let pairs: Vec<(usize, usize)> = (0..n_genes)
            .flat_map(|i| ((i + 1)..n_genes).map(move |j| (i, j)))
            .collect();

        let results: Vec<(usize, usize, f32, f32)> = pairs
            .par_iter()
            .map(|&(i, j)| {
                let x = counts_mat.row(i);
                let y = counts_mat.row(j);

                let lc = local_cov_pair(x, y, self.neighbours, &self.weights) * 2.0;

                // Use the minimum of the two Z-scores (more conservative)
                let eg = 0.0;

                let stdg_xy = eg2s[i].sqrt();
                let z_xy = (lc - eg) / stdg_xy;

                let stdg_yx = eg2s[j].sqrt();
                let z_yx = (lc - eg) / stdg_yx;

                let z = if z_xy.abs() < z_yx.abs() { z_xy } else { z_yx };

                (i, j, lc, z)
            })
            .collect();

        if verbose {
            println!(
                "Computed {} pairs in {:.2?}",
                n_pairs,
                start_pairs.elapsed()
            );
            println!("Building matrices...");
        }

        // generate symmetric matrices
        let mut lc_mat = Mat::zeros(n_genes, n_genes);
        let mut z_mat = Mat::zeros(n_genes, n_genes);

        for (i, j, lc, z) in results {
            lc_mat[(i, j)] = lc;
            lc_mat[(j, i)] = lc;
            z_mat[(i, j)] = z;
            z_mat[(j, i)] = z;
        }

        let lc_maxs = compute_local_cov_pairs_max(&self.node_degrees, &counts_mat);
        for i in 0..n_genes {
            for j in 0..n_genes {
                lc_mat[(i, j)] /= lc_maxs[(i, j)];
            }
        }

        if verbose {
            println!("Done!");
        }

        Ok(HotSpotPairRes {
            cor: lc_mat,
            z_scores: z_mat,
        })
    }

    /// Compute pairwise gene correlations (streaming version)
    ///
    /// Calculates local spatial correlations between all pairs of specified
    /// genes. Loads all gene data into memory for faster computation. Due to
    /// the nature of the problem, the function will calculate the correlation
    /// matrices in two passes with heavy nesting.
    ///
    /// ### Params
    ///
    /// * `gene_indices` - Indices of genes to analyse
    /// * `model` - Statistical model to use ("danb", "bernoulli", or "normal")
    /// * `verbose` - Whether to print progress information
    ///
    /// ### Returns
    ///
    /// Result containing HotSpotPairRes with correlation and Z-score matrices
    pub fn compute_gene_cor_streaming(
        &mut self,
        gene_indices: &[usize],
        model: &str,
        verbose: bool,
    ) -> Result<HotSpotPairRes, String> {
        const GENE_BATCH_SIZE: usize = 500;

        let gex_model =
            parse_gex_model(model).ok_or_else(|| format!("Invalid model type: {}", model))?;

        self.populate_umi_counts();

        let cell_set: IndexSet<u32> = self.cells_to_keep.iter().map(|&x| x as u32).collect();
        let reader = ParallelSparseReader::new(&self.f_path_gene).unwrap();

        let n_genes = gene_indices.len();
        let n_batches = n_genes.div_ceil(GENE_BATCH_SIZE);

        let mut lc_mat = Mat::zeros(n_genes, n_genes);
        let mut z_mat = Mat::zeros(n_genes, n_genes);

        let mut seen_pairs = HashSet::new();

        if verbose {
            println!("Processing {} genes in {} batches", n_genes, n_batches);
        }

        for batch_i in 0..n_batches {
            let start_batch_i = Instant::now();

            let start_i = batch_i * GENE_BATCH_SIZE;
            let end_i = ((batch_i + 1) * GENE_BATCH_SIZE).min(n_genes);
            let batch_i_indices = &gene_indices[start_i..end_i];

            if verbose {
                println!(
                    "\nProcessing batch {} / {} (genes {}-{})",
                    batch_i + 1,
                    n_batches,
                    start_i,
                    end_i - 1
                );
            }

            let mut batch_i_chunks = reader.read_gene_parallel(batch_i_indices);
            batch_i_chunks.par_iter_mut().for_each(|chunk| {
                chunk.filter_selected_cells(&cell_set);
            });

            let batch_i_centered: Vec<Vec<f32>> = batch_i_chunks
                .par_iter()
                .map(|gene| {
                    create_centered_counts_gene(
                        gene,
                        self.umi_counts.as_ref().unwrap(),
                        self.n_cells,
                        &gex_model,
                    )
                })
                .collect();

            let batch_i_eg2: Vec<f32> = batch_i_centered
                .par_iter()
                .map(|counts| conditional_eg2(counts, self.neighbours, &self.weights))
                .collect();

            let end_batch_i = start_batch_i.elapsed();

            if verbose {
                println!("Computed batch {} in: {:.2?}", batch_i + 1, end_batch_i);
            }

            let remaining_batches = n_batches - batch_i;
            if verbose {
                println!(
                    "  Computing pairs with {} remaining batches",
                    remaining_batches
                );
            }

            let start_remaining_batches = Instant::now();

            for batch_j in batch_i..n_batches {
                let start_j = batch_j * GENE_BATCH_SIZE;
                let end_j = ((batch_j + 1) * GENE_BATCH_SIZE).min(n_genes);

                if verbose {
                    println!("    Batch pair ({}, {})", batch_i + 1, batch_j + 1);
                }

                let (batch_j_centered, batch_j_eg2) = if batch_i == batch_j {
                    (batch_i_centered.clone(), batch_i_eg2.clone())
                } else {
                    let batch_j_indices = &gene_indices[start_j..end_j];
                    let mut batch_j_chunks = reader.read_gene_parallel(batch_j_indices);
                    batch_j_chunks.par_iter_mut().for_each(|chunk| {
                        chunk.filter_selected_cells(&cell_set);
                    });

                    let centered: Vec<Vec<f32>> = batch_j_chunks
                        .par_iter()
                        .map(|gene| {
                            create_centered_counts_gene(
                                gene,
                                self.umi_counts.as_ref().unwrap(),
                                self.n_cells,
                                &gex_model,
                            )
                        })
                        .collect();

                    let eg2: Vec<f32> = centered
                        .par_iter()
                        .map(|counts| conditional_eg2(counts, self.neighbours, &self.weights))
                        .collect();

                    (centered, eg2)
                };

                let pairs: Vec<(usize, usize)> = if batch_i == batch_j {
                    (0..batch_i_centered.len())
                        .flat_map(|i| ((i + 1)..batch_i_centered.len()).map(move |j| (i, j)))
                        .collect()
                } else {
                    (0..batch_i_centered.len())
                        .flat_map(|i| (0..batch_j_centered.len()).map(move |j| (i, j)))
                        .collect()
                };

                let results: Vec<(usize, usize, f32, f32, f32)> = pairs
                    .par_iter()
                    .map(|&(local_i, local_j)| {
                        let x = &batch_i_centered[local_i];
                        let y = &batch_j_centered[local_j];

                        let lc = local_cov_pair_vec(x, y, self.neighbours, &self.weights) * 2.0;

                        let max_i = compute_local_cov_max(&self.node_degrees, x);
                        let max_j = compute_local_cov_max(&self.node_degrees, y);
                        let lc_max = (max_i + max_j) / 2.0;

                        let eg = 0.0;
                        let stdg_xy = batch_i_eg2[local_i].sqrt();
                        let z_xy = (lc - eg) / stdg_xy;

                        let stdg_yx = batch_j_eg2[local_j].sqrt();
                        let z_yx = (lc - eg) / stdg_yx;

                        let z = if z_xy.abs() < z_yx.abs() { z_xy } else { z_yx };

                        let global_i = start_i + local_i;
                        let global_j = start_j + local_j;

                        (global_i, global_j, lc, z, lc_max)
                    })
                    .collect();

                for (i, j, lc, z, lc_max) in results {
                    if !seen_pairs.insert((i, j)) {
                        println!("DUPLICATE: ({}, {})", i, j);
                    }
                    let normalised_lc = if lc_max > 0.0 { lc / lc_max } else { 0.0 };

                    lc_mat[(i, j)] = normalised_lc;
                    lc_mat[(j, i)] = normalised_lc;
                    z_mat[(i, j)] = z;
                    z_mat[(j, i)] = z;
                }
            }

            let end_remaining_batches = start_remaining_batches.elapsed();

            if verbose {
                println!(
                    "Calculated all batches for batch {} in: {:.2?}",
                    batch_i + 1,
                    end_remaining_batches
                );
            }
        }

        if verbose {
            println!("Done!");
        }

        Ok(HotSpotPairRes {
            cor: lc_mat,
            z_scores: z_mat,
        })
    }

    /// Helper function to get the UMI counts per cell
    ///
    /// Reads and caches the total UMI count per cell for use in statistical
    /// models.
    fn populate_umi_counts(&mut self) {
        let reader = ParallelSparseReader::new(&self.f_path_cell).unwrap();
        let lib_sizes = reader.read_cell_library_sizes(self.cells_to_keep);
        let umi_counts = lib_sizes.iter().map(|x| *x as f32).collect::<Vec<f32>>();
        self.umi_counts = Some(umi_counts);
    }
}
