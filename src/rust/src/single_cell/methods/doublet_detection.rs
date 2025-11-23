use extendr_api::List;
use faer::{concat, Mat, MatRef};
use rand::prelude::*;
use rayon::prelude::*;
use rustc_hash::{FxHashMap, FxHashSet};
use std::time::Instant;

use crate::core::base::stats::hypergeom_pval;
use crate::core::data::sparse_io::*;
use crate::core::graph::community_detection::*;
use crate::core::graph::knn::*;
use crate::single_cell::methods::scrublet::*;
use crate::single_cell::processing::*;
use crate::single_cell::sc_knn_snn::*;

///////////
// Types //
///////////

/// Type alias for Boost PCA results
///
/// (Same as ScrubletPcaRes, but f--k DRY)
///
/// ### Fields
///
/// * `0` - PCA scores
/// * `1` - PCA loadings
/// * `2` - Gene means
/// * `3` - Gene standard deviations
type BoostPcaRes = (Mat<f32>, Mat<f32>, Vec<f32>, Vec<f32>);

////////////
// Params //
////////////

////////////////////////
// Params and results //
////////////////////////

/// Structure to store the Boost parameters
///
/// ### Fields
///
/// **General parameters:**
///
/// * `log_transform` - Shall the counts be log-transformed
/// * `mean_center` - Shall the data be mean-centred
/// * `normalise_variance` - Shall the data be variance normalised
/// * `target_size` - Optional target size. If not provided, will default to
///   the mean library size of the cells.
///
/// **HVG Detection:**
///
/// * `min_gene_var_pctl` - Percentile threshold for highly variable genes.
/// * `hvg_method` - Method for HVG selection. One of `"vst"`, `"mvb"`, or
///   `"dispersion"`.
/// * `loess_span` - Span parameter for loess fitting in VST method.
/// * `clip_max` - Optional maximum value for clipping in variance
///   stabilisation.
///
/// **Doublet Generation:**
///
/// * `boost_rate` - Number of doublets to simulate relative to the number of
///   observed cells (e.g., 1.5 simulates 1.5x as many doublets).
/// * `replace` - Whether to use replacement when sampling cell pairs.
///
/// **PCA:**
///
/// * `no_pcs` - Number of principal components to use for embedding.
/// * `random_svd` - Whether to use randomised SVD (faster) vs exact SVD.
///
/// **Clustering and Iteration:**
///
/// * `resolution` - Resolution parameter for Louvain clustering.
/// * `n_iters` - Number of boosting iterations to perform.
///
/// **Doublet Calling:**
///
/// * `p_thresh` - P-value threshold for doublet calling via voting.
/// * `voter_thresh` - Threshold for majority voting (0-1).
///
/// **kNN Graph:**
///
/// * `k` - Number of nearest neighbours for the kNN graph. If 0 (default),
///   automatically calculated as round(0.5 * sqrt(n_cells)).
/// * `knn_method` - Method for approximate nearest neighbour search.
///   One of `"annoy"` or `"hnsw"`.
/// * `dist_metric` - Distance metric to use. One of `"euclidean"` or
///   `"cosine"`.
/// * `search_budget` - Search budget for Annoy (higher = more accurate but
///   slower).
/// * `n_trees` - Number of trees for Annoy index generation.
/// * `nn_max_iter` - Maximum iterations for the NNDescent kNN search
/// * `rho` - Sampling rate for NNDescent.
/// * `delta` - Early termination criterium for NNDescent
#[derive(Clone, Debug)]
pub struct BoostParams {
    // general params
    pub log_transform: bool,
    pub mean_center: bool,
    pub normalise_variance: bool,
    pub target_size: Option<f32>,
    // hvg detection
    pub min_gene_var_pctl: f32,
    pub hvg_method: String,
    pub loess_span: f64,
    pub clip_max: Option<f32>,
    // doublet generation
    pub boost_rate: f32,
    pub replace: bool,
    // pca
    pub no_pcs: usize,
    pub random_svd: bool,
    // clustering
    pub resolution: f32,
    pub louvain_iters: usize,
    // iterations
    pub n_iters: usize,
    // doublet calling
    pub p_thresh: f32,
    pub voter_thresh: f32,
    // knn
    pub knn_params: KnnParams,
}

impl BoostParams {
    /// Generate BoostParams from an R list
    ///
    /// Should values not be found within the List, the parameters will default
    /// to sensible defaults based on the Boost algorithm.
    ///
    /// ### Params
    ///
    /// * `r_list` - The list with the Boost parameters.
    ///
    /// ### Returns
    ///
    /// The `BoostParams` with all parameters set.
    pub fn from_r_list(r_list: List) -> Self {
        let knn_params = KnnParams::from_r_list(r_list.clone());

        let params_list = r_list.into_hashmap();

        // norm parameters
        let log_transform = params_list
            .get("log_transform")
            .and_then(|v| v.as_bool())
            .unwrap_or(false);

        let mean_center = params_list
            .get("mean_center")
            .and_then(|v| v.as_bool())
            .unwrap_or(true);

        let normalise_variance = params_list
            .get("normalise_variance")
            .and_then(|v| v.as_bool())
            .unwrap_or(true);

        let target_size = params_list
            .get("target_size")
            .and_then(|v| v.as_real())
            .map(|x| x as f32);

        // hvg
        let min_gene_var_pctl = params_list
            .get("min_gene_var_pctl")
            .and_then(|v| v.as_real())
            .unwrap_or(0.85) as f32;

        let hvg_method = std::string::String::from(
            params_list
                .get("hvg_method")
                .and_then(|v| v.as_str())
                .unwrap_or("vst"),
        );

        let loess_span = params_list
            .get("loess_span")
            .and_then(|v| v.as_real())
            .unwrap_or(0.3);

        let clip_max = params_list
            .get("clip_max")
            .and_then(|v| v.as_real())
            .map(|x| x as f32);

        // doublet detection params
        let boost_rate = params_list
            .get("boost_rate")
            .and_then(|v| v.as_real())
            .unwrap_or(0.25) as f32;

        let replace = params_list
            .get("replace")
            .and_then(|v| v.as_bool())
            .unwrap_or(false);

        let resolution = params_list
            .get("resolution")
            .and_then(|v| v.as_real())
            .unwrap_or(1.0) as f32;

        let louvain_iters = params_list
            .get("louvain_iters")
            .and_then(|v| v.as_integer())
            .unwrap_or(10) as usize;

        let n_iters = params_list
            .get("n_iters")
            .and_then(|v| v.as_integer())
            .unwrap_or(10) as usize;

        // pca
        let no_pcs = params_list
            .get("no_pcs")
            .and_then(|v| v.as_integer())
            .unwrap_or(30) as usize;

        let random_svd = params_list
            .get("random_svd")
            .and_then(|v| v.as_bool())
            .unwrap_or(true);

        let p_thresh = params_list
            .get("p_thresh")
            .and_then(|v| v.as_real())
            .unwrap_or(1e-7) as f32;

        let voter_thresh = params_list
            .get("voter_thresh")
            .and_then(|v| v.as_real())
            .unwrap_or(0.9) as f32;

        Self {
            log_transform,
            mean_center,
            normalise_variance,
            target_size,
            min_gene_var_pctl,
            hvg_method,
            loess_span,
            clip_max,
            boost_rate,
            replace,
            no_pcs,
            random_svd,
            resolution,
            louvain_iters,
            n_iters,
            p_thresh,
            voter_thresh,
            knn_params,
        }
    }
}

/////////////
// Results //
/////////////

/// Result structure for Boost doublet detection
///
/// Contains predictions, scores, and voting statistics from the Boost
/// algorithm.
///
/// ### Fields
///
/// * `predicted_doublets` - Boolean vector indicating which observed cells are
///   predicted as doublets (true = doublet, false = singlet).
/// * `doublet_scores` - Doublet scores for each observed cell, typically
///   averaged across iterations. Higher scores indicate higher likelihood of
///   being a doublet.
/// * `voting_average` - Average voting fraction across iterations. Indicates
///   the consensus across boosting iterations for each cell. Only meaningful
///   when n_iters > 1.
#[derive(Clone, Debug)]
pub struct BoostResult {
    pub predicted_doublets: Vec<bool>,
    pub doublet_scores: Vec<f32>,
    pub voting_average: Vec<f32>,
}

/////////////
// Helpers //
/////////////

/// PCA for Boost (wrapper around pca_scrublet)
///
/// Computes PCA using the same methodology as Scrublet for consistency.
/// Reuses the existing `pca_scrublet()` function by constructing compatible
/// ScrubletParams.
///
/// ### Params
///
/// * `f_path_gene` - Path to the gene-based binary file.
/// * `cell_indices` - Slice of cell indices to include.
/// * `gene_indices` - Slice of gene indices (HVG).
/// * `library_sizes` - Library sizes for each cell (HVG only).
/// * `target_size` - Target normalisation size.
/// * `params` - Boost parameters.
/// * `seed` - Seed for randomised SVD.
/// * `verbose` - Controls verbosity.
///
/// ### Returns
///
/// Tuple of (PCA scores, loadings, gene means, gene standard deviations).
#[allow(clippy::too_many_arguments)]
fn pca_boost(
    f_path_gene: &str,
    cell_indices: &[usize],
    gene_indices: &[usize],
    library_sizes: &[usize],
    target_size: f32,
    params: &BoostParams,
    seed: usize,
    verbose: bool,
) -> BoostPcaRes {
    pca_scrublet(
        f_path_gene,
        cell_indices,
        gene_indices,
        library_sizes,
        target_size,
        // Create ScrubletParams from BoostParams for compatibility
        &ScrubletParams {
            log_transform: params.log_transform,
            mean_center: params.mean_center,
            normalise_variance: params.normalise_variance,
            target_size: params.target_size,
            min_gene_var_pctl: params.min_gene_var_pctl,
            hvg_method: params.hvg_method.clone(),
            loess_span: params.loess_span,
            clip_max: params.clip_max,
            sim_doublet_ratio: 0.0, // Not used in PCA
            expected_doublet_rate: 0.0,
            stdev_doublet_rate: 0.0,
            no_pcs: params.no_pcs,
            random_svd: params.random_svd,
            n_bins: 0,
            manual_threshold: None,
            knn_params: KnnParams::new(),
        },
        seed,
        verbose,
    )
}

/// Score communities based on synthetic doublet enrichment
///
/// Calculates enrichment of simulated doublets within each community using
/// hypergeometric test. Communities with high synthetic doublet enrichment
/// receive higher scores.
///
/// ### Params
///
/// * `orig_communities` - Community assignments for observed cells.
/// * `synth_communities` - Community assignments for simulated doublets.
///
/// ### Returns
///
/// Tuple of (enrichment scores, log p-values) for each observed cell.
fn score_communities(
    orig_communities: &[usize],
    synth_communities: &[usize],
) -> (Vec<f32>, Vec<f32>) {
    let mut synth_counts: FxHashMap<usize, usize> = FxHashMap::default();
    let mut orig_counts: FxHashMap<usize, usize> = FxHashMap::default();

    for &c in synth_communities {
        *synth_counts.entry(c).or_insert(0) += 1;
    }
    for &c in orig_communities {
        *orig_counts.entry(c).or_insert(0) += 1;
    }

    // Pre-compute per-community
    let unique_comms: FxHashSet<usize> = orig_communities.iter().copied().collect();
    let mut comm_scores: FxHashMap<usize, f32> = FxHashMap::default();
    let mut comm_log_p: FxHashMap<usize, f32> = FxHashMap::default();

    for &c in &unique_comms {
        let n_synth = *synth_counts.get(&c).unwrap_or(&0);
        let n_orig = *orig_counts.get(&c).unwrap_or(&0);

        comm_scores.insert(c, n_synth as f32 / (n_synth + n_orig) as f32);

        let k = n_synth + n_orig;
        let log_p =
            hypergeom_pval(n_synth, synth_communities.len(), orig_communities.len(), k).ln() as f32;
        comm_log_p.insert(c, log_p);
    }

    let scores: Vec<f32> = orig_communities.iter().map(|&c| comm_scores[&c]).collect();
    let log_p_values: Vec<f32> = orig_communities.iter().map(|&c| comm_log_p[&c]).collect();

    (scores, log_p_values)
}

/// Predict doublets via voting across iterations
///
/// Uses majority voting across iterations to call doublets. A cell is called
/// a doublet if it exceeds the voter threshold in the fraction of iterations
/// where it was significant.
///
/// ### Params
///
/// * `all_log_p_values` - Log p-values for each cell across all iterations.
/// * `p_thresh` - P-value threshold for significance in each iteration.
/// * `voter_thresh` - Fraction threshold for majority voting (0-1).
///
/// ### Returns
///
/// Tuple of (doublet predictions, average voting fraction) for each cell.
fn predict_voting(
    all_log_p_values: &[Vec<f32>],
    p_thresh: f32,
    voter_thresh: f32,
) -> (Vec<bool>, Vec<f32>) {
    let n_cells = all_log_p_values[0].len();
    let n_iters = all_log_p_values.len();
    let log_p_thresh = p_thresh.ln();

    let voting_avg: Vec<f32> = (0..n_cells)
        .map(|cell_idx| {
            let votes: usize = all_log_p_values
                .iter()
                .filter(|iter_vals| iter_vals[cell_idx] <= log_p_thresh)
                .count();
            votes as f32 / n_iters as f32
        })
        .collect();

    let labels: Vec<bool> = voting_avg.iter().map(|&avg| avg >= voter_thresh).collect();

    (labels, voting_avg)
}

/// Find score cutoff using largest gap heuristic
///
/// For single iteration mode, identifies the threshold between singlet and
/// doublet scores by finding the largest gap in the sorted score distribution.
///
/// ### Params
///
/// * `scores` - Community enrichment scores.
///
/// ### Returns
///
/// Threshold score value.
fn find_score_cutoff(scores: &[f32]) -> f32 {
    let mut sorted = scores.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());

    if sorted.len() < 2 {
        return sorted[0];
    }

    let mut max_gap = 0.0f32;
    let mut gap_idx = 0;

    for i in 1..sorted.len() {
        let gap = sorted[i] - sorted[i - 1];
        if gap > max_gap {
            max_gap = gap;
            gap_idx = i;
        }
    }

    sorted[gap_idx]
}

//////////
// Main //
//////////

/// Structure for Boost doublet detection algorithm
///
/// Implements the iterative Boost algorithm for doublet detection using
/// community detection and voting across multiple iterations.
///
/// ### Fields
///
/// * `f_path_gene` - Path to the binarised file in CSC format.
/// * `f_path_cell` - Path to the binarised file in CSR format.
/// * `params` - The Boost parameters
/// * `n_cells` - Number of observed cells.
/// * `n_cells_sim` - Number of simulated cells (varies per iteration).
/// * `cells_to_keep` - Indices of cells to keep/include in this analysis.
/// * `hvg_library_sizes` - Library sizes for HVGs only, computed once at start.
#[derive(Clone, Debug)]
pub struct BoostClassifier {
    f_path_gene: String,
    f_path_cell: String,
    params: BoostParams,
    n_cells: usize,
    n_cells_sim: usize,
    cells_to_keep: Vec<usize>,
    hvg_library_sizes: Vec<usize>,
}

impl BoostClassifier {
    /// Generate a new instance
    ///
    /// ### Params
    ///
    /// * `f_path_gene` - Path to the binarised file in CSC format.
    /// * `f_path_cell` - Path to the binarised file in CSR format.
    /// * `params` - The Boost parameters to use.
    /// * `cell_indices` - Slice of indices indicating which cells to keep/use.
    pub fn new(
        f_path_gene: &str,
        f_path_cell: &str,
        params: BoostParams,
        cell_indices: &[usize],
    ) -> Self {
        let n_cells = cell_indices.len();

        BoostClassifier {
            f_path_gene: f_path_gene.to_string(),
            f_path_cell: f_path_cell.to_string(),
            params,
            n_cells,
            n_cells_sim: 0,
            cells_to_keep: cell_indices.to_vec(),
            hvg_library_sizes: Vec::new(),
        }
    }

    /// Main function to run Boost
    ///
    /// Executes the full Boost algorithm across all iterations, identifying
    /// highly variable genes, simulating doublets, running PCA, building kNN
    /// graphs, clustering, and calling doublets via voting or score cutoff.
    ///
    /// ### Params
    ///
    /// * `streaming` - Shall the data be streamed. Reduces memory pressure
    ///   during HVG detection.
    /// * `seed` - Seed for reproducibility.
    /// * `verbose` - Controls verbosity of the function.
    ///
    /// ### Returns
    ///
    /// A `BoostResult` containing predictions and scores.
    pub fn run_boost(&mut self, streaming: bool, seed: usize, verbose: bool) -> BoostResult {
        let start_all = Instant::now();

        if verbose {
            println!("Identifying highly variable genes...");
        }
        let start_hvg = Instant::now();
        let hvg_genes = self.get_hvg(streaming, verbose);
        let end_hvg = start_hvg.elapsed();
        if verbose {
            println!(
                "Using {} highly variable genes. Done in {:.2?}",
                hvg_genes.len(),
                end_hvg
            );
        }

        let cell_reader = ParallelSparseReader::new(&self.f_path_cell).unwrap();
        let hvg_set: FxHashSet<usize> = hvg_genes.iter().copied().collect();

        let hvg_library_sizes: Vec<usize> = self
            .cells_to_keep
            .par_iter()
            .map(|&cell_idx| {
                let chunk = cell_reader.read_cell(cell_idx);
                chunk
                    .indices
                    .iter()
                    .enumerate()
                    .filter(|(_, &gene_idx)| hvg_set.contains(&(gene_idx as usize)))
                    .map(|(i, _)| chunk.data_raw[i] as usize)
                    .sum()
            })
            .collect();

        let target_size = self.params.target_size.unwrap_or_else(|| {
            let sum = hvg_library_sizes.iter().sum::<usize>() as f32;
            sum / hvg_library_sizes.len() as f32
        });

        self.hvg_library_sizes = hvg_library_sizes;
        self.n_cells_sim = (self.n_cells as f32 * self.params.boost_rate) as usize;

        let start_iters = Instant::now();

        let iter_results: Vec<(Vec<f32>, Vec<f32>)> = (0..self.params.n_iters)
            .map(|iter| {
                if verbose {
                    println!(
                        " == Running iteration {} of {} == ",
                        iter + 1,
                        self.params.n_iters
                    );
                }
                self.one_iteration(target_size, &hvg_genes, seed + iter, verbose)
            })
            .collect();

        let end_iters = start_iters.elapsed();
        if verbose {
            println!(
                "Completed {} iterations in {:.2?}",
                self.params.n_iters, end_iters
            );
        }

        // Aggregate results
        let all_scores: Vec<Vec<f32>> = iter_results.iter().map(|(s, _)| s.clone()).collect();
        let all_log_p_values: Vec<Vec<f32>> = iter_results.iter().map(|(_, p)| p.clone()).collect();

        let result = if self.params.n_iters > 1 {
            let (labels, voting_avg) = predict_voting(
                &all_log_p_values,
                self.params.p_thresh,
                self.params.voter_thresh,
            );

            // Average scores across iterations
            let avg_scores: Vec<f32> = (0..self.n_cells)
                .map(|i| {
                    all_scores.iter().map(|iter| iter[i]).sum::<f32>() / self.params.n_iters as f32
                })
                .collect();

            BoostResult {
                predicted_doublets: labels,
                doublet_scores: avg_scores,
                voting_average: voting_avg,
            }
        } else {
            // Single iteration - use score cutoff
            let scores = &all_scores[0];
            let cutoff = find_score_cutoff(scores);

            if verbose {
                println!("Score cutoff: {:.4}", cutoff);
            }

            let labels: Vec<bool> = scores.iter().map(|&s| s >= cutoff).collect();
            let voting_avg = vec![0.0; self.n_cells]; // Not meaningful for single iteration

            BoostResult {
                predicted_doublets: labels,
                doublet_scores: scores.clone(),
                voting_average: voting_avg,
            }
        };

        let end_all = start_all.elapsed();
        if verbose {
            let n_doublets = result.predicted_doublets.iter().filter(|&&d| d).count();
            println!(
                "Detected {} doublets ({:.1}%)",
                n_doublets,
                100.0 * n_doublets as f32 / self.n_cells as f32
            );
            println!("Total runtime: {:.2?}", end_all);
        }

        result
    }

    /// Execute a single Boost iteration
    ///
    /// Performs one complete cycle: simulating doublets, running PCA, building
    /// kNN graph, clustering, and scoring communities based on simulated
    /// doublet enrichment.
    ///
    /// ### Params
    ///
    /// * `target_size` - Target library size for normalisation.
    /// * `hvg_genes` - Indices of highly variable genes.
    /// * `seed` - Seed for reproducibility in this iteration.
    /// * `verbose` - Controls verbosity.
    ///
    /// ### Returns
    ///
    /// Tuple of (community scores, log p-values) for each observed cell.
    pub fn one_iteration(
        &self,
        target_size: f32,
        hvg_genes: &[usize],
        seed: usize,
        verbose: bool,
    ) -> (Vec<f32>, Vec<f32>) {
        // Simulate doublets
        let sim_chunks = self.simulate_doublets(target_size, hvg_genes, seed);

        // Run PCA
        let combined_pca = self.run_pca(&sim_chunks, hvg_genes, target_size, verbose, seed);

        // Build kNN graph
        let knn = self
            .build_combined_knn(combined_pca.as_ref(), seed, verbose)
            .unwrap();

        // Cluster
        let start_graph = Instant::now();

        let graph = knn_to_sparse_graph(&knn);

        let end_graph = start_graph.elapsed();

        if verbose {
            println!("Transformed kNN graph. Done in {:.2?}", end_graph);
        }

        let start_cluster = Instant::now();

        let communities = louvain_sparse_graph(
            &graph,
            self.params.resolution,
            self.params.louvain_iters,
            seed,
        );

        let end_cluster = start_cluster.elapsed();

        if verbose {
            println!(
                "Generated communities via Louvain clustering. Done in {:.2?}",
                end_cluster
            );
        }

        // Split communities
        let orig_communities = communities[..self.n_cells].to_vec();
        let synth_communities = communities[self.n_cells..].to_vec();

        // Score
        let (scores, log_p_values) = score_communities(&orig_communities, &synth_communities);

        (scores, log_p_values)
    }

    /// Get the indices of the highly variable genes
    ///
    /// Identify the HVGs for subsequent PCA.
    ///
    /// ### Params
    ///
    /// * `streaming` - Shall the data be loaded in a streaming fashion. Reduces
    ///   memory pressure
    /// * `verbose` - Controls verbosity
    ///
    /// ### Returns
    ///
    /// Vector of indices of the HVG.
    fn get_hvg(&self, streaming: bool, verbose: bool) -> Vec<usize> {
        // Same as Scrublet - reuse your existing code
        let hvg_type = parse_hvg_method(&self.params.hvg_method)
            .ok_or_else(|| format!("Invalid HVG method: {}", &self.params.hvg_method))
            .unwrap();

        let hvg_res: HvgRes = if streaming {
            match hvg_type {
                HvgMethod::Vst => get_hvg_vst_streaming(
                    &self.f_path_gene,
                    &self.cells_to_keep,
                    self.params.loess_span,
                    self.params.clip_max,
                    verbose,
                ),
                HvgMethod::MeanVarBin => get_hvg_mvb_streaming(),
                HvgMethod::Dispersion => get_hvg_dispersion_streaming(),
            }
        } else {
            match hvg_type {
                HvgMethod::Vst => get_hvg_vst(
                    &self.f_path_gene,
                    &self.cells_to_keep,
                    self.params.loess_span,
                    self.params.clip_max,
                    verbose,
                ),
                HvgMethod::MeanVarBin => get_hvg_mvb(),
                HvgMethod::Dispersion => get_hvg_dispersion(),
            }
        };

        let n_genes = hvg_res.mean.len() as f32;
        let n_genes_to_take = (n_genes * (1.0 - self.params.min_gene_var_pctl)).ceil() as usize;

        let mut indices: Vec<(usize, f64)> = hvg_res
            .var_std
            .iter()
            .enumerate()
            .map(|(i, &v)| (i, v))
            .collect();
        indices.sort_unstable_by(|a, b| b.1.partial_cmp(&a.1).unwrap());
        indices.truncate(n_genes_to_take);

        let mut result: Vec<usize> = indices.into_iter().map(|(i, _)| i).collect();
        result.sort_unstable();
        result
    }

    /// Simulate doublets
    ///
    /// Generate artifical doublets based on the real data. The algorithm
    /// will only construct doublets from the
    ///
    /// ### Params
    ///
    /// * `hvg_genes` - Indices of the highly variable genes.
    /// * `target_size` - The library target size. Needs to be the same as the
    ///   original one (in single cell typicall `1e4`).
    /// * `seed` - Seed for reproducibility.
    ///
    /// ### Returns
    ///
    /// A `Vec<CsrCellChunk>` that contains the artificially created doublets.
    fn simulate_doublets(
        &self,
        target_size: f32,
        hvg_genes: &[usize],
        seed: usize,
    ) -> Vec<CsrCellChunk> {
        let n_sim_doublets = (self.n_cells as f32 * self.params.boost_rate) as usize;
        let mut rng = StdRng::seed_from_u64(seed as u64);

        let pairs: Vec<(usize, usize)> = if self.params.replace {
            (0..n_sim_doublets)
                .map(|_| {
                    let i = rng.random_range(0..self.cells_to_keep.len());
                    let j = rng.random_range(0..self.cells_to_keep.len());
                    (i, j)
                })
                .collect()
        } else {
            let mut available: Vec<usize> = (0..self.cells_to_keep.len()).collect();
            available.shuffle(&mut rng);
            available
                .chunks(2)
                .take(n_sim_doublets)
                .filter_map(|chunk| {
                    if chunk.len() == 2 {
                        Some((chunk[0], chunk[1]))
                    } else {
                        None
                    }
                })
                .collect()
        };

        let hvg_set: FxHashSet<usize> = hvg_genes.iter().copied().collect();
        let gene_to_hvg_idx: FxHashMap<usize, u16> = hvg_genes
            .iter()
            .enumerate()
            .map(|(hvg_idx, &orig_idx)| (orig_idx, hvg_idx as u16))
            .collect();

        let reader = ParallelSparseReader::new(&self.f_path_cell).unwrap();

        pairs
            .par_iter()
            .enumerate()
            .map(|(doublet_idx, &(pos_i, pos_j))| {
                let cell_idx_i = self.cells_to_keep[pos_i];
                let cell_idx_j = self.cells_to_keep[pos_j];

                let cell1 = reader.read_cell(cell_idx_i);
                let cell2 = reader.read_cell(cell_idx_j);

                let hvg_combined_lib_size =
                    self.hvg_library_sizes[pos_i] + self.hvg_library_sizes[pos_j];

                let mut doublet = CsrCellChunk::add_cells_scrublet(
                    &cell1,
                    &cell2,
                    &hvg_set,
                    hvg_combined_lib_size,
                    target_size,
                    self.params.log_transform,
                    doublet_idx,
                );

                for idx in doublet.indices.iter_mut() {
                    *idx = gene_to_hvg_idx[&(*idx as usize)];
                }

                doublet
            })
            .collect()
    }

    /// Run PCA
    ///
    /// Runs PCA prior to kNN construction.
    ///
    /// ### Params
    ///
    /// * `sim_chunks` - Slice of `CsrCellChunk` containing the simulated
    ///   doublets.
    /// * `hvg_genes` - The indices of the highly variable genes.
    /// * `verbose` - Controls verbosity of the function.
    /// * `seed` - Seed for reproducibility. Relevant when using randomised
    ///   SVD.
    ///
    /// ### Returns
    ///
    /// The PCA scores with the top rows representing the actual data and
    /// the bottom rows the simulated data.
    fn run_pca(
        &self,
        sim_chunks: &[CsrCellChunk],
        hvg_genes: &[usize],
        target_size: f32,
        verbose: bool,
        seed: usize,
    ) -> Mat<f32> {
        let pca_res: BoostPcaRes = pca_boost(
            &self.f_path_gene,
            &self.cells_to_keep,
            hvg_genes,
            &self.hvg_library_sizes,
            target_size,
            &self.params,
            seed,
            verbose,
        );

        let scaled_sim = scale_cell_chunks_with_stats(
            sim_chunks,
            &pca_res.2,
            &pca_res.3,
            self.params.mean_center,
            self.params.normalise_variance,
            hvg_genes.len(),
        );

        let pca_sim = &scaled_sim * pca_res.1;
        concat![[pca_res.0], [pca_sim]]
    }

    /// Generate the kNN graph
    ///
    /// ### Params
    ///
    /// * `embd` - The embedding matrix to use for the generation of the kNN
    ///   graph. Usually the PCA of observed and simulated doublet cells.
    /// * `seed` - Seed for reproducibility. Relevant when using randomised
    /// * `verbose` - Controls verbosity of the function.
    ///
    /// ### Returns
    ///
    /// The kNN graph as a `Vec<Vec<usize>>`.
    fn build_combined_knn(
        &self,
        embd: MatRef<f32>,
        seed: usize,
        verbose: bool,
    ) -> Result<Vec<Vec<usize>>, String> {
        let knn_method = parse_knn_method(&self.params.knn_params.knn_method).ok_or_else(|| {
            format!(
                "Invalid KNN search method: {}",
                &self.params.knn_params.knn_method
            )
        })?;

        let k_adj = self.calculate_k_adj();

        if verbose {
            println!("Using {} neighbours in the kNN generation.", k_adj);
        }

        let knn = match knn_method {
            KnnSearch::Hnsw => generate_knn_hnsw(
                embd.as_ref(),
                &self.params.knn_params.ann_dist,
                k_adj,
                self.params.knn_params.m,
                self.params.knn_params.ef_construction,
                self.params.knn_params.ef_search,
                seed,
                verbose,
            ),
            KnnSearch::Annoy => generate_knn_annoy(
                embd.as_ref(),
                &self.params.knn_params.ann_dist,
                k_adj,
                self.params.knn_params.n_tree,
                self.params.knn_params.search_budget,
                seed,
                verbose,
            ),
            KnnSearch::NNDescent => generate_knn_nndescent(
                embd.as_ref(),
                &self.params.knn_params.ann_dist,
                k_adj,
                self.params.knn_params.max_iter,
                self.params.knn_params.delta,
                self.params.knn_params.rho,
                seed,
                verbose,
            ),
        };

        Ok(knn)
    }

    /// Calculates the adjusted k based on number actual cells and simulated
    /// cells
    fn calculate_k_adj(&self) -> usize {
        let k = if self.params.knn_params.k == 0 {
            ((self.n_cells as f32).sqrt() * 0.5).round() as usize
        } else {
            self.params.knn_params.k
        };

        let r = self.n_cells_sim as f32 / self.n_cells as f32;
        (k as f32 * (1.0 + r)).round() as usize
    }
}
