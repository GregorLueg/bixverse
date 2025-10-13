#![allow(unused_imports)]
#![allow(dead_code)]

use faer::{concat, Mat, MatRef};
use half::f16;
use indexmap::IndexSet;
use rand::prelude::*;
use rayon::prelude::*;
use rustc_hash::{FxHashMap, FxHashSet};
use std::path::Path;
use std::time::Instant;

use crate::core::base::pca_svd::{randomised_svd_f32, RandomSvdResults};
use crate::core::data::sparse_io::*;
use crate::core::data::sparse_structures::*;
use crate::single_cell::processing::*;
use crate::single_cell::sc_knn_snn::*;
use crate::utils::general::array_max_min;
use crate::utils::traits::F16;

///////////
// Types //
///////////

/// Type alias for Scrublet PCA results
///
/// ### Fields
///
/// * `0` - PCA scores
/// * `1` - PCA loadings
/// * `2` - Gene means
/// * `3` - Gene standard deviations
type ScrubletPcaRes = (Mat<f32>, Mat<f32>, Vec<f32>, Vec<f32>);

/// Type alias for Scrublet Doublet Scores
///
/// ### Fields
///
/// * `0` - Scores actual cells
/// * `1` - Errors actual cells
/// * `2` - Scores simulated cells
/// * `3` - Errors simulated cells
type ScrubletDoubletScores = (Vec<f32>, Vec<f32>, Vec<f32>, Vec<f32>);

////////////////////////
// Params and results //
////////////////////////

/// Structure that stores the Scrublet parameters
#[derive(Clone, Debug)]
pub struct ScrubletParams {
    // hvg detection
    pub min_gene_var_pctl: f32,
    pub hvg_method: String,
    pub loess_span: f64,
    pub clip_max: Option<f32>,

    // doublet generation
    pub sim_doublet_ratio: f32,
    pub expected_doublet_rate: f32,
    pub stdev_doublet_rate: f32,

    // pca
    pub no_pcs: usize,
    pub random_svd: bool,

    // knn
    pub k: usize,
    pub knn_method: String,
    pub dist_metric: String,
    pub search_budget: usize,
    pub n_trees: usize,
}

impl ScrubletParams {
    /// Generate Scrublet parameters with sensible defaults
    pub fn new() {
        // to be written
    }
}

/// Result structure for Scrublet
#[derive(Clone, Debug)]
pub struct ScrubletResult {
    pub predicted_doublets: Vec<bool>,
    pub doublet_scores_obs: Vec<f32>,
    pub doublet_scores_sim: Vec<f32>,
    pub doublet_errors_obs: Vec<f32>,
    pub z_scores: Vec<f32>,
    pub threshold: f32,
    pub detected_doublet_rate: f32,
    pub detectable_doublet_fraction: f32,
    pub overall_doublet_rate: f32,
}

/////////////
// Helpers //
/////////////

impl CsrCellChunk {
    /// Add two cells together for Scrublet doublet simulation
    ///
    /// This combines two cells by:
    /// 1. Filtering both to HVG genes only
    /// 2. Adding their raw counts
    /// 3. Using the combined library size (from ALL genes) for normalization
    /// 4. Normalizing to target size
    ///
    /// ### Params
    ///
    /// * `cell1` - First cell
    /// * `cell2` - Second cell  
    /// * `hvg_indices` - HashSet of HVG gene indices to keep
    /// * `target_size` - Target normalization size (e.g., 1e6 for CPM)
    /// * `doublet_index` - Index to assign to the new doublet
    ///
    /// ### Returns
    ///
    /// New `CsrCellChunk` representing the doublet
    pub fn add_cells_scrublet(
        cell1: &CsrCellChunk,
        cell2: &CsrCellChunk,
        hvg_indices: &FxHashSet<usize>,
        target_size: f32,
        doublet_index: usize,
    ) -> Self {
        // combined library size
        let combined_lib_size = cell1.library_size + cell2.library_size;
        let mut gene_counts: FxHashMap<u16, u32> = FxHashMap::default();

        // Add cell counts 1 and 2
        for i in 0..cell1.indices.len() {
            let gene_idx = cell1.indices[i] as usize;
            if hvg_indices.contains(&gene_idx) {
                *gene_counts.entry(cell1.indices[i]).or_insert(0) += cell1.data_raw[i] as u32;
            }
        }

        for i in 0..cell2.indices.len() {
            let gene_idx = cell2.indices[i] as usize;
            if hvg_indices.contains(&gene_idx) {
                *gene_counts.entry(cell2.indices[i]).or_insert(0) += cell2.data_raw[i] as u32;
            }
        }

        let mut gene_vec: Vec<(u16, u32)> = gene_counts.into_iter().collect();
        gene_vec.sort_unstable_by_key(|&(gene, _)| gene);

        let mut data_raw = Vec::with_capacity(gene_vec.len());
        let mut data_norm = Vec::with_capacity(gene_vec.len());
        let mut indices = Vec::with_capacity(gene_vec.len());

        let norm_factor = target_size / combined_lib_size as f32;

        for (gene, count) in gene_vec {
            // Clip to u16::MAX if needed
            let count_u16 = count.min(u16::MAX as u32) as u16;
            data_raw.push(count_u16);

            // Normalize: (count / combined_lib_size) * target_size, then ln(x+1)
            let normalized = (count as f32 * norm_factor).ln_1p();
            data_norm.push(F16::from(f16::from_f32(normalized)));

            indices.push(gene);
        }

        Self {
            data_raw,
            data_norm,
            library_size: combined_lib_size,
            indices,
            original_index: doublet_index,
            to_keep: true,
        }
    }
}

/// Calculate PCA and returns gene stats for downstream usage
///
/// ### Params
///
/// * `f_path` - Path to the gene-based binary file.
/// * `cell_indices` - Slice of indices for the cells.
/// * `gene_indices` - Slice of indices for the genes.
/// * `no_pcs` - Number of principal components to calculate
/// * `random_svd` - Shall randomised singular value decompostion be used. This
///   has the advantage of speed-ups, but loses precision.
/// * `seed` - Seed for randomised SVD.
/// * `verbose` - Boolean. Controls verbosity.
///
/// ### Return
///
/// A tuple of the samples projected on thePC space, gene loadings and singular
/// values.
pub fn pca_on_sc_with_stats(
    f_path: &str,
    cell_indices: &[usize],
    gene_indices: &[usize],
    no_pcs: usize,
    random_svd: bool,
    seed: usize,
    verbose: bool,
) -> ScrubletPcaRes {
    let start_total = Instant::now();
    let cell_set: IndexSet<u32> = cell_indices.iter().map(|&x| x as u32).collect();

    let start_reading = Instant::now();
    let reader = ParallelSparseReader::new(f_path).unwrap();
    let mut gene_chunks: Vec<CscGeneChunk> = reader.read_gene_parallel(gene_indices);
    let end_reading = start_reading.elapsed();

    if verbose {
        println!("Loaded in data : {:.2?}", end_reading);
    }

    let start_scaling = Instant::now();
    gene_chunks.par_iter_mut().for_each(|chunk| {
        chunk.filter_selected_cells(&cell_set);
    });

    let scaled_and_stats: Vec<(Vec<f32>, f32, f32)> = gene_chunks
        .par_iter()
        .map(|chunk| scale_csc_chunk(chunk, cell_indices.len()))
        .collect();

    let num_genes = scaled_and_stats.len();

    // Extract scaled data and statistics
    let scaled_data = Mat::from_fn(cell_indices.len(), num_genes, |row, col| {
        scaled_and_stats[col].0[row]
    });

    let means: Vec<f32> = scaled_and_stats.iter().map(|(_, m, _)| *m).collect();
    let stds: Vec<f32> = scaled_and_stats.iter().map(|(_, _, s)| *s).collect();

    let end_scaling = start_scaling.elapsed();
    if verbose {
        println!("Finished scaling : {:.2?}", end_scaling);
    }

    let start_svd = Instant::now();
    let (scores, loadings) = if random_svd {
        let res: RandomSvdResults<f32> =
            randomised_svd_f32(scaled_data.as_ref(), no_pcs, seed, Some(100_usize), None);
        let loadings = res.v.submatrix(0, 0, num_genes, no_pcs).to_owned();
        let scores = &scaled_data * &loadings;
        (scores, loadings)
    } else {
        let res = scaled_data.thin_svd().unwrap();
        let loadings = res.V().submatrix(0, 0, num_genes, no_pcs).to_owned();
        let scores = &scaled_data * &loadings;
        (scores, loadings)
    };

    let end_svd = start_svd.elapsed();
    if verbose {
        println!("Finished PCA calculations : {:.2?}", end_svd);
    }

    let end_total = start_total.elapsed();
    if verbose {
        println!("Total run time PCA detection: {:.2?}", end_total);
    }

    (scores, loadings, means, stds)
}

/// Scale Vec<CsrCellChunk> using pre-calculated gene means and stds
///
/// This ensures simulated doublets are scaled using the SAME statistics
/// as the observed cells (critical for proper PCA projection)
///
/// ### Params
///
/// * `chunks` - Vector of cell chunks (simulated doublets)
/// * `gene_means` - Mean for each gene (from observed data)
/// * `gene_stds` - Std dev for each gene (from observed data)
/// * `n_genes` - Total number of genes
///
/// ### Returns
///
/// Dense matrix (cells x genes) with z-scored values
pub fn scale_cell_chunks_with_stats(
    chunks: &[CsrCellChunk],
    gene_means: &[f32],
    gene_stds: &[f32],
    n_genes: usize,
) -> Mat<f32> {
    let n_cells = chunks.len();
    let mut scaled = Mat::<f32>::zeros(n_cells, n_genes);

    for (cell_idx, chunk) in chunks.iter().enumerate() {
        // First, fill zeros with -mean/std (z-score of 0)
        for gene in 0..n_genes {
            let z_score = -gene_means[gene] / gene_stds[gene];
            *scaled.get_mut(cell_idx, gene) = z_score;
        }

        // Then overwrite with actual values
        for i in 0..chunk.indices.len() {
            let gene = chunk.indices[i] as usize;
            let val = chunk.data_norm[i].to_f32();
            let z_score = (val - gene_means[gene]) / gene_stds[gene];
            *scaled.get_mut(cell_idx, gene) = z_score;
        }
    }

    scaled
}

pub fn find_threshold_min(scores: &[f32], n_bins: usize) -> f32 {
    let (min_score, max_score) = array_max_min(scores);

    if (max_score - min_score).abs() < 1e-6 {
        return (min_score + max_score) / 2.0;
    }

    let bin_width = (max_score - min_score) / n_bins as f32;

    let mut hist = vec![0usize; n_bins];

    for &score in scores {
        let bin = ((score - min_score) / bin_width).floor() as usize;
        let bin = bin.min(n_bins - 1);
        hist[bin] += 1;
    }

    let smoothed = moving_average(&hist, 3);

    let max_count = *smoothed.iter().max().unwrap_or(&1);
    let threshold_count = (max_count as f32 * 0.1) as usize;

    let mut found_first_peak = false;
    let mut min_idx = 0;
    let mut min_val = usize::MAX;

    for (i, &smoothed_i) in smoothed.iter().enumerate() {
        if !found_first_peak {
            // Look for first peak
            if smoothed_i > threshold_count {
                found_first_peak = true;
            }
        } else {
            // After first peak, find minimum
            if smoothed_i < min_val {
                min_idx = i;
                min_val = smoothed_i;
            }

            // Stop if we hit another peak (prevents going too far right)
            if smoothed_i > (min_val as f32 * 1.5) as usize {
                break;
            }
        }
    }

    // if no valley is found use median of simulated scores
    if !found_first_peak || min_val == usize::MAX {
        let mut sorted = scores.to_vec();
        sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
        return sorted[sorted.len() / 2];
    }

    min_score + (min_idx as f32 + 0.5) * bin_width
}

// Moving average function for the binning
fn moving_average(data: &[usize], window: usize) -> Vec<usize> {
    let half_window = window / 2;
    data.iter()
        .enumerate()
        .map(|(i, _)| {
            let start = i.saturating_sub(half_window);
            let end = (i + half_window + 1).min(data.len());
            let sum: usize = data[start..end].iter().sum();
            sum / (end - start)
        })
        .collect()
}

////////////////////
// Main structure //
////////////////////

#[derive(Clone, Debug)]
pub struct Scrublet {
    // paths
    f_path_gene: String,
    f_path_cell: String,
    // parameters
    params: ScrubletParams,
    // other
    n_cells: usize,
    n_cells_sim: usize,
    cells_to_keep: Vec<usize>,
}

impl Scrublet {
    /// Generate a new instance
    pub fn new(
        f_path_gene: &str,
        f_path_cells: &str,
        params: ScrubletParams,
        n_cells: usize,
        cell_indices: &[usize],
    ) -> Self {
        Scrublet {
            f_path_gene: f_path_gene.to_string(),
            f_path_cell: f_path_cells.to_string(),
            params,
            n_cells,
            n_cells_sim: 0,
            cells_to_keep: cell_indices.to_vec(),
        }
    }

    pub fn run_scrublet(
        &mut self,
        streaming: bool,
        manual_threshold: Option<f32>,
        n_bins: usize,
        seed: usize,
        verbose: bool,
    ) -> ScrubletResult {
        // HVG
        if verbose {
            println!("Identifying highly variable genes...");
        }
        let start_all = Instant::now();

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

        // simulate doublets
        if verbose {
            println!("Simulating doublets...");
        }
        let start_doublet_gen = Instant::now();

        let sim_chunks = self.simulate_doublets(&hvg_genes, seed);

        let end_doublet_gen = start_doublet_gen.elapsed();

        if verbose {
            println!(
                "Simulated {} doublets. Done in {:.2?}",
                sim_chunks.len(),
                end_doublet_gen
            );
        }

        // PCA
        if verbose {
            println!("Running PCA...");
        }

        let start_pca = Instant::now();

        let combined_pca = self.run_pca(&sim_chunks, &hvg_genes, verbose, seed);

        let end_pca = start_pca.elapsed();

        if verbose {
            println!("Done with PCA in {:.2?}", end_pca);
        }

        // KNN
        if verbose {
            println!("Building kNN graph...");
        }

        let start_knn = Instant::now();

        let knn_indices = self
            .build_combined_knn(combined_pca, seed, verbose)
            .expect("Failed to build kNN graph");

        let end_knn = start_knn.elapsed();

        if verbose {
            println!("Done with KNN generation in {:.2?}", end_knn);
        }

        if verbose {
            println!("Calculating doublet scores...");
        }
        let start_doublets = Instant::now();

        let doublet_scores: ScrubletDoubletScores = self.calculate_doublet_scores(&knn_indices);

        let res = self.call_doublets(doublet_scores, manual_threshold, n_bins, verbose);

        let end_doublets = start_doublets.elapsed();

        if verbose {
            println!("Done with doublet scoring and calling {:.2?}", end_doublets);
        }

        let end_all = start_all.elapsed();

        if verbose {
            println!("Finished Scrublet {:.2?}", end_all);
        }

        res
    }

    /// Get the indices of the highly variable genes
    fn get_hvg(&self, streaming: bool, verbose: bool) -> Vec<usize> {
        // identify highly variable genes
        let hvg_type = get_hvg_method(&self.params.hvg_method)
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
        let n_genes_to_take = (n_genes * self.params.min_gene_var_pctl).ceil() as usize;

        let mut indices: Vec<(usize, f64)> = hvg_res
            .var_std
            .iter()
            .enumerate()
            .map(|(i, &v)| (i, v))
            .collect();
        indices.sort_unstable_by(|a, b| b.1.partial_cmp(&a.1).unwrap());
        indices.truncate(n_genes_to_take);
        indices.into_iter().map(|(i, _)| i).collect()
    }

    /// Generate a vector of doublets
    fn simulate_doublets(&mut self, hvg_genes: &[usize], seed: usize) -> Vec<CsrCellChunk> {
        let n_sim_doublets = (self.n_cells as f32 * self.params.sim_doublet_ratio) as usize;
        self.n_cells_sim = n_sim_doublets;
        let mut rng = StdRng::seed_from_u64(seed as u64);

        let pairs: Vec<(usize, usize)> = (0..n_sim_doublets)
            .map(|_| {
                let i = self.cells_to_keep[rng.random_range(0..self.cells_to_keep.len())];
                let j = self.cells_to_keep[rng.random_range(0..self.cells_to_keep.len())];
                (i, j)
            })
            .collect();

        let hvg_set: FxHashSet<usize> = hvg_genes.iter().copied().collect();

        let reader = ParallelSparseReader::new(&self.f_path_cell).unwrap();

        let doublets: Vec<CsrCellChunk> = pairs
            .par_iter()
            .enumerate()
            .map(|(doublet_idx, &(i, j))| {
                let cell1 = reader.read_cell(i);
                let cell2 = reader.read_cell(j);

                CsrCellChunk::add_cells_scrublet(
                    &cell1,
                    &cell2,
                    &hvg_set,
                    1e6, // CPM normalization
                    doublet_idx,
                )
            })
            .collect();

        doublets
    }

    /// Run the PCA
    fn run_pca(
        &self,
        sim_chunks: &[CsrCellChunk],
        hvg_genes: &[usize],
        verbose: bool,
        seed: usize,
    ) -> Mat<f32> {
        let pca_res: ScrubletPcaRes = pca_on_sc_with_stats(
            &self.f_path_gene,
            &self.cells_to_keep,
            hvg_genes,
            self.params.no_pcs,
            self.params.random_svd,
            seed,
            verbose,
        );

        let scaled_sim =
            scale_cell_chunks_with_stats(sim_chunks, &pca_res.2, &pca_res.3, hvg_genes.len());

        let pca_sim = &scaled_sim * pca_res.1;

        concat![[pca_res.0], [pca_sim]]
    }

    /// Generate the kNN graph
    fn build_combined_knn(
        &self,
        embd: Mat<f32>,
        seed: usize,
        verbose: bool,
    ) -> Result<Vec<Vec<usize>>, String> {
        let knn_method = get_knn_method(&self.params.knn_method)
            .ok_or_else(|| format!("Invalid KNN search method: {}", &self.params.knn_method))?;

        let r = self.n_cells_sim as f32 / self.n_cells as f32;
        let k_adj = (self.params.k as f32 * (1.0 + r)).round() as usize;

        let knn = match knn_method {
            KnnSearch::Hnsw => generate_knn_hnsw(
                embd.as_ref(),
                &self.params.dist_metric,
                k_adj,
                seed,
                verbose,
            ),
            KnnSearch::Annoy => generate_knn_annoy(
                embd.as_ref(),
                &self.params.dist_metric,
                k_adj,
                self.params.n_trees,
                self.params.search_budget,
                seed,
                verbose,
            ),
        };

        Ok(knn)
    }

    /// Calculate the doublet scores
    fn calculate_doublet_scores(&self, knn_indices: &[Vec<usize>]) -> ScrubletDoubletScores {
        let n_obs = self.n_cells;
        let n_sim = self.n_cells_sim;

        let r = n_sim as f32 / n_obs as f32;
        let rho = self.params.expected_doublet_rate;
        let se_rho = self.params.stdev_doublet_rate;

        let k_adj = (self.params.k as f32 * (1.0 + r)).round() as usize;
        let n_adj = k_adj as f32;

        let scores_errors: Vec<(f32, f32)> = knn_indices
            .par_iter()
            .map(|neighbours| {
                let n_sim_neigh = neighbours.iter().filter(|&&idx| idx >= n_obs).count() as f32;
                let q = (n_sim_neigh + 1.0) / (n_adj + 2.0);
                let denominator = 1.0 - rho - q * (1.0 - rho - rho / r);
                let score = if denominator.abs() > 1e-10 {
                    (q * rho / r) / denominator
                } else {
                    0.0
                };

                let se_q = (q * (1.0 - q) / (n_adj + 3.0)).sqrt();
                let factor = q * rho / r / (denominator * denominator);
                let se_score = factor
                    * ((se_q / q * (1.0 - rho)).powi(2) + (se_rho / rho * (1.0 - q)).powi(2))
                        .sqrt();

                (score.max(0.0), se_score.max(1e-10))
            })
            .collect();

        let (scores_obs, errors_obs): (Vec<f32>, Vec<f32>) =
            scores_errors[..n_obs].iter().copied().unzip();

        let (scores_sim, errors_sim): (Vec<f32>, Vec<f32>) =
            scores_errors[n_obs..].iter().copied().unzip();

        (scores_obs, errors_obs, scores_sim, errors_sim)
    }

    fn call_doublets(
        &self,
        doublet_scores: ScrubletDoubletScores,
        manual_threshold: Option<f32>,
        n_bins: usize,
        verbose: bool,
    ) -> ScrubletResult {
        let threshold = manual_threshold.unwrap_or_else(|| {
            let t = find_threshold_min(&doublet_scores.2, n_bins);
            if verbose {
                println!("Automatically set threshold at doublet score = {:.4}", t);
            }
            t
        });

        // Call doublets
        let predicted_doublets: Vec<bool> = doublet_scores
            .0
            .iter()
            .map(|&score| score > threshold)
            .collect();

        // Calculate z-scores
        let z_scores: Vec<f32> = doublet_scores
            .0
            .iter()
            .zip(doublet_scores.1.iter())
            .map(|(&score, &error)| (score - threshold) / error)
            .collect();

        // Calculate statistics
        let n_detected = predicted_doublets.iter().filter(|&&x| x).count();
        let detected_doublet_rate = n_detected as f32 / doublet_scores.0.len() as f32;

        let n_detectable = doublet_scores.2.iter().filter(|&&s| s > threshold).count();
        let detectable_doublet_fraction = n_detectable as f32 / doublet_scores.2.len() as f32;

        let overall_doublet_rate = if detectable_doublet_fraction > 0.01 {
            detected_doublet_rate / detectable_doublet_fraction
        } else {
            0.0
        };

        if verbose {
            println!(
                "Detected doublet rate = {:.1}%",
                100.0 * detected_doublet_rate
            );
            println!(
                "Estimated detectable doublet fraction = {:.1}%",
                100.0 * detectable_doublet_fraction
            );
            println!("Overall doublet rate:");
            println!(
                "  Expected  = {:.1}%",
                100.0 * self.params.expected_doublet_rate
            );
            println!("  Estimated = {:.1}%", 100.0 * overall_doublet_rate);
        }

        ScrubletResult {
            predicted_doublets,
            doublet_scores_obs: doublet_scores.0,
            doublet_scores_sim: doublet_scores.2,
            doublet_errors_obs: doublet_scores.1,
            z_scores,
            threshold,
            detected_doublet_rate,
            detectable_doublet_fraction,
            overall_doublet_rate,
        }
    }
}
