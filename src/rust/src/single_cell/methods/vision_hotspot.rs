#![allow(dead_code)]

use extendr_api::{Conversions, List};
use rayon::prelude::*;
use std::time::Instant;

use crate::core::base::stats::z_scores_to_pval;
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
pub fn parse_model(s: &str) -> Option<GexModel> {
    match s.to_lowercase().as_str() {
        "danp" => Some(GexModel::DephAdjustNegBinom),
        "bernoulli" => Some(GexModel::Bernoulli),
        "normal" => Some(GexModel::Normal),
        _ => None,
    }
}

/// Structure for the gene results
///
/// ### Fields
///
/// To be written
#[derive(Debug, Clone)]
pub struct HotSpotGeneRes {
    pub gene_idx: usize,
    pub c: f32,
    pub z: f32,
}

/// Compute momentum weights
///
/// ### Params
///
/// ### Return
fn compute_moments_weights(
    mu: &[f32],
    x2: &[f32],
    neighbors: &[Vec<usize>],
    weights: &[Vec<f32>],
) -> (f32, f32) {
    let n = neighbors.len();

    let mu_sq: Vec<f32> = mu.iter().map(|&m| m * m).collect();

    let mut eg = 0_f32;
    let mut t1 = vec![0_f32; n];
    let mut t2 = vec![0_f32; n];

    for i in 0..n {
        let mu_i = mu[i];
        let mu_i_sq = mu_sq[i];

        for (k, &j) in neighbors[i].iter().enumerate() {
            let wij = weights[i][k];
            let mu_j = mu[j];

            eg += wij * mu_i * mu_j;

            t1[i] += wij * mu_j;
            let wij_sq = wij * wij;
            t2[i] += wij_sq * mu[j] * mu_j;
            t1[j] += wij * mu_i;
            t2[j] += wij_sq * mu_i_sq;
        }
    }

    let mut eg2 = 0_f32;

    for i in 0..n {
        eg2 += (x2[i] - mu_sq[i]) * (t1[i] * t1[i] - t2[i]);
    }

    for i in 0..n {
        let x2_i = x2[i];
        let mu_sq_i = mu_sq[i];

        for (k, &j) in neighbors[i].iter().enumerate() {
            let wij = weights[i][k];
            eg2 += wij * wij * (x2_i * x2[j] - mu_sq_i * mu_sq[j]);
        }
    }

    eg2 += eg * eg;

    (eg, eg2)
}

fn local_cov_weights(vals: &[f32], neighbors: &[Vec<usize>], weights: &[Vec<f32>]) -> f32 {
    let mut out = 0.0;

    for i in 0..vals.len() {
        let xi = vals[i];
        if xi == 0.0 {
            continue;
        }

        for (k, &j) in neighbors[i].iter().enumerate() {
            let xj = vals[j];
            let wij = weights[i][k];

            if xj != 0.0 && wij != 0.0 {
                out += xi * xj * wij;
            }
        }
    }

    out
}

fn compute_local_cov_max(node_degrees: &[f32], vals: &[f32]) -> f32 {
    let mut tot = 0.0;
    for i in 0..node_degrees.len() {
        tot += node_degrees[i] * vals[i] * vals[i];
    }
    tot / 2.0
}

fn center_values(vals: &mut [f32], mu: &[f32], var: &[f32]) {
    for i in 0..vals.len() {
        vals[i] = (vals[i] - mu[i]) / var[i].sqrt();
    }
}

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
    let vv = sum_sq * (n / (n - 1.0));

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

//////////
// Main //
//////////

#[derive(Clone, Debug)]
pub struct Hotspot<'a> {
    f_path_gene: String,
    f_path_cell: String,
    neighbors: &'a [Vec<usize>],
    weights: &'a [Vec<f32>],
    cells_to_keep: &'a [usize],
    node_degrees: Vec<f32>,
    umi_counts: Option<Vec<f32>>,
    wtot2: f32,
    n_cells: usize,
}

impl<'a> Hotspot<'a> {
    /// To be written
    pub fn new(
        f_path_gene: String,
        f_path_cell: String,
        cells_to_keep: &'a [usize],
        neighbors: &'a [Vec<usize>],
        weights: &'a [Vec<f32>],
    ) -> Self {
        let n_cells = neighbors.len();
        let node_degrees = weights
            .par_iter()
            .map(|w| w.iter().sum::<f32>())
            .collect::<Vec<f32>>();

        let wtot2: f32 = weights.iter().flatten().map(|&w| w * w).sum();

        Self {
            f_path_gene,
            f_path_cell,
            neighbors,
            weights,
            cells_to_keep,
            node_degrees,
            umi_counts: None,
            wtot2,
            n_cells,
        }
    }

    /// Compute a single gene
    ///
    /// ### Params
    ///
    /// ### Returns
    pub fn compute_single_gene(
        &self,
        gene: &CscGeneChunk,
        gex_model: &GexModel,
        centered: bool,
    ) -> HotSpotGeneRes {
        assert!(
            self.umi_counts.is_some(),
            "The internal UMI counts need to be populated"
        );

        let (mu, var, x2) = match gex_model {
            GexModel::DephAdjustNegBinom => {
                danb_model(gene, self.umi_counts.as_ref().unwrap(), self.n_cells)
            }
            GexModel::Bernoulli => panic!("Bernoulli model not implemented yet for HotSpot"),
            GexModel::Normal => panic!("Normal model not implemented yet for HotSpot"),
        };

        let mut vals = vec![0_f32; self.n_cells];
        for (&idx, &val) in gene.indices.iter().zip(&gene.data_raw) {
            vals[idx as usize] = val as f32;
        }

        if centered {
            center_values(&mut vals, &mu, &var);
        }

        let g = local_cov_weights(&vals, self.neighbors, self.weights);

        let (eg, eg2) = if centered {
            (0.0, self.wtot2)
        } else {
            compute_moments_weights(&mu, &x2, self.neighbors, self.weights)
        };

        let std_g = (eg2 - eg * eg).sqrt();
        let z = (g - eg) / std_g;

        let g_max = compute_local_cov_max(&self.node_degrees, &vals);
        let c = (g - eg) / g_max;

        HotSpotGeneRes {
            gene_idx: gene.original_index,
            c,
            z,
        }
    }

    /// Helper function to get the UMI counts per cell
    ///
    /// Add the umi counts for the cells to keep
    fn populate_umi_counts(&mut self) {
        let reader = ParallelSparseReader::new(&self.f_path_cell).unwrap();
        let lib_sizes = reader.read_cell_library_sizes(self.cells_to_keep);
        let umi_counts = lib_sizes.iter().map(|x| *x as f32).collect::<Vec<f32>>();
        self.umi_counts = Some(umi_counts);
    }
}
