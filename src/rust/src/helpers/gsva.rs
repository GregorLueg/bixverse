use extendr_api::prelude::*;

use faer::{Mat, MatRef};
use once_cell::sync::Lazy;
use rayon::prelude::*;
use rustc_hash::FxHashMap;
use statrs::distribution::{DiscreteCDF, Poisson};
use std::sync::Mutex;
use std::time::Instant;

use crate::utils::general::standard_deviation;
use crate::utils::r_rust_interface::NamedMatrix;

/////////////
// Globals //
/////////////

const SIGMA_FACTOR: f64 = 4.0;
const POISSON_BANDWIDTH: f64 = 0.5;

// Cache for Poisson distributions to avoid repeated creation
static POISSON_CACHE: Lazy<std::sync::Mutex<FxHashMap<u64, Poisson>>> =
    Lazy::new(|| std::sync::Mutex::new(FxHashMap::default()));

/////////////
// Helpers //
/////////////

// Fast approximate normal CDF using Abramowitz & Stegun approximation
#[inline]
fn fast_normal_cdf(x: f64, mean: f64, std_dev: f64) -> f64 {
    let z = (x - mean) / std_dev;

    if z < -8.0 {
        return 0.0;
    }
    if z > 8.0 {
        return 1.0;
    }

    let t = 1.0 / (1.0 + 0.2316419 * z.abs());
    let poly = t
        * (0.319381530
            + t * (-0.356563782 + t * (1.781477937 + t * (-1.821255978 + t * 1.330274429))));
    let cdf = 1.0 - 0.39894228 * (-0.5 * z * z).exp() * poly;

    if z >= 0.0 {
        cdf
    } else {
        1.0 - cdf
    }
}

// Get cached Poisson or create new one
#[inline]
fn get_poisson_cdf(lambda: f64, k: u64) -> f64 {
    let lambda_key = (lambda * 1000.0) as u64; // Cache with 0.001 precision

    let mut cache = POISSON_CACHE.lock().unwrap();
    let poisson = cache
        .entry(lambda_key)
        .or_insert_with(|| Poisson::new(lambda).unwrap());
    poisson.cdf(k)
}

/// Get gene set indices for GSVA
///
/// ### Params
///
/// * `exp` - `NamedMatrix` with rows representing the genes and columns
///           representing the samples
/// * `gs_list` - R list that contains the different gene sets
///
/// ### Returns
///
/// A vector of vectors with the index positions as usizes
pub fn get_gsva_gs_indices(exp: &NamedMatrix, gs_list: List) -> Result<Vec<Vec<usize>>> {
    let mut gs_indices: Vec<Vec<usize>> = Vec::with_capacity(gs_list.len());

    for i in 0..gs_list.len() {
        let list_elem = gs_list.elt(i).unwrap();
        let elem = list_elem.as_str_vector().unwrap();
        let gs_indices_i = exp.get_row_indices(Some(&elem[..]));
        gs_indices.push(gs_indices_i);
    }

    Ok(gs_indices)
}

/// Calculate the kernel density on a per column (gene) basis
fn row_kernel_density(
    density_row: &[f64],
    test_row: &[f64],
    use_gaussian: bool,
    bandwidth: f64,
) -> Vec<f64> {
    let mut results = Vec::with_capacity(test_row.len());
    let density_len_inv = 1.0 / density_row.len() as f64;

    if use_gaussian {
        // Pre-compute for vectorization
        for &test_val in test_row {
            let left_tail: f64 = density_row
                .iter()
                .map(|&density_val| fast_normal_cdf(test_val, density_val, bandwidth))
                .sum::<f64>()
                * density_len_inv;

            let log_odds = -((1.0 - left_tail) / left_tail.max(1e-15)).ln();
            results.push(log_odds);
        }
    } else {
        for &test_val in test_row {
            let test_val_u64 = test_val.max(0.0) as u64;
            let left_tail: f64 = density_row
                .iter()
                .map(|&density_val| {
                    let lambda = density_val + POISSON_BANDWIDTH;
                    if lambda > 0.0 {
                        get_poisson_cdf(lambda, test_val_u64)
                    } else {
                        0.0
                    }
                })
                .sum::<f64>()
                * density_len_inv;

            let log_odds = -((1.0 - left_tail) / left_tail.max(1e-15)).ln();
            results.push(log_odds);
        }
    }
    results
}

/// Matrix kernel density estimation
/// Matches the C function `matrix_d()` in kernel_estimation.c
pub fn matrix_kernel_density(
    density_matrix: &MatRef<f64>,
    test_matrix: &MatRef<f64>,
    use_gaussian: bool,
) -> Mat<f64> {
    let (n_genes, _) = density_matrix.shape();
    let (_, n_test_samples) = test_matrix.shape();

    let bandwidths: Vec<f64> = if use_gaussian {
        (0..n_genes)
            .map(|gene_idx| {
                let row_data: Vec<f64> = density_matrix.row(gene_idx).iter().copied().collect();
                (standard_deviation(&row_data) / SIGMA_FACTOR).max(0.001)
            })
            .collect()
    } else {
        vec![0.0; n_genes]
    };

    let results: Vec<Vec<f64>> = (0..n_genes)
        .into_par_iter()
        .map(|gene_idx| {
            let density_data: Vec<f64> = density_matrix.row(gene_idx).iter().copied().collect();
            let test_data: Vec<f64> = test_matrix.row(gene_idx).iter().copied().collect();
            row_kernel_density(
                &density_data,
                &test_data,
                use_gaussian,
                bandwidths[gene_idx],
            )
        })
        .collect();

    let mut result = Mat::zeros(n_genes, n_test_samples);
    for (gene_idx, row_data) in results.into_iter().enumerate() {
        for (sample_idx, val) in row_data.into_iter().enumerate() {
            result[(gene_idx, sample_idx)] = val;
        }
    }

    result
}

/// Order and rank statistic calculation
/// Matches the C function `order_rankstat()` in utils.c
pub fn order_rankstat(values: &[f64]) -> (Vec<usize>, Vec<f64>) {
    let n = values.len();
    let mut indices: Vec<usize> = (0..n).collect();

    // Sort indices by values in descending order
    // Matches C code: qsort(ord, n, sizeof(int), indirect_dbl_cmp_dec)
    indices.sort_by(|&a, &b| values[b].partial_cmp(&values[a]).unwrap());

    let mut rank_stats = vec![0.0; n];

    // Calculate symmetric rank statistic: |n - rank - n/2|
    // Matches C code: rst[ord[i]-1] = fabs(((double) n) - ((double) i) - (((double) n) / 2.0))
    for (rank, &original_idx) in indices.iter().enumerate() {
        rank_stats[original_idx] = ((n as f64) - (rank as f64) - (n as f64) / 2.0).abs();
    }

    (indices, rank_stats)
}

/// GSVA random walk calculation
/// Matches the C function `gsva_rnd_walk()` in ks_test.c
pub fn gsva_random_walk(
    gene_set_indices: &[usize],
    decreasing_order_indices: &[usize],
    symmetric_rank_stats: &[f64],
    tau: f64,
) -> (f64, f64) {
    let n = decreasing_order_indices.len();

    // Create rank lookup map for O(1) access instead of O(n) position() calls
    let rank_lookup: FxHashMap<usize, usize> = decreasing_order_indices
        .iter()
        .enumerate()
        .map(|(rank, &gene_idx)| (gene_idx, rank))
        .collect();

    // Get ranks of genes in gene set using HashMap lookup
    let mut gene_set_ranks = Vec::with_capacity(gene_set_indices.len());
    for &gene_idx in gene_set_indices {
        if let Some(&rank) = rank_lookup.get(&gene_idx) {
            gene_set_ranks.push(rank);
        }
    }

    if gene_set_ranks.is_empty() {
        return (f64::NAN, f64::NAN);
    }

    // Initialize step CDF arrays
    let mut step_cdf_in_geneset = vec![0.0; n];
    let mut step_cdf_out_geneset = vec![1; n];

    // Fill in gene set positions
    for &rank in &gene_set_ranks {
        let gene_idx = decreasing_order_indices[rank];
        step_cdf_in_geneset[rank] = if tau == 1.0 {
            symmetric_rank_stats[gene_idx]
        } else {
            symmetric_rank_stats[gene_idx].powf(tau)
        };
        step_cdf_out_geneset[rank] = 0;
    }

    // Convert to cumulative sums
    for i in 1..n {
        step_cdf_in_geneset[i] += step_cdf_in_geneset[i - 1];
        step_cdf_out_geneset[i] += step_cdf_out_geneset[i - 1];
    }

    let total_in = step_cdf_in_geneset[n - 1];
    let total_out = step_cdf_out_geneset[n - 1] as f64;

    if total_in <= 0.0 || total_out <= 0.0 {
        return (f64::NAN, f64::NAN);
    }

    let mut max_pos = 0.0;
    let mut max_neg = 0.0;
    let total_in_inv = 1.0 / total_in;
    let total_out_inv = 1.0 / total_out;

    // Calculate random walk statistic
    for i in 0..n {
        let walk_stat = step_cdf_in_geneset[i] * total_in_inv
            - (step_cdf_out_geneset[i] as f64) * total_out_inv;

        if walk_stat > max_pos {
            max_pos = walk_stat;
        }
        if walk_stat < max_neg {
            max_neg = walk_stat;
        }
    }

    (max_pos, max_neg)
}

/// Calculate GSVA enrichment scores for gene sets
/// Matches the C function `gsva_score_genesets_R()` in ks_test.c
pub fn gsva_score_genesets(
    gene_sets: &[Vec<usize>],
    decreasing_order_indices: &[usize],
    symmetric_rank_stats: &[f64],
    tau: f64,
    max_diff: bool,
    abs_rank: bool,
) -> Vec<f64> {
    let mut scores = Vec::with_capacity(gene_sets.len());

    for gene_set in gene_sets {
        let (walk_stat_pos, walk_stat_neg) = gsva_random_walk(
            gene_set,
            decreasing_order_indices,
            symmetric_rank_stats,
            tau,
        );

        let score = if walk_stat_pos.is_nan() || walk_stat_neg.is_nan() {
            f64::NAN
        } else if max_diff {
            // Normalized ES: walkstatpos + walkstatneg (or walkstatpos - walkstatneg if abs_rank)
            if abs_rank {
                walk_stat_pos - walk_stat_neg
            } else {
                walk_stat_pos + walk_stat_neg
            }
        } else {
            // Classical ES: maximum deviation from zero
            if walk_stat_pos > walk_stat_neg.abs() {
                walk_stat_pos
            } else {
                walk_stat_neg
            }
        };

        scores.push(score);
    }

    scores
}

/// Complete GSVA pipeline
/// Follows the algorithm described in the paper and implemented in the C code
pub fn gsva(
    expression_matrix: &MatRef<f64>,
    gene_sets: &[Vec<usize>],
    use_gaussian: bool,
    tau: f64,
    max_diff: bool,
    abs_rank: bool,
) -> Mat<f64> {
    let start_total = Instant::now();
    let (_, n_samples) = expression_matrix.shape();
    let mut result = Mat::zeros(gene_sets.len(), n_samples);

    println!(
        "Starting GSVA with {} samples and {} gene sets",
        n_samples,
        gene_sets.len()
    );

    // Step 1: Kernel density estimation for each gene across all samples
    // This creates log-odds transformed KCDF values
    let start_kcdf = Instant::now();
    let kcdf_matrix = matrix_kernel_density(
        expression_matrix, // density reference
        expression_matrix, // test data (same as density)
        use_gaussian,
    );
    let kcdf_time = start_kcdf.elapsed();
    println!("Step 1 - Kernel density estimation: {:.2?}", kcdf_time);

    // Parallel processing over samples
    let result_mutex = Mutex::new(&mut result);
    let start_parallel = Instant::now();

    (0..n_samples).into_par_iter().for_each(|sample_idx| {
        let sample_start = Instant::now();

        // Extract KCDF values for this sample - use column view
        let extract_start = Instant::now();
        let sample_kcdf: Vec<f64> = kcdf_matrix.col(sample_idx).iter().copied().collect();
        let extract_time = extract_start.elapsed();

        // Step 2: Ranking and symmetric rank statistic
        let rank_start = Instant::now();
        let (decreasing_order, rank_stats) = order_rankstat(&sample_kcdf);
        let rank_time = rank_start.elapsed();

        // Step 3 & 4: Calculate enrichment scores for each gene set
        let score_start = Instant::now();
        let scores = gsva_score_genesets(
            gene_sets,
            &decreasing_order,
            &rank_stats,
            tau,
            max_diff,
            abs_rank,
        );
        let score_time = score_start.elapsed();

        // Store results with lock
        let store_start = Instant::now();
        {
            let mut result_guard = result_mutex.lock().unwrap();
            for (gene_set_idx, &score) in scores.iter().enumerate() {
                result_guard[(gene_set_idx, sample_idx)] = score;
            }
        }
        let store_time = store_start.elapsed();

        let sample_total = sample_start.elapsed();

        // Print timing for every 100th sample or first/last few
        if sample_idx < 5 || sample_idx % 100 == 0 || sample_idx >= n_samples - 5 {
            println!(
                "Sample {}: total={:.2?}, extract={:.2?}, rank={:.2?}, score={:.2?}, store={:.2?}",
                sample_idx, sample_total, extract_time, rank_time, score_time, store_time
            );
        }
    });

    let parallel_time = start_parallel.elapsed();
    let total_time = start_total.elapsed();

    println!("Step 2-4 - Parallel processing: {:.2?}", parallel_time);
    println!("Total GSVA time: {:.2?}", total_time);
    println!(
        "Time breakdown: KCDF={:.1?} ({:.1}%), Parallel={:.1?} ({:.1}%)",
        kcdf_time,
        kcdf_time.as_secs_f64() / total_time.as_secs_f64() * 100.0,
        parallel_time,
        parallel_time.as_secs_f64() / total_time.as_secs_f64() * 100.0
    );

    result
}
