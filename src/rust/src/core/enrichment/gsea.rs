use extendr_api::List;

use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use rayon::prelude::*;

use crate::utils::general::{array_max, array_min, cumsum};

//////////////////
// Type aliases //
//////////////////

////////////////
// Structures //
////////////////

/////////////
// Results //
/////////////

/// Structure to store GSEA stats
///
/// ### Fields
///
/// * `es` - Enrichment score
/// * `leading_edge` - The index positions of the leading edge genes
/// * `top` - Top points for plotting purposes
/// * `bottom` - Bottom points for plotting purposes
#[derive(Clone, Debug)]
pub struct GseaStats {
    pub es: f64,
    pub leading_edge: Vec<i32>,
    pub top: Vec<f64>,
    pub bottom: Vec<f64>,
}

/// Structure for final GSEA results from any algorithm
///
/// ### Fields
///
/// * `es` - Enrichment score
/// * `nes` - Normalised enrichment score
/// * `pvals` - The p-value for this pathway based on permutations
/// * `n_more_extreme` - The number of permutations with higher/lower ES
/// * `le_zero` - The number of permutation ES ≤ 0
/// * `ge_zero` - The number of permutation ES ≥ 0
/// * `size` - The size of the pathway
#[derive(Clone, Debug)]
pub struct GseaResults<'a> {
    pub es: &'a [f64],
    pub nes: Vec<Option<f64>>,
    pub pvals: Vec<f64>,
    pub n_more_extreme: Vec<usize>,
    pub le_zero: Vec<usize>,
    pub ge_zero: Vec<usize>,
    pub size: &'a [usize],
}

/// Structure for results from the different GSEA permutation methods
///
/// ### Fields
///
/// * `le_es` - The number of permutations less than the score
/// * `ge_es` - The number of permutations greater than the score
/// * `le_zero` - The number of permutation ES ≤ 0
/// * `ge_zero` - The number of permutation ES ≥ 0
/// * `le_zero_sum` - The sum of permuted enrichment scores that were ≤ 0
/// * `ge_zero_sum` - The sum of permuted enrichment scores that were ≥ 0
#[derive(Clone, Debug)]
pub struct GseaBatchResults {
    pub le_es: Vec<usize>,
    pub ge_es: Vec<usize>,
    pub le_zero: Vec<usize>,
    pub ge_zero: Vec<usize>,
    pub le_zero_sum: Vec<f64>,
    pub ge_zero_sum: Vec<f64>,
}

////////////
// Params //
////////////

/// Structure to store GSEA params
///
/// ### Fields
///
/// * `gsea_param` - The GSEA parameter, usually 1.0
/// * `max_size` - The maximum size of the allowed pathways
/// * `min_size` - The minimum size of the allowed pathways
#[derive(Clone, Debug)]
pub struct GseaParams {
    pub gsea_param: f64,
    pub max_size: usize,
    pub min_size: usize,
}

/// Prepare GSEA parameters from R list input
///
/// ### Params
///
/// * `r_list` - R list containing parameter values
///
/// ### Returns
///
/// `GseaParams` struct with parsed parameters (defaults: gsea_param=1.0,
/// min_size=5, max_size=500)
pub fn prepare_gsea_params(r_list: List) -> GseaParams {
    let gsea_params = r_list.into_hashmap();

    let gsea_param = gsea_params
        .get("gsea_param")
        .and_then(|v| v.as_real())
        .unwrap_or(1.0);

    let min_size = gsea_params
        .get("min_size")
        .and_then(|v| v.as_integer())
        .unwrap_or(5) as usize;

    let max_size = gsea_params
        .get("max_size")
        .and_then(|v| v.as_integer())
        .unwrap_or(500) as usize;

    GseaParams {
        gsea_param,
        max_size,
        min_size,
    }
}

//////////////////
// Segment tree //
//////////////////

/// Structure from the fgsea simple algorithm implementing a segment tree data structure
///
/// ### Fields
///
/// * `t` - Tree array for storing segment values
/// * `b` - Block array for storing block-level aggregates
/// * `n` - Total size of the tree (padded to power of 2)
/// * `k2` - Number of blocks in the structure
/// * `log_k` - Log base 2 of block size for bit operations
/// * `block_mask` - Bitmask for extracting block positions
#[derive(Clone, Debug)]
pub struct SegmentTree<T> {
    t: Vec<T>,
    b: Vec<T>,
    n: usize,
    k2: usize,
    log_k: usize,
    block_mask: usize,
}

impl<T> SegmentTree<T>
where
    T: Copy + std::ops::AddAssign + Default + std::ops::Add<Output = T>,
{
    /// Create a new segment tree with size n
    ///
    /// ### Params
    ///
    /// * `n_` - Size of the tree
    ///
    /// ### Returns
    ///
    /// Initialised structure
    pub fn new(n_: usize) -> Self {
        let mut k = 1;
        let mut log_k = 0;

        while k * k < n_ {
            k <<= 1;
            log_k += 1;
        }

        let k2 = (n_ - 1) / k + 1;
        let n = k * k;
        let block_mask = k - 1;
        let t = vec![T::default(); n];
        let b = vec![T::default(); k2];

        SegmentTree {
            t,
            b,
            n,
            k2,
            log_k,
            block_mask,
        }
    }

    /// Increment the value at position p by delta
    ///
    /// ### Params
    ///
    /// * `p` - Position to increment (mutable)
    /// * `delta` - Value to add
    ///
    /// ### Note
    ///
    /// p NEEDS to be mutable for the algorithm to work correctly
    pub fn increment(&mut self, mut p: usize, delta: T) {
        let block_end = p - (p & self.block_mask) + self.block_mask + 1;
        while p < block_end && p < self.n {
            self.t[p] += delta;
            p += 1;
        }
        let mut p1 = p >> self.log_k; // Now uses the modified p
        while p1 < self.k2 {
            self.b[p1] += delta;
            p1 += 1;
        }
    }

    /// Calculate the sum in range 0 to r
    ///
    /// ### Params
    ///
    /// * `r` - Right bound of range
    ///
    /// ### Returns
    ///
    /// Sum in range [0, r)
    pub fn query_r(&self, mut r: usize) -> T {
        if r == 0 {
            return T::default();
        }
        r -= 1;
        self.t[r] + self.b[r >> self.log_k]
    }
}

//////////////
// ES Ruler //
//////////////

/// Calculate the ES and leading edge genes
///
/// ### Params
///
/// * `stats` - Gene statistics
/// * `gs_idx` - Gene set indices
/// * `gsea_param` - GSEA parameter for weighting
/// * `return_leading_edge` - Whether to return leading edge genes
/// * `return_all_extreme` - Whether to return the all points for plotting
/// * `one_indexed` - Whether indices are one-based
///
/// ### Returns
///
/// Tuple of (gene statistic, leading edge genes)
pub fn calc_gsea_stats(
    stats: &[f64],
    gs_idx: &[i32],
    gsea_param: f64,
    return_leading_edge: bool,
    return_all_extreme: bool,
    one_indexed: bool,
) -> GseaStats {
    let n = stats.len();
    let m = gs_idx.len();
    let mut r_adj = Vec::with_capacity(m);
    for &i in gs_idx {
        let idx = if one_indexed { i - 1 } else { i } as usize;
        r_adj.push(stats[idx].abs().powf(gsea_param));
    }
    let nr: f64 = r_adj.iter().sum();
    let r_cum_sum: Vec<f64> = if nr == 0.0 {
        gs_idx
            .iter()
            .enumerate()
            .map(|(i, _)| i as f64 / r_adj.len() as f64)
            .collect()
    } else {
        cumsum(&r_adj).iter().map(|x| x / nr).collect()
    };

    let top_tmp: Vec<f64> = gs_idx
        .iter()
        .enumerate()
        .map(|(i, x)| (*x as f64 - (i + 1) as f64) / (n as f64 - m as f64))
        .collect();
    let tops: Vec<f64> = r_cum_sum
        .iter()
        .zip(top_tmp.iter())
        .map(|(x1, x2)| x1 - x2)
        .collect();
    let bottoms: Vec<f64> = if nr == 0.0 {
        tops.iter().map(|x| x - (1.0 / m as f64)).collect()
    } else {
        tops.iter()
            .zip(r_adj.iter())
            .map(|(top, adj)| top - (adj / nr))
            .collect()
    };
    let max_p = array_max(&tops);
    let min_p = array_min(&bottoms);

    let gene_stat = if max_p == -min_p {
        0.0
    } else if max_p > -min_p {
        max_p
    } else {
        min_p
    };
    // leading edges
    let leading_edge = if return_leading_edge {
        if max_p > -min_p {
            let max_idx = bottoms
                .iter()
                .enumerate()
                .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
                .map(|(idx, _)| idx)
                .unwrap_or(0);

            gs_idx.iter().take(max_idx + 1).cloned().collect()
        } else if max_p < -min_p {
            let min_idx = bottoms
                .iter()
                .enumerate()
                .min_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
                .map(|(idx, _)| idx)
                .unwrap_or(0);
            gs_idx.iter().skip(min_idx).cloned().rev().collect()
        } else {
            Vec::new()
        }
    } else {
        Vec::new()
    };
    // extreme points
    if return_all_extreme {
        GseaStats {
            es: gene_stat,
            leading_edge,
            top: tops,
            bottom: bottoms,
        }
    } else {
        GseaStats {
            es: gene_stat,
            leading_edge,
            top: Vec::new(),
            bottom: Vec::new(),
        }
    }
}

/// Convert order indices to ranks (from fgsea C++ implementation)
///
/// ### Params
///
/// * `order` - Ordering indices
///
/// ### Returns
///
/// Rank vector
pub fn ranks_from_order(order: &[usize]) -> Vec<i32> {
    let mut res = vec![0; order.len()];
    for (i, _) in order.iter().enumerate() {
        let idx = order[i];
        res[idx] = i as i32;
    }
    res
}

/// Returns a vector of indices that would sort the input slice
/// Implementation of the C++ code
///
/// ### Params
///
/// * `x` - Input slice to get ordering for
///
/// ### Returns
///
/// Vector of sorting indice
pub fn fgsea_order<T>(x: &[T]) -> Vec<usize>
where
    T: PartialOrd,
{
    let mut res: Vec<usize> = (0..x.len()).collect();
    res.sort_by(|&i, &j| x[i].partial_cmp(&x[j]).unwrap_or(std::cmp::Ordering::Equal));
    res
}

/// Extract values at specified indices from a vector
///
/// ### Params
///
/// * `from` - Source vector
/// * `indices` - Indices to extract (1-indexed)
///
/// # Returns
///
/// Option containing extracted values or None if invalid index
fn subvector(from: &[f64], indices: &[usize]) -> Option<Vec<f64>> {
    let mut result = Vec::with_capacity(indices.len());
    for &idx in indices {
        if idx > 0 && idx <= from.len() {
            result.push(from[idx - 1]);
        } else {
            return None;
        }
    }
    Some(result)
}

/// Generate random gene set indices using parallel processing
///
/// ### Params
///
/// * `iter_number` - Number of iterations/permutations
/// * `max_len` - Maximum length of each sample
/// * `universe_length` - Total number of genes
/// * `seed` - Random seed
/// * `one_indexed` - Whether to use 1-based indexing
///
/// # Returns
///
/// Vector of random index vectors
pub fn create_random_gs_indices(
    iter_number: usize,
    max_len: usize,
    universe_length: usize,
    seed: u64,
    one_indexed: bool,
) -> Vec<Vec<usize>> {
    // Parallelize creation of random samples
    (0..iter_number)
        .into_par_iter()
        .map(|i| {
            // Create a unique seed for each iteration
            let iter_seed = seed.wrapping_add(i as u64);
            let mut rng = StdRng::seed_from_u64(iter_seed);

            // Fisher-Yates shuffle for efficient sampling without replacement
            // Also deal with potential issues in indexing here
            let adjusted_universe = universe_length - 1;
            let actual_len = std::cmp::min(max_len, adjusted_universe);

            let mut indices: Vec<usize> = (0..adjusted_universe).collect();

            // Perform partial Fisher-Yates shuffle (we only need actual_len elements)
            for i in 0..actual_len {
                let j = rng.random_range(i..indices.len());
                indices.swap(i, j);
            }

            indices.truncate(actual_len);

            // Deal with R's 1-index, if need be
            if one_indexed {
                indices.iter_mut().for_each(|x| *x += 1);
            }

            indices
        })
        .collect()
}

/// Transform batch results into final GSEA results (es, nes, pval and size)
///
/// ### Params
///
/// * `pathway_scores` - Enrichment scores for pathways
/// * `pathway_sizes` - Sizes of pathways
/// * `gsea_res` - Batch results from permutations
///
/// ### Returns
///
/// Final GSEA results structure
pub fn calculate_nes_es_pval<'a>(
    pathway_scores: &'a [f64],
    pathway_sizes: &'a [usize],
    gsea_res: &GseaBatchResults,
) -> GseaResults<'a> {
    let le_zero_mean: Vec<f64> = gsea_res
        .le_zero_sum
        .iter()
        .zip(gsea_res.le_zero.iter())
        .map(|(a, b)| a / *b as f64)
        .collect();

    let ge_zero_mean: Vec<f64> = gsea_res
        .ge_zero_sum
        .iter()
        .zip(gsea_res.ge_zero.iter())
        .map(|(a, b)| a / *b as f64)
        .collect();

    let nes: Vec<Option<f64>> = pathway_scores
        .iter()
        .zip(ge_zero_mean.iter().zip(le_zero_mean.iter()))
        .map(|(&score, (&ge_mean, &le_mean))| {
            if (score > 0.0 && ge_mean != 0.0) || (score < 0.0 && le_mean != 0.0) {
                Some(score / if score > 0.0 { ge_mean } else { le_mean.abs() })
            } else {
                None
            }
        })
        .collect();

    let pvals: Vec<f64> = gsea_res
        .le_es
        .iter()
        .zip(gsea_res.le_zero.iter())
        .zip(gsea_res.ge_es.iter())
        .zip(gsea_res.ge_zero.iter())
        .map(|(((le_es, le_zero), ge_es), ge_zero)| {
            ((1 + le_es) as f64 / (1 + le_zero) as f64)
                .min((1 + ge_es) as f64 / (1 + ge_zero) as f64)
        })
        .collect();

    let n_more_extreme: Vec<usize> = pathway_scores
        .iter()
        .zip(gsea_res.ge_es.iter())
        .zip(gsea_res.le_es.iter())
        .map(|((es, ge_es), le_es)| if es > &0.0 { *ge_es } else { *le_es })
        .collect();

    GseaResults {
        es: pathway_scores,
        nes,
        pvals,
        n_more_extreme,
        le_zero: gsea_res.le_zero.clone(),
        ge_zero: gsea_res.ge_zero.clone(),
        size: pathway_sizes,
    }
}

////////////////////
// Main functions //
////////////////////

////////////////////////////////////////////
// Classical Gene Set enrichment analysis //
////////////////////////////////////////////

//////////////////
// FGSEA simple //
//////////////////

/// Square root heuristic approximation from fgsea paper with convex hull updates
///
/// Selected stats needs to be one-indexed!!!
///
/// ### Params
///
/// * `stats` - Gene statistics
/// * `selected_stats` - Selected gene indices (one-indexed)
/// * `selected_order` - Order of selected genes
/// * `gsea_param` - GSEA parameter for weighting
/// * `rev` - Whether to reverse direction
///
/// # Returns
///
/// Vector of enrichment scores
fn gsea_stats_sq(
    stats: &[f64],
    selected_stats: &[usize],
    selected_order: &[usize],
    gsea_param: f64,
    rev: bool,
) -> Vec<f64> {
    let n = stats.len() as i32;
    let k = selected_stats.len();
    let mut nr = 0.0;

    // Create results vector with 0
    let mut res = vec![0.0; k];

    // Create the SegmentTree like in the fgsea algorithm
    let mut xs = SegmentTree::<i32>::new(k + 1);
    let mut ys = SegmentTree::<f64>::new(k + 1);

    let mut selected_ranks = ranks_from_order(selected_order);

    // Initialize x-coordinates based on gene positions
    if !rev {
        // Forward direction with pre-ordered stats
        let mut prev: i32 = -1;
        for (i, _) in (0..k).enumerate() {
            let j = selected_order[i];
            let t = (selected_stats[j] - 1) as i32;
            assert!(t - prev >= 1);
            xs.increment(i, t - prev);
            prev = t;
        }
        xs.increment(k, n - 1 - prev);
    } else {
        // Reverse direction with pre-ordered stats
        let mut prev = n;
        for i in 0..k {
            selected_ranks[i] = k as i32 - 1 - selected_ranks[i];
            let j = selected_order[k - 1 - i];
            let t = (selected_stats[j] - 1) as i32;
            assert!(prev - t >= 1);
            xs.increment(i, prev - t);
            prev = t;
        }
        xs.increment(k, prev); // They had a random -0 here which clippy did not like
    }

    let mut st_prev = vec![0; k + 2];
    let mut st_next = vec![0; k + 2];

    let k1 = std::cmp::max(1, (k as f64 + 1.0).sqrt() as usize);
    let k2 = (k + 1) / k1 + 1;

    let mut block_summit = vec![0; k2];
    let mut block_start = vec![0; k2];
    let mut block_end = vec![0; k2];

    // Initialize blocks
    for i in (0..=k + 1).step_by(k1) {
        let block = i / k1;
        block_start[block] = i;
        block_end[block] = std::cmp::min(i + k1 - 1, k + 1);

        // Set up linked list for convex hull <- crazy stuff happening here...
        for j in 1..=block_end[block] - i {
            st_prev[i + j] = i + j - 1;
            st_next[i + j - 1] = i + j;
        }
        st_prev[i] = i;
        st_next[block_end[block]] = block_end[block];
        block_summit[block] = block_end[block];
    }

    // Calculate statistic epsilon for handling small values
    let mut stat_eps: f64 = 1e-5;
    for (i, _) in (0..k).enumerate() {
        let t = selected_stats[i];
        let xx = stats[t].abs();
        if xx > 0.0 {
            stat_eps = stat_eps.min(xx);
        }
    }
    stat_eps /= 1024.0;

    // Main algorithm loop - process each gene
    for i in 0..k {
        let t = selected_stats[i] - 1;
        let t_rank = selected_ranks[i] as usize;

        // Calculate adjusted statistic value with gene weight
        let adj_stat = stats[t].abs().max(stat_eps).powf(gsea_param);

        // Update coordinates
        xs.increment(t_rank, -1);
        ys.increment(t_rank, adj_stat);
        nr += adj_stat;

        let m = i + 1;

        // Calculate block indices
        let cur_block = (t_rank + 1) / k1;
        let bs = block_start[cur_block];
        let be = block_end[cur_block];

        // Recalculate convex hull for invalidated block
        let mut cur_top = t_rank.max(bs);

        for j in t_rank + 1..=be {
            let c = j;

            let xc = xs.query_r(c) as f64;
            let yc = ys.query_r(c);

            let mut b = cur_top;

            let mut xb = xs.query_r(b) as f64;
            let mut yb = ys.query_r(b);

            while st_prev[cur_top] != cur_top {
                let a = st_prev[cur_top];

                let xa = xs.query_r(a) as f64;
                let ya = ys.query_r(a);

                // Calculate cross product to determine turn direction
                let pr = (xb - xa) * (yc - yb) - (yb - ya) * (xc - xb);

                // Handle numerical precision issues
                let pr = if yc - ya < 1e-13 { 0.0 } else { pr };

                if pr <= 0.0 {
                    // Right turn - keep current point
                    break;
                }
                // Left turn - pop point from hull
                cur_top = a;
                st_next[b] = usize::MAX; // Mark as removed (-1 equivalent)
                b = a;
                xb = xa;
                yb = ya;
            }

            // Add point to convex hull
            st_prev[c] = cur_top;
            st_next[cur_top] = c;
            cur_top = c;

            if st_next[c] != usize::MAX {
                break;
            }
        }

        // Calculate coefficient for distance
        let coef = (n - m as i32) as f64 / nr;

        // Initialize max distance tracking
        let mut max_p: f64 = 0.0;

        // Update block summit
        block_summit[cur_block] = cur_top.max(block_summit[cur_block]);

        // Find farthest point across all blocks
        for (block, _) in (0..k2).collect::<std::vec::Vec<usize>>().iter().enumerate() {
            let mut cur_summit = block_summit[block];

            let mut cur_dist = ys.query_r(cur_summit) * coef - xs.query_r(cur_summit) as f64;

            // Find optimal point in this block
            loop {
                let next_summit = st_prev[cur_summit];
                let next_dist = ys.query_r(next_summit) * coef - xs.query_r(next_summit) as f64;

                if next_dist <= cur_dist {
                    break;
                }
                cur_dist = next_dist;
                cur_summit = next_summit;
            }

            block_summit[block] = cur_summit;
            max_p = max_p.max(if cur_dist.is_sign_negative() {
                0.0
            } else {
                cur_dist
            });
        }

        // Normalize and store result
        max_p /= (n - m as i32) as f64;

        res[i] = max_p;
    }

    res
}

/// Calculate cumulative enrichment scores via block-wise approximation
///
/// ### Params
///
/// * `stats` - Gene statistics
/// * `selected_stats` - Selected gene indices
/// * `gsea_param` - GSEA parameter for weighting
///
/// ### Returns
///
/// Vector of cumulative enrichment scores
pub fn calc_gsea_stat_cumulative(
    stats: &[f64],
    selected_stats: &[usize],
    gsea_param: f64,
) -> Vec<f64> {
    let selected_order = fgsea_order(selected_stats);

    // Calculate positive scores
    let res = gsea_stats_sq(stats, selected_stats, &selected_order, gsea_param, false);

    // Calculate negative scores
    let res_down = gsea_stats_sq(stats, selected_stats, &selected_order, gsea_param, true);

    // Combine results (take max magnitude with sign)
    res.into_iter()
        .zip(res_down)
        .map(|(up, down)| {
            if up == down {
                0.0
            } else if up < down {
                -down
            } else {
                up
            }
        })
        .collect()
}

/// Create permutations for the fgsea simple method
///
/// ### Params
///
/// * `stats` - Gene statistics
/// * `gsea_param` - GSEA parameter
/// * `iters` - Number of iterations
/// * `max_len` - Maximum pathway length
/// * `universe_length` - Total number of genes
/// * `seed` - Random seed
/// * `one_indexed` - Whether to use 1-based indexing
///
/// ### Returns
///
/// Vector of permutation enrichment score vectors
pub fn create_perm_es_simple(
    stats: &[f64],
    gsea_param: f64,
    iters: usize,
    max_len: usize,
    universe_length: usize,
    seed: u64,
    one_indexed: bool,
) -> Vec<Vec<f64>> {
    let shared_perm = create_random_gs_indices(iters, max_len, universe_length, seed, one_indexed);

    let chunk_size = std::cmp::max(1, iters / (rayon::current_num_threads() * 4));

    let rand_es: Vec<Vec<f64>> = shared_perm
        .into_par_iter()
        .chunks(chunk_size)
        .flat_map(|chunk| {
            let mut local_results = Vec::with_capacity(chunk.len());

            for i in chunk {
                let x = calc_gsea_stat_cumulative(stats, &i, gsea_param);
                local_results.push(x)
            }

            local_results
        })
        .collect();

    rand_es
}

/// Abstraction wrapper to be used in different parts of the package
///
/// ### Params
///
/// * `pathway_scores` - Pathway enrichment scores
/// * `pathway_sizes` - Pathway sizes
/// * `shared_perm` - Shared permutation results
///
/// ### Returns
///
/// Batch results from permutation analysis
pub fn calc_gsea_stats_wrapper(
    pathway_scores: &[f64],
    pathway_sizes: &[usize],
    shared_perm: &Vec<Vec<f64>>,
) -> GseaBatchResults {
    let m = pathway_scores.len();

    let mut le_es = vec![0; m];
    let mut ge_es = vec![0; m];
    let mut le_zero = vec![0; m];
    let mut ge_zero = vec![0; m];
    let mut le_zero_sum = vec![0.0; m];
    let mut ge_zero_sum = vec![0.0; m];

    // Process each permutation
    for rand_es_i in shared_perm {
        // Get the subvector once
        let rand_es_p = subvector(rand_es_i, pathway_sizes).unwrap();

        // Process all stats in a single pass for each pathway
        for i in 0..m {
            // Compare with pathway score
            if rand_es_p[i] <= pathway_scores[i] {
                le_es[i] += 1;
            } else {
                ge_es[i] += 1;
            }

            // Compare with zero
            if rand_es_p[i] <= 0.0 {
                le_zero[i] += 1;
                le_zero_sum[i] += rand_es_p[i]; // Already negative or zero
            } else {
                ge_zero[i] += 1;
                ge_zero_sum[i] += rand_es_p[i]; // Already positive
            }
        }
    }

    GseaBatchResults {
        le_es,
        ge_es,
        le_zero,
        ge_zero,
        le_zero_sum,
        ge_zero_sum,
    }
}
