use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use rayon::prelude::*;
use std::collections::HashSet;

use crate::utils_rust::array_max;

////////////////
// Structures //
////////////////

/// Structure for results from the Batch GSEA algorithm
#[derive(Clone, Debug)]
pub struct GseaBatchResults {
    pub le_es: Vec<usize>,
    pub ge_es: Vec<usize>,
    pub le_zero: Vec<usize>,
    pub ge_zero: Vec<usize>,
    pub le_zero_sum: Vec<f64>,
    pub ge_zero_sum: Vec<f64>,
}

/// Structure from the fgsea algorithm
#[derive(Clone, Debug)]
pub struct SegmentTree<T> {
    t: Vec<T>,
    b: Vec<T>,
    n: usize,
    k2: usize,
    log_k: usize,
    block_mask: usize,
}

/////////////////////
// Implementations //
/////////////////////

impl<T> SegmentTree<T>
where
    T: Copy + std::ops::AddAssign + Default + std::ops::Add<Output = T>,
{
    /// Create a new tree with size n
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
    /// p NEEDS to be mutable...
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
    pub fn query_r(&self, mut r: usize) -> T {
        if r == 0 {
            return T::default();
        }
        r -= 1;
        self.t[r] + self.b[r >> self.log_k]
    }
}

//////////////////////
// Helper functions //
//////////////////////

/// The ranks from order function from the fgsea C++ implementation
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
pub fn fgsea_order<T>(x: &[T]) -> Vec<usize>
where
    T: PartialOrd,
{
    let mut res: Vec<usize> = (0..x.len()).collect();
    res.sort_by(|&i, &j| x[i].partial_cmp(&x[j]).unwrap_or(std::cmp::Ordering::Equal));
    res
}

/// Subvector function to extract the relevant values
fn subvector(from: &[f64], indices: &[usize]) -> Result<Vec<f64>, String> {
    indices
        .iter()
        .map(|&idx| {
            if idx > 0 && idx <= from.len() {
                Ok(from[idx - 1])
            } else {
                Err(format!("Index out of bounds: {}", idx))
            }
        })
        .collect()
}

/// Get the indices of a given gene set. Assumes that gene universe
/// is based on the sorted gene level statistics
pub fn get_gene_set_indices(gene_universe: &[String], gs_genes: &[String]) -> Vec<usize> {
    let gs_set: HashSet<&String> = gs_genes.iter().collect();
    gene_universe
        .iter()
        .enumerate()
        .filter_map(|(a, b)| if gs_set.contains(b) { Some(a) } else { None })
        .collect()
}

/// Generate random gene set indices.
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
            let adjusted_universe = universe_length - 1;
            let actual_len = std::cmp::min(max_len, adjusted_universe);

            // Create array of indices
            let mut indices: Vec<usize> = (0..adjusted_universe).collect();

            // Perform partial Fisher-Yates shuffle (we only need actual_len elements)
            for i in 0..actual_len {
                let j = rng.random_range(i..indices.len());
                indices.swap(i, j);
            }

            // Take only the shuffled portion and adjust indexing if needed
            indices.truncate(actual_len);

            // Deal with R's 1-index, if need be
            if one_indexed {
                indices.iter_mut().for_each(|x| *x += 1);
            }

            indices
        })
        .collect()
}

////////////////////
// Main functions //
////////////////////

////////////////////////////////////////////
// Classical Gene Set enrichment analysis //
////////////////////////////////////////////

/// Calculate the Enrichment score. The function assumes
/// that stats is sorted and pathway are the index positions
/// of the genes in the pathway
pub fn calculate_es(stats: &[f64], pathway: &[usize]) -> f64 {
    // stats and p must be sorted
    let no_genes = stats.len();
    let p_total = pathway.len();
    let mut nr = 0.0;
    for p in pathway {
        nr += stats[*p].abs()
    }
    let weight_hit = 1.0 / nr;
    let weight_miss = 1.0 / (no_genes - p_total) as f64;
    let mut running_sum = 0.0;
    let mut max_run_sum = 0.0;
    let mut min_run_sum = 0.0;
    for (i, x) in stats.iter().enumerate() {
        if pathway.contains(&i) {
            running_sum += x.abs() * weight_hit
        } else {
            running_sum -= weight_miss
        }
        max_run_sum = f64::max(running_sum, max_run_sum);
        min_run_sum = f64::min(running_sum, min_run_sum);
    }
    if max_run_sum > min_run_sum.abs() {
        max_run_sum
    } else {
        min_run_sum
    }
}

/// Calculate the p-value given a permutation es vector and the normalised enrichment score.
/// This is the approach to calculate the p-values in the traditional GSEA.
pub fn calculate_pval(actual_es: f64, perm_es: &[f64], n_perm: usize, pos: bool) -> f64 {
    let count = if pos {
        perm_es.iter().filter(|val| val >= &&actual_es).count()
    } else {
        perm_es.iter().filter(|val| val <= &&actual_es).count()
    };
    (count as f64 + 1.0) / (n_perm as f64 + 1.0) * 2.0
}

/// Calculate the normalised enrichment score
pub fn calculate_nes(actual_es: f64, perm_es: &[f64], pos: bool) -> f64 {
    let filtered: Vec<f64> = if pos {
        perm_es.iter().filter(|x| x >= &&0.0).copied().collect()
    } else {
        perm_es.iter().filter(|x| x <= &&0.0).copied().collect()
    };
    let sum: f64 = filtered.iter().sum();
    let mean = (sum / filtered.len() as f64).abs();
    actual_es / mean
}

//////////////////
// FGSEA simple //
//////////////////

/// Crazy approximation from the fgsea paper leveraging a square root heuristic
/// and convex hull updates translated into Rust
/// Selected stats needs to be one-indexed!!!
pub fn gsea_stats_sq(
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
        for (i, _) in (0..k).collect::<std::vec::Vec<usize>>().iter().enumerate() {
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
    for (i, _) in (0..k).collect::<std::vec::Vec<usize>>().iter().enumerate() {
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

/// Calculate the cumulative enrichment scores via the
/// block-wise approximation
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
    let mut final_res = vec![0.0; res.len()];
    for i in 0..res.len() {
        if res[i] == res_down[i] {
            final_res[i] = 0.0;
        } else if res[i] < res_down[i] {
            final_res[i] = -res_down[i];
        } else {
            final_res[i] = res[i];
        }
    }

    final_res
}

pub fn calc_gsea_stat_cumulative_batch(
    stats: &[f64],
    pathway_scores: &[f64],
    pathway_sizes: &[usize],
    iters: usize,
    gsea_param: f64,
    seed: u64,
) -> Result<GseaBatchResults, String> {
    let n = stats.len();
    let k = array_max(pathway_sizes);

    let m = pathway_scores.len();

    let shared_perm = create_random_gs_indices(iters, k, n, seed, true);

    let rand_es: Vec<Vec<f64>> = shared_perm
        .par_iter()
        .map(|selected_genes_random| {
            calc_gsea_stat_cumulative(stats, selected_genes_random, gsea_param)
        })
        .collect();

    let mut le_es = vec![0; m];
    let mut ge_es = vec![0; m];
    let mut le_zero = vec![0; m];
    let mut ge_zero = vec![0; m];
    let mut le_zero_sum = vec![0.0; m];
    let mut ge_zero_sum = vec![0.0; m];

    for rand_es_i in rand_es {
        let rand_es_p = subvector(&rand_es_i, pathway_sizes)?;

        let aux: Vec<bool> = rand_es_p
            .iter()
            .zip(pathway_scores.iter())
            .map(|(a, b)| a <= b)
            .collect();
        let diff: Vec<usize> = aux.iter().map(|&x| if x { 1 } else { 0 }).collect();
        for i in 0..m {
            le_es[i] += diff[i];
        }

        let diff: Vec<usize> = aux.iter().map(|&x| if x { 0 } else { 1 }).collect();
        for i in 0..m {
            ge_es[i] += diff[i];
        }

        let aux: Vec<bool> = rand_es_p.iter().map(|&a| a <= 0.0).collect();

        let diff: Vec<usize> = aux.iter().map(|&x| if x { 1 } else { 0 }).collect();
        for i in 0..m {
            le_zero[i] += diff[i];
        }

        let diff: Vec<usize> = aux.iter().map(|&x| if x { 0 } else { 1 }).collect();
        for i in 0..m {
            ge_zero[i] += diff[i];
        }

        for i in 0..m {
            le_zero_sum[i] += rand_es_p[i].min(0.0);
        }

        for i in 0..m {
            ge_zero_sum[i] += rand_es_p[i].max(0.0);
        }
    }

    Ok(GseaBatchResults {
        le_es,
        ge_es,
        le_zero,
        ge_zero,
        le_zero_sum,
        ge_zero_sum,
    })
}

//////////////////////
// FGSEA multilevel //
//////////////////////

///////////////
// Dead code //
///////////////

#[allow(dead_code)]
/// Helper function to extract the maximum absolute ES while preserving the sign
pub fn get_max_es(es_vec: &[f64]) -> f64 {
    let (max_val, _) = es_vec
        .iter()
        .map(|&x| (x, x.abs()))
        .max_by(|(_, abs_a), (_, abs_b)| {
            abs_a
                .partial_cmp(abs_b)
                .unwrap_or(std::cmp::Ordering::Equal)
        })
        .unwrap();

    max_val
}
