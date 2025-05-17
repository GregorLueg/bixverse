use rand::distr::Uniform;
use rand::rngs::StdRng;
use rand::{seq::SliceRandom, Rng, SeedableRng};
use rand_distr::Distribution;
use rayon::prelude::*;
use statrs::function::gamma::digamma;
use std::collections::HashMap;

use crate::utils_rust::{array_max, array_min, cumsum, unique};

//////////////////
// Type aliases //
//////////////////

////////////////
// Structures //
////////////////

/// Structure for final GSEA results from any algorithm
#[derive(Clone, Debug)]
pub struct GseaResults<'a> {
    pub es: &'a [f64],
    pub nes: Vec<Option<f64>>,
    pub pvals: Vec<f64>,
    pub n_more_extreme: Vec<usize>,
    pub size: &'a [usize],
}

/// Structure for results from the different GSEA permutation methods
#[derive(Clone, Debug)]
pub struct GseaBatchResults {
    pub le_es: Vec<usize>,
    pub ge_es: Vec<usize>,
    pub le_zero: Vec<usize>,
    pub ge_zero: Vec<usize>,
    pub le_zero_sum: Vec<f64>,
    pub ge_zero_sum: Vec<f64>,
}

/// Structure for results from the different GSEA permutation methods
#[derive(Clone, Debug)]
pub struct GseaMultiLevelresults {
    pub pvals: Vec<f64>,
    pub is_cp_ge_half: Vec<bool>,
}

/// Structure from the fgsea simple algorithm
#[derive(Clone, Debug)]
pub struct SegmentTree<T> {
    t: Vec<T>,
    b: Vec<T>,
    n: usize,
    k2: usize,
    log_k: usize,
    block_mask: usize,
}

/// Structure for fgsea multi-level
/// Stores and manages chunked samples for efficient processing

#[derive(Clone, Debug)]
struct SampleChunks {
    chunk_sum: Vec<f64>,
    chunks: Vec<Vec<i32>>,
}

/// Further structure for fgsea multi-level
/// This is EsRuler implementation

#[derive(Clone, Debug)]
struct EsRuler {
    ranks: Vec<f64>,                  // Gene rankings
    sample_size: usize,               // Number of random samples
    pathway_size: usize,              // Size of gene set
    current_samples: Vec<Vec<usize>>, // Current set of samples
    enrichment_scores: Vec<f64>,      // Calculated ES values
    prob_corrector: Vec<usize>,       // Correction factors for p-values
    chunks_number: i32,               // Number of chunks for optimization
    chunk_last_element: Vec<i32>,     // Last element index in each chunk
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

impl SampleChunks {
    /// Creates new SampleChunks with specified number of chunks
    fn new(chunks_number: i32) -> Self {
        Self {
            chunk_sum: vec![0.0; chunks_number as usize],
            chunks: vec![Vec::new(); chunks_number as usize],
        }
    }
}

impl EsRuler {
    /// Creates a new ES ruler
    fn new(inp_ranks: &[f64], inp_sample_size: usize, inp_pathway_size: usize) -> Self {
        let mut current_samples: Vec<Vec<usize>> = Vec::with_capacity(inp_sample_size);
        current_samples.resize_with(inp_sample_size, Vec::new);

        Self {
            ranks: inp_ranks.to_vec(),
            sample_size: inp_sample_size,
            pathway_size: inp_pathway_size,
            current_samples,
            enrichment_scores: Vec::new(),
            prob_corrector: Vec::new(),
            chunks_number: 0,
            chunk_last_element: Vec::new(),
        }
    }

    /// Removes samples with low ES and duplicates samples with high ES
    /// This is to drive the sampling process to higher and higher ES
    fn duplicate_samples(&mut self) {
        let mut stats: Vec<(f64, usize)> = vec![(0.0, 0); self.sample_size];
        let mut pos_es_indxs: Vec<usize> = Vec::new();
        let mut total_pos_es_count = 0;

        // Calculate ES for each sample
        let loop_iter = 0..self.sample_size; // Avoids clippy complaints
        for sample_id in loop_iter {
            let sample_es_pos = calc_positive_es(&self.ranks, &self.current_samples[sample_id]);
            let sample_es = calc_es(&self.ranks, &self.current_samples[sample_id]);
            if sample_es > 0.0 {
                total_pos_es_count += 1;
                pos_es_indxs.push(sample_id);
            }
            stats[sample_id] = (sample_es_pos, sample_id);
        }

        stats.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

        // Store lower half ES values and track positive ES counts
        let loop_iter = 0..(self.sample_size) / 2;
        for sample_id in loop_iter {
            self.enrichment_scores.push(stats[sample_id].0);
            if pos_es_indxs.contains(&stats[sample_id].1) {
                total_pos_es_count -= 1;
            }
            self.prob_corrector.push(total_pos_es_count);
        }

        // Create new samples by duplicating high-ES samples
        let mut new_sets = Vec::new();
        for sample_id in 0..(self.sample_size - 2) / 2 {
            for _ in 0..2 {
                new_sets.push(
                    self.current_samples[stats[self.sample_size - 1 - sample_id].1 as usize]
                        .clone(),
                );
            }
        }
        new_sets.push(self.current_samples[stats[self.sample_size >> 1].1 as usize].clone());

        // Replace current samples with new sets <- Interesting command to swap stuff in memory
        std::mem::swap(&mut self.current_samples, &mut new_sets);
        self.sample_size = self.current_samples.len(); // VERY importand to add...
    }

    /// Attempts to improve a sample by swapping genes in/out
    /// Returns number of successful perturbations
    #[allow(unused_assignments)]
    fn perturbate(
        &self,
        ranks: &[f64],
        k: i32,
        sample_chunks: &mut SampleChunks,
        bound: f64,
        rng: &mut StdRng,
    ) -> i32 {
        let pert_prmtr = 0.1; // Controls perturbation intensity
        let n = ranks.len() as i32;
        let uid_n = Uniform::new(0, n).unwrap();
        let uid_k = Uniform::new(0, k).unwrap();

        // Calculate initial rank sum
        let mut ns = 0.0;
        for i in 0..sample_chunks.chunks.len() {
            for &pos in &sample_chunks.chunks[i] {
                ns += ranks[pos as usize];
            }
        }

        let q1 = 1.0 / (n as f64 - k as f64);
        let iters = std::cmp::max(1, (k as f64 * pert_prmtr) as i32);
        let mut moves = 0;

        let mut cand_val = -1;
        let mut has_cand = false;
        let mut cand_x = 0;
        let mut cand_y = 0.0;

        // Try multiple perturbations
        for _ in 0..iters {
            // Select random position to replace
            let old_ind = uid_k.sample(rng);
            let mut old_chunk_ind = 0;
            let mut old_ind_in_chunk = 0;
            let old_val;

            // Find chunk containing the position
            {
                let mut tmp = old_ind;
                while old_chunk_ind < sample_chunks.chunks.len()
                    && sample_chunks.chunks[old_chunk_ind].len() as i32 <= tmp
                {
                    tmp -= sample_chunks.chunks[old_chunk_ind].len() as i32;
                    old_chunk_ind += 1;
                }
                old_ind_in_chunk = tmp as usize;
                old_val = sample_chunks.chunks[old_chunk_ind][old_ind_in_chunk];
            }

            // Select random new value
            let new_val = uid_n.sample(rng);

            // Find insertion position
            let new_chunk_ind = match self.chunk_last_element.binary_search(&new_val) {
                Ok(idx) => idx,
                Err(idx) => idx,
            };

            let new_ind_in_chunk = match sample_chunks.chunks[new_chunk_ind].binary_search(&new_val)
            {
                Ok(idx) => idx,
                Err(idx) => idx,
            };

            // Skip if new value already in sample
            if new_ind_in_chunk < sample_chunks.chunks[new_chunk_ind].len()
                && sample_chunks.chunks[new_chunk_ind][new_ind_in_chunk] == new_val
            {
                if new_val == old_val {
                    moves += 1;
                }
                continue;
            }

            // Remove old and insert new value
            sample_chunks.chunks[old_chunk_ind].remove(old_ind_in_chunk);
            let adjust = if old_chunk_ind == new_chunk_ind && old_ind_in_chunk < new_ind_in_chunk {
                1
            } else {
                0
            };
            sample_chunks.chunks[new_chunk_ind].insert(new_ind_in_chunk - adjust, new_val);

            // Update sums
            ns = ns - ranks[old_val as usize] + ranks[new_val as usize];
            sample_chunks.chunk_sum[old_chunk_ind] -= ranks[old_val as usize];
            sample_chunks.chunk_sum[new_chunk_ind] += ranks[new_val as usize];

            // Update candidate tracking
            if has_cand {
                // This is how I understand the clippy complain, but really not sure about this
                // To revisit if stuff breaks
                match old_val.cmp(&cand_val) {
                    std::cmp::Ordering::Equal => {
                        has_cand = false;
                    }
                    std::cmp::Ordering::Less => {
                        cand_x += 1;
                        cand_y -= ranks[old_val as usize];
                    }
                    std::cmp::Ordering::Greater => {
                        // No action required for this case
                    }
                }

                if new_val < cand_val {
                    cand_x -= 1;
                    cand_y += ranks[new_val as usize];
                }
            }

            let q2 = 1.0 / ns;

            // Quick validity check with candidate
            if has_cand && -q1 * cand_x as f64 + q2 * cand_y > bound {
                moves += 1;
                continue;
            }

            // Check if new configuration is valid by finding a leading edge
            let mut cur_x = 0;
            let mut cur_y = 0.0;
            let mut ok = false;
            let mut last = -1;

            // Check each chunk
            for i in 0..sample_chunks.chunks.len() {
                if q2 * (cur_y + sample_chunks.chunk_sum[i]) - q1 * (cur_x as f64) < bound {
                    // Skip this chunk - not promising
                    cur_y += sample_chunks.chunk_sum[i];
                    cur_x += self.chunk_last_element[i]
                        - last
                        - 1
                        - sample_chunks.chunks[i].len() as i32;
                    last = self.chunk_last_element[i] - 1;
                } else {
                    // Check each position in chunk
                    for &pos in &sample_chunks.chunks[i] {
                        cur_y += ranks[pos as usize];
                        cur_x += pos - last - 1;
                        if q2 * cur_y - q1 * cur_x as f64 > bound {
                            // Found valid leading edge
                            ok = true;
                            has_cand = true;
                            cand_x = cur_x;
                            cand_y = cur_y;
                            cand_val = pos;
                            break;
                        }
                        last = pos;
                    }
                    if ok {
                        break;
                    }
                    cur_x += self.chunk_last_element[i] - 1 - last;
                    last = self.chunk_last_element[i] - 1;
                }
            }

            // If not valid, revert changes
            if !ok {
                ns = ns - ranks[new_val as usize] + ranks[old_val as usize];
                sample_chunks.chunk_sum[old_chunk_ind] += ranks[old_val as usize];
                sample_chunks.chunk_sum[new_chunk_ind] -= ranks[new_val as usize];

                sample_chunks.chunks[new_chunk_ind].remove(new_ind_in_chunk - adjust);
                sample_chunks.chunks[old_chunk_ind].insert(old_ind_in_chunk, old_val);

                // Revert candidate tracking
                if has_cand {
                    if new_val == cand_val {
                        has_cand = false;
                    } else if old_val < cand_val {
                        cand_x -= 1;
                        cand_y += ranks[old_val as usize];
                    }

                    if new_val < cand_val {
                        cand_x += 1;
                        cand_y -= ranks[new_val as usize];
                    }
                }
            } else {
                moves += 1; // Count successful move
            }
        }

        moves
    }

    /// Extends the ES distribution to include the target ES value
    /// Uses an adaptive sampling approach to explore higher ES values
    fn extend(&mut self, es: f64, seed: u64, eps: f64) {
        let mut rng = StdRng::seed_from_u64(seed);

        for sample_id in 0..self.sample_size {
            self.current_samples[sample_id] =
                combination(0, self.ranks.len() - 1, self.pathway_size, &mut rng);
            self.current_samples[sample_id].sort();
            let _ = calc_es(&self.ranks, &self.current_samples[sample_id]);
        }

        // Setup chunk processing structure for optimization
        self.chunks_number = std::cmp::max(1, (self.pathway_size as f64).sqrt() as i32);
        self.chunk_last_element = vec![0; self.chunks_number as usize];
        self.chunk_last_element[self.chunks_number as usize - 1] = self.ranks.len() as i32;
        let mut tmp: Vec<i32> = vec![0; self.sample_size];
        let mut samples_chunks = vec![SampleChunks::new(self.chunks_number); self.sample_size];

        // Initial sample processing
        self.duplicate_samples();
        while self.enrichment_scores.last().unwrap_or(&0.0) <= &(es - 1e-10) {
            // Update chunk boundaries
            for i in 0..self.chunks_number - 1 {
                let pos = (self.pathway_size as i32 + i) / self.chunks_number;
                let loop_iter = 0..self.sample_size;

                for j in loop_iter {
                    tmp[j] = self.current_samples[j][pos as usize] as i32;
                }

                nth_element(&mut tmp, self.sample_size / 2);
                self.chunk_last_element[i as usize] = tmp[self.sample_size / 2];
            }

            // Reorganize samples into chunks
            let loop_iter = 0..self.sample_size;
            for i in loop_iter {
                for j in 0..self.chunks_number as usize {
                    samples_chunks[i].chunk_sum[j] = 0.0;
                    samples_chunks[i].chunks[j].clear();
                }

                let mut cnt = 0;
                for &pos in &self.current_samples[i] {
                    while cnt < self.chunk_last_element.len()
                        && self.chunk_last_element[cnt] <= pos as i32
                    {
                        cnt += 1;
                    }
                    samples_chunks[i].chunks[cnt].push(pos as i32);
                    samples_chunks[i].chunk_sum[cnt] += self.ranks[pos];
                }
            }

            // Perturbate samples to improve ES
            let mut moves = 0;
            while moves < (self.sample_size * self.pathway_size) as i32 {
                let loop_iter = 0..self.sample_size;
                for sample_id in loop_iter {
                    moves += self.perturbate(
                        &self.ranks,
                        self.pathway_size as i32,
                        &mut samples_chunks[sample_id],
                        *self.enrichment_scores.last().unwrap_or(&0.0),
                        &mut rng,
                    );
                }
            }

            // Reconstruct samples from chunks
            let loop_iter = 0..self.sample_size;
            for i in loop_iter {
                self.current_samples[i].clear();
                for j in 0..self.chunks_number as usize {
                    self.current_samples[i]
                        .extend(samples_chunks[i].chunks[j].iter().map(|&x| x as usize));
                }
            }

            // Check progress and stop conditions
            let prev_top_score = *self.enrichment_scores.last().unwrap_or(&0.0);
            self.duplicate_samples();
            if self.enrichment_scores.last().unwrap_or(&0.0) <= &prev_top_score {
                // No improvement, stop iterating
                break;
            }

            // Stop if precision requirement met
            if eps != 0.0 {
                let k = self.enrichment_scores.len() / ((self.sample_size + 1) / 2);
                if k as f64 > -f64::log2(0.5 * eps) {
                    break;
                }
            }
        }
    }

    fn get_pval(&self, es: f64, _eps: f64, sign: bool) -> (f64, bool) {
        let half_size = (self.sample_size + 1) / 2;
        let it_index;
        let mut good_error = true;

        // Find position in the ES distribution
        if es >= *self.enrichment_scores.last().unwrap_or(&0.0) {
            it_index = self.enrichment_scores.len() - 1;
            if es > self.enrichment_scores[it_index] + 1e-10 {
                good_error = false; // Beyond measured range
            }
        } else {
            // Binary search for position
            it_index = match self
                .enrichment_scores
                .binary_search_by(|probe| probe.partial_cmp(&es).unwrap())
            {
                Ok(index) => index,
                Err(index) => index, // Insert position if not found exactly
            };
        }

        let indx = if it_index > 0 { it_index } else { 0 };
        let k = indx / half_size;
        let remainder = self.sample_size - (indx % half_size);

        let adj_log = beta_mean_log(half_size, self.sample_size);
        let adj_log_pval = k as f64 * adj_log + beta_mean_log(remainder + 1, self.sample_size);

        if sign {
            (0.0_f64.max(1.0_f64.min(adj_log_pval.exp())), good_error)
        } else {
            // Apply correction for multiple testing
            let correction = calc_log_correction(&self.prob_corrector, indx, self.sample_size);
            let res_log = adj_log_pval + correction.0;
            (
                0.0_f64.max(1.0_f64.min(res_log.exp())),
                good_error && correction.1,
            )
        }
    }
}

//////////////////////
// Helper functions //
//////////////////////

/// Calculate the enrichment score (based on the fgsea C++ implementation)
fn calc_es(ranks: &[f64], pathway_indices: &[usize]) -> f64 {
    let mut ns = 0.0;
    for p in pathway_indices {
        ns += ranks[*p]
    }
    let n = ranks.len();
    let k = pathway_indices.len();
    let mut res: f64 = 0.0;
    let mut cur: f64 = 0.0;
    let q1 = 1.0 / (n - k) as f64;
    let q2 = 1.0 / ns;
    let mut last: i64 = -1;
    for p in pathway_indices {
        cur -= q1 * (*p as i64 - last - 1) as f64;
        if cur.abs() > res.abs() {
            res = cur;
        }
        cur += q2 * ranks[*p];
        if cur.abs() > res.abs() {
            res = cur;
        }
        last = *p as i64;
    }
    res
}

/// Calculate the positive enrichment score (based on the fgsea C++ implementation)
fn calc_positive_es(ranks: &[f64], pathway_indices: &[usize]) -> f64 {
    let mut ns = 0.0;
    for p in pathway_indices {
        ns += ranks[*p]
    }
    let n = ranks.len();
    let k = pathway_indices.len();
    let mut res: f64 = 0.0;
    let mut cur: f64 = 0.0;
    let q1 = 1.0 / (n - k) as f64;
    let q2 = 1.0 / ns;
    let mut last: i64 = -1;
    for p in pathway_indices {
        cur += q2 * ranks[*p] - q1 * (*p as i64 - last - 1) as f64;
        res = res.max(cur);
        last = *p as i64;
    }
    res
}

/// Calculate the ES and leading edge genes
pub fn calc_gsea_stats(
    stats: &[f64],
    gs_idx: &[i32],
    gsea_param: f64,
    return_leading_edge: bool,
    one_indexed: bool,
) -> (f64, Vec<i32>) {
    let n = stats.len();
    let m = gs_idx.len();
    let r_adj: Vec<f64> = gs_idx
        .iter()
        .map(|i| {
            let idx = if one_indexed { *i - 1 } else { *i } as usize;
            stats[idx].abs().powf(gsea_param)
        })
        .collect();
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

    (gene_stat, leading_edge)
}

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

/// Transform the batch results into final GSEA results,
/// i.e., es, nes, pval and size
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
        size: pathway_sizes,
    }
}

/// Generate k random numbers from [a, b] (inclusive range)
/// Uses different algorithms for sparse (k < n/2) and dense (k >= n/2) sampling
fn combination(a: usize, b: usize, k: usize, rng: &mut impl Rng) -> Vec<usize> {
    let n = b - a + 1;
    if k > n {
        panic!("k cannot be greater than range size n");
    }

    let mut v = Vec::with_capacity(k);
    let mut used = vec![false; n];

    if (k as f64) < (n as f64) * 0.5 {
        // Sparse sampling: randomly select k unique elements
        let range = Uniform::new_inclusive(a, b).unwrap();

        for _ in 0..k {
            // Try up to 100 times to find an unused element
            for _ in 0..100 {
                let x = range.sample(rng);
                if !used[x - a] {
                    v.push(x);
                    used[x - a] = true;
                    break;
                }
            }
        }
    } else {
        // Dense sampling: sample n-k positions to exclude
        for r in (n - k)..n {
            let range = Uniform::new(0, r as i32).unwrap();
            let x = range.sample(rng) as usize;

            if !used[x] {
                v.push(a + x);
                used[x] = true;
            } else {
                v.push(a + r);
                used[r] = true;
            }
        }

        // Shuffle the resulting vector
        v.shuffle(rng);
    }

    v
}

/// Rearranges array so nth element is in its sorted position
fn nth_element(arr: &mut [i32], n: usize) {
    if arr.is_empty() || n >= arr.len() {
        return;
    }
    arr.select_nth_unstable(n);
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

/// Calculate once for each size the permutation-based enrichment scores.
fn create_perm_es(
    stats: &[f64],
    gene_set_sizes: &[usize],
    shared_perms: &[Vec<usize>],
) -> HashMap<usize, Vec<f64>> {
    let mut shared_perm_es = HashMap::with_capacity(gene_set_sizes.len());
    for size in gene_set_sizes {
        let perm_es: Vec<f64> = shared_perms
            .into_par_iter()
            .map(|perm| calculate_es(stats, &perm[..*size]))
            .collect();
        shared_perm_es.insert(*size, perm_es);
    }
    shared_perm_es
}

/// Calculate the permutations in the 'traditional' way.
pub fn calc_gsea_stat_traditional_batch(
    stats: &[f64],
    pathway_scores: &[f64],
    pathway_sizes: &[usize],
    iters: usize,
    // gsea_param: f64,
    seed: u64,
) -> GseaBatchResults {
    let n = stats.len();
    let k = array_max(pathway_sizes);
    let k_unique = unique(pathway_sizes);

    let m = pathway_scores.len();

    let shared_perm = create_random_gs_indices(iters, k, n, seed, false);

    let shared_perm_es = create_perm_es(stats, &k_unique, &shared_perm);

    let mut le_es = Vec::with_capacity(m);
    let mut ge_es = Vec::with_capacity(m);
    let mut le_zero = Vec::with_capacity(m);
    let mut ge_zero = Vec::with_capacity(m);
    let mut le_zero_sum = Vec::with_capacity(m);
    let mut ge_zero_sum = Vec::with_capacity(m);

    for (es, size) in pathway_scores.iter().zip(pathway_sizes.iter()) {
        let random_es = shared_perm_es.get(size).unwrap();
        let le_es_i = random_es.iter().filter(|x| x <= &es).count();
        let ge_es_i = iters - le_es_i;
        let le_zero_i = random_es.iter().filter(|x| x <= &&0.0).count();
        let ge_zero_i = iters - le_zero_i;
        let le_zero_sum_i: f64 = random_es.iter().map(|x| x.min(0.0)).sum();
        let ge_zero_sum_i: f64 = random_es.iter().map(|x| x.max(0.0)).sum();
        le_es.push(le_es_i);
        ge_es.push(ge_es_i);
        le_zero.push(le_zero_i);
        ge_zero.push(ge_zero_i);
        le_zero_sum.push(le_zero_sum_i);
        ge_zero_sum.push(ge_zero_sum_i);
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

//////////////////
// FGSEA simple //
//////////////////

/// Crazy approximation from the fgsea paper leveraging a square root heuristic
/// and convex hull updates translated into Rust
/// Selected stats needs to be one-indexed!!!
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

/// Create the permutations for the fgsea simple method
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

    let rand_es: Vec<Vec<f64>> = shared_perm
        .par_iter()
        .map(|selected_genes_random| {
            calc_gsea_stat_cumulative(stats, selected_genes_random, gsea_param)
        })
        .collect();

    rand_es
}

/// Abstraction wrapper to be used in different parts of the package
pub fn calc_gsea_stats_wrapper(
    pathway_scores: &[f64],
    pathway_sizes: &[usize],
    shared_perm: &Vec<Vec<f64>>,
) -> Result<GseaBatchResults, String> {
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
        let rand_es_p = subvector(rand_es_i, pathway_sizes)?;

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

    Ok(GseaBatchResults {
        le_es,
        ge_es,
        le_zero,
        ge_zero,
        le_zero_sum,
        ge_zero_sum,
    })
}

/// Calculate random scores batch-wise
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

    let shared_perm = create_perm_es_simple(stats, gsea_param, iters, k, n, seed, true);

    calc_gsea_stats_wrapper(pathway_scores, pathway_sizes, &shared_perm)
}

//////////////////////
// FGSEA multilevel //
//////////////////////

/// Helper to calculate the beta mean log
fn beta_mean_log(a: usize, b: usize) -> f64 {
    digamma(a as f64) - digamma((b + 1) as f64)
}

/// Calculates the log corrections
fn calc_log_correction(
    prob_corrector: &[usize],
    prob_corr_idx: usize,
    sample_size: usize,
) -> (f64, bool) {
    let mut result = 0.0;
    let half_size = (sample_size + 1) / 2;
    let remainder = sample_size - (prob_corr_idx % half_size);
    let cond_prob = beta_mean_log(prob_corrector[prob_corr_idx] + 1, remainder);
    result += cond_prob;
    if cond_prob.exp() >= 0.5 {
        (result, true)
    } else {
        (result, false)
    }
}

/// Function to do the multi-level magic in fgsea
pub fn fgsea_multilevel(
    enrichment_scores: &[f64],
    ranks: &[f64],
    pathway_size: usize,
    sample_size: usize,
    seed: u64,
    eps: f64,
    sign: bool,
) -> GseaMultiLevelresults {
    let pos_ranks: Vec<f64> = ranks.iter().map(|&r| r.abs()).collect();
    let mut neg_ranks = pos_ranks.clone();
    neg_ranks.reverse();

    // Initialise the EsRulers
    let mut es_ruler_pos = EsRuler::new(&pos_ranks, sample_size, pathway_size);
    let mut es_ruler_neg = EsRuler::new(&neg_ranks, sample_size, pathway_size);

    let max_es = enrichment_scores
        .iter()
        .fold(f64::NEG_INFINITY, |a, &b| a.max(b));
    let min_es = enrichment_scores
        .iter()
        .fold(f64::INFINITY, |a, &b| a.min(b));

    if max_es >= 0.0 {
        es_ruler_pos.extend(max_es.abs(), seed, eps);
    }

    if min_es < 0.0 {
        es_ruler_neg.extend(min_es.abs(), seed, eps);
    }

    let mut pval_res = Vec::with_capacity(enrichment_scores.len());
    let mut is_cp_ge_half = Vec::with_capacity(enrichment_scores.len());

    for &current_es in enrichment_scores {
        let res_pair = if current_es >= 0.0 {
            es_ruler_pos.get_pval(current_es.abs(), eps, sign)
        } else {
            es_ruler_neg.get_pval(current_es.abs(), eps, sign)
        };

        pval_res.push(res_pair.0);
        is_cp_ge_half.push(res_pair.1);
    }

    GseaMultiLevelresults {
        pvals: pval_res,
        is_cp_ge_half,
    }
}
