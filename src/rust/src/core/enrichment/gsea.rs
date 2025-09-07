use extendr_api::List;

use rand::distr::Uniform;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use rand_distr::Distribution;
use rayon::prelude::*;
use rustc_hash::{FxBuildHasher, FxHashMap, FxHashSet};
use statrs::distribution::{Beta, ContinuousCDF};
use statrs::function::gamma::digamma;

use crate::core::base::stats::trigamma;
use crate::utils::general::{array_max, array_min, cumsum, unique};

//////////////////
// Type aliases //
//////////////////

/// Type alias for multi level error calculations.
///
/// ### Fields
///
/// * `0` - Simple error for the multi level fgsea
/// * `1` - Multi error for the multi level fgsea
pub type MultiLevelErrRes = (Vec<f64>, Vec<f64>);

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

///////////////////
// Sample Chunks //
///////////////////

/// Structure for fgsea multi-level
/// Stores and manages chunked samples for efficient processing
///
/// ### Fields
///
/// * `chunk_sum` - Sum of rank values in each chunk for fast computation
/// * `chunks` - Vector of chunks, each containing gene indices for that chunk
#[derive(Clone, Debug)]
struct SampleChunks {
    chunk_sum: Vec<f64>,
    chunks: Vec<Vec<i32>>,
}

impl SampleChunks {
    /// Creates new SampleChunks with specified number of chunks
    ///
    /// ### Params
    ///
    /// * `chunks_number` - Number of chunks to create
    ///
    /// ### Returns
    ///
    /// Initialised structure
    fn new(chunks_number: usize) -> Self {
        Self {
            chunk_sum: vec![0.0; chunks_number],
            chunks: vec![Vec::new(); chunks_number],
        }
    }
}

//////////////
// ES Ruler //
//////////////

/// This is EsRuler implementation for adaptive enrichment score sampling
///
/// Further structure for fgsea multi-level
///
/// ### Fields
///
/// * `ranks` - Gene ranks used for ES calculation
/// * `sample_size` - Current sample size (may change during duplication)
/// * `original_sample_size` - Original sample size for p-value calculation
/// * `pathway_size` - Number of genes in the pathway
/// * `current_samples` - Current sample sets being processed
/// * `enrichment_scores` - Calculated enrichment scores from samples
/// * `prob_corrector` - Probability correction factors for p-value adjustment
/// * `chunks_number` - Number of chunks used for optimization
/// * `chunk_last_element` - Last element index in each chunk
#[derive(Clone, Debug)]
struct EsRuler {
    ranks: Vec<f64>,
    sample_size: usize,
    original_sample_size: usize,
    pathway_size: usize,
    current_samples: Vec<Vec<usize>>,
    enrichment_scores: Vec<f64>,
    prob_corrector: Vec<usize>,
    chunks_number: i32,
    chunk_last_element: Vec<i32>,
}

impl EsRuler {
    /// Creates a new ES ruler for adaptive sampling
    ///
    /// ### Params
    ///
    /// * `inp_ranks` - Input gene ranks
    /// * `inp_sample_size` - Sample size
    /// * `inp_pathway_size` - Pathway size
    ///
    /// ### Returns
    ///
    /// Initialised structure
    fn new(inp_ranks: &[f64], inp_sample_size: usize, inp_pathway_size: usize) -> Self {
        let mut current_samples: Vec<Vec<usize>> = Vec::with_capacity(inp_sample_size);
        current_samples.resize_with(inp_sample_size, Vec::new);

        Self {
            ranks: inp_ranks.to_vec(),
            sample_size: inp_sample_size,
            original_sample_size: inp_sample_size,
            pathway_size: inp_pathway_size,
            current_samples,
            enrichment_scores: Vec::new(),
            prob_corrector: Vec::new(),
            chunks_number: 0,
            chunk_last_element: Vec::new(),
        }
    }

    /// Removes samples with low ES and duplicates samples with high ES
    ///
    /// This drives the sampling process toward higher and higher ES values
    fn duplicate_samples(&mut self) {
        let mut stats: Vec<(f64, usize)> = vec![(0.0, 0); self.sample_size];
        let mut pos_es_indxs: FxHashSet<usize> = FxHashSet::default();
        let mut total_pos_es_count: i32 = 0;

        for (sample_id, stat) in stats.iter_mut().enumerate().take(self.sample_size) {
            let sample_es_pos = calc_positive_es(&self.ranks, &self.current_samples[sample_id]);
            let sample_es = calc_es(&self.ranks, &self.current_samples[sample_id]);
            if sample_es > 0.0 {
                total_pos_es_count += 1;
                pos_es_indxs.insert(sample_id);
            }
            *stat = (sample_es_pos, sample_id);
        }

        stats.sort_unstable_by(|a, b| {
            a.0.partial_cmp(&b.0)
                .unwrap_or(std::cmp::Ordering::Equal)
                .then_with(|| a.1.cmp(&b.1))
        });

        // Memory pre-allocation for speed
        let half_size = self.sample_size / 2;
        self.enrichment_scores.reserve(half_size);
        self.prob_corrector.reserve(half_size);

        let mut sample_id = 0;
        while 2 * sample_id < self.sample_size {
            self.enrichment_scores.push(stats[sample_id].0);
            if pos_es_indxs.contains(&stats[sample_id].1) {
                total_pos_es_count -= 1;
            }
            self.prob_corrector.push(total_pos_es_count as usize);
            sample_id += 1;
        }

        let mut new_sets = Vec::with_capacity(self.sample_size);
        let mut sample_id = 0;
        while 2 * sample_id < self.sample_size - 2 {
            for _ in 0..2 {
                new_sets
                    .push(self.current_samples[stats[self.sample_size - 1 - sample_id].1].clone());
            }
            sample_id += 1;
        }
        new_sets.push(self.current_samples[stats[self.sample_size >> 1].1].clone());

        std::mem::swap(&mut self.current_samples, &mut new_sets);

        // We need to update this here, to avoid potential out of index errors.
        self.sample_size = self.current_samples.len();
    }

    /// Attempts to improve a sample by swapping genes in/out using perturbation
    ///
    /// ### Params
    /// * `ranks` - Gene ranks
    /// * `k` - Number of genes in sample
    /// * `sample_chunks` - Sample chunks to modify
    /// * `bound` - ES boundary threshold
    /// * `rng` - Random number generator
    ///
    /// ### Returns
    ///
    /// Number of successful perturbations
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
        let mut ns: f64 = sample_chunks.chunk_sum.iter().sum();

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
                    && sample_chunks.chunks[old_chunk_ind].len() <= tmp as usize
                {
                    tmp -= sample_chunks.chunks[old_chunk_ind].len() as i32;
                    old_chunk_ind += 1;
                }
                old_ind_in_chunk = tmp;
                old_val = sample_chunks.chunks[old_chunk_ind][old_ind_in_chunk as usize];
            }

            // Select random new value
            let new_val = uid_n.sample(rng);

            // Find insertion position
            let new_chunk_ind = match self.chunk_last_element.binary_search(&(new_val)) {
                Ok(idx) => idx + 1,
                Err(idx) => idx,
            };

            let new_ind_in_chunk =
                match sample_chunks.chunks[new_chunk_ind].binary_search(&(new_val)) {
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
            sample_chunks.chunks[old_chunk_ind].remove(old_ind_in_chunk as usize);
            let adjust =
                if old_chunk_ind == new_chunk_ind && old_ind_in_chunk < new_ind_in_chunk as i32 {
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
                sample_chunks.chunks[old_chunk_ind].insert(old_ind_in_chunk as usize, old_val);

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
    ///
    /// Uses an adaptive sampling approach to explore higher ES values
    ///
    /// ### Params
    ///
    /// * `es` - Target enrichment score
    /// * `seed` - Random seed
    /// * `eps` - Precision parameter (0.0 for no precision requirement)
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
        let mut samples_chunks =
            vec![SampleChunks::new(self.chunks_number as usize); self.sample_size];

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

    /// Calculate the p-value for a given enrichment score
    ///
    /// ### Params
    ///
    /// * `es` - Enrichment score
    /// * `sign` - Whether to consider sign in calculation
    ///
    /// # Returns
    ///
    /// Tuple of (p-value, error quality flag)
    fn get_pval(&self, es: f64, sign: bool) -> (f64, bool) {
        let half_size = (self.original_sample_size + 1) / 2; // Use original sample_size here!
        let it_index;
        let mut good_error = true;

        if es >= *self.enrichment_scores.last().unwrap_or(&0.0) {
            it_index = self.enrichment_scores.len() - 1;
            if es > self.enrichment_scores[it_index] + 1e-10 {
                good_error = false;
            }
        } else {
            it_index = match self.enrichment_scores.binary_search_by(|probe| {
                probe.partial_cmp(&es).unwrap_or(std::cmp::Ordering::Equal)
            }) {
                Ok(index) => index,
                Err(index) => index,
            };
        }

        let indx = if it_index > 0 { it_index } else { 0 };
        let k = indx / half_size;
        let remainder = self.original_sample_size - (indx % half_size); // Use original sample_size here!

        let adj_log = beta_mean_log(half_size, self.original_sample_size); // Use original sample_size here!
        let adj_log_pval =
            k as f64 * adj_log + beta_mean_log(remainder + 1, self.original_sample_size); // Use original sample_size here!

        if sign {
            (0.0_f64.max(1.0_f64.min(adj_log_pval.exp())), good_error)
        } else {
            let correction =
                calc_log_correction(&self.prob_corrector, indx, self.original_sample_size); // Use original sample_size here!
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
///
/// ### Params
///
/// * `ranks` - Gene ranks array
/// * `pathway_indices` - Indices of genes in the pathway
///
/// ### Returns
///
/// Enrichment score value
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
///
/// ### Params
///
/// * `ranks` - Gene ranks array
/// * `pathway_indices` - Indices of genes in the pathway
///
/// # Returns
///
/// Positive enrichment score value
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

/// Generate k random numbers from [a, b] inclusive range using Fisher-Yates shuffle
///
/// ### Params
///
/// * `a` - Range start
/// * `b` - Range end
/// * `k` - Number of elements to select
/// * `rng` - Random number generator
///
/// ### Returns
///
/// Sorted vector of k random elements
///
/// ### Panics
///
/// If k > range size (b - a + 1)
fn combination(a: usize, b: usize, k: usize, rng: &mut impl Rng) -> Vec<usize> {
    let n = b - a + 1;
    if k > n {
        panic!("k cannot be greater than range size n");
    }

    // Create a vector with all possible values
    let mut indices: Vec<usize> = (a..=b).collect();

    // Use Fisher-Yates shuffle to get k random elements
    // This is equivalent to partial_shuffle in C++
    for i in 0..k {
        let j = rng.random_range(i..n);
        indices.swap(i, j);
    }

    // Take first k elements and sort them
    let mut result: Vec<usize> = indices.into_iter().take(k).collect();
    result.sort_unstable();

    result
}

/// Rearranges array so nth element is in its sorted position
///
/// ### Params
///
/// * `arr` - Array to rearrange
/// * `n` - Target position
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

/// Calculate the Enrichment score assuming stats is sorted and pathway contains index positions
///
/// ### Params
///
/// * `stats` - Sorted gene statistics
/// * `pathway` - Index positions of genes in the pathway
///
/// ### Returns
///
/// Enrichment score
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

/// Calculate once for each size the permutation-based enrichment scores
///
/// ### Params
///
/// * `stats` - Gene statistics
/// * `gene_set_sizes` - Unique gene set sizes
/// * `shared_perms` - Shared permutations
///
/// ### Returns
///
/// HashMap mapping sizes to permutation scores
fn create_perm_es(
    stats: &[f64],
    gene_set_sizes: &[usize],
    shared_perms: &[Vec<usize>],
) -> FxHashMap<usize, Vec<f64>> {
    let mut shared_perm_es =
        FxHashMap::with_capacity_and_hasher(gene_set_sizes.len(), FxBuildHasher);
    for size in gene_set_sizes {
        let perm_es: Vec<f64> = shared_perms
            .into_par_iter()
            .map(|perm| calculate_es(stats, &perm[..*size]))
            .collect();
        shared_perm_es.insert(*size, perm_es);
    }
    shared_perm_es
}

/// Calculate the permutations in the 'traditional' way
///
/// ### Params
///
/// * `stats` - Gene statistics
/// * `pathway_scores` - Pathway enrichment scores
/// * `pathway_sizes` - Pathway sizes
/// * `iters` - Number of iterations
/// * `seed` - Random seed
///
/// ### Returns
///
/// Batch results from traditional, permutation-based method
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

/// Calculate random scores batch-wise using the fgsea simple method
///
/// ### Params
///
/// * `stats` - Gene statistics
/// * `pathway_scores` - Pathway enrichment scores
/// * `pathway_sizes` - Pathway sizes
/// * `iters` - Number of iterations
/// * `gsea_param` - GSEA parameter
/// * `seed` - Random seed
///
/// ### Returns
///
/// Batch results from fgsea simple method
pub fn calc_gsea_stat_cumulative_batch(
    stats: &[f64],
    pathway_scores: &[f64],
    pathway_sizes: &[usize],
    iters: usize,
    gsea_param: f64,
    seed: u64,
) -> GseaBatchResults {
    let n = stats.len();
    let k = array_max(pathway_sizes);

    let shared_perm = create_perm_es_simple(stats, gsea_param, iters, k, n, seed, true);

    calc_gsea_stats_wrapper(pathway_scores, pathway_sizes, &shared_perm)
}

//////////////////////
// FGSEA multilevel //
//////////////////////

/// Helper to calculate the beta mean log
///
/// ### Params
///
/// * `a` - First beta parameter
/// * `b` - Second beta parameter
///
/// ### Returns
///
/// Beta mean log value
fn beta_mean_log(a: usize, b: usize) -> f64 {
    digamma(a as f64) - digamma((b + 1) as f64)
}

/// Calculates the log corrections for p-value adjustment
///
/// ### Params
///
/// * `prob_corrector` - Probability correction vector
/// * `prob_corr_idx` - Correction index
/// * `sample_size` - Sample size
///
/// ### Returns
///
/// Tuple of (log correction, validity flag)
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

/// Calculates multilevel error for a given p-value and sample size
///
/// ### Params
///
/// * `pval` - P-value
/// * `sample_size` - Sample size
///
/// ### Returns
///
/// Multilevel error estimate
fn multilevel_error(pval: &f64, sample_size: &f64) -> f64 {
    let floor_term = (-pval.log2() + 1.0).floor();
    let trigamma_diff = trigamma((sample_size + 1.0) / 2.0) - trigamma(sample_size + 1.0);

    (floor_term * trigamma_diff).sqrt() / f64::ln(2.0)
}

/// Function to do the multi-level magic in fgsea
///
/// ### Params
///
/// * `enrichment_score` - Target enrichment score
/// * `ranks` - Gene ranks
/// * `pathway_size` - Size of pathway
/// * `sample_size` - Sample size
/// * `seed` - Random seed
/// * `eps` - Precision parameter (0.0 for no precision requirement)
/// * `sign` - Whether to consider sign in calculation
///
/// ### Returns
///
/// Tuple of (p-value, error quality flag)
pub fn fgsea_multilevel_helper(
    enrichment_score: f64,
    ranks: &[f64],
    pathway_size: usize,
    sample_size: usize,
    seed: u64,
    eps: f64,
    sign: bool,
) -> (f64, bool) {
    let ranks: Vec<f64> = if enrichment_score >= 0.0 {
        ranks.iter().map(|&r| r.abs()).collect()
    } else {
        // Only create what we need
        let mut neg_ranks: Vec<f64> = ranks.iter().map(|&r| r.abs()).collect();
        neg_ranks.reverse();
        neg_ranks
    };

    let mut es_ruler = EsRuler::new(&ranks, sample_size, pathway_size);
    es_ruler.extend(enrichment_score.abs(), seed, eps);
    es_ruler.get_pval(enrichment_score.abs(), sign)
}

/// Calculates the simple and multi error estimates
///
/// ### Params
///
/// * `n_more_extreme` - Number of more extreme permutations
/// * `nperm` - Total permutations
/// * `sample_size` - Sample size
///
/// ### Returns
///
/// Tuple of (simple errors, multi errors)
pub fn calc_simple_and_multi_error(
    n_more_extreme: &[usize],
    nperm: usize,
    sample_size: usize,
) -> MultiLevelErrRes {
    let no_tests = n_more_extreme.len();
    let n_more_extreme_f64: Vec<f64> = n_more_extreme.iter().map(|x| *x as f64).collect();
    let nperm_f64 = nperm as f64;
    let sample_size_f64 = sample_size as f64;

    let mut left_border = Vec::with_capacity(no_tests);
    let mut right_border = Vec::with_capacity(no_tests);
    let mut crude_est = Vec::with_capacity(no_tests);
    for n in &n_more_extreme_f64 {
        // Left border
        if n > &0.0 {
            let beta = Beta::new(*n, nperm_f64 - n + 1.0).unwrap();
            left_border.push(beta.inverse_cdf(0.025).log2());
        } else {
            // Need to deal with special case of n = 0; this then becomes neg infinity
            left_border.push(f64::NEG_INFINITY);
        }
        // Right broder
        let beta = Beta::new(n + 1.0, nperm_f64 - n).unwrap();
        right_border.push(beta.inverse_cdf(1.0 - 0.025).log2());
        // Crude
        crude_est.push(((n + 1.0) / (nperm_f64 + 1.0)).log2());
    }

    let mut simple_err = Vec::with_capacity(no_tests);
    for i in 0..no_tests {
        simple_err.push(
            0.5 * f64::max(
                crude_est[i] - left_border[i],
                right_border[i] - crude_est[i],
            ),
        );
    }

    let multi_err: Vec<f64> = n_more_extreme_f64
        .iter()
        .map(|n| {
            let pval = (*n + 1.0) / (nperm_f64 + 1.0);
            multilevel_error(&pval, &sample_size_f64)
        })
        .collect();

    (simple_err, multi_err)
}
