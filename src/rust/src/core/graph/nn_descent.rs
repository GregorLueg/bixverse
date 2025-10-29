use faer::MatRef;
use rand::rngs::SmallRng;
use rand::{Rng, SeedableRng};
use rayon::prelude::*;
use std::cell::RefCell;
use std::cmp::Reverse;
use std::collections::BinaryHeap;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::time::Instant;

use crate::core::graph::annoy::*;
use crate::core::graph::knn::{parse_ann_dist, AnnDist};

////////////////
// Neighbours //
////////////////

/// Packed neighbours structure for efficient memory layout
///
/// ### Fields
///
/// * `pid` - Neighbour's point ID
/// * `distance` - Distance to this neighbour
/// * `is_new` - 1 or 0 to identify if new. Using `u32` for better memory
///   alignment (avoids padding bytes in the struct)
#[repr(C)]
#[derive(Clone, Copy, Debug)]
struct Neighbour {
    pid: u32,
    dist: f32,
    is_new: u32,
}

impl Neighbour {
    /// Generate a new Neighbour
    ///
    /// ### Params
    ///
    /// * `pid` - Current point ID
    /// * `dist` - Distance
    /// * `is_new` - Boolean indicating if the neighbour is new.
    #[inline(always)]
    fn new(pid: usize, dist: f32, is_new: bool) -> Self {
        Self {
            pid: pid as u32,
            dist,
            is_new: is_new as u32,
        }
    }

    /// Getter for is True
    ///
    /// ### Returns
    ///
    /// Is this a new neighbour.
    #[inline(always)]
    fn is_new(&self) -> bool {
        self.is_new != 0
    }

    /// Getter for point ID
    ///
    /// ### Returns
    ///
    /// The point identifier.
    #[inline(always)]
    fn pid(&self) -> usize {
        self.pid as usize
    }
}

/// Wrapper for f32 that implements Ord for use in BinaryHeap
///
/// Faster than the sorts on full vectors and allows to keep data on heap
#[derive(Clone, Copy, Debug)]
struct OrderedFloat(f32);

impl PartialEq for OrderedFloat {
    fn eq(&self, other: &Self) -> bool {
        self.0 == other.0
    }
}

impl Eq for OrderedFloat {}

impl PartialOrd for OrderedFloat {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for OrderedFloat {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.0
            .partial_cmp(&other.0)
            .unwrap_or(std::cmp::Ordering::Equal)
    }
}

//////////////////////////
// Thread-local buffers //
//////////////////////////

thread_local! {
    /// HEAP_BUFFER -> Heap storage for `update_neighbours()`
    static HEAP_BUFFER: RefCell<BinaryHeap<Reverse<(OrderedFloat, usize, bool)>>> =
        const { RefCell::new(BinaryHeap::new()) };
    /// PID_SET -> SETS for the point IDs for `update_neighbours()`
    static PID_SET: RefCell<Vec<bool>> = const { RefCell::new(Vec::new()) };
    /// CANDIDATE_SET -> Candidate sets for `local_join()`
    static CANDIDATE_SET: RefCell<Vec<bool>> = const {RefCell::new(Vec::new())};
    /// SAMPLE_INDICES -> Sample indices for `local_join()`
    static SAMPLE_INDICES: RefCell<Vec<usize>> = const{RefCell::new(Vec::new())};
}

////////////////////
// Main algorithm //
////////////////////

/// NN-Descent graph builder
///
/// ### Fields
///
/// * `vectors_flat` - Flat structure of the initial data for better cache
///   locality
/// * `dim` - Initial feature dimensions
/// * `n` - Number of samples
/// * `metric` - Which distance metric to use. One of Euclidean or Cosine
/// * `norms` - Normalised values for each data point. Will be pre-computed
///   once if distance metric is Cosine.
///
/// ### Algorithm Overview
///
/// NN-Descent iteratively improves a k-NN graph through "local joins":
///
/// 1. For each vertex, examine neighbours of its neighbours
/// 2. Compute distances between these candidates
/// 3. Update the k-NN lists if improvements found
/// 4. Mark updated neighbours as "new" to guide next iteration
/// 5. Repeat until convergence
pub struct NNDescent {
    vectors_flat: Vec<f32>,
    dim: usize,
    n: usize,
    metric: AnnDist,
    norms: Vec<f32>,
}

impl NNDescent {
    /// Build the kNN graph with NN-Descent
    ///
    /// ### Params
    ///
    /// * `mat` - The embedding matrix with samples x features
    /// * `k` - Number of neighbours to look for
    /// * `max_iter` - Maximum number of iterations to run for
    /// * `delta` - Tolerance parameter. Should the proportion of changes in a
    ///   given iteration fall below that value, the algorithm stops.
    /// * `rho` - Sampling rate. Will be adaptively reduced in each iteration.
    /// * `seed` - Random seed for reproduction
    /// * `verbose` - Controls verbosity of the function
    ///
    /// ### Returns
    ///
    /// A nested vector of the updated nearest neighbours with their distances.
    #[allow(clippy::too_many_arguments)]
    pub fn build(
        mat: MatRef<f32>,
        k: usize,
        dist_metric: &str,
        max_iter: usize,
        delta: f32,
        rho: f32,
        seed: usize,
        verbose: bool,
    ) -> Vec<Vec<(usize, f32)>> {
        let metric = parse_ann_dist(dist_metric).unwrap_or(AnnDist::Cosine);
        let n = mat.nrows();
        let n_features = mat.ncols();

        // flatten matrix into contiguous memory for faster access
        let mut vectors_flat = Vec::with_capacity(n * n_features);
        for i in 0..n {
            vectors_flat.extend(mat.row(i).iter().copied());
        }

        // compute norms for cosine distance
        let norms = if metric == AnnDist::Cosine {
            (0..n)
                .map(|i| {
                    let start = i * n_features;
                    let end = start + n_features;
                    vectors_flat[start..end]
                        .iter()
                        .map(|x| x * x)
                        .sum::<f32>()
                        .sqrt()
                })
                .collect()
        } else {
            Vec::new()
        };

        let builder = NNDescent {
            vectors_flat,
            dim: n_features,
            n,
            metric,
            norms,
        };

        // initialise with Annoy for better-than-random starting point
        // 32 trees is what they use in the original paper/repo

        let start_initial_index = Instant::now();

        let annoy_index = AnnoyIndex::new(mat, 32, seed);

        let end_initial_index = start_initial_index.elapsed();

        if verbose {
            println!("Generated initial Annoy index: {:.2?}", end_initial_index);
        }

        builder.run(k, max_iter, &annoy_index, delta, rho, seed, verbose)
    }

    /// Run the underlying algorithm
    ///
    /// ### Params
    ///
    /// * `k` - Number of neighbours to search for
    /// * `max_iter` - Maximum number of iterations for the algorithm.
    /// * `annoy_index` - Annoy index for the initial fast initialisation of
    ///   the graph.
    /// * `delta` - Tolerance parameter. Should the proportion of changes in a
    ///   given iteration fall below that value, the algorithm stops.
    /// * `rho` - Sampling rate. Will be adaptively reduced in each iteration.
    /// * `seed` - Random seed for reproduction
    /// * `verbose` - Controls verbosity of the function
    ///
    /// ### Returns
    ///
    /// A nested vector of the updated nearest neighbours with their distances.
    ///
    /// ### Algorithm Details
    ///
    /// Each iteration consists of:
    /// 1. Candidate Generation: *For each vertex, find candidates via
    ///    local_join*
    /// 2. Bidirectional Distribution: *Candidates are both forward (i→j) and
    ///    reverse (j→i), ensuring the graph remains undirected*
    /// 3. Parallel Update: *Each vertex updates its k-NN list independently*
    /// 4. Convergence Check: *Stop if update rate falls below delta*
    #[allow(clippy::too_many_arguments)]
    fn run(
        &self,
        k: usize,
        max_iter: usize,
        annoy_index: &AnnoyIndex,
        delta: f32,
        rho: f32,
        seed: usize,
        verbose: bool,
    ) -> Vec<Vec<(usize, f32)>> {
        if verbose {
            println!("Initialising NN-Descent with {} cells - k={}", self.n, k);
        }

        let start_total = Instant::now();

        let mut graph = self.initialise_with_annoy(k, annoy_index);

        for iter in 0..max_iter {
            let updates = AtomicUsize::new(0);

            let current_rho = if iter == 0 {
                rho
            } else {
                (rho * 0.8_f32.powi(iter as i32 - 1)).max(0.3)
            };

            // phase 1: parallel candidate collection
            let all_candidates: Vec<(usize, Vec<(usize, f32)>)> = (0..self.n)
                .into_par_iter()
                .map(|i| {
                    let candidates = self.local_join(i, &graph, k, current_rho, seed + iter);
                    (i, candidates)
                })
                .collect();

            // phase 2: bidirectional candidate distribution
            let mut forward_candidates: Vec<Vec<(usize, f32)>> =
                vec![Vec::with_capacity(k * 5); self.n];
            let mut reverse_candidates: Vec<Vec<(usize, f32)>> =
                vec![Vec::with_capacity(k * 5); self.n];

            let chunk_size = all_candidates.len().div_ceil(rayon::current_num_threads());
            let reverse_chunks: Vec<Vec<Vec<(usize, f32)>>> = all_candidates
                .par_chunks(chunk_size)
                .map(|chunk| {
                    let mut local_reverse = vec![Vec::new(); self.n];
                    for (i, candidates) in chunk {
                        for &(j, dist) in candidates {
                            local_reverse[j].push((*i, dist));
                        }
                    }
                    local_reverse
                })
                .collect();

            for (i, candidates) in &all_candidates {
                forward_candidates[*i].extend(candidates.iter().copied());
            }

            for local_reverse in reverse_chunks {
                for (j, edges) in local_reverse.into_iter().enumerate() {
                    reverse_candidates[j].extend(edges);
                }
            }

            // phase 3: parallel graph update
            let new_graph: Vec<Vec<Neighbour>> = (0..self.n)
                .into_par_iter()
                .map(|i| {
                    let mut combined = Vec::with_capacity(
                        forward_candidates[i].len() + reverse_candidates[i].len(),
                    );
                    combined.extend_from_slice(&forward_candidates[i]);
                    combined.extend_from_slice(&reverse_candidates[i]);

                    self.update_neighbours(i, &graph[i], &combined, k, &updates)
                })
                .collect();

            graph = new_graph;

            // phase 4: convergence check
            let update_count = updates.load(Ordering::Relaxed);
            let update_rate = update_count as f32 / self.n as f32;

            if verbose {
                println!(
                    "Iteration {}: {} updates ({:.2}% of nodes), rho={:.3}",
                    iter + 1,
                    update_count,
                    update_rate * 100.0,
                    current_rho
                );
            }

            // early termination if convergence achieved
            if update_rate < delta {
                if verbose {
                    println!(
                        "Converged after {} iterations (update rate {:.4} < {:.4})",
                        iter + 1,
                        update_rate,
                        delta
                    );
                }
                break;
            }
        }

        let end_total = start_total.elapsed();

        if verbose {
            println!("Total run-time for NNDescent: {:.2?}", end_total);
        }

        graph
            .into_iter()
            .map(|neighbours| neighbours.into_iter().map(|n| (n.pid(), n.dist)).collect())
            .collect()
    }

    /// Initialise a first set of neighbours with annoy
    ///
    /// ### Params
    ///
    /// * `k` - Number of neighbours to sample
    /// * `annoy_index` - The Annoy index.
    ///
    /// ### Return
    ///
    /// A nested Vec of `Neighbour` structures.
    fn initialise_with_annoy(&self, k: usize, annoy_index: &AnnoyIndex) -> Vec<Vec<Neighbour>> {
        (0..self.n)
            .into_par_iter()
            .map(|i| {
                let query_vec = &self.vectors_flat[i * self.dim..(i + 1) * self.dim];

                // Use Annoy to get initial k neighbors
                // search_k controls quality - higher = better but slower
                let search_k = (k * 3).min(100);
                let (indices, distances) =
                    annoy_index.query(query_vec, &self.metric, k, Some(search_k));

                indices
                    .into_iter()
                    .zip(distances)
                    .filter(|(idx, _)| *idx != i) // exclude self
                    .map(|(idx, dist)| Neighbour::new(idx, dist, true))
                    .collect()
            })
            .collect()
    }

    /// Local join: find candidates from neighbours of neighbours
    ///
    ///
    /// ### Params
    ///
    /// * `node` - Node to check
    /// * `graph` - Current kNN graph
    /// * `k` - Number of neighbours
    /// * `rho` - Sampling rate for old neighbours
    /// * `seed` - Random seed for sampling
    ///
    /// ### Returns
    ///
    /// Vec of potential candidates as tuple `(index, dist)`.
    ///
    /// ### Algorithm Details
    ///
    /// 1. Separate new vs old neighbours
    /// 2. Sample old neighbours according to rho
    /// 3. Explore neighbours-of-neighbours (NN's NN)
    /// 4. Only compute distances if they might beat worst current distance
    /// 5. Return promising candidates
    fn local_join(
        &self,
        node: usize,
        graph: &[Vec<Neighbour>],
        k: usize,
        rho: f32,
        seed: usize,
    ) -> Vec<(usize, f32)> {
        let mut rng = SmallRng::seed_from_u64((seed as u64).wrapping_mul((node + 1) as u64));

        let worst_current_dist = graph[node].last().map(|n| n.dist).unwrap_or(f32::INFINITY);

        let mut new_neighbours = Vec::new();
        let mut old_neighbours = Vec::new();

        for n in &graph[node] {
            if n.is_new() {
                new_neighbours.push(n.pid());
            } else {
                old_neighbours.push(n.pid());
            }
        }

        let n_old_sample =
            ((old_neighbours.len() as f32 * rho).ceil() as usize).min(old_neighbours.len());
        if old_neighbours.len() > n_old_sample {
            for i in 0..n_old_sample {
                let j = rng.random_range(i..old_neighbours.len());
                old_neighbours.swap(i, j);
            }
            old_neighbours.truncate(n_old_sample);
        }

        let max_per_neighbour = k.min(10);
        let max_per_old = (k / 2).min(5);
        let max_candidates = k * 4;

        CANDIDATE_SET.with(|set_cell| {
            SAMPLE_INDICES.with(|indices_cell| {
                let mut candidate_set = set_cell.borrow_mut();
                let mut sample_indices = indices_cell.borrow_mut();

                if candidate_set.len() < self.n {
                    candidate_set.resize(self.n, false);
                }

                let mut candidate_ids = Vec::with_capacity(max_candidates);

                // Process new neighbours
                for &new_nb in &new_neighbours {
                    let neighbour_list = &graph[new_nb];
                    let take = max_per_neighbour.min(neighbour_list.len());

                    if neighbour_list.len() <= take {
                        for nn in neighbour_list {
                            let pid = nn.pid();
                            if pid != node && !candidate_set[pid] {
                                candidate_set[pid] = true;
                                candidate_ids.push(pid);
                                if candidate_ids.len() >= max_candidates {
                                    break;
                                }
                            }
                        }
                    } else {
                        sample_indices.clear();
                        sample_indices.extend(0..neighbour_list.len());

                        for i in 0..take {
                            let j = rng.random_range(i..neighbour_list.len());
                            sample_indices.swap(i, j);
                        }

                        for &idx in &sample_indices[..take] {
                            if candidate_ids.len() >= max_candidates {
                                break;
                            }
                            let pid = neighbour_list[idx].pid();
                            if pid != node && !candidate_set[pid] {
                                candidate_set[pid] = true;
                                candidate_ids.push(pid);
                            }
                        }
                    }
                }

                // Process old neighbours
                for &old_nb in &old_neighbours {
                    if candidate_ids.len() >= max_candidates {
                        break;
                    }
                    let neighbour_list = &graph[old_nb];
                    let take = max_per_old.min(neighbour_list.len());

                    if neighbour_list.len() <= take {
                        for nn in neighbour_list {
                            let pid = nn.pid();
                            if pid != node && !candidate_set[pid] {
                                candidate_set[pid] = true;
                                candidate_ids.push(pid);
                                if candidate_ids.len() >= max_candidates {
                                    break;
                                }
                            }
                        }
                    } else {
                        sample_indices.clear();
                        sample_indices.extend(0..neighbour_list.len());

                        for i in 0..take {
                            let j = rng.random_range(i..neighbour_list.len());
                            sample_indices.swap(i, j);
                        }

                        for &idx in &sample_indices[..take] {
                            if candidate_ids.len() >= max_candidates {
                                break;
                            }
                            let pid = neighbour_list[idx].pid();
                            if pid != node && !candidate_set[pid] {
                                candidate_set[pid] = true;
                                candidate_ids.push(pid);
                            }
                        }
                    }
                }

                let result = match self.metric {
                    AnnDist::Euclidean => candidate_ids
                        .iter()
                        .filter_map(|&c| {
                            let dist = unsafe { self.euclidean_distance(node, c) };
                            if dist < worst_current_dist {
                                Some((c, dist))
                            } else {
                                None
                            }
                        })
                        .collect(),
                    AnnDist::Cosine => candidate_ids
                        .iter()
                        .filter_map(|&c| {
                            let dist = unsafe { self.cosine_distance(node, c) };
                            if dist < worst_current_dist {
                                Some((c, dist))
                            } else {
                                None
                            }
                        })
                        .collect(),
                };

                // Clear candidate set
                for &pid in &candidate_ids {
                    candidate_set[pid] = false;
                }

                result
            })
        })
    }

    /// Update the neighbours with the improvements
    ///
    /// ### Params
    ///
    /// * `node` - Current node index
    /// * `current` - Current best neighbours
    /// * `candidates` - Potential new neighbours
    /// * `k` - Number of neighbours to find
    /// * `updates` - Borrewed AtomicUsize to check if an update happened
    ///
    /// ### Returns
    ///
    /// Vec of updates `Neigbour`s.
    fn update_neighbours(
        &self,
        node: usize,
        current: &[Neighbour],
        candidates: &[(usize, f32)],
        k: usize,
        updates: &AtomicUsize,
    ) -> Vec<Neighbour> {
        if candidates.is_empty() {
            return current
                .iter()
                .map(|n| Neighbour::new(n.pid(), n.dist, false))
                .collect();
        }

        HEAP_BUFFER.with(|heap_cell| {
            PID_SET.with(|set_cell| {
                let mut heap = heap_cell.borrow_mut();
                let mut pid_set = set_cell.borrow_mut();

                heap.clear();

                // ensure pid_set is large enough and clear it
                if pid_set.len() < self.n {
                    pid_set.resize(self.n, false);
                }

                // build max-heap and mark existing PIDs
                for n in current {
                    let pid = n.pid();
                    if pid != node {
                        heap.push(Reverse((OrderedFloat(n.dist), pid, false)));
                        pid_set[pid] = true;
                    }
                }

                // process candidates
                for &(pid, dist) in candidates {
                    if pid == node || pid_set[pid] {
                        continue;
                    }

                    if heap.len() < k {
                        heap.push(Reverse((OrderedFloat(dist), pid, true)));
                        pid_set[pid] = true;
                    } else if let Some(&Reverse((OrderedFloat(worst_dist), _, _))) = heap.peek() {
                        if dist < worst_dist {
                            if let Some(Reverse((_, old_pid, _))) = heap.pop() {
                                pid_set[old_pid] = false;
                            }
                            heap.push(Reverse((OrderedFloat(dist), pid, true)));
                            pid_set[pid] = true;
                        }
                    }
                }

                // convert heap to sorted vector
                let mut result: Vec<_> = heap
                    .drain()
                    .map(|Reverse((OrderedFloat(d), p, is_new))| {
                        pid_set[p] = false;
                        (d, p, is_new)
                    })
                    .collect();
                result.sort_unstable_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

                // check for changes
                let changed = result.len() != current.len()
                    || result
                        .iter()
                        .zip(current.iter())
                        .any(|(a, b)| a.1 != b.pid() || (a.0 - b.dist).abs() > 1e-6);

                if changed {
                    updates.fetch_add(1, Ordering::Relaxed);
                }

                result
                    .into_iter()
                    .map(|(dist, pid, is_new)| Neighbour::new(pid, dist, is_new && changed))
                    .collect()
            })
        })
    }

    /// Fast distance calculation with unsafe pointer arithmetic (Euclidean)
    ///
    /// ### Implementation note
    ///
    /// Manual loop unrolling (processing 4 elements at a time) helps the
    /// compiler generate better SIMD instructions. Modern CPUs can compute
    /// 4 float operations in parallel, so this can be ~4x faster than a
    /// naive loop.
    ///
    /// ### Params
    ///
    /// * `i` - Sample index i
    /// * `j` - Sample index j
    ///
    /// ### Returns
    ///
    /// The Euclidean distance between the two samples
    #[inline(always)]
    unsafe fn euclidean_distance(&self, i: usize, j: usize) -> f32 {
        let ptr_i = self.vectors_flat.as_ptr().add(i * self.dim);
        let ptr_j = self.vectors_flat.as_ptr().add(j * self.dim);

        let mut sum = 0_f32;
        let mut k = 0;

        // Process 4 elements at a time (helps with SIMD)
        while k + 4 <= self.dim {
            let d0 = *ptr_i.add(k) - *ptr_j.add(k);
            let d1 = *ptr_i.add(k + 1) - *ptr_j.add(k + 1);
            let d2 = *ptr_i.add(k + 2) - *ptr_j.add(k + 2);
            let d3 = *ptr_i.add(k + 3) - *ptr_j.add(k + 3);

            sum += d0 * d0 + d1 * d1 + d2 * d2 + d3 * d3;
            k += 4;
        }

        // Handle remaining elements
        while k < self.dim {
            let diff = *ptr_i.add(k) - *ptr_j.add(k);
            sum += diff * diff;
            k += 1;
        }

        sum
    }

    /// Fast distance calculation with unsafe pointer arithmetic (Cosine)
    ///
    /// ### Cosine Distance
    ///
    /// cosine_dist(u, v) = 1 - (u·v) / (||u|| ||v||)
    /// We pre-compute norms during initialisation to avoid repeated sqrt calls.
    ///
    /// ### Params
    ///
    /// * `i` - Sample index i
    /// * `j` - Sample index j
    ///
    /// ### Returns
    ///
    /// The Cosine distance between the two samples
    #[inline(always)]
    unsafe fn cosine_distance(&self, i: usize, j: usize) -> f32 {
        let ptr_i = self.vectors_flat.as_ptr().add(i * self.dim);
        let ptr_j = self.vectors_flat.as_ptr().add(j * self.dim);

        let mut dot = 0_f32;
        let mut k = 0;

        // Process 4 elements at a time (SIMD-friendly)
        while k + 4 <= self.dim {
            dot += *ptr_i.add(k) * *ptr_j.add(k)
                + *ptr_i.add(k + 1) * *ptr_j.add(k + 1)
                + *ptr_i.add(k + 2) * *ptr_j.add(k + 2)
                + *ptr_i.add(k + 3) * *ptr_j.add(k + 3);
            k += 4;
        }

        // Handle remaining elements
        while k < self.dim {
            dot += *ptr_i.add(k) * *ptr_j.add(k);
            k += 1;
        }

        1_f32 - (dot / (*self.norms.get_unchecked(i) * *self.norms.get_unchecked(j)))
    }
}
