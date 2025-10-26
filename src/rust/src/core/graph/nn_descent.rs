use faer::MatRef;
use rand::rngs::SmallRng;
use rand::{Rng, SeedableRng};
use rayon::prelude::*;
use rustc_hash::FxHashSet;
use std::cell::RefCell;
use std::sync::atomic::{AtomicUsize, Ordering};

use crate::core::graph::annoy::*;
use crate::core::graph::knn::{parse_ann_dist, AnnDist};

////////////////
// Neighbours //
////////////////

/// Packed neighbours
///
/// ### Fields
///
/// * `pid` - Id
/// * `distance` - Distance
/// * `is_new` - 1 or 0 to identify if new. Using `u32` for better memory
///   alignment
#[repr(C)]
#[derive(Clone, Copy, Debug)]
struct Neighbour {
    pid: u32,
    dist: f32,
    is_new: u32,
}

impl Neighbour {
    /// Generate a new neighbour
    ///
    /// ### Params
    ///
    /// * `pid` - Pair id
    /// * `dist` - Distance
    /// * `is_new` - Is this a new neighbour
    ///
    /// ### Returns
    ///
    /// Returns the initialised self
    #[inline(always)]
    fn new(pid: usize, dist: f32, is_new: bool) -> Self {
        Self {
            pid: pid as u32,
            dist,
            is_new: is_new as u32,
        }
    }

    /// Getter to see if the cell is new.
    ///
    /// ### Returns
    ///
    /// Returns true if the neighbour is new
    #[inline(always)]
    fn is_new(&self) -> bool {
        self.is_new != 0
    }

    /// Getter for the pair id.
    ///
    /// ### Returns
    ///
    /// Returns the pair id
    #[inline(always)]
    fn pid(&self) -> usize {
        self.pid as usize
    }
}

//////////////////////////
// Thread-local buffers //
//////////////////////////

thread_local! {
    static CANDIDATE_SET: RefCell<FxHashSet<usize>> = RefCell::new(FxHashSet::default());
    static MERGED_BUFFER: RefCell<Vec<(usize, f32, bool)>> = const { RefCell::new(Vec::new()) };
}

////////////////////
// Main algorithm //
////////////////////

/// NN-Descent graph builder
pub struct NNDescent {
    vectors_flat: Vec<f32>,
    dim: usize,
    n: usize,
    metric: AnnDist,
    norms: Vec<f32>,
}

impl NNDescent {
    /// Build the kNN graph with NNDescent
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

        let mut vectors_flat = Vec::with_capacity(n * n_features);
        for i in 0..n {
            vectors_flat.extend(mat.row(i).iter().copied());
        }

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

        let annoy_index = AnnoyIndex::new(mat, 20, seed);

        builder.run(k, max_iter, &annoy_index, delta, rho, seed, verbose)
    }

    /// Run the underlying algorithm
    ///
    /// ### Params
    ///
    /// * `k` - Number of neighbours to search for
    /// * `max_iter` - Maximum number of iterations for the algorithm.
    /// * `annoy_index` - The (small) initial Annoy index for better than
    ///   random initialisation.
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

        let mut graph = self.initialise_with_annoy(k, annoy_index);

        for iter in 0..max_iter {
            let updates = AtomicUsize::new(0);

            let current_rho = if iter == 0 {
                1.0
            } else {
                (rho * 0.8_f32.powi(iter as i32 - 1)).max(0.3)
            };

            // Parallel candidate collection
            let all_candidates: Vec<(usize, Vec<(usize, f32)>)> = (0..self.n)
                .into_par_iter()
                .map(|i| {
                    let candidates = self.local_join(i, &graph, k, current_rho, seed + iter);
                    (i, candidates)
                })
                .collect();

            // Pre-allocate with capacity
            let mut forward_candidates: Vec<Vec<(usize, f32)>> =
                vec![Vec::with_capacity(k * 5); self.n];
            let mut reverse_candidates: Vec<Vec<(usize, f32)>> =
                vec![Vec::with_capacity(k * 5); self.n];

            // Fast sequential distribution
            for (i, candidates) in all_candidates {
                for (j, dist) in candidates {
                    forward_candidates[i].push((j, dist));
                    reverse_candidates[j].push((i, dist));
                }
            }

            // Parallel graph update
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

            let update_count = updates.load(Ordering::Relaxed);
            let update_rate = update_count as f32 / self.n as f32;

            if verbose {
                println!(
                    "Iteration {}: {:.2}% nodes updated ({} / {})",
                    iter + 1,
                    update_rate * 100.0,
                    update_count,
                    self.n
                );
            }

            if update_rate < delta {
                if verbose {
                    println!("Converged after {} iterations", iter + 1);
                }
                break;
            }
        }

        graph
            .into_iter()
            .map(|neighbours| {
                neighbours
                    .into_iter()
                    .filter(|n| n.pid() < self.n)
                    .map(|n| (n.pid(), n.dist))
                    .collect()
            })
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
                let search_k = (k * 3).min(k * 10);
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
    /// ### Params
    ///
    /// * `node` - Node to check
    /// * `graph` - Current kNN graph
    /// * `k` - Number of neighbours
    /// * `rho` - Sampling rate
    /// * `seed` - Random seed for sampling
    ///
    /// ### Returns
    ///
    /// Returns a Vec of potential candidates as tuple `(index, dist)`.
    fn local_join(
        &self,
        node: usize,
        graph: &[Vec<Neighbour>],
        k: usize,
        rho: f32,
        seed: usize,
    ) -> Vec<(usize, f32)> {
        CANDIDATE_SET.with(|cell| {
            let mut candidates = cell.borrow_mut();
            candidates.clear();

            let mut rng = SmallRng::seed_from_u64((seed as u64).wrapping_mul((node + 1) as u64));

            let mut new_neighbours = Vec::new();
            let mut old_neighbours = Vec::new();

            for n in &graph[node] {
                if n.is_new() {
                    new_neighbours.push(n.pid());
                } else {
                    old_neighbours.push(n.pid());
                }
            }

            // Sample old neighbours
            let n_old_sample =
                ((old_neighbours.len() as f32 * rho).ceil() as usize).min(old_neighbours.len());
            if old_neighbours.len() > n_old_sample {
                for i in 0..n_old_sample {
                    let j = rng.random_range(i..old_neighbours.len());
                    old_neighbours.swap(i, j);
                }
                old_neighbours.truncate(n_old_sample);
            }

            // Reduced exploration limits for speed
            let max_per_neighbour = k.min(10);
            let max_per_old = (k / 2).min(5);

            // Explore new neighbours
            for &new_nb in &new_neighbours {
                let neighbour_list = &graph[new_nb];
                let take = max_per_neighbour.min(neighbour_list.len());

                if neighbour_list.len() <= take {
                    for nn in neighbour_list {
                        let pid = nn.pid();
                        if pid != node {
                            candidates.insert(pid);
                        }
                    }
                } else {
                    // Sample randomly
                    for _ in 0..take {
                        let idx = rng.random_range(0..neighbour_list.len());
                        let pid = neighbour_list[idx].pid();
                        if pid != node {
                            candidates.insert(pid);
                        }
                    }
                }
            }

            // Explore sampled old neighbours
            for &old_nb in &old_neighbours {
                let neighbour_list = &graph[old_nb];
                let take = max_per_old.min(neighbour_list.len());

                if neighbour_list.len() <= take {
                    for nn in neighbour_list {
                        let pid = nn.pid();
                        if pid != node {
                            candidates.insert(pid);
                        }
                    }
                } else {
                    for _ in 0..take {
                        let idx = rng.random_range(0..neighbour_list.len());
                        let pid = neighbour_list[idx].pid();
                        if pid != node {
                            candidates.insert(pid);
                        }
                    }
                }
            }

            // Limit total candidates to prevent excessive distance computations
            let max_candidates = k * 4;
            let candidate_vec: Vec<usize> = if candidates.len() > max_candidates {
                candidates.iter().copied().take(max_candidates).collect()
            } else {
                candidates.iter().copied().collect()
            };

            match self.metric {
                AnnDist::Euclidean => candidate_vec
                    .into_iter()
                    .map(|c| (c, unsafe { self.euclidean_distance(node, c) }))
                    .collect(),
                AnnDist::Cosine => candidate_vec
                    .into_iter()
                    .map(|c| (c, unsafe { self.cosine_distance(node, c) }))
                    .collect(),
            }
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
    /// * `updates` - Borrowed AtomicUsize to track if an update happened
    ///
    /// ### Returns
    ///
    /// Vec of updated `Neighbour`s.
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

        MERGED_BUFFER.with(|cell| {
            let mut merged = cell.borrow_mut();
            merged.clear();

            for n in current {
                if n.pid() != node {
                    merged.push((n.pid(), n.dist, false));
                }
            }

            for &(pid, dist) in candidates {
                if pid != node {
                    merged.push((pid, dist, true));
                }
            }

            merged.sort_unstable_by(|a, b| {
                a.1.partial_cmp(&b.1).unwrap().then_with(|| a.0.cmp(&b.0))
            });
            merged.dedup_by_key(|x| x.0);
            merged.truncate(k);

            let changed = merged.len() != current.len()
                || merged
                    .iter()
                    .zip(current.iter())
                    .any(|(a, b)| a.0 != b.pid() || (a.1 - b.dist).abs() > 1e-6);

            if changed {
                updates.fetch_add(1, Ordering::Relaxed);
            }

            merged
                .iter()
                .map(|&(pid, dist, is_new)| Neighbour::new(pid, dist, is_new && changed))
                .collect()
        })
    }

    /// Fast distance calculation with unsafe pointer arithmetic (Euclidean)
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

        while k + 4 <= self.dim {
            let d0 = *ptr_i.add(k) - *ptr_j.add(k);
            let d1 = *ptr_i.add(k + 1) - *ptr_j.add(k + 1);
            let d2 = *ptr_i.add(k + 2) - *ptr_j.add(k + 2);
            let d3 = *ptr_i.add(k + 3) - *ptr_j.add(k + 3);

            sum += d0 * d0 + d1 * d1 + d2 * d2 + d3 * d3;
            k += 4;
        }

        while k < self.dim {
            let diff = *ptr_i.add(k) - *ptr_j.add(k);
            sum += diff * diff;
            k += 1;
        }

        sum
    }

    /// Fast distance calculation with unsafe pointer arithmetic (Cosine)
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

        while k + 4 <= self.dim {
            dot += *ptr_i.add(k) * *ptr_j.add(k)
                + *ptr_i.add(k + 1) * *ptr_j.add(k + 1)
                + *ptr_i.add(k + 2) * *ptr_j.add(k + 2)
                + *ptr_i.add(k + 3) * *ptr_j.add(k + 3);
            k += 4;
        }

        while k < self.dim {
            dot += *ptr_i.add(k) * *ptr_j.add(k);
            k += 1;
        }

        1_f32 - (dot / (*self.norms.get_unchecked(i) * *self.norms.get_unchecked(j)))
    }
}
