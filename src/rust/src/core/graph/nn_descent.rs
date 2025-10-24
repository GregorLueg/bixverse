#![allow(dead_code)]

use faer::MatRef;
use rand::rngs::SmallRng;
use rand::{Rng, SeedableRng};
use rayon::prelude::*;
use rustc_hash::FxHashSet;
use std::sync::atomic::{AtomicUsize, Ordering};

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
/// * `is_new` - 1 or 0 to identify if new. Using using `u32` for better memory
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
    ///
    #[inline(always)]
    fn is_new(&self) -> bool {
        self.is_new != 0
    }

    /// Getter to see if the cell is new.
    ///
    /// ### Returns
    ///
    /// Returns the pair id
    #[inline(always)]
    fn pid(&self) -> usize {
        self.pid as usize
    }
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

        // flatten for better cache locality
        let mut vectors_flat = Vec::with_capacity(n * n_features);
        for i in 0..n {
            vectors_flat.extend(mat.row(i).iter().copied());
        }

        let builder = NNDescent {
            vectors_flat,
            dim: n_features,
            n,
            metric,
        };

        builder.run(k, max_iter, delta, rho, seed, verbose)
    }

    /// Run the underlying algorithm
    ///
    /// ### Params
    ///
    /// * `k` - Number of neighbours to search for
    /// * `max_iter` - Maximum number of iterations for the algorithm.
    /// * `delta` - Tolerance parameter. Should the proportion of changes in a
    ///   given iteration fall below that value, the algorithm stops.
    /// * `rho` - Sampling rate. Will be adaptively reduced in each iteration.
    /// * `seed` - Random seed for reproduction
    /// * `verbose` - Controls verbosity of the function
    ///
    /// ### Returns
    ///
    /// A nested vector of the updated nearest neighbours with their distances.
    fn run(
        &self,
        k: usize,
        max_iter: usize,
        delta: f32,
        rho: f32,
        seed: usize,
        verbose: bool,
    ) -> Vec<Vec<(usize, f32)>> {
        if verbose {
            println!("Initialising NN-Descent with {} cells - k={}", self.n, k);
        }

        // initialise graph
        let mut graph = self.initialise_random(k, seed);

        for iter in 0..max_iter {
            let updates = AtomicUsize::new(0);

            // adaptive sampling; reduce as iterations progress
            let current_rho = if iter < 2 {
                rho
            } else {
                (rho * 0.5_f32.powi((iter - 1) as i32)).max(0.2)
            };

            let candidates: Vec<Vec<(usize, f32)>> = (0..self.n)
                .into_par_iter()
                .map(|i| self.local_join(i, &graph, k, current_rho, seed + iter))
                .collect();

            let new_graph: Vec<Vec<Neighbour>> = (0..self.n)
                .into_par_iter()
                .map(|i| {
                    let improved =
                        self.update_neighbours(i, &graph[i], &candidates[i], k, &updates);
                    improved
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

        // Convert to output format
        graph
            .into_iter()
            .map(|neighbors| {
                neighbors
                    .into_iter()
                    .filter(|n| n.pid() < self.n)
                    .map(|n| (n.pid(), n.dist))
                    .collect()
            })
            .collect()
    }

    /// Initialise a random set of neigbours
    ///
    /// ### Params
    ///
    /// * `k` - Number of neighbours to sample
    /// * `seed` - Seed for reproducibility
    ///
    /// ### Return
    ///
    /// A nested Vec of `Neighbour` structures.
    fn initialise_random(&self, k: usize, seed: usize) -> Vec<Vec<Neighbour>> {
        (0..self.n)
            .into_par_iter()
            .map(|i| {
                let mut rng = SmallRng::seed_from_u64((seed as u64).wrapping_mul((i + 1) as u64));
                let mut neighbours: Vec<Neighbour> = Vec::with_capacity(k);

                let mut sampled = FxHashSet::default();
                sampled.insert(i);

                while neighbours.len() < k && sampled.len() < self.n {
                    let j = rng.random_range(0..self.n);
                    if sampled.insert(j) {
                        let dist = unsafe { self.distance(i, j) };
                        neighbours.push(Neighbour::new(j, dist, true));
                    }
                }

                neighbours.sort_unstable_by(|a, b| a.dist.partial_cmp(&b.dist).unwrap());
                neighbours
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
    /// * `seed`
    ///
    /// ### Params
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
        let mut candidates = FxHashSet::default();
        let mut rng = SmallRng::seed_from_u64((seed as u64).wrapping_mul((node + 1) as u64));

        let new_neighbours: Vec<usize> = graph[node]
            .iter()
            .filter(|n| n.is_new())
            .map(|n| n.pid())
            .collect();

        if new_neighbours.is_empty() {
            return Vec::new();
        }

        for &neighbour in &new_neighbours {
            for nn in &graph[neighbour] {
                if nn.is_new() {
                    let pid = nn.pid();
                    if pid != node {
                        candidates.insert(pid);
                    }
                }
            }
        }

        // sample based on rho
        let sample_size = ((candidates.len() as f32 * rho).ceil() as usize).min(k * 5);
        let mut candidate_vec: Vec<usize> = candidates.into_iter().collect();

        if candidate_vec.len() > sample_size {
            // Fisher-Yates shuffle...
            for i in 0..sample_size {
                let j = rng.random_range(i..candidate_vec.len());
                candidate_vec.swap(i, j);
            }
            candidate_vec.truncate(sample_size);
        }

        candidate_vec
            .into_iter()
            .map(|c| {
                let dist = unsafe { self.distance(node, c) };
                (c, dist)
            })
            .collect()
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
            // Mark all as old
            return current
                .iter()
                .map(|n| Neighbour::new(n.pid(), n.dist, false))
                .collect();
        }

        let mut merged = Vec::with_capacity(current.len() + candidates.len());

        // add current ones
        for n in current {
            if n.pid() != node {
                merged.push((n.pid(), n.dist, false));
            }
        }

        // add new ones
        for &(pid, dist) in candidates {
            if pid != node {
                merged.push((pid, dist, true));
            }
        }

        merged.sort_unstable_by(|a, b| a.1.partial_cmp(&b.1).unwrap().then_with(|| a.0.cmp(&b.0)));
        merged.dedup_by_key(|x| x.0);
        merged.truncate(k);

        // check changes
        let changed = merged.len() != current.len()
            || merged
                .iter()
                .zip(current.iter())
                .any(|(a, b)| a.0 != b.pid() || (a.1 - b.dist).abs() > 1e-6);

        if changed {
            updates.fetch_add(1, Ordering::Relaxed);
        }

        // mark recent additions as new, keep old ones as old
        merged
            .into_iter()
            .map(|(pid, dist, is_new)| Neighbour::new(pid, dist, is_new))
            .collect()
    }

    /// Fast distance calculation with unsafe pointer arithmetic
    ///
    /// ### Params
    ///
    /// * `i` - Sample index i
    /// * `j` - Sample index j
    ///
    /// ### Returns
    ///
    /// The distance between the two samples
    #[inline(always)]
    unsafe fn distance(&self, i: usize, j: usize) -> f32 {
        let ptr_i = self.vectors_flat.as_ptr().add(i * self.dim);
        let ptr_j = self.vectors_flat.as_ptr().add(j * self.dim);

        match self.metric {
            AnnDist::Euclidean => {
                let mut sum = 0_f32;
                for k in 0..self.dim {
                    let diff = *ptr_i.add(k) - *ptr_j.add(k);
                    sum += diff * diff;
                }
                // no square rooting needed
                sum
            }
            AnnDist::Cosine => {
                let mut dot = 0_f32;
                let mut norm_a = 0_f32;
                let mut norm_b = 0_f32;

                for k in 0..self.dim {
                    let a = *ptr_i.add(k);
                    let b = *ptr_j.add(k);
                    dot += a * b;
                    norm_a += a * a;
                    norm_b += b * b;
                }

                1_f32 - (dot / (norm_a.sqrt() * norm_b.sqrt()))
            }
        }
    }
}
