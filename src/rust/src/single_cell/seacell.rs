#![allow(dead_code)]

use faer::MatRef;
use rand::prelude::*;
use rand::rngs::StdRng;
use rustc_hash::FxHashSet;

use crate::core::data::sparse_structures::*;

///////////
// Enums //
///////////

/// SNN similarity method
#[derive(Clone, Copy)]
pub enum SeaCellGraphGen {
    /// Only intersecting nearest neigbhbours will be considered
    Intersection,
    /// The union of nearest neighbours will be considered
    Union,
}

/// Helper function to parse the SEACell graph generation
///
/// ### Params
///
/// * `s` - Type of graph to build
///
/// ### Returns
///
/// Option of the SeaCellGraphGen
pub fn parse_seacell_graph(s: &str) -> Option<SeaCellGraphGen> {
    match s.to_lowercase().as_str() {
        "intersection" => Some(SeaCellGraphGen::Intersection),
        "union" => Some(SeaCellGraphGen::Union),
        _ => None,
    }
}

////////////
// Params //
////////////

pub struct SEACellsParams {
    pub n_sea_cells: usize,
    pub max_fw_iters: usize,
    pub convergence_epsilon: f32,
    pub max_iter: usize,
    pub min_iter: usize,
}

//////////
// Main //
//////////

pub struct SEACells<'a> {
    n_cells: usize,
    k: usize,
    kernel_mat: Option<CompressedSparseData<f32>>,
    k_squared: Option<CompressedSparseData<f32>>,
    a: Option<CompressedSparseData<f32>>,
    b: Option<CompressedSparseData<f32>>,
    archetypes: Option<Vec<usize>>,
    rss_history: Vec<f32>,
    convergence_threshold: Option<f32>,
    params: &'a SEACellsParams,
}

impl<'a> SEACells<'a> {
    /// Generate new instance
    pub fn new(n_cells: usize, n_seacells: usize, params: &'a SEACellsParams) -> Self {
        Self {
            n_cells,
            k: n_seacells,
            kernel_mat: None,
            k_squared: None,
            a: None,
            b: None,
            archetypes: None,
            convergence_threshold: None,
            rss_history: Vec::new(),
            params,
        }
    }

    /// Generate the kernel matrix
    pub fn construct_kernel_mat(
        &mut self,
        pca: MatRef<f32>,
        knn_indices: &[Vec<usize>],
        knn_distances: &[Vec<f32>],
        graph_construction: &SeaCellGraphGen,
        verbose: bool,
    ) {
        let n = pca.nrows();
        let k = knn_indices[0].len();

        if verbose {
            println!("Computing adaptive bandwidth RBF kernel...");
        }

        let median_idx = k / 2;
        let median_dist = knn_distances
            .iter()
            .map(|d| d[median_idx])
            .collect::<Vec<f32>>();

        // Use HashSet for sparse graph
        let mut edges = FxHashSet::default();
        for (i, neighbours) in knn_indices.iter().enumerate() {
            for &j in neighbours {
                edges.insert((i, j));
            }
        }

        match graph_construction {
            SeaCellGraphGen::Union => {
                let to_add: Vec<_> = edges
                    .iter()
                    .filter_map(|&(i, j)| (!edges.contains(&(j, i))).then_some((j, i)))
                    .collect();
                edges.extend(to_add);
            }
            SeaCellGraphGen::Intersection => {
                let to_keep: FxHashSet<_> = edges
                    .iter()
                    .copied()
                    .filter(|&(i, j)| edges.contains(&(j, i)))
                    .collect();
                edges = to_keep;
            }
        }

        // compute RBF...
        let mut rows: Vec<usize> = Vec::new();
        let mut cols: Vec<usize> = Vec::new();
        let mut vals: Vec<f32> = Vec::new();

        for &(i, j) in &edges {
            let mut dist_square = 0_f32;
            for dim in 0..pca.ncols() {
                let diff = pca.get(i, dim) - pca.get(j, dim);
                dist_square += diff * diff;
            }
            let sigma_prod = median_dist[i] * median_dist[j];
            let val = (-dist_square / sigma_prod).exp();
            if val > 1e-8 {
                rows.push(i);
                cols.push(j);
                vals.push(val);
            }
        }

        if verbose {
            println!("Built kernel with {} non-zeros", vals.len());
        }

        let kernel = coo_to_csr(&rows, &cols, &vals, (n, n));

        let k_squared = csr_matmul_csr(&kernel, &kernel.transform());

        self.kernel_mat = Some(kernel);
        self.k_squared = Some(k_squared);
    }

    /// Fit
    pub fn fit(&mut self, seed: usize, verbose: bool) {
        assert!(
            self.kernel_mat.is_some(),
            "Must construct kernel matrix first"
        );

        self.initialise_archetypes_greedy(verbose);
        self.initialise_matrices(verbose, seed as u64);

        let a = self.a.as_ref().unwrap();
        let b = self.b.as_ref().unwrap();

        let initial_rss = self.compute_rss(a, b);
        self.rss_history.push(initial_rss);
        self.convergence_threshold = Some(self.params.convergence_epsilon * initial_rss);

        if verbose {
            println!("Initial RSS: {:.6}", initial_rss);
            println!(
                "Convergence threshold: {:.6}",
                self.convergence_threshold.unwrap()
            );
        }

        let mut converged = false;
        let mut n_iter = 0;

        while (!converged && n_iter < self.params.max_iter) || n_iter < self.params.min_iter {
            n_iter += 1;

            if verbose && (n_iter == 1 || n_iter % 10 == 0) {
                println!("Starting iteration {}...", n_iter);
            }

            // Update A and B
            let b_current = self.b.clone().unwrap();
            let a_current = self.a.clone().unwrap();

            let a_new = self.update_a_mat(&b_current, &a_current);
            let b_new = self.update_b_mat(&a_new, &b_current);

            let rss = self.compute_rss(&a_new, &b_new);
            self.rss_history.push(rss);

            self.a = Some(a_new);
            self.b = Some(b_new);

            if verbose && (n_iter == 1 || n_iter % 10 == 0) {
                println!("  RSS: {:.6}", rss);
            }

            // Check convergence
            if n_iter > 1 {
                let rss_diff = (self.rss_history[n_iter - 1] - self.rss_history[n_iter]).abs();
                if rss_diff < self.convergence_threshold.unwrap() {
                    if verbose {
                        println!("Converged after {} iterations!", n_iter);
                    }
                    converged = true;
                }
            }
        }

        if !converged && verbose {
            println!(
                "Warning: Algorithm did not converge after {} iterations",
                self.params.max_iter
            );
        }
    }

    /// Initialise the archetypes in a greedy fashion
    fn initialise_archetypes_greedy(&mut self, verbose: bool) {
        let k_mat = self.k_squared.as_ref().unwrap();
        let n = k_mat.shape.0;
        let k = self.k;

        if verbose {
            println!("Initialising {} archetypes via greedy CSSP...", k);
        }

        // Convert to CSC for efficient column access
        let k_mat_csc = k_mat.transform();

        let mut f = vec![0_f32; n];
        let mut g = vec![0_f32; n];

        for i in 0..n {
            let row_start = k_mat.indptr[i];
            let row_end = k_mat.indptr[i + 1];

            for idx in row_start..row_end {
                let j = k_mat.indices[idx];
                let val = k_mat.data[idx];

                f[j] += val * val;

                if i == j {
                    g[i] = val;
                }
            }
        }

        let mut omega: Vec<Vec<f32>> = vec![vec![0_f32; n]; k];
        let mut centers: Vec<usize> = Vec::with_capacity(k);

        for iter in 0..k {
            let mut best_idx = 0;
            let mut best_score = f32::MIN;

            for i in 0..n {
                if g[i] > 1e-15 {
                    let score = f[i] / g[i];
                    if score > best_score {
                        best_score = score;
                        best_idx = i;
                    }
                }
            }

            centers.push(best_idx);

            // Extract column from CSC (much faster!)
            let mut k_col_p = vec![0_f32; n];
            let col_start = k_mat_csc.indptr[best_idx];
            let col_end = k_mat_csc.indptr[best_idx + 1];
            for idx in col_start..col_end {
                k_col_p[k_mat_csc.indices[idx]] = k_mat_csc.data[idx];
            }

            let mut delta = k_col_p.clone();
            for i in 0..n {
                let mut omega_sum = 0.0f32;
                for r in 0..iter {
                    omega_sum += omega[r][best_idx] * omega[r][i];
                }
                delta[i] -= omega_sum;
            }

            delta[best_idx] = delta[best_idx].max(0.0);

            let delta_p_sqrt = delta[best_idx].sqrt().max(1e-6);
            let mut omega_new = vec![0.0f32; n];
            for i in 0..n {
                omega_new[i] = delta[i] / delta_p_sqrt;
            }

            let omega_sq_norm: f32 = omega_new.iter().map(|&x| x * x).sum();

            let mut k_omega_new = vec![0.0f32; n];
            for i in 0..n {
                let row_start = k_mat.indptr[i];
                let row_end = k_mat.indptr[i + 1];
                for idx in row_start..row_end {
                    k_omega_new[i] += k_mat.data[idx] * omega_new[k_mat.indices[idx]];
                }
            }

            let dots: Vec<f32> = (0..iter)
                .map(|r| omega[r].iter().zip(&omega_new).map(|(a, b)| a * b).sum())
                .collect();

            for i in 0..n {
                let omega_hadamard = omega_new[i] * omega_new[i];
                let term1 = omega_sq_norm * omega_hadamard;

                let pl: f32 = (0..iter).map(|r| dots[r] * omega[r][i]).sum();
                let term2 = omega_new[i] * (k_omega_new[i] - pl);

                f[i] += -2.0 * term2 + term1;
                g[i] += omega_hadamard;
            }

            omega[iter] = omega_new;

            if verbose && (iter + 1) % 10 == 0 {
                println!("  Selected {} / {} archetypes", iter + 1, k);
            }
        }

        self.archetypes = Some(centers);
    }

    /// Initialise A and B matrices
    fn initialise_matrices(&mut self, verbose: bool, seed: u64) {
        let archetypes = self.archetypes.as_ref().unwrap();
        let k = archetypes.len();
        let n = self.n_cells;

        if verbose {
            println!("Initialising A and B matrices...");
        }

        let mut b_rows = Vec::new();
        let mut b_cols = Vec::new();
        let mut b_vals = Vec::new();

        for (col, &row) in archetypes.iter().enumerate() {
            b_rows.push(row);
            b_cols.push(col);
            b_vals.push(1_f32);
        }

        let b = coo_to_csr(&b_rows, &b_cols, &b_vals, (n, k));

        let archetypes_per_cell = (k as f32 * 0.25).ceil() as usize;
        let mut rng = StdRng::seed_from_u64(seed);

        let mut a_rows = Vec::new();
        let mut a_cols = Vec::new();
        let mut a_vals = Vec::new();

        for cell in 0..n {
            for _ in 0..archetypes_per_cell {
                let archetype = rng.random_range(0..k);
                a_rows.push(archetype);
                a_cols.push(cell);
                a_vals.push(rng.random::<f32>());
            }
        }

        let mut a = coo_to_csr(&a_rows, &a_cols, &a_vals, (k, n));

        normalise_csr_columns_l1(&mut a);

        a = self.update_a_mat(&b, &a);

        self.a = Some(a);
        self.b = Some(b);
    }

    /// Update the A matrix using Frank-Wolfe
    fn update_a_mat(
        &self,
        b: &CompressedSparseData<f32>,
        a_prev: &CompressedSparseData<f32>,
    ) -> CompressedSparseData<f32> {
        let k_mat = self.k_squared.as_ref().unwrap();

        let k_b = csr_matmul_csr(k_mat, b);
        let t2 = k_b.transform();
        let t1 = csr_matmul_csr(&t2, b);

        let mut a = a_prev.clone();
        let n = a.shape.1;
        let k = a.shape.0;

        for t in 0..self.params.max_fw_iters {
            let t1_a = csr_matmul_csr(&t1, &a);
            let g_mat = sparse_subtract_csr(&t1_a, &t2);
            let g_scaled = sparse_scalar_multiply_csr(&g_mat, 2.0);

            // Convert to CSC for O(sparsity) column access instead of O(k × sparsity)
            let g_csc = g_scaled.transform();

            let mut e_rows = Vec::with_capacity(n);
            let mut e_cols = Vec::with_capacity(n);
            let mut e_vals = Vec::with_capacity(n);

            for col in 0..n {
                let col_start = g_csc.indptr[col];
                let col_end = g_csc.indptr[col + 1];

                let mut min_val = 0.0f32; // Implicit zeros
                let mut min_idx = 0;

                // Check sparse entries
                for idx in col_start..col_end {
                    let row = g_csc.indices[idx];
                    let val = g_csc.data[idx];
                    if val < min_val {
                        min_val = val;
                        min_idx = row;
                    }
                }

                // If min is still 0.0 and we have sparse entries, find first missing row
                if min_val >= 0.0 && col_end > col_start {
                    let mut present = vec![false; k];
                    for idx in col_start..col_end {
                        present[g_csc.indices[idx]] = true;
                    }
                    if let Some(first_missing) = (0..k).find(|&r| !present[r]) {
                        min_idx = first_missing;
                    }
                }

                e_rows.push(min_idx);
                e_cols.push(col);
                e_vals.push(1.0f32);
            }

            let e = coo_to_csr(&e_rows, &e_cols, &e_vals, (k, n));

            let step_size = 2.0 / (t as f32 + 2.0);
            let e_minus_a = sparse_subtract_csr(&e, &a);
            let update = sparse_scalar_multiply_csr(&e_minus_a, step_size);
            a = sparse_add_csr(&a, &update);
        }

        a
    }

    /// Update the B matrix using Frank-Wolfe
    fn update_b_mat(
        &self,
        a: &CompressedSparseData<f32>,
        b_prev: &CompressedSparseData<f32>,
    ) -> CompressedSparseData<f32> {
        let k_mat = self.k_squared.as_ref().unwrap();

        let a_t = a.transform();
        let t1 = csr_matmul_csr(a, &a_t);
        let t2 = csr_matmul_csr(k_mat, &a_t);

        let mut b = b_prev.clone();
        let n = b.shape.0;
        let k = b.shape.1;

        for t in 0..self.params.max_fw_iters {
            let k_b = csr_matmul_csr(k_mat, &b);
            let k_b_t1 = csr_matmul_csr(&k_b, &t1);
            let g_mat = sparse_subtract_csr(&k_b_t1, &t2);
            let g_scaled = sparse_scalar_multiply_csr(&g_mat, 2.0);

            // Convert to CSC for efficient column access
            let g_csc = g_scaled.transform();

            let mut e_rows = Vec::with_capacity(k);
            let mut e_cols = Vec::with_capacity(k);
            let mut e_vals = Vec::with_capacity(k);

            for col in 0..k {
                let col_start = g_csc.indptr[col];
                let col_end = g_csc.indptr[col + 1];

                let mut min_val = 0.0f32; // Implicit zeros
                let mut min_idx = 0;

                // Check sparse entries
                for idx in col_start..col_end {
                    let row = g_csc.indices[idx];
                    let val = g_csc.data[idx];
                    if val < min_val {
                        min_val = val;
                        min_idx = row;
                    }
                }

                // If min is still 0.0 and we have sparse entries, find first missing row
                if min_val >= 0.0 && col_end > col_start {
                    let mut present = vec![false; n];
                    for idx in col_start..col_end {
                        present[g_csc.indices[idx]] = true;
                    }
                    if let Some(first_missing) = (0..n).find(|&r| !present[r]) {
                        min_idx = first_missing;
                    }
                }

                e_rows.push(min_idx);
                e_cols.push(col);
                e_vals.push(1.0f32);
            }

            let e = coo_to_csr(&e_rows, &e_cols, &e_vals, (n, k));

            let step_size = 2.0 / (t as f32 + 2.0);
            let e_minus_b = sparse_subtract_csr(&e, &b);
            let update = sparse_scalar_multiply_csr(&e_minus_b, step_size);
            b = sparse_add_csr(&b, &update);
        }

        b
    }

    /// Compute RSS
    ///
    /// Formula: ```||K - K @ B @ A||_F^2```
    fn compute_rss(&self, a: &CompressedSparseData<f32>, b: &CompressedSparseData<f32>) -> f32 {
        let k_mat = self.kernel_mat.as_ref().unwrap();

        // Reconstruction: K @ B @ A
        let k_b = csr_matmul_csr(k_mat, b);
        let reconstruction = csr_matmul_csr(&k_b, a);

        // Difference: K - reconstruction
        let diff = sparse_subtract_csr(k_mat, &reconstruction);

        // Frobenius norm
        frobenius_norm(&diff)
    }
}
