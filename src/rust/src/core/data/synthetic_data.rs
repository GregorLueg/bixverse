use extendr_api::*;
use faer::{Mat, MatRef};
use rand::prelude::*;
use rand_distr::weighted::WeightedAliasIndex;
use rand_distr::{Beta, Binomial, Distribution, Gamma, Normal, Poisson};
use rayon::prelude::*;

use crate::core::data::sparse_structures::{CompressedSparseData, CompressedSparseFormat};

///////////
// Enums //
///////////

/// Enum for the sparsity types
#[derive(Clone, Debug)]
pub enum SparsityFunction {
    /// Use the `Logistic` type function for sparisification
    Logistic,
    /// Use the `PowerDecay` type function for sparisification
    PowerDecay,
}

/////////////
// Helpers //
/////////////

/// Function to parse the sparsification type
///
/// ### Params
///
/// * `s` - String to parse
///
/// ### Returns
///
/// Option<SparsityFunction>
pub fn parse_sparsification(s: &str) -> Option<SparsityFunction> {
    match s.to_lowercase().as_str() {
        "log" => Some(SparsityFunction::Logistic),
        "powerdecay" => Some(SparsityFunction::PowerDecay),
        _ => None,
    }
}

////////////////
// Structures //
////////////////

/// Structure for synthetic RNAseq data
///
/// ### Fields
///
/// * `count_matrix` - The synthetic counts.
/// * `gene_modules` - A vector indicating to which gene module a gene belongs.
///   `0` indicating background.
#[derive(Clone, Debug)]
pub struct SyntheticRnaSeqData {
    pub count_matrix: Mat<f64>,
    pub gene_modules: Vec<usize>,
}

/////////////////////////
// Synthetic bulk data //
/////////////////////////

/// Generate synthetic bulk RNAseq data with optional correlation structure
///
/// ### Params
///
/// * `num_samples` - Number of samples in the data.
/// * `num_genes` - Number of genes in the data.
/// * `seed` - Seed for reproducibility purposes.
/// * `add_modules` - Shall correlation structure be added to the data.
/// * `module_sizes` - A vector indicating the size of the modules. Caution:
///   sum of `module_sizes` needs to be ≤ num_genes.
///
/// ### Returns
///
/// The `SyntheticRnaSeqData` data.
pub fn generate_bulk_rnaseq(
    num_samples: usize,
    num_genes: usize,
    seed: u64,
    add_modules: bool,
    module_sizes: Option<Vec<usize>>,
) -> SyntheticRnaSeqData {
    let module_sizes = module_sizes.unwrap_or(vec![300, 250, 200, 300, 500]);

    if add_modules {
        let total_module_genes: usize = module_sizes.iter().sum();
        assert!(
            total_module_genes <= num_genes,
            "Module sizes exceed total number of genes"
        );
    }

    let mut rng = StdRng::seed_from_u64(seed);

    // pre-compute gene level statistics
    let gamma = Gamma::new(5.0, 10.0).unwrap();
    let mean_exp: Vec<f64> = (0..num_genes).map(|_| gamma.sample(&mut rng)).collect();
    let dispersion: Vec<f64> = mean_exp
        .iter()
        .map(|&mean| 1.0 / (0.2 + mean * 0.3))
        .collect();

    // pre-compute base distributions for background genes
    let base_distributions: Vec<Gamma<f64>> = (0..num_genes)
        .map(|i| {
            let r = 1.0 / dispersion[i];
            let p = r / (r + mean_exp[i]);
            let scale = (1.0 - p) / p;
            Gamma::new(r, scale).unwrap()
        })
        .collect();

    let mut gene_modules = vec![0; num_genes];

    // pre-compute also the constant
    let inv_num_samples = 1.0 / num_samples as f64;

    let gene_data: Vec<Vec<f64>> = if !add_modules {
        // simple case without modules - parallelize gene processing
        (0..num_genes)
            .into_par_iter()
            .map(|i| {
                let mut local_rng = StdRng::seed_from_u64(seed.wrapping_add(i as u64));
                let mut gene_counts = Vec::with_capacity(num_samples);

                for _ in 0..num_samples {
                    let lambda = base_distributions[i].sample(&mut local_rng);
                    let poisson_dist = Poisson::new(lambda).unwrap();
                    gene_counts.push(poisson_dist.sample(&mut local_rng));
                }
                gene_counts
            })
            .collect()
    } else {
        // assign gene modules
        let mut gene_idx = 0;
        for (module_id, &size) in module_sizes.iter().enumerate() {
            for _ in 0..size {
                if gene_idx < num_genes {
                    gene_modules[gene_idx] = module_id + 1; // +1 since 0 indicates background
                    gene_idx += 1;
                }
            }
        }

        // generate module factors
        let num_modules = module_sizes.len();
        let mut module_factors: Mat<f64> = Mat::zeros(num_modules, num_samples);
        let normal = Normal::new(1.0, 0.7).unwrap();

        for i in 0..num_modules {
            for j in 0..num_samples {
                module_factors[(i, j)] = normal.sample(&mut rng);
            }
        }

        // pre-compute correlation strengths for all genes
        let beta_dist = Beta::new(5.0, 2.0).unwrap();
        let correlation_strengths: Vec<f64> = (0..num_genes)
            .map(|i| {
                if gene_modules[i] > 0 {
                    beta_dist.sample(&mut rng)
                } else {
                    0.0
                }
            })
            .collect();

        // extract module factors as vectors for easier parallel access
        let module_factor_vecs: Vec<Vec<f64>> = (0..num_modules)
            .map(|i| (0..num_samples).map(|j| module_factors[(i, j)]).collect())
            .collect();

        // parallelize gene processing
        (0..num_genes)
            .into_par_iter()
            .map(|i| {
                let mut local_rng = StdRng::seed_from_u64(seed.wrapping_add((i as u64) * 1000));
                let mut gene_counts = Vec::with_capacity(num_samples);

                if gene_modules[i] > 0 {
                    let module_id = gene_modules[i] - 1;
                    let correlation_strength = correlation_strengths[i];

                    let mut module_signal = Vec::with_capacity(num_samples);
                    for &factor in &module_factor_vecs[module_id] {
                        module_signal.push((correlation_strength * factor).exp());
                    }

                    let signal_mean = module_signal.iter().sum::<f64>() * inv_num_samples;
                    let scale_factor = mean_exp[i] / signal_mean;

                    #[allow(clippy::needless_range_loop)]
                    for j in 0..num_samples {
                        let scaled_signal = module_signal[j] * scale_factor;

                        let r = 1.0 / dispersion[i];
                        let p = r / (r + scaled_signal);
                        let scale = (1.0 - p) / p;
                        let gamma_dist = Gamma::new(r, scale).unwrap();
                        let lambda = gamma_dist.sample(&mut local_rng);
                        let poisson_dist = Poisson::new(lambda).unwrap();
                        gene_counts.push(poisson_dist.sample(&mut local_rng));
                    }
                } else {
                    // background genes
                    for _ in 0..num_samples {
                        let lambda = base_distributions[i].sample(&mut local_rng);
                        let poisson_dist = Poisson::new(lambda).unwrap();
                        gene_counts.push(poisson_dist.sample(&mut local_rng));
                    }
                }
                gene_counts
            })
            .collect()
    };
    // move results to matrix
    let mut count_matrix: Mat<f64> = Mat::zeros(num_genes, num_samples);

    for (i, gene_counts) in gene_data.into_iter().enumerate() {
        for (j, count) in gene_counts.into_iter().enumerate() {
            count_matrix[(i, j)] = count;
        }
    }

    SyntheticRnaSeqData {
        count_matrix,
        gene_modules,
    }
}

////////////////////
// Sparsification //
////////////////////

/// Uses a logistic function that plateaus at high expression, preserving some
/// variance through partial dropout (binomial thinning). Better for maintaining
/// mean-variance relationships in synthetic single-cell data.
///
/// **Dropout probability:**
///
/// `P(dropout) = clamp(1 / (1 + exp(shape * (ln(exp+1) - ln(midpoint+1)))), 0.3, 0.8) * (1 - global_sparsity) + global_sparsity`
///
/// **Characteristics:**
/// - Plateaus at ~30% dropout for high expression genes
/// - Partial dropout preserves count structure via binomial thinning
/// - Good for preserving variance-mean relationships
///
/// ### Params
///
/// * `original_counts` - The original count matrix
/// * `dropout_midpoint` - Expression level where dropout probability = 50%
/// * `dropout_shape` - Steepness of logistic curve (higher = steeper)
/// * `global_sparsity` - Additional dropout applied to all expression levels
/// * `seed` - Random seed for reproducibility
///
/// ### Returns
///
/// Count matrix with logistic dropouts and partial dropout applied
pub fn simulate_dropouts_logistic(
    original_counts: &MatRef<f64>,
    dropout_midpoint: f64,
    dropout_shape: f64,
    global_sparsity: f64,
    seed: u64,
) -> Mat<f64> {
    let (n_genes, n_samples) = original_counts.shape();

    // pre-calculate these ones
    let midpoint_ln1p = dropout_midpoint.ln_1p();
    let one_minus_global_sparsity = 1.0 - global_sparsity;

    // parallelise over the genes
    let gene_data: Vec<Vec<f64>> = (0..n_genes)
        .into_par_iter()
        .map(|i| {
            let mut local_rng = StdRng::seed_from_u64(seed.wrapping_add(i as u64));
            let mut gene_row = Vec::with_capacity(n_samples);

            let random_vals: Vec<f64> = (0..n_samples).map(|_| local_rng.random::<f64>()).collect();

            #[allow(clippy::needless_range_loop)]
            for j in 0..n_samples {
                let exp_val = *original_counts.get(i, j);

                // inlining this function makes it faster
                let exp_ln1p = exp_val.ln_1p();
                let logit = dropout_shape * (exp_ln1p - midpoint_ln1p);
                let prob = (1.0 / (1.0 + logit.exp())).clamp(0.3, 0.8);
                let final_prob = (prob * one_minus_global_sparsity + global_sparsity).min(1.0);

                gene_row.push(if random_vals[j] > final_prob {
                    exp_val
                } else if random_vals[j] > final_prob * 0.5 {
                    let retention_prob = local_rng.random_range(0.3..0.7);
                    let binomial = Binomial::new(exp_val as u64, retention_prob).unwrap();
                    binomial.sample(&mut local_rng) as f64
                } else {
                    0.0 // True zero
                });
            }
            gene_row
        })
        .collect();

    // add data to the matrix
    let mut sparse_mat = Mat::zeros(n_genes, n_samples);
    for (i, gene_row) in gene_data.into_iter().enumerate() {
        for (j, val) in gene_row.into_iter().enumerate() {
            sparse_mat[(i, j)] = val;
        }
    }
    sparse_mat
}

/// Simulate single cell dropouts using power decay function
///
/// Uses power decay that continues affecting high-expression genes without
/// plateau. More aggressive on highly expressed genes than logistic function.
/// Uses complete dropout only.
///
/// **Dropout probability:**
///
/// `P(dropout) = (midpoint / (exp + midpoint))^power * scale_factor * (1 - global_sparsity) + global_sparsity`
///
/// **Characteristics:**
/// - No plateau - high expression genes get substantial dropout
/// - Complete dropout only (no partial dropout)
/// - More uniform dropout across expression range
///
/// ### Params
///
/// * `original_counts` - The original count matrix
/// * `dropout_midpoint` - Reference expression level for decay calculation
/// * `power_factor` - Controls decay steepness (0.2-0.5 typical range)
/// * `global_sparsity` - Additional dropout applied to all expression levels
/// * `seed` - Random seed for reproducibility
///
/// ### Returns
///
/// Count matrix with power decay dropouts applied
pub fn simulate_dropouts_power_decay(
    original_counts: &MatRef<f64>,
    dropout_midpoint: f64,
    power_factor: f64,
    global_sparsity: f64,
    seed: u64,
) -> Mat<f64> {
    let (n_genes, n_samples) = original_counts.shape();

    let one_minus_global_sparsity = 1.0 - global_sparsity;
    let max_dropout = 0.9;
    // Cap maximum dropout at 90%
    // The function is very aggressive otherwise

    // parallelise over the genes
    let gene_data: Vec<Vec<f64>> = (0..n_genes)
        .into_par_iter()
        .map(|i| {
            let mut local_rng = StdRng::seed_from_u64(seed.wrapping_add(i as u64));
            let mut gene_row = Vec::with_capacity(n_samples);

            let random_vals: Vec<f64> = (0..n_samples).map(|_| local_rng.random::<f64>()).collect();

            #[allow(clippy::needless_range_loop)]
            for j in 0..n_samples {
                let exp_val = *original_counts.get(i, j);

                let base_ratio = dropout_midpoint / (exp_val + dropout_midpoint);
                let expr_dropout = base_ratio.powf(power_factor);
                let scaled_dropout = expr_dropout * 0.3; // scale down overall dropout
                let final_prob =
                    (scaled_dropout * one_minus_global_sparsity + global_sparsity).min(max_dropout);

                gene_row.push(if random_vals[j] < final_prob {
                    0.0
                } else {
                    exp_val
                });
            }
            gene_row
        })
        .collect();

    // add data to the matrix
    let mut sparse_mat = Mat::zeros(n_genes, n_samples);
    for (i, gene_row) in gene_data.into_iter().enumerate() {
        for (j, val) in gene_row.into_iter().enumerate() {
            sparse_mat[(i, j)] = val;
        }
    }
    sparse_mat
}

/////////////////////////
// Random sparse data ///
/////////////////////////

/// Create weighted sparse data resembling single cell counts in CSC
///
/// ### Params
///
/// * `nrow` - Number of rows (cells)
/// * `ncol` - Number of columns (genes)
/// * `n_cells` - Total no of cells
/// * `no_genes_exp` - Tuple representing the minimum number and the maximum
///   number of genes expressed per cell
/// * `max_exp` - Maximum expression a given gene can reach. Expression values
///   will be between `1..max_exp`
/// * `seed` - Seed for reproducibility purposes
///
/// ### Returns
///
/// The `CscData` type with the synthetic data.
pub fn create_sparse_csc_data(
    nrow: usize,
    ncol: usize,
    genes_per_cell: (usize, usize),
    max_exp: i32,
    seed: usize,
) -> CompressedSparseData<i32> {
    let weights: Vec<f64> = (1..=ncol).map(|i| 1.0 / i as f64).collect();
    let alias = WeightedAliasIndex::new(weights).unwrap();

    // Imitate what's going on the the CSR
    let mut gene_data: Vec<Vec<(usize, i32)>> = vec![Vec::new(); ncol];

    for cell_idx in 0..nrow {
        let mut rng = StdRng::seed_from_u64(seed as u64 + cell_idx as u64);
        let no_genes_expressed = rng.random_range(genes_per_cell.0..=genes_per_cell.1);

        let mut temp_vec = Vec::with_capacity(genes_per_cell.1);

        for _ in 0..no_genes_expressed {
            let gene_idx = alias.sample(&mut rng);
            let count = rng.random_range(1..=max_exp);
            temp_vec.push((gene_idx, count));
        }

        temp_vec.sort_unstable_by_key(|(gene_idx, _)| *gene_idx);

        for (gene_idx, count) in temp_vec {
            gene_data[gene_idx].push((cell_idx, count));
        }
    }

    // generate the CSC structure
    let estimated_total: usize = gene_data.iter().map(|v| v.len()).sum();
    let mut indptr = Vec::with_capacity(ncol + 1);
    let mut indices = Vec::with_capacity(estimated_total);
    let mut data = Vec::with_capacity(estimated_total);
    indptr.push(0);

    #[allow(clippy::needless_range_loop)]
    for gene_idx in 0..ncol {
        // Sort cells for this gene
        gene_data[gene_idx].sort_unstable_by_key(|(cell_idx, _)| *cell_idx);

        // Add ALL data for this gene (including duplicates from same cell)
        for (cell_idx, count) in &gene_data[gene_idx] {
            indices.push(*cell_idx);
            data.push(*count);
        }

        indptr.push(indices.len());
    }

    CompressedSparseData {
        data,
        indices,
        indptr,
        cs_type: CompressedSparseFormat::Csc,
        data_2: None::<Vec<i32>>,
        shape: (nrow, ncol),
    }
}

/// Create weighted sparse data resembling single cell counts in CSR format
///
/// ### Params
///
/// * `nrow` - Number of rows (cells)
/// * `ncol` - Number of columns (genes)
/// * `no_genes_exp` - Tuple representing the min and max number of genes expressed
///   per cell
/// * `max_exp` - Maximum expression a given gene can reach. Expression values will
///   be between `1..max_exp`
/// * `seed` - Seed for reproducibility purposes
///
/// ### Returns
///
/// The `CsrData` type with the synthetic data (cells as rows, genes as columns).
pub fn create_sparse_csr_data(
    nrow: usize,
    ncol: usize,
    no_genes_exp: (usize, usize),
    max_exp: i32,
    seed: usize,
) -> CompressedSparseData<i32> {
    let weights: Vec<f64> = (1..=ncol).map(|i| 1.0 / i as f64).collect();
    let alias = WeightedAliasIndex::new(weights).unwrap();

    let avg_genes = (no_genes_exp.0 + no_genes_exp.1) / 2;
    let estimated_total = ncol * avg_genes;

    let mut indptr = Vec::with_capacity(ncol + 1);
    let mut indices = Vec::with_capacity(estimated_total);
    let mut data = Vec::with_capacity(estimated_total);
    indptr.push(0);

    let mut temp_vec = Vec::with_capacity(no_genes_exp.1);

    for cell_idx in 0..nrow {
        let mut rng = StdRng::seed_from_u64(seed as u64 + cell_idx as u64);
        let no_genes_expressed = rng.random_range(no_genes_exp.0..=no_genes_exp.1);

        temp_vec.clear();

        for _ in 0..no_genes_expressed {
            let gene_idx = alias.sample(&mut rng);
            let count = rng.random_range(1..=max_exp);
            temp_vec.push((gene_idx, count));
        }

        // Sort by gene index
        temp_vec.sort_unstable_by_key(|(gene_idx, _)| *gene_idx);

        for (gene_idx, count) in temp_vec.iter() {
            indices.push(*gene_idx);
            data.push(*count);
        }

        indptr.push(indices.len());
    }

    CompressedSparseData {
        data,
        indices,
        indptr,
        cs_type: CompressedSparseFormat::Csr,
        data_2: None::<Vec<i32>>,
        shape: (nrow, ncol),
    }
}

///////////////////////////
// Specific sparse data ///
///////////////////////////

/// Structure to keep the CellTypeConfig
///
/// ### Fields
///
/// * `marker_genes` - Which indices are the marker genes for this specific
///   cell type
/// * `marker_exp_range` - Range of expression of the marker genes with
///   `(min, max)`.
/// * `markers_per_cell` - Number of markers per cell.
#[derive(Clone, Debug)]
pub struct CellTypeConfig {
    pub marker_genes: Vec<usize>,
    pub marker_exp_range: (u32, u32),
    pub markers_per_cell: (usize, usize),
}

impl CellTypeConfig {
    /// Generate the CellTypeConfig from an R list
    ///
    /// If values are not found, will use default values
    ///
    /// ### Params
    ///
    /// * `r_list` - The R list containing the parameters. Should have the
    ///   elements `"marker_genes"`, `"marker_exp_range"`, `"markers_per_cell"`.
    ///
    /// ### Returns
    ///
    /// The `CellTypeConfig` based on the R list.
    pub fn from_r_list(r_list: List) -> Self {
        let map = r_list.into_hashmap();

        let marker_genes = map
            .get("marker_genes")
            .and_then(|v| v.as_integer_vector())
            .map(|v| v.iter().map(|x| *x as usize).collect())
            .unwrap_or_default();

        let marker_exp_range = map
            .get("marker_exp_range")
            .and_then(|v| v.as_integer_vector())
            .and_then(|v| {
                if v.len() >= 2 {
                    Some((v[0] as u32, v[1] as u32))
                } else {
                    None
                }
            })
            .unwrap_or((10, 50));

        let markers_per_cell = map
            .get("markers_per_cell")
            .and_then(|v| v.as_integer_vector())
            .and_then(|v| {
                if v.len() >= 2 {
                    Some((v[0] as usize, v[1] as usize))
                } else {
                    None
                }
            })
            .unwrap_or((2, 4));

        CellTypeConfig {
            marker_genes,
            marker_exp_range,
            markers_per_cell,
        }
    }
}

/// Helper function to create synthetic data with specific cell types
///
/// ### Params
///
/// * `nrow` - Number of rows (cells).
/// * `ncol` - Number of columns (genes).
/// * `cell_type_configs` - A vector of cell type configurations.
/// * `background_genes_exp` - The number of background genes to express per
///   cell type.
/// * `background_exp_range` - The range of the expression in the background
///   genes.
/// * `seed` - Integer for reproducibility purposes
///
/// ### Returns
///
/// A tuple with `(csr data, indices of cell types)`
pub fn create_celltype_sparse_csr_data(
    nrow: usize,
    ncol: usize,
    cell_type_configs: Vec<CellTypeConfig>,
    background_genes_exp: (usize, usize),
    background_exp_range: (u32, u32),
    seed: usize,
) -> (CompressedSparseData<u32>, Vec<usize>) {
    let mut indptr = Vec::with_capacity(nrow + 1);
    let mut indices = Vec::with_capacity(nrow * 100);
    let mut data = Vec::with_capacity(nrow * 100);
    let mut cell_type_labels = Vec::with_capacity(nrow);
    indptr.push(0);
    let n_cell_types = cell_type_configs.len();
    let mut temp_vec = Vec::with_capacity(ncol);

    // Create weights for background genes - inverse distribution makes genes sparse
    let weights: Vec<f64> = (1..=ncol).map(|i| 1.0 / i as f64).collect();
    let alias = WeightedAliasIndex::new(weights).unwrap();

    for cell_idx in 0..nrow {
        let mut rng = StdRng::seed_from_u64(seed as u64 + cell_idx as u64);
        let cell_type = cell_idx % n_cell_types;
        cell_type_labels.push(cell_type);
        let config = &cell_type_configs[cell_type];
        temp_vec.clear();

        // sample marker genes for this cell type
        let n_markers = rng.random_range(config.markers_per_cell.0..=config.markers_per_cell.1);
        let marker_indices: Vec<usize> = config
            .marker_genes
            .choose_multiple(&mut rng, n_markers.min(config.marker_genes.len()))
            .copied()
            .collect();

        for &gene_idx in &marker_indices {
            let count = rng.random_range(config.marker_exp_range.0..=config.marker_exp_range.1);
            temp_vec.push((gene_idx, count));
        }

        // sample background genes with weighted distribution
        let n_background = rng.random_range(background_genes_exp.0..=background_genes_exp.1);
        for _ in 0..n_background {
            let gene_idx = alias.sample(&mut rng);
            if !marker_indices.contains(&gene_idx) {
                let count = rng.random_range(background_exp_range.0..=background_exp_range.1);
                temp_vec.push((gene_idx, count));
            }
        }

        // sort and deduplicate by summing counts
        temp_vec.sort_unstable_by_key(|(gene_idx, _)| *gene_idx);
        let mut i = 0;
        while i < temp_vec.len() {
            let gene_idx = temp_vec[i].0;
            let mut total_count = temp_vec[i].1;
            let mut j = i + 1;
            while j < temp_vec.len() && temp_vec[j].0 == gene_idx {
                total_count += temp_vec[j].1;
                j += 1;
            }
            indices.push(gene_idx);
            data.push(total_count);
            i = j;
        }
        indptr.push(indices.len());
    }

    let csr = CompressedSparseData {
        data,
        indices,
        indptr,
        cs_type: CompressedSparseFormat::Csr,
        data_2: None::<Vec<u32>>,
        shape: (nrow, ncol),
    };
    (csr, cell_type_labels)
}
