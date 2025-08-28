use faer::{Mat, MatRef};
use rand::prelude::*;
use rand_distr::{Beta, Distribution, Gamma, Normal, Poisson};

////////////////
// Structures //
////////////////

/// Contains synthetic RNAseq data
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

/// Generate synthetic bulk RNAseq data with correlation structure
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

    // mean expression and dispersion parameter
    let gamma = Gamma::new(5.0, 10.0).unwrap();
    let mean_exp: Vec<f64> = (0..num_genes).map(|_| gamma.sample(&mut rng)).collect();
    let dispersion: Vec<f64> = mean_exp.iter().map(|&mean| 1.0 / (0.5 + mean)).collect();

    // generate the matrix
    let mut count_matrix: Mat<f64> = Mat::zeros(num_genes, num_samples);

    // assign the gene modules
    let mut gene_modules = vec![0; num_genes];

    if !add_modules {
        // easy version with no module assignment
        for i in 0..num_genes {
            let r = 1_f64 / dispersion[i];
            let p = r / (r + mean_exp[i]);

            // negative binomials were misbehaving, using this mathematical
            // equivalence trick with gamma, poisson
            let scale = (1.0 - p) / p;
            let gamma_dist = Gamma::new(r, scale).unwrap();

            for j in 0..num_samples {
                let lambda = gamma_dist.sample(&mut rng);
                let poisson_dist = Poisson::new(lambda).unwrap();
                count_matrix[(i, j)] = poisson_dist.sample(&mut rng);
            }
        }
    } else {
        let mut gene_idx = 0;
        for (module_id, &size) in module_sizes.iter().enumerate() {
            for _ in 0..size {
                if gene_idx < num_genes {
                    gene_modules[gene_idx] = module_id;
                    gene_idx += 1;
                }
            }
        }

        let num_modules = module_sizes.len();
        let mut module_factors: Mat<f64> = Mat::zeros(num_modules, num_samples);
        let normal = Normal::new(1.0, 0.7).unwrap();

        for i in 0..num_modules {
            for j in 0..num_samples {
                module_factors[(i, j)] = normal.sample(&mut rng);
            }
        }

        // add correlation structure
        let beta_dist = Beta::new(5.0, 2.0).unwrap();

        for i in 0..num_genes {
            if gene_modules[i] > 0 {
                let module_id = gene_modules[i] - 1;
                let correlation_strength = beta_dist.sample(&mut rng);

                let mut module_signal = vec![0.0; num_samples];
                for j in 0..num_samples {
                    module_signal[j] =
                        (correlation_strength * module_factors[(module_id, j)]).exp();
                }

                let signal_mean = module_signal.iter().sum::<f64>() / num_samples as f64;
                let scale_factor = mean_exp[i] / signal_mean;

                for j in 0..num_samples {
                    module_signal[j] *= scale_factor;

                    // generate counts with module signal using same gamma-poisson trick
                    let r = 1.0 / dispersion[i];
                    let p = r / (r + module_signal[j]);
                    let scale = (1.0 - p) / p;
                    let gamma_dist = Gamma::new(r, scale).unwrap();
                    let lambda = gamma_dist.sample(&mut rng);
                    let poisson_dist = Poisson::new(lambda).unwrap();
                    count_matrix[(i, j)] = poisson_dist.sample(&mut rng) as f64;
                }
            } else {
                // background genes
                let r = 1.0 / dispersion[i];
                let p = r / (r + mean_exp[i]);
                let scale = (1.0 - p) / p;
                let gamma_dist = Gamma::new(r, scale).unwrap();

                for j in 0..num_samples {
                    let lambda = gamma_dist.sample(&mut rng);
                    let poisson_dist = Poisson::new(lambda).unwrap();
                    count_matrix[(i, j)] = poisson_dist.sample(&mut rng);
                }
            }
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

/// Helper function to calculate the drop out probability
///
/// ### Params
#[inline(always)]
fn calculate_dropout_prob(exp: f64, midpoint: f64, shape: f64, global_sparsity: f64) -> f64 {
    let prob = 1_f64 / (1_f64 + (shape * (exp.ln_1p() - midpoint.ln_1p())).exp());
    let final_prob = prob + global_sparsity * (1_f64 - prob);

    final_prob.min(1_f64)
}

pub fn simulate_dropouts(
    original_counts: &MatRef<f64>,
    dropout_midpoint: f64,
    dropout_shape: f64,
    global_sparsity: f64,
    seed: u64,
) -> Mat<f64> {
    let mut rng = StdRng::seed_from_u64(seed);
    let (n_genes, n_samples) = original_counts.shape();

    let mut sparse_mat: Mat<f64> = Mat::zeros(n_genes, n_samples);

    for i in 0..n_genes {
        for j in 0..n_samples {
            let exp_val = original_counts.get(i, j);
            let dropout_prob =
                calculate_dropout_prob(*exp_val, dropout_midpoint, dropout_shape, global_sparsity);

            if rng.random::<f64>() > dropout_prob {
                sparse_mat[(i, j)] = *exp_val
            }
        }
    }

    sparse_mat
}
