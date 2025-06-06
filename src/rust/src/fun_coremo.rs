use extendr_api::prelude::*;

use faer::Mat;
use rand::prelude::*;
use rayon::prelude::*;
use std::collections::HashMap;

use crate::utils_r_rust::{faer_to_r_matrix, r_matrix_to_faer};
use crate::utils_stats::*;

/// Calculates the TOM over an affinity matrix
///
/// @description Calculates the topological overlap measure for a given affinity matrix
/// x. Has the option to calculate the signed and unsigned version.
///
/// @param x Numerical matrix. Affinity matrix.
/// @param tom_type String. One of `c("v1", "v2")` - pending on choice, a different
/// normalisation method will be used.
/// @param signed Boolean. Shall the signed TOM be calculated. If set to `FALSE`, values
/// should be ≥ 0.
///
/// @return Returns the TOM matrix.
///
/// @export
#[extendr]
fn rs_tom(
    x: RMatrix<f64>,
    tom_type: &str,
    signed: bool,
) -> extendr_api::Result<extendr_api::RArray<f64, [usize; 2]>> {
    let x = r_matrix_to_faer(&x);
    let tom_version = parse_tom_types(tom_type).ok_or_else(|| {
        extendr_api::Error::Other(format!("Invalid TOM version type: {}", tom_type))
    })?;

    let tom_mat: Mat<f64> = match tom_version {
        TomType::Version1 => calc_tom(x, signed),
        TomType::Version2 => calc_tom_v2(x, signed),
    };

    let res = faer_to_r_matrix(tom_mat.as_ref());

    Ok(res)
}

/// Helper function to assess CoReMo cluster quality
///
/// @description This function assesses the quality of the clusters
/// with a given cut `k`. Returns the median R2 (cor^2) and the median absolute
/// deviation (MAD) of the clusters. Large clusters (≥1000) are subsampled
/// to a random set of 1000 genes.
///
/// @param cluster_genes A list. Contains the cluster and their respective genes.
/// @param cor_mat Numerical matrix. Contains the correlation coefficients.
/// @param row_names String vector. The row names (or column names) of the
/// correlation matrix.
/// @param seed Integer. Random seed for the sub sampling of genes.
///
/// @return A list containing:
///  \itemize{
///   \item r2med - median R2 of the cluster.
///   \item r2mad - median absolute deviation of the R2 in the cluster.
///   \item size - size of the cluster.
/// }
#[extendr]
fn rs_coremo_quality(
    cluster_genes: List,
    cor_mat: RMatrix<f64>,
    row_names: Vec<String>,
    seed: u64,
) -> extendr_api::Result<List> {
    // R matrix to faer
    let cor_mat = r_matrix_to_faer(&cor_mat);

    // Faster look-ups
    let gene_map: HashMap<&str, usize> = row_names
        .iter()
        .enumerate()
        .map(|(i, gene)| (gene.as_str(), i))
        .collect();

    // Get the cluster indices
    let mut all_indices = Vec::with_capacity(cluster_genes.len());

    for i in 0..cluster_genes.len() {
        let element_i = cluster_genes.elt(i)?;
        if let Some(internal_vals) = element_i.as_string_vector() {
            let mut indices = Vec::with_capacity(internal_vals.len());
            for gene in &internal_vals {
                if let Some(&idx) = gene_map.get(gene.as_str()) {
                    indices.push(idx);
                }
            }
            all_indices.push(indices)
        }
    }

    let all_sizes: Vec<usize> = all_indices.iter().map(|vec| vec.len()).collect();

    // Iterate through the clusters and get the data
    let res: Vec<(f64, f64)> = all_indices
        .par_iter()
        .map(|index_vec| {
            let n = index_vec.len();

            let indices: Vec<usize> = if n > 1000 {
                let mut rng = StdRng::seed_from_u64(seed);
                index_vec.choose_multiple(&mut rng, 1000).cloned().collect()
            } else {
                index_vec.to_vec()
            };

            let indices_diagonal: Vec<(usize, &[usize])> = indices
                .iter()
                .enumerate()
                .map(|(i, first)| (*first, &index_vec[i + 1..]))
                .take_while(|(_, rest)| !rest.is_empty())
                .collect();

            let expected_size = index_vec.len() * (index_vec.len() - 1) / 2;
            let mut vals: Vec<f64> = Vec::with_capacity(expected_size);

            for (r_idx, col_idx) in indices_diagonal {
                for col in col_idx {
                    vals.push(cor_mat[(r_idx, *col)].powi(2))
                }
            }

            let r2_med = median(&vals);
            let r2_mad = mad(&vals);

            (r2_med, r2_mad)
        })
        .collect();

    let mut r2_med_vec: Vec<f64> = Vec::with_capacity(res.len());
    let mut r2_mad_vec: Vec<f64> = Vec::with_capacity(res.len());

    for (med, mad) in res {
        r2_med_vec.push(med);
        r2_mad_vec.push(mad);
    }

    Ok(list!(
        r2med = r2_med_vec,
        r2mad = r2_mad_vec,
        size = all_sizes
    ))
}

extendr_module! {
    mod fun_coremo;
    fn rs_tom;
    fn rs_coremo_quality;
}
