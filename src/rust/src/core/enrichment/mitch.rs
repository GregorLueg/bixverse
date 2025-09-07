use extendr_api::*;
use faer::{Mat, MatRef};
use rayon::prelude::*;
use rustc_hash::FxHashMap;
use statrs::statistics::Statistics;

use crate::core::base::stats::*;
use crate::core::base::utils::*;

//////////////////
// Type aliases //
//////////////////

/// Type alias for processed MitchPathways
///
/// ### Slots
///
/// * `0` - Name of the pathway
/// * `1` - Index positions of the pathway relative to the input matrix
type MitchPathways = (Vec<String>, Vec<Vec<usize>>);

////////////////////
// Result structs //
////////////////////

#[derive(Clone, Debug)]
pub struct MitchResult<'a> {
    pub pathway_name: &'a str,
    pub pathway_size: usize,
    pub manova_pval: f64,
    pub scores: Vec<f64>,
    pub anova_pvals: Vec<f64>,
    pub s_dist: f64,
    pub mysd: f64,
}

///////////////
// Functions //
///////////////

pub fn manova_mitch(x: MatRef<f64>, group1_indices: &[usize]) -> ManovaResult {
    let (n, p) = x.shape();

    assert!(
        group1_indices.iter().all(|&idx| idx < n),
        "All indices must be less than the number of rows"
    );

    let mut is_group1 = vec![false; n];
    for &idx in group1_indices {
        is_group1[idx] = true;
    }

    // Collect group0 indices
    let group0_indices: Vec<usize> = (0..n).filter(|&i| !is_group1[i]).collect();

    let n1 = group1_indices.len();
    let n0 = group0_indices.len();

    let mut x0 = Mat::zeros(n0, p);
    let mut x1 = Mat::zeros(n1, p);

    for (new_i, &orig_i) in group0_indices.iter().enumerate() {
        for j in 0..p {
            x0[(new_i, j)] = x[(orig_i, j)];
        }
    }

    for (new_i, &orig_i) in group1_indices.iter().enumerate() {
        for j in 0..p {
            x1[(new_i, j)] = x[(orig_i, j)];
        }
    }

    let mean_0 = col_means(x0.as_ref());
    let mean_1 = col_means(x1.as_ref());
    let mean_overall = col_means(x);

    // scale the matrices
    let x_centered = scale_matrix_col(&x.as_ref(), false);
    let x0_centered = scale_matrix_col(&x0.as_ref(), false);
    let x1_centered = scale_matrix_col(&x1.as_ref(), false);

    let sscp_within =
        x0_centered.transpose() * &x0_centered + x1_centered.transpose() * &x1_centered;
    let sscp_total = x_centered.transpose() * &x_centered;
    let sscp_between = &sscp_total - &sscp_within;

    ManovaResult {
        sscp_between,
        sscp_within,
        sscp_total,
        df_between: 1,
        df_within: n - 2,
        df_total: n - 1,
        n_vars: p,
        group_means: vec![mean_0, mean_1],
        overall_mean: mean_overall,
    }
}

pub fn prepare_mitch_pathways(
    row_names: &[String],
    pathway_list: List,
    min_size: usize,
) -> Result<MitchPathways> {
    let gene_map: FxHashMap<&str, usize> = row_names
        .iter()
        .enumerate()
        .map(|(i, gene)| (gene.as_str(), i))
        .collect();

    let list_names: Vec<String> = pathway_list
        .names()
        .unwrap()
        .map(|s| s.to_string())
        .collect();

    let mut filtered_pathways = Vec::new();
    let mut filtered_names = Vec::new();

    #[allow(clippy::needless_range_loop)]
    for i in 0..pathway_list.len() {
        let element = pathway_list.elt(i)?;
        if let Some(internal_vals) = element.as_string_vector() {
            let mut indices = Vec::with_capacity(internal_vals.len());

            for gene in &internal_vals {
                if let Some(&idx) = gene_map.get(gene.as_str()) {
                    indices.push(idx);
                }
            }

            if indices.len() >= min_size {
                indices.sort_unstable();
                filtered_pathways.push(indices);
                filtered_names.push(list_names[i].clone());
            }
        }
    }

    Ok((filtered_names, filtered_pathways))
}

pub fn mitch_rank(mat: &MatRef<f64>) -> Mat<f64> {
    let mut ranked_mat = Mat::zeros(mat.nrows(), mat.ncols());

    // parallel ranking directly into the matrix
    ranked_mat
        .par_col_iter_mut()
        .enumerate()
        .for_each(|(col_idx, mut col)| {
            let original_col: Vec<f64> = mat.col(col_idx).iter().copied().collect();
            let mut zeros = 0_f64;
            let mut neg = 0_f64;
            for x in &original_col {
                if *x == 0_f64 {
                    zeros += 1_f64;
                } else if *x < 0_f64 {
                    neg += 1_f64;
                }
            }
            let adj = neg + (zeros / 2.0);

            let ranks = rank_vector(&original_col);
            let ranks = ranks.iter().map(|x| *x - adj).collect::<Vec<f64>>();

            for (row_idx, &rank) in ranks.iter().enumerate() {
                col[row_idx] = rank;
            }
        });

    ranked_mat
}

pub fn process_mitch_pathway<'a>(
    ranked_mat: MatRef<f64>,
    pathway_name: &'a str,
    pathway_indices: &[usize],
) -> MitchResult<'a> {
    let nrow = ranked_mat.nrows() as f64;
    let manova_res: ManovaResult = manova_mitch(ranked_mat, pathway_indices);
    let sum_manova: ManovaSummary = ManovaSummary::from_manova_res(&manova_res);
    let sum_anova = summary_aov(&manova_res);

    let p_manova = sum_manova.p_val_pillai;
    let mut p_aovs = Vec::with_capacity(sum_anova.len());

    for sum in sum_anova {
        p_aovs.push(sum.p_val);
    }

    let scores = &manova_res.group_means[0]
        .iter()
        .zip(&manova_res.group_means[1])
        .map(|(mean_1, mean_0)| (2.0 * (mean_1 - mean_0)) / nrow)
        .collect::<Vec<f64>>();

    let sd_scores = scores.std_dev();
    let hypotenuse = scores.iter().map(|x| x.powi(2)).sum::<f64>().sqrt();

    MitchResult {
        pathway_name,
        pathway_size: pathway_indices.len(),
        manova_pval: p_manova,
        scores: scores.clone(),
        anova_pvals: p_aovs,
        s_dist: hypotenuse,
        mysd: sd_scores,
    }
}
