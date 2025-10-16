use faer::ColRef;

/// Calculate the mutual information between two column references
///
/// The columns need to be binned
///
/// ### Params
///
/// * `col_i` - Column reference to the first column to compare
/// * `col_j` - Column reference to the second column to compare
/// * `n_bins` - Optional number of bins. If not provided, will default to
///   `sqrt(nrows)`.
///
/// ### Returns
///
/// The mutual information between the two columns
pub fn calculate_mi(col_i: ColRef<usize>, col_j: ColRef<usize>, n_bins: Option<usize>) -> f64 {
    let n_rows = col_i.nrows();
    let n_bins = n_bins.unwrap_or_else(|| (n_rows as f64).sqrt() as usize);

    let mut joint_counts = vec![vec![0usize; n_bins]; n_bins];
    let mut marginal_i = vec![0usize; n_bins];
    let mut marginal_j = vec![0usize; n_bins];

    for i in 0..n_rows {
        let bin_i = col_i[i];
        let bin_j = col_j[i];

        joint_counts[bin_i][bin_j] += 1;
        marginal_i[bin_i] += 1;
        marginal_j[bin_j] += 1;
    }

    let n = n_rows as f64;
    let mut mi = 0.0;

    
    for i in 0..n_bins {
        for j in 0..n_bins {
            let joint_prob = joint_counts[i][j] as f64 / n;
            if joint_prob > 0.0 {
                let marginal_prob_i = marginal_i[i] as f64 / n;
                let marginal_prob_j = marginal_j[j] as f64 / n;

                mi += joint_prob * (joint_prob / (marginal_prob_i * marginal_prob_j)).ln();
            }
        }
    }

    mi
}

/// Calculates the joint entropy between two column references
///
/// ### Params
///
/// * `col_i` - Column reference to the first column to compare
/// * `col_j` - Column reference to the second column to compare
/// * `n_bins` - Optional number of bins. If not provided, will default to
///   `sqrt(nrows)`.
///
/// ### Returns
///
/// The joint entropy
pub fn calculate_joint_entropy(
    col_i: ColRef<usize>,
    col_j: ColRef<usize>,
    n_bins: Option<usize>,
) -> f64 {
    let n_rows = col_i.nrows();
    let n_bins = n_bins.unwrap_or_else(|| (n_rows as f64).sqrt() as usize);

    let mut joint_counts = vec![vec![0usize; n_bins]; n_bins];

    for i in 0..n_rows {
        joint_counts[col_i[i]][col_j[i]] += 1;
    }

    let n = n_rows as f64;
    let mut joint_entropy = 0.0;

    
    for i in 0..n_bins {
        for j in 0..n_bins {
            let joint_prob = joint_counts[i][j] as f64 / n;
            if joint_prob > 0.0 {
                joint_entropy -= joint_prob * joint_prob.ln();
            }
        }
    }

    joint_entropy
}

/// Calculates the entropy of a column reference
///
/// The column needs to be binned
///
/// ### Params
///
/// * `col` - Column reference for which to calculate entropy
/// * `n_bins` - Optional number of bins. If not provided, will default to
///   `sqrt(nrows)`.
///
/// ### Returns
///
/// The entropy of the column
pub fn calculate_entropy(col: ColRef<usize>, n_bins: Option<usize>) -> f64 {
    let n_rows = col.nrows();
    let n_bins = n_bins.unwrap_or_else(|| (n_rows as f64).sqrt() as usize);

    let mut counts = vec![0usize; n_bins];

    for i in 0..n_rows {
        counts[col[i]] += 1;
    }

    let n = n_rows as f64;
    let mut entropy = 0.0;

    for &count in &counts {
        if count > 0 {
            let prob = count as f64 / n;
            entropy -= prob * prob.ln();
        }
    }

    entropy
}
