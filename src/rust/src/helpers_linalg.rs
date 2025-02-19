use faer::Mat;
use rayon::iter::*;
use crate::utils_rust::*;
use crate::utils_stats::z_scores_to_pval;

//////////////////////////////
// ENUMS, TYPES, STRUCTURES //
//////////////////////////////

/// Enum for the ICA types
#[derive(Debug)]
pub enum IcaType {
  Exp,
  LogCosh,
}

/// Type alias of the ICA results
type IcaRes = (faer::Mat<f64>, f64);

/// Structure for DiffCor results
pub struct DiffCorRes {
  pub r_a: Vec<f64>,
  pub r_b: Vec<f64>,
  pub z_score: Vec<f64>,
  pub p_vals: Vec<f64>
}

//////////////////////////////
// SCALING, COVAR, COR, PCA //
//////////////////////////////

/// Scale a matrix by its mean (column wise)
pub fn scale_matrix_col(
  mat: &Mat<f64>,
  scale_sd: bool
) -> Mat<f64>{
  let n_rows = mat.nrows();
  let ones = Mat::from_fn(n_rows, 1, |_, _| 1.0);
  let means = (ones.transpose() * mat) / n_rows as f64;
  let centered = mat - &ones * &means;

  if !scale_sd {
    return centered;
  }

  let squared_diff = Mat::from_fn(
    centered.nrows(),
    centered.ncols(),
  |i, j| {
      let x = centered.get(i, j);
      x.powi(2)
    }
  );

  let sum_squared_diff = ones.transpose() * squared_diff;
  let variances = sum_squared_diff / (n_rows as f64 - 1.0);

  let standard_dev = Mat::from_fn(
    variances.nrows(),
    variances.ncols(),
    |i, j| {
      let x = variances.get(i, j);
      let x = x.sqrt();
      if x < 1e-10 { 1.0 } else { x }
    }
  );

  let mut scaled = centered.clone();
  let n_cols = mat.ncols();

  for j in 0..n_cols {
    for i in 0..n_rows {
      scaled[(i, j)] = centered[(i, j)] / standard_dev[(0, j)];
    }
  }
  
  scaled
}

/// Calculate the column-wise co-variance
pub fn column_covariance(mat: &Mat<f64>) -> Mat<f64> {
  let n_rows = mat.nrows();
  let centered = scale_matrix_col(mat, false);
  let covariance = (centered.transpose() * &centered) / (n_rows - 1) as f64;
    
  covariance
}

/// Calculate the column-wise correlation. Option to use spearman.
pub fn column_correlation(
  mat: &Mat<f64>,
  spearman: bool
) -> Mat<f64>{
  let mat = if spearman {
    let ranked_vecs: Vec<Vec<f64>> = mat
      .par_col_iter()
      .map(|x_i| { 
        let x_i: Vec<f64> = x_i.iter().copied().collect();
        rank_vector(&x_i)
      }).collect();

    nested_vector_to_faer_mat(ranked_vecs)
  } else {
    mat.cloned()
  };

  let scaled = scale_matrix_col(&mat, true);

  let nrow = scaled.nrows() as f64;

  let cor = scaled.transpose() * &scaled / (nrow - 1_f64);

  cor
}

/// Calculate differential calculations
pub fn calculate_diff_correlation(
  mat_a: &Mat<f64>,
  mat_b: &Mat<f64>,
  no_sample_a: usize,
  no_sample_b: usize,
  spearman: bool,
) -> DiffCorRes {
  let mut cors_a: Vec<f64> = Vec::new();
  let mut cors_b: Vec<f64> = Vec::new();

  let upper_triangle_indices = upper_triangle_indices(
    mat_a.ncols(),
    1
  );

  for (&r, &c) in upper_triangle_indices.0.iter().zip(upper_triangle_indices.1.iter()) {
    cors_a.push(*mat_a.get(r, c));
    cors_b.push(*mat_b.get(r, c));
  }

  cors_a
    .par_iter_mut()
    .for_each(|x| *x = x.atanh());
  cors_b
    .par_iter_mut()
    .for_each(|x| *x = x.atanh());

  // Constant will depend on if Spearman or Pearson  
  let constant = if spearman { 1.06 } else { 1.0 } ;
  let denominator = ((constant / (no_sample_a as f64 - 3.0)) + (constant / (no_sample_b as f64 - 3.0))).sqrt();

  let z_scores: Vec<f64> = cors_a
    .par_iter()
    .zip(
      cors_b.par_iter()
    )
    .map(|(a, b)| {
      (a - b) / denominator
    })
    .collect();

  let p_values = z_scores_to_pval(&z_scores);

  DiffCorRes{
    r_a: cors_a,
    r_b: cors_b,
    z_score: z_scores,
    p_vals: p_values
  }
}

/// Get the eigenvalues and vectors from a symmetric matrix
pub fn get_top_eigenvalues(
  matrix: &Mat<f64>, 
  top_n: usize
) -> Vec<(f64, Vec<f64>)> {
  // Ensure the matrix is square
  assert!(matrix.nrows() == matrix.ncols(), "Matrix must be square");

  let eigendecomp = matrix.eigen_from_real().unwrap();

  let s = eigendecomp.S();
  let u = eigendecomp.U();

  // Extract the real part of the eigenvalues and vectors
  let mut eigenpairs = s
    .column_vector()
    .iter()
    .zip(u.col_iter())
    .map(|(l, v)| {
      let l_real = l.re;
      let v_real = v
        .iter()
        .map(|v_i| {v_i.re})
        .collect::<Vec<f64>>();
      (l_real, v_real)
    })
    .collect::<Vec<(f64, Vec<f64>)>>();

  // Sort and return Top N
  eigenpairs.sort_by(|a, b| b.0.abs().total_cmp(&a.0.abs()));

  let res: Vec<(f64, Vec<f64>)> = eigenpairs
    .into_iter()
    .take(top_n)
    .collect();

  res
}

/////////
// ICA //
/////////

/// Whiten a matrix. This is needed pre-processing for ICA.
pub fn prepare_whitening(
  x: Mat<f64>
) -> (faer::Mat<f64>, faer::Mat<f64>) {
  let n = x.nrows();  

  let centered = scale_matrix_col(&x, false);

  let centered = centered.transpose();

  let v =  centered * centered.transpose() / n as f64;

  // SVD
  let svd_res = v.svd().unwrap();

  // Get d
  let s = svd_res.S();
  let s = s
    .column_vector()
    .iter()
    .map(|x| {1_f64 / x.sqrt()})
    .collect::<Vec<_>>();

  let d = faer_diagonal_from_vec(s);

  // Get k
  let u = svd_res.U();

  let u_transpose = u.transpose();

  let k = d * u_transpose;

  (centered.cloned(), k)
}


/// Update the mixing matrix for ICA
pub fn update_mix_mat(
  w: Mat<f64>
) -> faer::Mat<f64> {
  // SVD
  let svd_res = w.svd().unwrap();

  let s = svd_res.S();
  let u = svd_res.U();
  let s = s
    .column_vector()
    .iter()
    .map(|x| {
      1_f64 / x
    })
    .collect::<Vec<_>>();
  let d = faer_diagonal_from_vec(s);

  u * d * u.transpose() * w
}


/// Parsing the ICA types
pub fn parse_ica_type(s: &str) -> Option<IcaType> {
  match s.to_lowercase().as_str() {
    "exp" => Some(IcaType::Exp),
    "logcosh" => Some(IcaType::LogCosh),
    _ => None,
  }
}

/// Fast ICA implementation based on logcosh.
pub fn fast_ica_logcosh(
  x: Mat<f64>,
  w_init: Mat<f64>,
  tol: f64,
  alpha: f64,
  maxit: usize,
  verbose: bool,
) -> IcaRes {
  let p = x.ncols();
  let mut w = update_mix_mat(w_init);
  let mut lim = vec![1000_f64; maxit];

  let mut it = 0;

  while it < maxit && lim[it] > tol {
    let wx: Mat<f64> = &w * &x;

    let gwx = Mat::from_fn(
      wx.nrows(),
      wx.ncols(),
    |i, j| {
        let x = wx.get(i, j);
        (alpha * x).tanh()
      }
    );

    let v1 = &gwx * x.transpose() / p as f64; 

    let gwx_2 = alpha * Mat::from_fn(
      gwx.nrows(),
      gwx.ncols(),
    |i, j| {
        let x = gwx.get(i, j);
        1_f64 - x.powi(2)
      }
    );

    let ones = Mat::from_fn(p, 1, |_, _| 1.0);
    let row_means = (&gwx_2 * &ones) * (1.0 / p as f64);

    let row_means_vec: Vec<f64> = row_means.as_ref().col_iter()
      .flat_map(|col| col.iter())
      .copied()
      .collect();

    let v2 = faer_diagonal_from_vec(row_means_vec) * w.clone();

    let w1 = update_mix_mat(v1 - v2);

    let w1_up = w1.clone() * w.transpose();

    let tol_it = w1_up.diagonal()
      .column_vector()
      .iter()
      .map(|x| (x.abs() - 1.0).abs())
      .fold(f64::NEG_INFINITY, f64::max);

    if it + 1 < maxit {
      lim[it + 1] = tol_it
    }

    if verbose {
      println!("Iteration: {:?}, tol: {:?}", it + 1, tol_it)
    }

    w = w1;

    it += 1;
  }

  let min_tol = array_f64_min(&lim);

  (w, min_tol)
}


/// Fast ICA implementation based on exp.
pub fn fast_ica_exp(
  x: Mat<f64>,
  w_init: Mat<f64>,
  tol: f64,
  maxit: usize,
  verbose: bool,
) -> IcaRes {
  let p = x.ncols();
  let mut w = update_mix_mat(w_init);
  let mut lim = vec![1000_f64; maxit];

  let mut it = 0;
  while it < maxit && lim[it] > tol {
    let wx: Mat<f64> = &w * &x;

    let gwx = Mat::from_fn(
      wx.nrows(),
      wx.ncols(),
    |i, j| {
        let x = wx.get(i, j);
        x * (-x.powi(2) / 2.0).exp()
      }
    );

    let v1 = &gwx * x.transpose() / p as f64;   

    let gwx_2 = Mat::from_fn(
      wx.nrows(),
      wx.ncols(),
    |i, j| {
        let x = wx.get(i, j);
        (1.0 - x.powi(2)) * (-x.powi(2) / 2.0).exp()
      }
    );

    let ones = Mat::from_fn(p, 1, |_, _| 1.0);
    let row_means = (&gwx_2 * &ones) * (1.0 / p as f64);

    let row_means_vec: Vec<f64> = row_means.as_ref().col_iter()
      .flat_map(|col| col.iter())
      .copied()
      .collect();

    let v2 = faer_diagonal_from_vec(row_means_vec) * w.clone();

    let w1 = update_mix_mat(v1 - v2);

    let w1_up = w1.clone() * w.transpose();
    
    let tol_it = w1_up.diagonal()
      .column_vector()
      .iter()
      .map(|x| (x.abs() - 1.0).abs())
      .fold(f64::NEG_INFINITY, f64::max);
    
    if it + 1 < maxit {
      lim[it + 1] = tol_it
    }

    if verbose {
      println!("Iteration: {:?}, tol: {:?}", it + 1, tol_it)
    }

    w = w1;

    it += 1;
  }

  let min_tol = array_f64_min(&lim);

  (w, min_tol)
}




