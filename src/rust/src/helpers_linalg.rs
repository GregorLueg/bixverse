use faer::Mat;

use crate::utils_rust::{faer_diagonal_from_vec, array_f64_min};

/// Enum for the ICA types
#[derive(Debug)]
pub enum IcaType {
  Exp,
  LogCosh,
}

/// Type alias of the ICA results
type IcaRes = (faer::Mat<f64>, f64);

/// Scale a matrix by its mean (column wise)
pub fn scale_matrix_col(mat: &Mat<f64>) -> Mat<f64>{
  let n_rows = mat.nrows();
  let ones = Mat::from_fn(n_rows, 1, |_, _| 1.0);
  let means = (ones.transpose() * mat) / n_rows as f64;
  mat - &ones * &means
}

/// Calculate the co-variance
pub fn column_covariance(mat: &Mat<f64>) -> Mat<f64> {
  let n_rows = mat.nrows();
  let centered = scale_matrix_col(mat);
  let covariance = (centered.transpose() * &centered) / (n_rows - 1) as f64;
    
  covariance
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


/// Whiten a matrix. This is needed pre-processing for ICA.
pub fn prepare_whitening(
  x: Mat<f64>
) -> (faer::Mat<f64>, faer::Mat<f64>) {
  let n = x.nrows();  

  let centered = scale_matrix_col(&x);

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




