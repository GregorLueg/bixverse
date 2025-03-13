use faer::{Mat,linalg::solvers::{Solve, PartialPivLu}};
use rayon::iter::*;
use rand::prelude::*;
use rand_distr::Normal;

use crate::utils_rust::*;
use crate::helpers_linalg::scale_matrix_col;

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
  w: &Mat<f64>
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

/// Generate a w_init matrix of size n_comp * n_comp given a random seed.
pub fn create_w_init(
  n_comp: usize,
  seed: u64,
) -> faer::Mat<f64> {
  let mut rng = StdRng::seed_from_u64(seed);
  let normal = Normal::new(0.0, 1.0).unwrap();
  let vec_size = n_comp.pow(2);
  let data: Vec<f64> = (0..vec_size)
    .map(|_| {normal.sample(&mut rng)})
    .collect();

  Mat::from_fn(
      n_comp, n_comp, |i, j| data[i + j * n_comp]
  )
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
  x: &Mat<f64>,
  w_init: &Mat<f64>,
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
    let wx: Mat<f64> = &w * x;

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

    let w1 = update_mix_mat(&(v1 - v2));

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
  x: &Mat<f64>,
  w_init: &Mat<f64>,
  tol: f64,
  maxit: usize,
  verbose: bool,
) -> IcaRes {
  let p = x.ncols();
  let mut w = update_mix_mat(w_init);
  let mut lim = vec![1000_f64; maxit];

  let mut it = 0;
  while it < maxit && lim[it] > tol {
    let wx: Mat<f64> = &w * x;

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

    let w1 = update_mix_mat(&(v1 - v2));

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


/// Iterate through a set of 
pub fn stabilised_ica_iters(
  x_whiten: Mat<f64>,
  k: Mat<f64>,
  no_comp: usize,
  no_iters: usize,
  maxit: usize,
  alpha: f64,
  tol: f64,
  ica_type: &str,
  random_seed: usize,
  verbose: bool,
) -> (Mat<f64>, Vec<bool>) {
  // -> (Mat<f64>, Vec<bool>) 
  // Generate the random w_inits
  let w_inits: Vec<Mat<f64>> = (0..no_iters)
    .map(|iter| {
      create_w_init(
        no_comp, (random_seed + iter) as u64
      )
    })
    .collect();
  let k_ncol = k.nrows();  
  let k_red = k.get(0..no_comp, 0..k_ncol);
  let x1 = k_red * x_whiten;

  let ica_type = parse_ica_type(ica_type).unwrap();
  let iter_res: Vec<(Mat<f64>, f64)>  = w_inits
    .par_iter()
    .map(|w_init| {
      match ica_type {
        IcaType::Exp => fast_ica_exp(&x1, w_init, tol, maxit, verbose),
        IcaType::LogCosh => fast_ica_logcosh(&x1, w_init, tol, alpha, maxit, verbose),
      }
    }).collect();

  let mut convergence = Vec::new();
  let mut a_matrices = Vec::new();

  for (a, final_tol) in iter_res {
    a_matrices.push(a);
    convergence.push(final_tol < tol);
  };

  let s_matrices: Vec<Mat<f64>> = a_matrices
    .par_iter()
    .map(|a| {
      let w = a * k_red;
      let to_solve = w.clone() * w.transpose();
      let identity = Mat::<f64>::identity(to_solve.nrows(), to_solve.ncols());
      let lu = PartialPivLu::new(to_solve.as_ref());
      let solved = lu.solve(&identity);
      w.transpose() * solved
    })
    .collect();

  let s_combined = colbind_matrices(s_matrices);

  (s_combined, convergence)
}