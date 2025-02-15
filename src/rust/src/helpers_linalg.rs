use faer::Mat;

/// Scale a matrix by its mean (column wise)
pub fn scale_matrix(mat: &Mat<f64>) -> Mat<f64>{
  let n_rows = mat.nrows();
  let ones = Mat::from_fn(n_rows, 1, |_, _| 1.0);
  let means = (ones.transpose() * mat) / n_rows as f64;
  mat - &ones * &means
}

/// Calculate the co-variance
pub fn column_covariance(mat: &Mat<f64>) -> Mat<f64> {
  let n_rows = mat.nrows();
    
  // Center the matrix by subtracting means from each column
  let centered = scale_matrix(mat);
    
  // Calculate covariance: (1/(n-1)) * (X - 1*means)^T * (X - 1*means)
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

  // Extract the real part of the eigen values and vectors
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


pub fn whiten_matrix(
  x: Mat<f64>
) -> faer::Mat<f64> {
  let n = x.nrows();  

  let centered = scale_matrix(&x);

  let centered = centered.transpose();

  let covariance =  centered * centered.transpose() / n as f64;

  // SVD
  let svd_res = covariance.svd().unwrap();

  // Get d
  let s = svd_res.S();
  let s = s
    .column_vector()
    .iter()
    .map(|x| {1_f64 / x.sqrt()})
    .collect::<Vec<_>>();
  let len_s = s.len();
  let d = Mat::from_fn(len_s, len_s, |row, col| if row == col {s[row]} else {0.0});

  // Get k
  let u = svd_res.U();
  let k = d * u.transpose();

  k * centered
}

pub fn fast_ica_logcosh(
  x: Mat<f64>,
  n_ica: usize,
  tol: f64,
  alpha: f64,
  maxit: usize,
) {

}