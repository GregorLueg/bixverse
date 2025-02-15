use faer::Mat;

pub fn column_covariance(mat: &Mat<f64>) -> Mat<f64> {
    let n_rows = mat.nrows();
    
    // Calculate column means: (1/n) * 1^T * X
    let ones = Mat::from_fn(n_rows, 1, |_, _| 1.0);
    let means = (ones.transpose() * mat) / n_rows as f64;
    
    // Center the matrix by subtracting means from each column
    let centered = mat - &ones * &means;
    
    // Calculate covariance: (1/(n-1)) * (X - 1*means)^T * (X - 1*means)
    let covariance = (centered.transpose() * &centered) / (n_rows - 1) as f64;
    
    covariance
}