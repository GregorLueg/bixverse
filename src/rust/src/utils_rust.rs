use faer::Mat;

/// Flatten a nested vector
pub fn flatten_vector<T>(
    vec: Vec<Vec<T>>
) -> Vec<T> {
    vec.into_iter().flatten().collect()
}

/// Get the maximum value from an f64 array. 
pub fn array_f64_max(
    arr: &Vec<f64>
) -> f64 {
    let mut max_val = arr[0];
    for number in arr{
        if *number > max_val {
            max_val = *number
        }
    }
    max_val
}

// /// Get the maximum value from an f64 array. 
// pub fn array_f64_min(
//     arr: &Vec<f64>
// ) -> f64 {
//     let mut min_val = arr[0];
//     for number in arr{
//         if *number < min_val {
//             min_val = *number
//         }
//     }
//     min_val
// }

/// Transform a nested vector into a faer matrix
pub fn nested_vector_to_faer_mat(
  nested_vec: Vec<Vec<f64>>
) -> faer::Mat<f64> {
  let nrow = nested_vec[0].len();
  let ncol = nested_vec.len();
  let data = flatten_vector(nested_vec);
  Mat::from_fn(nrow, ncol, |i, j| data[i + j * nrow])
}


/// Create a diagonal matrix with the vector values in the diagonal
/// and the rest being 0's
pub fn faer_diagonal_from_vec(
    vec: Vec<f64>
) -> Mat<f64> {
    let len = vec.len();
    Mat::from_fn(len, len, |row, col| if row == col {vec[row]} else {0.0})
}