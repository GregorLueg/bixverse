use faer::Mat;

//////////////////
// VECTOR STUFF //
//////////////////

/// Flatten a nested vector
pub fn flatten_vector<T>(
  vec: Vec<Vec<T>>
) -> Vec<T> {
  vec.into_iter().flatten().collect()
}

/// Get the maximum value from an f64 array. 
pub fn array_f64_max(
  arr: &[f64]
) -> f64 {
  let mut max_val = arr[0];
  for number in arr{
    if *number > max_val {
      max_val = *number
    }
  }
  max_val
}

/// Get the minimum value from an f64 array. 
pub fn array_f64_min(
  arr: &[f64]
) -> f64 {
  let mut min_val = arr[0];
  for number in arr{
    if *number < min_val {
      min_val = *number
    }
  }
  min_val
}

// Get the mean value from an f64 array
// pub fn array_f64_mean(
//   x: &[f64]
// ) -> f64 {
//   let len_x = x.len();
//   let sum_x: f64 = x
//     .iter()
//     .sum();
//   sum_x / len_x as f64
// }

// Get the variance from an f64 array
// pub fn array_f64_var(
//   x: &[f64]
// ) -> f64 {
//   let mean_a = array_f64_mean(x);
//   let var: f64 = x
//     .iter()
//     .map(|x| {
//       (x - mean_a).powi(2)
//     })
//     .sum();

//   var / (x.len() - 1) as f64
// }

/// Generate the rank of a vector with tie correction.
pub fn rank_vector(
  vec: &[f64]
) -> Vec<f64> {
  let mut vec_index: Vec<(f64, usize)> = vec
    .iter()
    .copied()
    .enumerate()
    .map(|(i, x)| {
      (x, i)
    })
    .collect();
    
  vec_index
    .sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));

  let mut ranks = vec![0.0; vec.len()];
  let mut i = 0;

  while i < vec_index.len() {
    let value = vec_index[i].0;
    let mut j = i + 1;

    // Tie correction
    while j < vec_index.len() && vec_index[j].0 == value {
      j += 1;
    }

    let rank = (i + j - 1) as f64 / 2.0 + 1.0;

    vec_index[i..j].iter().for_each(|&(_, original_index)| {
      ranks[original_index] = rank;
    });

    i = j;
  }

  ranks  
}

//////////////////
// MATRIX STUFF //
//////////////////

/// Transform a nested vector into a faer matrix
pub fn nested_vector_to_faer_mat(
  nested_vec: Vec<Vec<f64>>
) -> faer::Mat<f64> {
  let nrow = nested_vec[0].len();
  let ncol = nested_vec.len();
  let data = flatten_vector(nested_vec);
  Mat::from_fn(nrow, ncol, |i, j| data[i + j * nrow])
}

/// Create a diagonal matrix with the vector values in the diagonal and the rest being 0's
pub fn faer_diagonal_from_vec(
  vec: Vec<f64>
) -> Mat<f64> {
  let len = vec.len();
  Mat::from_fn(
    len,
    len, 
    |row, col| if row == col { vec[row] } else { 0.0 }
  )
}

/// Get the index positions of the upper triangle of a symmetric matrix
pub fn upper_triangle_indices(
  n_dim: usize,
  offset: usize
) -> (Vec<usize>, Vec<usize>) {
  let mut row_indices: Vec<usize> = Vec::new();
  let mut col_indices: Vec<usize> = Vec::new();

  for row in 0..n_dim{
    let start_col = std::cmp::max(row + offset, 0) as usize;
    if start_col < n_dim {
      for col in start_col..n_dim {
        row_indices.push(row);
        col_indices.push(col);
      }
    }
  }

  (row_indices, col_indices)
}



