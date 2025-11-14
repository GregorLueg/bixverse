// use extendr_api::*;

// /// Extract the kNN data from an R sc_knn
// ///
// /// ### Params
// ///
// /// * `knn_data` - The `sc_knn` class from which to extract the data
// fn get_knn_data_from_r(knn_data: Robj) -> (Vec<Vec<usize>>, Vec<Vec<f32>>) {
//     let indices: RMatrix<i32> = knn_data.dollar("indices").unwrap().as_matrix().unwrap();
//     let dist: RMatrix<f64> = knn_data.dollar("dist").unwrap().as_matrix().unwrap();

//     let nrow = indices.nrows();
//     let ncol = indices.ncols();

//     let indices_data = indices.data();
//     let dist_data = dist.data();

//     let knn_indices: Vec<Vec<usize>> = (0..ncol)
//         .map(|j| {
//             (0..nrow)
//                 .map(|i| indices_data[i + j * nrow] as usize)
//                 .collect::<Vec<usize>>()
//         })
//         .collect();

//     let knn_dist: Vec<Vec<f32>> = (0..ncol)
//         .map(|j| {
//             (0..nrow)
//                 .map(|i| dist_data[i + j * nrow] as f32)
//                 .collect::<Vec<f32>>()
//         })
//         .collect();

//     (knn_indices, knn_dist)
// }
