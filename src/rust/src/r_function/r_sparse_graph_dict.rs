use extendr_api::prelude::*;

use crate::helpers::sparse_graph_dict::*;
use crate::helpers::structs_sparse::SparseColumnMatrix;
use crate::utils::r_rust_interface::{faer_to_r_matrix, r_matrix_to_faer, sparse_matrix_to_list};

/// Generate a sparse dictionary with DGRDL
///
/// @description This is the Rust implementation of dual graph regularised
/// dictionary learning in the implementation of Pan, et al., Cell Systems,
/// 2022.
///
/// @param x Numerical matrix. Rows = samples, columns = features.
/// @param dgrdl_params A list with the parameters for the algorithm. Expects
/// the following items.
/// \itemize{
///   \item sparsity - Sparsity constraint (max non-zero coefficients per signal).
///   \item dict_size - Size of the dictionary.
///   \item alpha - Float. Sample context regularisation weight. The higher the stronger
///   the regularisation.
///   \item beta - Float. Feature context regularisation weight. The higher the stronger
///   the regularisation.
///   \item max_iter - Integer. Maximum iteration for the algorithm.
///   \item k_neighbours - Integer. Number of k neighbours for the sample and feature
///   Laplacian matrix for the regularisation
///   \item admm_iter Integer. Number of iterations for using alternating direction
///   method of multipliers (ADMM).
///   \item rho Float. ADMM step size.
/// }
/// @param seed Integer. Seed for the initialisation of the algorithm.
/// @param verbose Boolean. Controls the verbosity of the function and reports timing
/// of individual steps.
///
/// @returns A list with the following elements:
///  \itemize{
///   \item dictionary - The dictionary of samples x dict_size.
///   \item coefficients - The feature loadings of size dict_size x features.
///   \item feature_laplacian - The KNN graph laplacian of the features in a
///   sparse format list.
///   \item sample_laplacian - The KNN graph laplacian of the samples in a
///   sparse format list.
/// }
///
/// @export
#[extendr]
fn rs_sparse_dict_dgrdl(x: RMatrix<f64>, dgrdl_params: List, seed: usize, verbose: bool) -> List {
    let x = r_matrix_to_faer(&x);

    let dgrdl_params = DgrdlParams::from_r_list(dgrdl_params);

    let dgrdl_object = Dgrdl::new(dgrdl_params);

    let res: DgrdlResults = dgrdl_object.fit(&x, seed, verbose);

    let feature_laplacian = SparseColumnMatrix::from_dense_matrix(res.feature_laplacian.as_ref());
    let sample_laplacian = SparseColumnMatrix::from_dense_matrix(res.sample_laplacian.as_ref());

    list!(
        dictionary = faer_to_r_matrix(res.dictionary.as_ref()),
        coefficients = faer_to_r_matrix(res.coefficients.as_ref()),
        feature_laplacian = sparse_matrix_to_list(feature_laplacian),
        sample_laplacian = sparse_matrix_to_list(sample_laplacian),
    )
}

extendr_module! {
    mod r_sparse_graph_dict;
    fn rs_sparse_dict_dgrdl;
}
