use extendr_api::prelude::*;

use crate::helpers::sparse_graph_dict::*;
use crate::helpers::structs_sparse::SparseColumnMatrix;
use crate::utils::r_rust_interface::{faer_to_r_matrix, r_matrix_to_faer, sparse_matrix_to_list};

/// @export
#[extendr]
fn rs_sparse_dict_dgrdl(dat: RMatrix<f64>, dgrdl_params: List, verbose: bool) -> List {
    let x = r_matrix_to_faer(&dat);

    let dgrdl_params = DgrdlParams::from_r_list(dgrdl_params);

    let dgrdl_object = Dgrdl::new(dgrdl_params);

    let res: DgrdlResults = dgrdl_object.fit(&x, verbose);

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
