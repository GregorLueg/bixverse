use bixverse_rs::prelude::*;
use bixverse_rs::single_cell::multi_modal::wnn::WnnParams;
use extendr_api::*;

use crate::single_cell::utils::knn_data_to_rust;

////////////////////
// extendr Module //
////////////////////

extendr_module! {
    mod r_sc_multimodal;
}

/////////
// WNN //
/////////

fn rs_wnn(
    modality_emb_one: RMatrix<f64>,
    modality_emb_two: RMatrix<f64>,
    knn_data_one: Nullable<List>,
    knn_data_two: Nullable<List>,
    wnn_params: List,
    seed: usize,
    verbose: usize,
) -> Result<()> {
    let verbosity = parse_verbosity_level(verbose);
    let wnn_params = WnnParams::from_r_list(wnn_params)?;

    let modality_emb_one = r_matrix_to_faer_fp32(&modality_emb_one);
    let modality_emb_two = r_matrix_to_faer_fp32(&modality_emb_two);

    let knn_one_provided = knn_data_one != Nullable::Null;
    let knn_two_provided = knn_data_two != Nullable::Null;

    // deal with the knn graphs
    // first modality
    let (indices_1, dist_1) = if knn_one_provided {
        if verbosity.normal_verbosity() {
            println!("Using the provided kNN data for modality one.")
        }

        let knn_data = knn_data_one
            .into_robj()
            .as_list()
            .ok_or_else(|| Error::Other("'knn_data' is not a list".into()))?;

        let (knn_indices, knn_dist, _, _) = knn_data_to_rust(knn_data)?;

        (knn_indices, knn_dist)
    } else {
        if verbosity.normal_verbosity() {
            println!("No kNN provided for modality one. Running the kNN generation.")
        }

        let (knn, dist) = generate_knn_with_dist(
            modality_emb_one.as_ref(),
            &wnn_params.knn_params,
            true,
            false,
            seed,
            verbosity.detailed_verbosity(),
        )
        .to_extendr()?;

        (knn, dist.unwrap())
    };

    // second modality
    let (indices_2, dist_2) = if knn_two_provided {
        if verbosity.normal_verbosity() {
            println!("Using the provided kNN data for modality two.")
        }

        let knn_data = knn_data_two
            .into_robj()
            .as_list()
            .ok_or_else(|| Error::Other("'knn_data' is not a list".into()))?;

        let (knn_indices, knn_dist, _, _) = knn_data_to_rust(knn_data)?;

        (knn_indices, knn_dist)
    } else {
        if verbosity.normal_verbosity() {
            println!("No kNN provided for modality one. Running the kNN generation.")
        }

        let (knn, dist) = generate_knn_with_dist(
            modality_emb_two.as_ref(),
            &wnn_params.knn_params,
            true,
            false,
            seed,
            verbosity.detailed_verbosity(),
        )
        .to_extendr()?;

        (knn, dist.unwrap())
    };

    Ok(())
}
