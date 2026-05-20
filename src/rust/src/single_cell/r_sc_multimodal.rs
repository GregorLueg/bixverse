use bixverse_rs::prelude::*;
use bixverse_rs::single_cell::multi_modal::wnn::*;
use bixverse_rs::single_cell::sc_batch_correction::batch_utils::cosine_normalise;
use extendr_api::*;
use faer::Mat;

////////////////////
// extendr Module //
////////////////////

extendr_module! {
    // module
    mod r_sc_multimodal;
    // functions
    fn rs_wnn;
}

/////////
// WNN //
/////////

/// Run the weighted nearest neighbour algorithm
///
/// @param modality_emb_one Numerical matrix of the first modality. For example
/// the PCA (or other embeddings) from the transcriptomics.
/// @param modality_emb_two Numerical matrix of the second modality. For example
/// the PCA (or other embeddings) from the ADT counts.
/// @param wnn_params Named list. The weighted nearest neighbour parameters.
/// @param seed Integer. For reproducibility purposes.
/// @param verbose Integer. `0L` - quiet; `1L` - normal verbosity; `2L` -
/// detailed verbosity.
///
/// @returns A list with
/// \itemize{
///   \item indices An integer matrix representing the indices of the
///  approximate nearest neighbours.
///  \item dist - An numerical matrix representing the distances to the nearest
///  neighbours.
///  \item modality_one_weights - The weights of the first modality.
///  \item modality_two_weights - The weights of the second modality.
/// }
///
/// @export
#[extendr]
fn rs_wnn(
    modality_emb_one: RMatrix<f64>,
    modality_emb_two: RMatrix<f64>,
    wnn_params: List,
    seed: usize,
    verbose: usize,
) -> Result<List> {
    let verbosity = parse_verbosity_level(verbose);
    let mut wnn_params = WnnParams::from_r_list(wnn_params)?;

    let k_min = wnn_params.knn_range;

    wnn_params.knn_params.k = k_min;

    let is_squared_dist = wnn_params.knn_params.ann_dist == "euclidean";

    let sqrt_distances = |x: &mut [Vec<f32>]| {
        for row in x.iter_mut() {
            for d in row.iter_mut() {
                *d = d.sqrt();
            }
        }
    };

    let modality_emb_one = r_matrix_to_faer_fp32(&modality_emb_one);
    let modality_emb_two = r_matrix_to_faer_fp32(&modality_emb_two);

    let modality_emb_one = cosine_normalise(&modality_emb_one);
    let modality_emb_two = cosine_normalise(&modality_emb_two);

    // deal with the knn graphs
    // first modality

    if verbosity.normal_verbosity() {
        println!(
            "Running the kNN generation for the first modality with {:?} neighbours.",
            k_min
        )
    }

    let (indices_1, dist_1) = generate_knn_with_dist(
        modality_emb_one.as_ref(),
        &wnn_params.knn_params,
        true,
        false,
        seed,
        verbosity.detailed_verbosity(),
    )
    .to_extendr()?;

    let mut dist_1 = dist_1.unwrap();

    if is_squared_dist {
        sqrt_distances(&mut dist_1);
    }

    // second modality
    if verbosity.normal_verbosity() {
        println!(
            "Running the kNN generation for the second modality with {:?} neighbours.",
            k_min
        )
    }

    let (indices_2, dist_2) = generate_knn_with_dist(
        modality_emb_two.as_ref(),
        &wnn_params.knn_params,
        true,
        false,
        seed,
        verbosity.detailed_verbosity(),
    )
    .to_extendr()?;

    let mut dist_2 = dist_2.unwrap();

    if is_squared_dist {
        sqrt_distances(&mut dist_2);
    }

    let wnn_modality_1 = ModalityInput::new(modality_emb_one.as_ref(), &indices_1, &dist_1);
    let wnn_modality_2 = ModalityInput::new(modality_emb_two.as_ref(), &indices_2, &dist_2);

    let wnn_res: WnnResult =
        compute_wnn([wnn_modality_1, wnn_modality_2], &wnn_params, verbose).to_extendr()?;

    let index_mat = Mat::from_fn(
        wnn_res.wnn_indices.len(),
        wnn_res.wnn_indices[0].len(),
        |i, j| wnn_res.wnn_indices[i][j] as i32,
    );
    let dist_mat = Mat::from_fn(
        wnn_res.wnn_indices.len(),
        wnn_res.wnn_indices[0].len(),
        |i, j| wnn_res.wnn_distances[i][j] as f64,
    );

    let modality_one_weights = wnn_res.modality_weights[0].clone().r_float_convert();
    let modality_two_weights = wnn_res.modality_weights[1].clone().r_float_convert();

    Ok(list!(
        indices = faer_to_r_matrix(index_mat.as_ref()),
        dist = faer_to_r_matrix(dist_mat.as_ref()),
        dist_metric = "kernelised pseudo-distance".to_string(),
        modality_one_weights = modality_one_weights,
        modality_two_weights = modality_two_weights
    ))
}
