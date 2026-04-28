use bixverse_rs::single_cell::sc_batch_correction::harmony_v2::harmony_v2;
use bixverse_rs::single_cell::sc_batch_correction::harmony_v2::HarmonyParamsV2;
use extendr_api::*;
use faer::Mat;

use bixverse_rs::prelude::*;
use bixverse_rs::single_cell::sc_batch_correction::bbknn::*;
use bixverse_rs::single_cell::sc_batch_correction::fast_mnn::*;
use bixverse_rs::single_cell::sc_batch_correction::harmony::*;
use bixverse_rs::single_cell::sc_processing::metrics::*;

////////////////////
// extendr Module //
////////////////////

extendr_module! {
    mod r_sc_batch_corr;
    // metrics
    fn rs_kbet;
    fn rs_batch_silhouette_width;
    fn rs_batch_lisi;
    // batch corrections
    fn rs_bbknn;
    fn rs_bbknn_filtering;
    fn rs_mnn;
    fn rs_harmony;
    fn rs_harmony_v2;
}

///////////////
// Functions //
///////////////

/// Calculate kBET type scores
///
/// @description
/// The function takes in a kNN matrix and a batch vector indicating which
/// cell belongs to which batch. The function will check for the neighbourhood
/// of each cell if the proportion of represented batches are different from
/// the overall batch proportions. Good mixing of batches would mean very
/// cells have significant differences; bad mixing a lot of the batches
/// have bad mixing.
///
/// @param knn_mat Integer matrix. The rows represent the cells and the
/// columns the neighbour indices.
/// @param batch_vector Integer vector. The integers indicate to which
/// batch a given cell belongs.
///
/// @return A list with the following items
/// \itemize{
///   \item pval - The p-values from the ChiSquare test
///   \item chi_square_stats - ChiSquare statistics
///   \item mean_chi_square - The mean ChiSquare value
///   \item median_chi_square - The median ChiSquare value
/// }
///
/// @export
#[extendr]
fn rs_kbet(knn_mat: RMatrix<i32>, batch_vector: Vec<i32>) -> List {
    let n_cells = knn_mat.nrows();
    let k_neighbours = knn_mat.ncols();

    let knn_matrix: Vec<Vec<usize>> = (0..n_cells)
        .map(|i| {
            (0..k_neighbours)
                .map(|j| knn_mat[[i, j]] as usize)
                .collect()
        })
        .collect();

    // Convert batch_vector to Vec<usize>
    let batches: Vec<usize> = batch_vector.r_int_convert();

    let kbet_res = kbet(&knn_matrix, &batches);

    list![
        pval = kbet_res.p_values,
        chi_square_stats = kbet_res.chi_square_stats,
        mean_chi_square = kbet_res.mean_chi_square,
        median_chi_square = kbet_res.median_chi_square
    ]
}

/// Calculate batch silhouette width from an embedding
///
/// @description
/// Computes the average silhouette width on batch labels using pairwise
/// distances in the embedding space. Values near 0 indicate good batch
/// mixing, values near 1 indicate batch separation.
///
/// @param embedding Numeric matrix. The embedding to assess (e.g. PCA or
/// corrected embedding). Rows are cells, columns are dimensions.
/// @param batch_vector Integer vector. The integers indicate to which
/// batch a given cell belongs.
/// @param max_cells Integer or NULL. If not NULL, subsample to this many
/// cells for performance. Defaults to 5000.
/// @param seed Integer. Seed for subsampling reproducibility.
///
/// @return A list with the following items
/// \itemize{
///   \item per_cell - Per-cell silhouette scores
///   \item mean_asw - Mean silhouette width
///   \item median_asw - Median silhouette width
/// }
///
/// @export
#[extendr]
fn rs_batch_silhouette_width(
    embedding: RMatrix<f64>,
    batch_vector: Vec<i32>,
    max_cells: Nullable<i32>,
    seed: i32,
) -> List {
    let embd = r_matrix_to_faer_fp32(&embedding);
    let batches: Vec<usize> = batch_vector.r_int_convert();
    let subsample = match max_cells {
        Nullable::NotNull(n) => Some(n as usize),
        Nullable::Null => None,
    };

    let res = batch_silhouette_width(embd.as_ref(), &batches, subsample, seed as usize);

    list![
        per_cell = res.per_cell,
        mean_asw = res.mean_asw,
        median_asw = res.median_asw
    ]
}

/// Calculate batch LISI scores
///
/// @description
/// Computes the Local Inverse Simpson's Index on batch labels using the
/// kNN graph. Measures the effective number of batches in each cell's
/// neighbourhood. Under perfect mixing LISI equals the number of batches,
/// under no mixing LISI equals 1.
///
/// @param knn_mat Integer matrix. The rows represent the cells and the
/// columns the neighbour indices.
/// @param batch_vector Integer vector. The integers indicate to which
/// batch a given cell belongs.
///
/// @return A list with the following items
/// \itemize{
///   \item per_cell - Per-cell LISI scores
///   \item mean_lisi - Mean LISI
///   \item median_lisi - Median LISI
/// }
///
/// @export
#[extendr]
fn rs_batch_lisi(knn_mat: RMatrix<i32>, batch_vector: Vec<i32>) -> List {
    let n_cells = knn_mat.nrows();
    let k = knn_mat.ncols();

    let knn_indices: Vec<Vec<usize>> = (0..n_cells)
        .map(|i| (0..k).map(|j| knn_mat[[i, j]] as usize).collect())
        .collect();

    let batches: Vec<usize> = batch_vector.r_int_convert();
    let res = batch_lisi(&knn_indices, &batches);

    list![
        per_cell = res.per_cell,
        mean_lisi = res.mean_lisi,
        median_lisi = res.median_lisi
    ]
}

///////////
// BBKNN //
///////////

/// BBKNN implementation in Rust
///
/// @description
/// This function implements the BBKNN algorithm from Polański, et al.
///
/// @param embd Numerical matrix. The embedding matrix to use to generate the
/// BBKNN parameters. Usually PCA. Rows represent cells.
/// @param batch_labels Integer vector. These represent to which batch a given
/// cell belongs.
/// @param bbknn_params List. Contains all of the BBKNN parameters.
/// @param seed Integer. Seed for reproducibility purposes.
/// @param verbose Boolean. Controls verbosity of the function.
///
/// @return A list of two lists representing the sparse matrix representation
/// of the distances and the connectivities.
///
/// @export
///
/// @references Polański, et al., Bioinformatics, 2020
#[extendr]
fn rs_bbknn(
    embd: RMatrix<f64>,
    batch_labels: Vec<i32>,
    bbknn_params: List,
    seed: usize,
    verbose: bool,
) -> List {
    let bbknn_params = BbknnParams::from_r_list(bbknn_params);
    let embd = r_matrix_to_faer_fp32(&embd);
    let batch_labels = batch_labels
        .iter()
        .map(|x| *x as usize)
        .collect::<Vec<usize>>();

    let (distances, connectivities) =
        bbknn(embd.as_ref(), &batch_labels, &bbknn_params, seed, verbose);

    list!(
        distances = sparse_data_to_list(distances),
        connectivities = sparse_data_to_list(connectivities)
    )
}

/// Reduce BBKNN results to Top X neighbours
///
/// @param indptr Integer vector. The index pointers of the underlying data.
/// @param indices Integer vector. The indices of the nearest neighbours.
/// @param data Numeric vector. The distances to the nearest neighbours.
/// @param no_neighbours_to_keep Integer. Number of nearest neighbours to keep.
///
/// @return A list with `indices` (integer matrix) and `dist` (numeric matrix),
/// each with shape (n_cells, no_neighbours_to_keep). Positions without
/// neighbours are filled with -1 (indices) or NaN (distances).
///
/// @export
#[extendr]
fn rs_bbknn_filtering(
    indptr: Vec<i32>,
    indices: Vec<i32>,
    data: Vec<f64>,
    no_neighbours_to_keep: usize,
) -> List {
    let nrow = indptr.len() - 1;
    let ncol = no_neighbours_to_keep;

    let mut idx_mat: Mat<f64> = Mat::from_fn(nrow, ncol, |_, _| f64::NAN);
    let mut dist_mat: Mat<f64> = Mat::from_fn(nrow, ncol, |_, _| f64::NAN);

    for i in 0..nrow {
        let start_i = indptr[i] as usize;
        let end_i = indptr[i + 1] as usize;
        let take = (end_i - start_i).min(no_neighbours_to_keep);

        for j in 0..take {
            idx_mat[(i, j)] = indices[start_i + j] as f64;
            dist_mat[(i, j)] = data[start_i + j];
        }
    }

    list!(
        indices = faer_to_r_matrix(idx_mat.as_ref()),
        dist = faer_to_r_matrix(dist_mat.as_ref())
    )
}

/////////////
// FastMNN //
/////////////

/// FastMNN batch correction in Rust
///
/// @description
/// This function implements the (fast) MNN algorithm from Haghverdi, et al.
/// Instead of working on the full matrix, it uses under the hood PCA and
/// generates a batch-aligned embedding space.
///
/// @param f_path_gene String. Path to the `counts_genes.bin` file.
/// @param cell_indices Integer. The cell indices to use. (0-indexed!)
/// @param gene_indices Integer. The gene indices to use. (0-indexed!) Ideally
/// these are batch-aware highly variable genes.
/// @param batch_indices Integer vector. These represent to which batch a given
/// cell belongs.
/// @param mnn_params List. Contains all of the fastMNN parameters.
/// @param precomputed_pca Optional PCA matrix. If you want to provide a
/// pre-computed matrix.
/// @param seed Integer. Seed for reproducibility purposes.
/// @param verbose Boolean. Controls verbosity of the function.
///
/// @return The batch-corrected embedding space.
///
/// @export
#[extendr]
#[allow(clippy::too_many_arguments)]
fn rs_mnn(
    f_path_gene: &str,
    cell_indices: Vec<i32>,
    gene_indices: Vec<i32>,
    batch_indices: Vec<i32>,
    precomputed_pca: Option<RMatrix<f64>>,
    mnn_params: List,
    verbose: bool,
    seed: usize,
) -> Result<RArray<f64, [usize; 2]>> {
    let cell_indices = cell_indices.r_int_convert();
    let gene_indices = gene_indices.r_int_convert();
    let batch_indices = batch_indices.r_int_convert();
    let mnn_params = FastMnnParams::from_r_list(mnn_params);

    let pre_computed_pca = precomputed_pca.map(|embd| r_matrix_to_faer_fp32(&embd));

    let corrected_embd = fast_mnn_main(
        f_path_gene,
        &cell_indices,
        &gene_indices,
        &batch_indices,
        pre_computed_pca,
        &mnn_params,
        verbose,
        seed,
    )
    .to_extendr()?;

    Ok(faer_to_r_matrix(corrected_embd.as_ref()))
}

/////////////
// Harmony //
/////////////

/// Harmony batch correction in Rust
///
/// @description
/// This function implements the Harmony algorithm from Korsunsky et al., 2019.
///
/// @param pca Numerical matrix, i.e., the PCA matrix you want to correct.
/// @param harmony_params List. The parameters for the Harmony algorithm.
/// @param batch_labels List. Each element in the list needs to be a 0-indexed
/// integer that represents the batch effects you wish to regress out.
/// @param seed Integer. Seed for reproducibility purposes.
/// @param verbose Boolean. Controls verbosity of the function.
///
/// @return The batch-corrected Harmony embedding space.
///
/// @export
#[extendr]
fn rs_harmony(
    pca: RMatrix<f64>,
    harmony_params: List,
    batch_labels: List,
    seed: usize,
    verbose: bool,
) -> extendr_api::Result<RArray<f64, [usize; 2]>> {
    let mut batch_indices: Vec<Vec<usize>> = Vec::new();

    for i in 0..batch_labels.len() {
        let batch_indices_i = batch_labels.elt(i)?;
        let batch_indices_i = batch_indices_i.as_integer_vector().unwrap();
        batch_indices.push(batch_indices_i.r_int_convert());
    }

    let harmony_params = HarmonyParams::from_r_list(harmony_params);

    let embd = r_matrix_to_faer_fp32(&pca);

    let res = harmony(
        embd.as_ref(),
        &batch_indices,
        &harmony_params,
        seed,
        verbose,
    );

    Ok(faer_to_r_matrix(res.as_ref()))
}

/// Harmony batch correction in Rust (version 2)
///
/// @description
/// This function implements the version 2 Harmony algorithm from Patikas, et
/// al., 2026.
///
/// @param pca Numerical matrix, i.e., the PCA matrix you want to correct.
/// @param harmony_params List. The parameters for the Harmony (v2) algorithm.
/// @param batch_labels List. Each element in the list needs to be a 0-indexed
/// integer that represents the batch effects you wish to regress out.
/// @param seed Integer. Seed for reproducibility purposes.
/// @param verbose Boolean. Controls verbosity of the function.
///
/// @return The batch-corrected Harmony (v2) embedding space.
///
/// @export
#[extendr]
fn rs_harmony_v2(
    pca: RMatrix<f64>,
    harmony_params: List,
    batch_labels: List,
    seed: usize,
    verbose: bool,
) -> extendr_api::Result<RArray<f64, [usize; 2]>> {
    let mut batch_indices: Vec<Vec<usize>> = Vec::new();

    for i in 0..batch_labels.len() {
        let batch_indices_i = batch_labels.elt(i)?;
        let batch_indices_i = batch_indices_i.as_integer_vector().unwrap();
        batch_indices.push(batch_indices_i.r_int_convert());
    }

    let harmony_params = HarmonyParamsV2::from_r_list(harmony_params);

    let embd = r_matrix_to_faer_fp32(&pca);

    let res = harmony_v2(
        embd.as_ref(),
        &batch_indices,
        &harmony_params,
        seed,
        verbose,
    );

    Ok(faer_to_r_matrix(res.as_ref()))
}
