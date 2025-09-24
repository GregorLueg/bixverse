use extendr_api::prelude::*;
use faer::Mat;
use rustc_hash::FxHashSet;
use std::time::Instant;

use crate::core::graph::leiden::leiden_clustering;
use crate::single_cell::processing::*;
use crate::single_cell::sc_knn_snn::*;
use crate::utils::r_rust_interface::{faer_to_r_matrix, r_matrix_to_faer_fp32};

//////////////////////////
// Gene set proportions //
//////////////////////////

/// Calculate the percentage of gene sets in the cells
///
/// @description
/// This function allows to calculate for example the proportion of
/// mitochondrial genes, or ribosomal genes in the cells for QC purposes.
///
/// @param f_path_cell String. Path to the `counts_cells.bin` file.
/// @param gene_set_idx Named list with integer(!) positions (1-indexed!) as
/// elements of the genes of interest.
///
/// @return A list with the percentages of counts per gene set group detected
/// in the cells.
///
/// @export
#[extendr]
fn rs_sc_get_gene_set_perc(f_path_cell: &str, gene_set_idx: List) -> extendr_api::Result<List> {
    let mut gene_set_indices = Vec::with_capacity(gene_set_idx.len());

    for i in 0..gene_set_idx.len() {
        let element = gene_set_idx.elt(i).unwrap();
        let indices_i: Vec<u16> = element
            .as_integer_vector()
            .unwrap()
            .iter()
            .map(|x| *x as u16)
            .collect();
        gene_set_indices.push(indices_i);
    }

    let res = get_gene_set_perc(f_path_cell, gene_set_indices);

    let mut result_list = List::new(gene_set_idx.len());
    if let Some(names) = gene_set_idx.names() {
        result_list.set_names(names).unwrap();
    }

    #[allow(clippy::needless_range_loop)]
    for i in 0..result_list.len() {
        let res_i = &res[i];
        result_list.set_elt(i, Robj::from(res_i)).unwrap();
    }

    Ok(result_list)
}

///////////////////////////
// Highly variable genes //
///////////////////////////

/// Enum for the different methods
enum HvgMethod {
    /// Variance stabilising transformation
    Vst,
    /// Binned version by average expression
    MeanVarBin,
    /// Simple dispersion
    Dispersion,
}

/// Helper function to parse the HVG
///
/// ### Params
///
/// * `s` - Type of HVG calculation to do
///
/// ### Returns
///
/// Option of the HvgMethod (some not yet implemented)
fn get_hvg_method(s: &str) -> Option<HvgMethod> {
    match s.to_lowercase().as_str() {
        "vst" => Some(HvgMethod::Vst),
        "meanvarbin" => Some(HvgMethod::MeanVarBin),
        "dispersion" => Some(HvgMethod::Dispersion),
        _ => None,
    }
}

/// Calculate the percentage of gene sets in the cells
///
/// @description
/// This function allows to calculate for example the proportion of
/// mitochondrial genes, or ribosomal genes in the cells for QC purposes.
///
/// @param f_path_gene String. Path to the `counts_genes.bin` file.
/// @param hvg_method String. Which HVG detection method to use. Options
/// are `c("vst", "meanvarbin", "dispersion")`. So far, only the first is
/// implemented.
/// @param cell_indices Integer positions (0-indexed!) that defines the cells
/// to keep.
/// @param loess_span Numeric. The span parameter for the loess function.
/// @param clip_max Optional clipping number. Defaults to `sqrt(no_cells)` if
/// not provided.
/// @param verbose Boolean. Controls verbosity of the function.
///
/// @return A list with the percentages of counts per gene set group detected
/// in the cells:
/// \itemize{
///   \item mean - The average expression of the gene.
///   \item var - The variance of the gene.
///   \item var_exp - The expected variance of the gene.
///   \item var_std - The standardised variance of the gene.
/// }
///
/// @export
#[extendr]
fn rs_sc_hvg(
    f_path_gene: &str,
    hvg_method: &str,
    cell_indices: Vec<i32>,
    loess_span: f64,
    clip_max: Option<f32>,
    verbose: bool,
) -> List {
    let cell_set: FxHashSet<u32> = cell_indices.iter().map(|x| *x as u32).collect();
    let hvg_type = get_hvg_method(hvg_method)
        .ok_or_else(|| format!("Invalid HVG method: {}", hvg_method))
        .unwrap();

    let hvg_res: HvgRes = match hvg_type {
        HvgMethod::Vst => get_hvg_vst(f_path_gene, &cell_set, loess_span, clip_max, verbose),
        HvgMethod::MeanVarBin => get_hvg_mvb(),
        HvgMethod::Dispersion => get_hvg_dispersion(),
    };

    list!(
        mean = hvg_res.mean,
        var = hvg_res.var,
        var_exp = hvg_res.var_exp,
        var_std = hvg_res.var_std
    )
}

/////////
// PCA //
/////////

/// Calculates PCA for single cell
///
/// @description
/// Helper function that will calculate the PCA for the specified highly
/// variable genes. Has the option to use randomised SVD for faster solving
/// of the PCA.
///
/// @param f_path_gene String. Path to the `counts_genes.bin` file.
/// @param no_pcs Integer. Number of PCs to calculate.
/// @param random_svd Boolean. Shall randomised SVD be used.
/// @param cell_indices Integer. The cell indices to use. (0-indexed!)
/// @param gene_indices Integer. The gene indices to use. (0-indexed!)
/// @param seed Integer. Random seed for the randomised SVD.
/// @param verbose Boolean. Controls verbosity of the function.
///
/// @returns A list with with the following items
/// \itemize{
///   \item scores - The samples projected on the PCA space.
///   \item loadings - The loadings of the features for the PCA.
/// }
///
/// @export
#[extendr]
fn rs_sc_pca(
    f_path_gene: &str,
    no_pcs: usize,
    random_svd: bool,
    cell_indices: Vec<i32>,
    gene_indices: Vec<i32>,
    seed: usize,
    verbose: bool,
) -> List {
    let cell_set: FxHashSet<u32> = cell_indices.iter().map(|x| *x as u32).collect();
    let gene_indices = gene_indices.iter().map(|x| *x as usize).collect();

    let res = pca_on_sc(
        f_path_gene,
        &cell_set,
        gene_indices,
        no_pcs,
        random_svd,
        seed,
        verbose,
    );

    let scores = Mat::from_fn(res.0.nrows(), res.0.ncols(), |i, j| {
        let val = res.0.get(i, j);
        *val as f64
    });
    let loadings = Mat::from_fn(res.1.nrows(), res.1.ncols(), |i, j| {
        let val = res.1.get(i, j);
        *val as f64
    });

    list!(
        scores = faer_to_r_matrix(scores.as_ref()),
        loadings = faer_to_r_matrix(loadings.as_ref())
    )
}

///////////////
// kNN / sNN //
///////////////

/// Enum for the different methods
enum KnnSearch {
    /// Annoy-based
    Annoy,
    /// Hierarchical Navigable Small World
    Hnsw,
}

/// Helper function to get the KNN method
///
/// ### Params
///
/// * `s` - Type of HVG calculation to do
///
/// ### Returns
///
/// Option of the HvgMethod (some not yet implemented)
fn get_knn_method(s: &str) -> Option<KnnSearch> {
    match s.to_lowercase().as_str() {
        "annoy" => Some(KnnSearch::Annoy),
        "hnsw" => Some(KnnSearch::Hnsw),
        _ => None,
    }
}

#[extendr]
fn rs_sc_knn(
    embd: RMatrix<f64>,
    no_neighbours: usize,
    seed: usize,
    n_trees: usize,
    search_budget: usize,
    verbose: bool,
    algorithm_type: String,
) -> extendr_api::RArray<i32, [usize; 2]> {
    let embd = r_matrix_to_faer_fp32(&embd);

    let start_knn = Instant::now();

    let knn_method = get_knn_method(&algorithm_type)
        .ok_or_else(|| format!("Invalid KNN search method: {}", algorithm_type))
        .unwrap();

    let knn = match knn_method {
        KnnSearch::Hnsw => generate_knn_hnsw(embd.as_ref(), no_neighbours, seed),
        KnnSearch::Annoy => {
            generate_knn_annoy(embd.as_ref(), no_neighbours, n_trees, search_budget, seed)
        }
    };

    let end_knn = start_knn.elapsed();

    if verbose {
        println!("KNN generation : {:.2?}", end_knn);
    }

    let index_mat = Mat::from_fn(embd.nrows(), no_neighbours, |i, j| knn[i][j] as i32);

    faer_to_r_matrix(index_mat.as_ref())
}

////////////
// Leiden //
////////////

#[extendr]
fn rs_leiden_clustering(
    from: Vec<i32>,
    to: Vec<i32>,
    weights: Vec<f64>,
    max_iterations: usize,
    res: f64,
    seed: Option<u64>,
) -> Vec<i32> {
    let from = from.iter().map(|x| *x as usize).collect::<Vec<usize>>();
    let to = to.iter().map(|x| *x as usize).collect::<Vec<usize>>();
    let weights = weights.iter().map(|x| *x as f32).collect::<Vec<f32>>();

    let leiden_communities = leiden_clustering(from, to, weights, max_iterations, res as f32, seed);

    leiden_communities.iter().map(|x| *x as i32).collect()
}

extendr_module! {
    mod r_preprocessing;
    fn rs_sc_get_gene_set_perc;
    fn rs_sc_hvg;
    fn rs_sc_pca;
    fn rs_sc_knn;
    fn rs_leiden_clustering;
}
