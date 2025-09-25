use extendr_api::prelude::*;
use faer::Mat;
use rustc_hash::FxHashSet;
use std::time::Instant;

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

/// Generates the kNN graph
///
/// @description
/// This function is a wrapper over the Rust-based generation of the approximate
/// nearest neighbours. You have two options to generate the kNNs. `"annoy"` or
/// `"hnsw"`.
///
/// @param embd Numerical matrix. The embedding matrix to use to generate the
/// kNN graph.
/// @param no_neighbours Integer. Number of neighbours to return
/// @param n_trees Integer. Number of trees to use for the `"annoy"` algorithm.
/// @param search_budget Integer. Search budget per tree for the `"annoy"`
/// algorithm.
/// @param verbose Boolean. Controls verbosity of the function and returns
/// how long certain operations took.
/// @param seed Integer. Seed for reproducibility purposes.
///
/// @return A integer matrix of N x k with N being the number of cells and k the
/// number of neighbours.
///
/// @export
#[extendr]
fn rs_sc_knn(
    embd: RMatrix<f64>,
    no_neighbours: usize,
    n_trees: usize,
    search_budget: usize,
    algorithm_type: String,
    verbose: bool,
    seed: usize,
) -> extendr_api::Result<extendr_api::RArray<i32, [usize; 2]>> {
    let embd = r_matrix_to_faer_fp32(&embd);

    let start_knn = Instant::now();

    let knn_method = get_knn_method(&algorithm_type)
        .ok_or_else(|| format!("Invalid KNN search method: {}", algorithm_type))?;

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

    Ok(faer_to_r_matrix(index_mat.as_ref()))
}

/// Generates the sNN graph for igraph
///
/// @description
/// This function takes a kNN matrix and generates the inputs for an SNN
/// graph based on it.
///
/// @param knn_mat Integer matrix. Rows represent cells and the columns
/// represent the neighbours.
/// @param snn_method String. Which method to use to calculate the similarity.
/// Choice of `c("jaccard", "rank")`.
/// @param pruning Float. Below which value for the Jaccard similarity to prune
/// the weight to 0.
///
/// @return A integer matrix of N x k with N being the number of cells and k the
/// number of neighbours.
///
/// @export
#[extendr]
fn rs_sc_snn(
    knn_mat: RMatrix<i32>,
    snn_method: String,
    limited_graph: bool,
    pruning: f64,
) -> extendr_api::Result<List> {
    let n_neighbours = knn_mat.ncols();
    let data = knn_mat
        .data()
        .iter()
        .map(|x| *x as usize)
        .collect::<Vec<usize>>();

    let snn_method = get_snn_similiarity_method(&snn_method)
        .ok_or_else(|| format!("Invalid SNN similarity method: {}", snn_method))?;

    let snn_data = if limited_graph {
        generate_snn_limited(
            &data,
            n_neighbours,
            knn_mat.nrows(),
            pruning as f32,
            snn_method,
        )
    } else {
        generate_snn_full(
            &data,
            n_neighbours,
            knn_mat.nrows(),
            pruning as f32,
            snn_method,
        )
    };

    Ok(list!(
        edges = snn_data.0.iter().map(|x| *x as i32).collect::<Vec<i32>>(),
        weights = snn_data.1
    ))
}

extendr_module! {
    mod r_preprocessing;
    fn rs_sc_get_gene_set_perc;
    fn rs_sc_hvg;
    fn rs_sc_pca;
    fn rs_sc_knn;
    fn rs_sc_snn;
}
