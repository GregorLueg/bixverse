use extendr_api::prelude::*;
use faer::Mat;
use std::time::Instant;

use crate::single_cell::methods::doublet_detection::*;
use crate::single_cell::methods::scrublet::*;
use crate::single_cell::processing::*;
use crate::single_cell::sc_knn_snn::*;
use crate::utils::r_rust_interface::{faer_to_r_matrix, r_matrix_to_faer_fp32};
use crate::utils::traits::*;

extendr_module! {
    mod r_sc_processing;
    fn rs_sc_scrublet;
    fn rs_sc_doublet_detection;
    fn rs_sc_get_top_genes_perc;
    fn rs_sc_get_gene_set_perc;
    fn rs_sc_hvg;
    fn rs_sc_hvg_batch_aware;
    fn rs_sc_pca;
    fn rs_sc_knn;
    fn rs_sc_snn;
}

///////////////////////
// Doublet detection //
///////////////////////

/// Scrublet Rust interface
///
/// @param f_path_gene String. Path to the `counts_genes.bin` file.
/// @param f_path_cell String. Path to the `counts_cells.bin` file.
/// @param cells_to_keep Integer vector. The indices (0-indexed!) of the cells
/// to include in this analysis.
/// @param scrublet_params List. Parameter list, see
/// [bixverse::params_scrublet()].
/// @param seed Integer. Seed for reproducibility purposes.
/// @param verbose Boolean. Controls verbosity
/// @param streaming Boolean. Shall the data be streamed for the HVG
/// calculations.
/// @param return_combined_pca Boolean. Shall the generated PCA be returned.
/// @param return_pairs Boolean. Shall the parents of the simulated cells
/// be returned.
///
/// @returns A list with
/// \itemize{
///  \item predicted_doublets - Boolean vector indicating which observed cells
///  predicted as doublets (TRUE = doublet, FALSE = singlet).
///  \item doublet_scores_obs - Numerical vector with the likelihood of being
///  a doublet for the observed cells.
///  \item doublet_scores_sim - Numerical vector with the likelihood of being
///  a doublet for the simulated cells.
///  \item doublet_errors_obs - Numerical vector with the standard errors of
///  the scores for the observed cells.
///  \item z_scores - Z-scores for the observed cells. Represents:
///  `score - threshold / error`.
///  \item threshold - Used threshold.
///  \item detected_doublet_rate - Fraction of cells that are called as
///  doublet.
///  \item detectable_doublet_fraction - Fraction of simulated doublets with
///  scores above the threshold.
///  \item overall_doublet_rate - Estimated overall doublet rate. Should roughly
///  match the expected doublet rate.
///  \item pca - Optional PCA embeddings across the original cells and simulated
///  doublets.
///  \item pair_1 - Optional integer vector representing the first parent of the
///  simulated doublets.
///  \item pair_2 -  Optional integer vector representing the second parent of
///  the simulated doublets.
/// }
///
/// @export
#[extendr]
#[allow(clippy::too_many_arguments)]
fn rs_sc_scrublet(
    f_path_gene: &str,
    f_path_cell: &str,
    cells_to_keep: Vec<i32>,
    scrublet_params: List,
    seed: usize,
    verbose: bool,
    streaming: bool,
    return_combined_pca: bool,
    return_pairs: bool,
) -> List {
    let scrublet_params = ScrubletParams::from_r_list(scrublet_params);
    let cells_to_keep = cells_to_keep.r_int_convert();
    let mut scrublet = Scrublet::new(f_path_gene, f_path_cell, scrublet_params, &cells_to_keep);
    let (scrublet_res, pca, pair_1, pair_2): FinalScrubletRes =
        scrublet.run_scrublet(streaming, seed, verbose, return_combined_pca, return_pairs);

    let pca_out = pca.map(|m| faer_to_r_matrix(m.as_ref()));
    let pair_1_out: Robj = match pair_1 {
        Some(p) => p.r_int_convert().into(),
        None => NULL.into(),
    };
    let pair_2_out: Robj = match pair_2 {
        Some(p) => p.r_int_convert().into(),
        None => NULL.into(),
    };

    list!(
        predicted_doublets = scrublet_res.predicted_doublets,
        doublet_scores_obs = scrublet_res.doublet_scores_obs.r_float_convert(),
        doublet_scores_sim = scrublet_res.doublet_scores_sim.r_float_convert(),
        doublet_errors_obs = scrublet_res.doublet_errors_obs.r_float_convert(),
        z_scores = scrublet_res.z_scores.r_float_convert(),
        threshold = scrublet_res.threshold as f64,
        detected_doublet_rate = scrublet_res.detected_doublet_rate as f64,
        detectable_doublet_fraction = scrublet_res.detectable_doublet_fraction as f64,
        overall_doublet_rate = scrublet_res.overall_doublet_rate as f64,
        pca = pca_out,
        pair_1 = pair_1_out,
        pair_2 = pair_2_out
    )
}

/// Detect Doublets via BoostClassifier (in Rust)
///
/// @param f_path_gene String. Path to the `counts_genes.bin` file.
/// @param f_path_cell String. Path to the `counts_cells.bin` file.
/// @param cells_to_keep Integer vector. The indices (0-indexed!) of the cells
/// to include in this analysis.
/// @param boost_params List. Parameter list, see
/// [bixverse::params_boost()].
/// @param seed Integer. Seed for reproducibility purposes.
/// @param verbose Boolean. Controls verbosity
/// @param streaming Boolean. Shall the data be streamed for the HVG
/// calculations.
///
/// @returns A list with
/// \itemize{
///  \item predicted_doublets - Boolean vector indicating which observed cells
///  predicted as doublets (TRUE = doublet, FALSE = singlet).
///  \item doublet_scores_obs - Numerical vector with the likelihood of being
///  a doublet for the observed cells.
///  \item voting_avg - Voting average across the different iterations.
/// }
///
/// @export
#[extendr]
fn rs_sc_doublet_detection(
    f_path_gene: &str,
    f_path_cell: &str,
    cells_to_keep: Vec<i32>,
    boost_params: List,
    seed: usize,
    streaming: bool,
    verbose: bool,
) -> List {
    let boost_params = BoostParams::from_r_list(boost_params);
    let cells_to_keep = cells_to_keep.r_int_convert();

    let mut boost_classifier =
        BoostClassifier::new(f_path_gene, f_path_cell, boost_params, &cells_to_keep);

    let boost_res: BoostResult = boost_classifier.run_boost(streaming, seed, verbose);

    list!(
        doublet = boost_res.predicted_doublets,
        doublet_score = boost_res.doublet_scores,
        voting_avg = boost_res.voting_average
    )
}

//////////////////////
// Cumulative genes //
//////////////////////

/// Calculates the cumulative proportion of the top X genes
///
/// @description
/// This calculates the cumulative proportion of the top X genes, for example
/// Top10, 50, 100. High values here indicate low complexity samples, i.e.,
/// bad quality.
///
/// @param f_path_cell String. Path to the `counts_cells.bin` file.
/// @param top_n_vals Integer. The different Top X to look for.
/// @param cell_indices Integer. The indices of the cells for which to calculate
/// the proportions. (0-indexed!)
/// @param streaming Boolean. Shall the data be worked on in chunks.
/// @param verbose Boolean. Controls verbosity of the function.
///
/// @return A list with the cumulative percentages of the Top X genes defined
/// as in `top_n_vals`.
///
/// @export
#[extendr]
fn rs_sc_get_top_genes_perc(
    f_path_cell: &str,
    top_n_vals: &[i32],
    cell_indices: &[i32],
    streaming: bool,
    verbose: bool,
) -> List {
    let cell_indices = cell_indices.r_int_convert();
    let top_n_vals = top_n_vals.r_int_convert();

    let res = if streaming {
        get_top_genes_perc_streaming(f_path_cell, &top_n_vals, &cell_indices, verbose)
    } else {
        get_top_genes_perc(f_path_cell, &top_n_vals, &cell_indices, verbose)
    };

    let mut result_list = List::new(top_n_vals.len());

    for i in 0..result_list.len() {
        let res_i = &res[i];
        result_list.set_elt(i, Robj::from(res_i)).unwrap();
    }

    result_list
}

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
/// @param gene_set_idx Named list with integer(!) positions (0-indexed!) as
/// elements of the genes of interest.
/// @param cell_indices Integer. The indices of the cells for which to calculate
/// the proportions. (0-indexed!)
/// @param streaming Boolean. Shall the data be worked on in chunks.
/// @param verbose Boolean. Controls verbosity of the function.
///
/// @return A list with the percentages of counts per gene set group detected
/// in the cells.
///
/// @export
#[extendr]
fn rs_sc_get_gene_set_perc(
    f_path_cell: &str,
    gene_set_idx: List,
    cell_indices: Vec<i32>,
    streaming: bool,
    verbose: bool,
) -> extendr_api::Result<List> {
    let mut gene_set_indices = Vec::with_capacity(gene_set_idx.len());

    let cell_indices = cell_indices.r_int_convert();

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

    let res = if streaming {
        get_gene_set_perc_streaming(f_path_cell, gene_set_indices, &cell_indices, verbose)
    } else {
        get_gene_set_perc(f_path_cell, gene_set_indices, &cell_indices, verbose)
    };

    let mut result_list = List::new(gene_set_idx.len());
    if let Some(names) = gene_set_idx.names() {
        result_list.set_names(names).unwrap();
    }

    for i in 0..result_list.len() {
        let res_i = &res[i];
        result_list.set_elt(i, Robj::from(res_i)).unwrap();
    }

    Ok(result_list)
}

///////////////////////////
// Highly variable genes //
///////////////////////////

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
/// @param streaming Boolean. Shall the genes be streamed in to reduce memory
/// pressure.
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
    streaming: bool,
    verbose: bool,
) -> List {
    let cell_set = cell_indices.r_int_convert();

    let hvg_type = get_hvg_method(hvg_method)
        .ok_or_else(|| format!("Invalid HVG method: {}", hvg_method))
        .unwrap();

    let hvg_res: HvgRes = if streaming {
        match hvg_type {
            HvgMethod::Vst => {
                get_hvg_vst_streaming(f_path_gene, &cell_set, loess_span, clip_max, verbose)
            }
            HvgMethod::MeanVarBin => get_hvg_mvb_streaming(),
            HvgMethod::Dispersion => get_hvg_dispersion_streaming(),
        }
    } else {
        match hvg_type {
            HvgMethod::Vst => get_hvg_vst(f_path_gene, &cell_set, loess_span, clip_max, verbose),
            HvgMethod::MeanVarBin => get_hvg_mvb(),
            HvgMethod::Dispersion => get_hvg_dispersion(),
        }
    };

    list!(
        mean = hvg_res.mean,
        var = hvg_res.var,
        var_exp = hvg_res.var_exp,
        var_std = hvg_res.var_std
    )
}

/// Calculate HVG per batch
///
/// @description
/// Batch-aware highly variable gene detection. Calculates HVG statistics
/// separately for each batch, allowing for downstream selection strategies
/// such as union of top genes per batch.
///
/// @param f_path_gene String. Path to the `counts_genes.bin` file.
/// @param hvg_method String. Which HVG detection method to use. Currently
/// only `"vst"` is implemented for batch-aware mode.
/// @param cell_indices Integer positions (0-indexed!) that defines the cells
/// to keep.
/// @param batch_labels Integer vector (0-indexed!) defining batch membership
/// for each cell. Must be same length as `cell_indices`.
/// @param loess_span Numeric. The span parameter for the loess function.
/// @param clip_max Optional clipping number. Defaults to `sqrt(no_cells)` per
/// batch if not provided.
/// @param streaming Boolean. Shall the genes be streamed in to reduce memory
/// pressure.
/// @param verbose Boolean. Controls verbosity of the function.
///
/// @return A list with HVG statistics concatenated across all batches:
/// \itemize{
/// \item mean - The average expression of each gene in each batch.
/// \item var - The variance of each gene in each batch.
/// \item var_exp - The expected variance of each gene in each batch.
/// \item var_std - The standardised variance of each gene in each batch.
/// \item batch - Batch index for each gene (length = n_genes * n_batches).
/// \item gene_idx - Gene index for each entry (0-indexed, length = n_genes * n_batches).
/// }
///
/// @export
#[allow(clippy::too_many_arguments)]
#[extendr]
fn rs_sc_hvg_batch_aware(
    f_path_gene: &str,
    hvg_method: &str,
    cell_indices: Vec<i32>,
    batch_labels: Vec<i32>,
    loess_span: f64,
    clip_max: Option<f32>,
    streaming: bool,
    verbose: bool,
) -> List {
    let cell_set = cell_indices.r_int_convert();
    let batch_set = batch_labels.r_int_convert();

    let hvg_type = get_hvg_method(hvg_method)
        .ok_or_else(|| format!("Invalid HVG method: {}", hvg_method))
        .unwrap();

    let hvg_results: Vec<HvgRes> = if streaming {
        match hvg_type {
            HvgMethod::Vst => get_hvg_vst_batch_aware_streaming(
                f_path_gene,
                &cell_set,
                &batch_set,
                loess_span,
                clip_max,
                verbose,
            ),
            HvgMethod::MeanVarBin => panic!("MeanVarBin not implemented for batch-aware mode"),
            HvgMethod::Dispersion => panic!("Dispersion not implemented for batch-aware mode"),
        }
    } else {
        match hvg_type {
            HvgMethod::Vst => get_hvg_vst_batch_aware(
                f_path_gene,
                &cell_set,
                &batch_set,
                loess_span,
                clip_max,
                verbose,
            ),
            HvgMethod::MeanVarBin => panic!("MeanVarBin not implemented for batch-aware mode"),
            HvgMethod::Dispersion => panic!("Dispersion not implemented for batch-aware mode"),
        }
    };

    // Flatten results for easy R manipulation
    let n_genes = hvg_results[0].mean.len();
    let n_batches = hvg_results.len();
    let total_len = n_genes * n_batches;

    let mut mean_flat = Vec::with_capacity(total_len);
    let mut var_flat = Vec::with_capacity(total_len);
    let mut var_exp_flat = Vec::with_capacity(total_len);
    let mut var_std_flat = Vec::with_capacity(total_len);
    let mut batch_idx = Vec::with_capacity(total_len);
    let mut gene_idx = Vec::with_capacity(total_len);

    for (batch, hvg_res) in hvg_results.into_iter().enumerate() {
        mean_flat.extend(hvg_res.mean);
        var_flat.extend(hvg_res.var);
        var_exp_flat.extend(hvg_res.var_exp);
        var_std_flat.extend(hvg_res.var_std);
        batch_idx.extend(vec![batch as i32; n_genes]);
        gene_idx.extend(0..n_genes as i32);
    }

    list!(
        mean = mean_flat,
        var = var_flat,
        var_exp = var_exp_flat,
        var_std = var_std_flat,
        batch = batch_idx,
        gene_idx = gene_idx
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
/// @param return_scaled Boolean. Shall the scaled data be returned.
/// @param verbose Boolean. Controls verbosity of the function.
///
/// @returns A list with with the following items
/// \itemize{
///   \item scores - The samples projected on the PCA space.
///   \item loadings - The loadings of the features for the PCA.
///   \item scaled - The scaled matrix if you set return_scaled to `TRUE`.
/// }
///
/// @export
#[allow(clippy::too_many_arguments)]
#[extendr]
fn rs_sc_pca(
    f_path_gene: &str,
    no_pcs: usize,
    random_svd: bool,
    cell_indices: Vec<i32>,
    gene_indices: Vec<i32>,
    seed: usize,
    return_scaled: bool,
    verbose: bool,
) -> List {
    let cell_set = cell_indices
        .iter()
        .map(|x| *x as usize)
        .collect::<Vec<usize>>();
    let gene_indices = gene_indices
        .iter()
        .map(|x| *x as usize)
        .collect::<Vec<usize>>();

    let res = pca_on_sc(
        f_path_gene,
        &cell_set,
        &gene_indices,
        no_pcs,
        random_svd,
        seed,
        return_scaled,
        verbose,
    );

    let singular_values_f64: Vec<f64> = res.2.iter().map(|&x| x as f64).collect();
    let scaled = res.3.map(|s| faer_to_r_matrix(s.as_ref()));

    list!(
        scores = faer_to_r_matrix(res.0.as_ref()),
        loadings = faer_to_r_matrix(res.1.as_ref()),
        singular_values = singular_values_f64,
        scaled = scaled
    )
}

///////////////
// kNN / sNN //
///////////////

/// Generates the kNN graph
///
/// @description
/// This function is a wrapper over the Rust-based generation of the approximate
/// nearest neighbours. You have several options to get the approximate nearest
/// neighbours:
///
/// - `"annoy"`: leverages binary trees to generate rapidly in a parallel manner
///   an index. Good compromise of index generation, querying speed.
/// - `"hnsw"`: uses a hierarchical navigatable small worlds index under the
///   hood. The index generation takes more long, but higher recall and ideal
///   for very large datasets due to subdued memory pressure.
/// - `"nndescent"`: an index-free approximate nearest neighbour algorithm
///   that is ideal for small, ephemeral kNN graphs.
///
/// @param embd Numerical matrix. The embedding matrix to use to generate the
/// kNN graph.
/// @param knn_params List. The kNN parameters defined by
/// [bixverse::params_sc_neighbours()].
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
    knn_params: List,
    verbose: bool,
    seed: usize,
) -> extendr_api::Result<extendr_api::RArray<i32, [usize; 2]>> {
    let embd = r_matrix_to_faer_fp32(&embd);

    let knn_params = KnnParams::from_r_list(knn_params);

    let start_knn = Instant::now();

    let knn_method = parse_knn_method(&knn_params.knn_method)
        .ok_or_else(|| format!("Invalid KNN search method: {}", knn_params.knn_method))?;

    let knn = match knn_method {
        KnnSearch::Hnsw => generate_knn_hnsw(
            embd.as_ref(),
            &knn_params.ann_dist,
            knn_params.k,
            seed,
            verbose,
        ),
        KnnSearch::Annoy => generate_knn_annoy(
            embd.as_ref(),
            &knn_params.ann_dist,
            knn_params.k,
            knn_params.n_tree,
            knn_params.search_budget,
            seed,
            verbose,
        ),
        KnnSearch::NNDescent => generate_knn_nndescent(
            embd.as_ref(),
            &knn_params.ann_dist,
            knn_params.k,
            knn_params.max_iter,
            knn_params.delta,
            knn_params.rho,
            seed,
            verbose,
        ),
    };

    let end_knn = start_knn.elapsed();

    if verbose {
        println!("KNN generation done : {:.2?}", end_knn);
    }

    let index_mat = Mat::from_fn(embd.nrows(), knn_params.k, |i, j| knn[i][j] as i32);

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
/// @param limited_graph Boolean. Shall the sNNs only be calculated between
/// direct neighbours in the graph, or between all possible combinations.
/// @param pruning Float. Below which value for the Jaccard similarity to prune
/// the weight to 0.
/// @param verbose Boolean. Controls verbosity of the function.
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
    verbose: bool,
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
            verbose,
        )
    } else {
        generate_snn_full(
            &data,
            n_neighbours,
            knn_mat.nrows(),
            pruning as f32,
            snn_method,
            verbose,
        )
    };

    Ok(list!(
        edges = snn_data.0.iter().map(|x| *x as i32).collect::<Vec<i32>>(),
        weights = snn_data.1
    ))
}
