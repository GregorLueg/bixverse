use extendr_api::*;
use faer::Mat;

use crate::core::data::sparse_io::ParallelSparseReader;
use crate::core::data::sparse_structures::CompressedSparseData;
use crate::single_cell::cell_aggregations::*;
use crate::single_cell::dge_aucs::*;
use crate::single_cell::sc_knn_snn::*;
use crate::single_cell::seacells::*;
use crate::utils::r_rust_interface::{faer_to_r_matrix, r_matrix_to_faer_fp32};
use crate::utils::traits::*;

//////////
// DGEs //
//////////

/// Calculate DGEs between cells based on Mann Whitney stats
///
/// @description
/// The function will take two sets of cell indices and calculate the
/// differential gene expression based on the Mann Whitney test between the
/// two groups.
///
/// @param f_path String. Path to the `counts_cells.bin` file.
/// @param cell_indices_1 Integer. Index positions (0-indexed) of the cells
/// of group 1.
/// @param cell_indices_2 Integer. Index positions (0-indexed) of the cells
/// of group 2.
/// @param min_prop Minimum proportion of expression in at least one of the
/// two groups to be tested.
/// @param alternative String. One of `c("twosided", "greater", "less")`. Null
/// hypothesis.
/// @param verbose Boolean. Controls verbosity of the function.
///
/// @return A list with the following elements
/// \itemize{
///   \item lfc - Log fold changes between the two groups.
///   \item prop1 - Proportion of cells expressing the gene in group 1.
///   \item prop2 - Proportion of cells expressing the gene in group 2.
///   \item z_scores - Z-scores based on the Mann Whitney statistic.
///   \item p_values - P-values of the Mann Whitney statistic.
///   \item fdr - False discovery rate after BH adjustment
///   \item genes_to_keep - Boolean indicating which genes were tested.
/// }
///
/// @export
#[extendr]
fn rs_calculate_dge_mann_whitney(
    f_path: String,
    cell_indices_1: &[i32],
    cell_indices_2: &[i32],
    min_prop: f64,
    alternative: String,
    verbose: bool,
) -> extendr_api::Result<List> {
    let cell_indices_1 = cell_indices_1
        .iter()
        .map(|x| *x as usize)
        .collect::<Vec<usize>>();
    let cell_indices_2 = cell_indices_2
        .iter()
        .map(|x| *x as usize)
        .collect::<Vec<usize>>();

    let dge_results: DgeMannWhitneyRes = calculate_dge_grps_mann_whitney(
        &f_path,
        &cell_indices_1,
        &cell_indices_2,
        min_prop as f32,
        &alternative,
        verbose,
    )?;

    Ok(list!(
        lfc = dge_results
            .lfc
            .iter()
            .map(|x| *x as f64)
            .collect::<Vec<f64>>(),
        prop1 = dge_results
            .prop1
            .iter()
            .map(|x| *x as f64)
            .collect::<Vec<f64>>(),
        prop2 = dge_results
            .prop2
            .iter()
            .map(|x| *x as f64)
            .collect::<Vec<f64>>(),
        z_scores = dge_results.z_scores,
        p_values = dge_results.p_vals,
        fdr = dge_results.fdr,
        genes_to_keep = dge_results.genes_to_keep
    ))
}

////////////////////////
// Pathway activities //
////////////////////////

/// Calculate AUCell in Rust
///
/// @description
/// The function will take in a list of gene set indices (0-indexed!) and
/// calculate an AUCell type statistic. Two options here: calculate this
/// with proper AUROC calculations (useful for marker gene expression) or
/// based on the Mann-Whitney statistic (useful for pathway activity
/// measurs). Data can be streamed in chunks of 50k cells per or loaded in
/// in one go.
///
/// @param f_path String. Path to the `counts_cells.bin` file.
/// @param gs_list List. List with the gene set indices (0-indexed!) of the
/// genes of interest.
/// @param cells_to_keep Integer. Vector of indices of the cells to keep.
/// @param auc_type String. One of `"wilcox"` or `"auroc"`, pending on
/// which statistic you wish to calculate.
/// @param streaming Boolean. Shall the data be streamed.
/// @param verbose Boolean. Controls verbosity of the function.
///
/// @return A matrix of gene set AUCs x cells.
///
/// @export
#[extendr]
fn rs_aucell(
    f_path: String,
    gs_list: List,
    cells_to_keep: Vec<i32>,
    auc_type: &str,
    streaming: bool,
    verbose: bool,
) -> extendr_api::Result<extendr_api::RArray<f64, [usize; 2]>> {
    let mut gs_indices: Vec<Vec<usize>> = Vec::with_capacity(gs_list.len());

    let cells_to_keep = cells_to_keep
        .iter()
        .map(|x| *x as usize)
        .collect::<Vec<usize>>();

    for i in 0..gs_list.len() {
        let r_obj = gs_list.elt(i).unwrap();
        let int = r_obj
            .as_integer_vector()
            .unwrap()
            .iter()
            .map(|x| *x as usize)
            .collect::<Vec<usize>>();
        gs_indices.push(int);
    }

    let res = if streaming {
        calculate_aucell_streaming(&f_path, &gs_indices, &cells_to_keep, auc_type, verbose)?
    } else {
        calculate_aucell(&f_path, &gs_indices, &cells_to_keep, auc_type, verbose)?
    };

    let auc_mat = Mat::from_fn(res.len(), res[0].len(), |i, j| res[i][j] as f64);

    Ok(faer_to_r_matrix(auc_mat.as_ref()))
}

////////////////
// Meta cells //
////////////////

/// Generate meta cells
///
/// @description This function implements the approach from Morabito, et al.
/// to generate meta cells. You can provide an already pre-computed kNN matrix
/// or an embedding to regenerate the kNN matrix with specified parameters in
/// the meta_cell_params. If `knn_mat` is provided, this one will be used. You
/// need to at least provide `knn_mat` or `embd`!
///
/// @param f_path String. Path to the `counts_cells.bin` file.
/// @param knn_mat Optional integer matrix. The kNN matrix you wish to use
/// for the generation of the meta cells. This function expects 0-indices!
/// @param embd Optional numerical matrix. The embedding matrix (for example
/// PCA embedding) you wish to use for the generation of the kNN graph that
/// is used subsequently for aggregation of the meta cells.
/// @param meta_cell_params A list containing the meta cell parameters.
/// @param target_size Numeric. Target library size for re-normalisation of
/// the meta cells. Typically `1e4`.
/// @param seed Integer. For reproducibility purposes.
/// @param verbose Boolean. Controls verbosity of the function.
/// @param return_aggregated Boolean. If TRUE, aggregates counts into meta cells.
///
/// @returns A list with the following elements:
/// \itemize{
///  \item assignments - A list containing assignment information with elements:
///    assignments (vector), metacells (list), unassigned (vector), n_metacells,
///    n_cells, n_unassigned
///  \item aggregated - If return_aggregated is TRUE, a list with indptr,
///    indices, raw_counts, norm_counts, nrow, ncol in sparse format. NULL
///    otherwise.
/// }
///
/// @export
#[extendr]
#[allow(clippy::too_many_arguments)]
fn rs_get_metacells(
    f_path: String,
    knn_mat: Option<RMatrix<i32>>,
    embd: Option<RMatrix<f64>>,
    meta_cell_params: List,
    target_size: f64,
    seed: usize,
    verbose: bool,
    return_aggregated: bool,
) -> extendr_api::Result<List> {
    let meta_cell_params = MetaCellParams::from_r_list(meta_cell_params);

    let knn_graph: Vec<Vec<usize>> = match (knn_mat, embd) {
        // first case - knn_mat provided
        (Some(knn_mat), _) => {
            if verbose {
                println!("Using provided kNN matrix");
            }
            let ncol = knn_mat.ncols();
            let nrow = knn_mat.nrows();
            let data = knn_mat.data();

            (0..nrow)
                .map(|j| {
                    (0..ncol)
                        .filter_map(|i| {
                            let val = data[j + i * nrow];
                            if val > 0 {
                                Some((val - 1) as usize)
                            } else {
                                None
                            }
                        })
                        .collect()
                })
                .collect()
        }
        // second case - knn needs to be re-calculated
        (None, Some(embd)) => {
            if verbose {
                println!("Calculating the kNN matrix from the provided data.");
            }

            let embd = r_matrix_to_faer_fp32(&embd);
            let knn_method = parse_knn_method(&meta_cell_params.knn_method).ok_or_else(|| {
                format!("Invalid KNN search method: {}", meta_cell_params.knn_method)
            })?;

            match knn_method {
                KnnSearch::Hnsw => generate_knn_hnsw(
                    embd.as_ref(),
                    &meta_cell_params.ann_dist,
                    meta_cell_params.k,
                    seed,
                    verbose,
                ),
                KnnSearch::Annoy => generate_knn_annoy(
                    embd.as_ref(),
                    &meta_cell_params.ann_dist,
                    meta_cell_params.k,
                    meta_cell_params.n_trees,
                    meta_cell_params.search_budget,
                    seed,
                    verbose,
                ),
                KnnSearch::NNDescent => generate_knn_nndescent(
                    embd.as_ref(),
                    &meta_cell_params.ann_dist,
                    meta_cell_params.k,
                    meta_cell_params.nn_max_iter,
                    meta_cell_params.delta,
                    meta_cell_params.rho,
                    seed,
                    verbose,
                ),
            }
        }
        // case three - nothing has been provided
        (None, None) => {
            return Err("Must provide either 'knn_mat' or 'embd' parameter".into());
        }
    };

    let nn_map = build_nn_map(&knn_graph);

    if verbose {
        println!("Identifying meta cells.");
    }

    let meta_cell_indices: Vec<&[usize]> = identify_meta_cells(
        &nn_map,
        meta_cell_params.max_shared,
        meta_cell_params.target_no_metacells,
        meta_cell_params.max_iter,
        seed,
        verbose,
    );

    let n_cells = nn_map.len();

    let assignments = metacells_to_assignments(&meta_cell_indices, n_cells);

    let assignment_list = assignments_to_r_list(&assignments, n_cells);

    if return_aggregated {
        if verbose {
            println!("Aggregating meta cells.");
        }

        let reader = ParallelSparseReader::new(&f_path).unwrap();
        let n_genes = reader.get_header().total_genes;

        let aggregated: CompressedSparseData<u32, f32> =
            aggregate_meta_cells(&reader, &meta_cell_indices, target_size as f32, n_genes);

        Ok(list!(
            assignments = assignment_list,
            aggregated = list!(
                indptr = aggregated.indptr.r_int_convert(),
                indices = aggregated.indices.r_int_convert(),
                raw_counts = aggregated.data.r_int_convert(),
                norm_counts = aggregated.data_2.unwrap().r_float_convert(),
                nrow = aggregated.shape.0,
                ncol = aggregated.shape.1
            )
        ))
    } else {
        Ok(list!(assignments = assignment_list, aggregated = NULL))
    }
}

/// Generate SEACells
///
/// @description This function implements the SEACells algorithm for generating
/// meta cells from Persad et al. An embedding matrix must be provided which is
/// used to construct the kNN graph and kernel matrix for the SEACells
/// algorithm. This version is highly memory and speed-optimised and will
/// truncate small values during matrix operations which can affect convergence.
///
/// @param f_path String. Path to the `counts_cells.bin` file.
/// @param embd Numerical matrix. The embedding matrix (for example PCA embedding)
/// used for the generation of the kNN graph and kernel matrix.
/// @param seacells_params A list containing the SEACells parameters.
/// @param target_size Numeric. Target library size for re-normalisation of
/// the meta cells. Typically `1e4`.
/// @param seed Integer. For reproducibility purposes.
/// @param verbose Boolean. Controls verbosity of the function.
/// @param return_aggregated Boolean. If TRUE, aggregates counts into meta cells.
///
/// @returns A list with the following elements:
/// \itemize{
///  \item assignments - A list containing assignment information with elements:
///    assignments (vector), metacells (list), unassigned (vector), n_metacells,
///    n_cells, n_unassigned
///  \item aggregated - If return_aggregated is TRUE, a list with indptr,
///    indices, raw_counts, norm_counts, nrow, ncol in sparse format. NULL
///    otherwise.
/// }
///
/// @export
///
/// @references Persad, et al., Nat. Biotechnol., 2023.
#[extendr]
fn rs_get_seacells(
    f_path: String,
    embd: RMatrix<f64>,
    seacells_params: List,
    target_size: f64,
    seed: usize,
    verbose: bool,
    return_aggregated: bool,
) -> extendr_api::Result<List> {
    let seacells_params = SEACellsParams::from_r_list(seacells_params);
    let embd = r_matrix_to_faer_fp32(&embd);
    let knn_method = parse_knn_method(&seacells_params.knn_method)
        .ok_or_else(|| format!("Invalid KNN search method: {}", seacells_params.knn_method))
        .unwrap();

    if verbose {
        println!(
            "Building the kNN graph with {} and {} neighbours",
            seacells_params.knn_method, seacells_params.k
        );
    }

    let (knn_indices, knn_dist) = match knn_method {
        KnnSearch::Annoy => {
            let annoy_index = build_annoy_index(embd.as_ref(), seacells_params.n_trees, seed);
            query_annoy_index(
                embd.as_ref(),
                &annoy_index,
                &seacells_params.ann_dist,
                seacells_params.k,
                seacells_params.search_budget,
                true,
                verbose,
            )
        }
        KnnSearch::Hnsw => {
            let hnsw_index = build_hnsw_index(embd.as_ref(), &seacells_params.ann_dist, seed);
            query_hnsw_index(
                embd.as_ref(),
                &hnsw_index,
                &seacells_params.ann_dist,
                seacells_params.k,
                true,
                verbose,
            )
        }
        KnnSearch::NNDescent => generate_knn_nndescent_with_dist(
            embd.as_ref(),
            &seacells_params.ann_dist,
            seacells_params.k,
            seacells_params.nn_max_iter,
            seacells_params.delta,
            seacells_params.rho,
            seed,
            verbose,
            true,
        ),
    };

    let mut seacell = SEACells::new(embd.nrows(), &seacells_params);

    println!("Number of neighbours: {}", knn_indices[0].len());

    seacell.construct_kernel_mat(embd.as_ref(), &knn_indices, &knn_dist.unwrap(), verbose);

    seacell.fit(seed, verbose);

    let assignments = seacell.get_hard_assignments();
    let k = seacells_params.n_sea_cells;

    let rss = seacell.get_rss_history();

    let assignments_opt: Vec<Option<usize>> = assignments.iter().map(|&x| Some(x)).collect();

    let assignment_list = assignments_to_r_list(&assignments_opt, embd.nrows());

    if return_aggregated {
        if verbose {
            println!("Aggregating meta cells.");
        }

        let meta_cell_indices = assignments_to_metacells(&assignments, k);
        let meta_cell_indices: Vec<&[usize]> =
            meta_cell_indices.iter().map(|v| v.as_slice()).collect();

        let reader = ParallelSparseReader::new(&f_path).unwrap();
        let n_genes = reader.get_header().total_genes;

        let aggregated: CompressedSparseData<u32, f32> =
            aggregate_meta_cells(&reader, &meta_cell_indices, target_size as f32, n_genes);

        Ok(list!(
            assignments = assignment_list,
            aggregated = list!(
                indptr = aggregated.indptr.r_int_convert(),
                indices = aggregated.indices.r_int_convert(),
                raw_counts = aggregated.data.r_int_convert(),
                norm_counts = aggregated.data_2.unwrap().r_float_convert(),
                nrow = aggregated.shape.0,
                ncol = aggregated.shape.1
            ),
            rss = rss.to_vec()
        ))
    } else {
        Ok(list!(
            assignments = assignment_list,
            aggregated = NULL,
            rss = rss.to_vec()
        ))
    }
}

extendr_module! {
    mod r_sc_analysis;
    fn rs_calculate_dge_mann_whitney;
    fn rs_aucell;
    fn rs_get_metacells;
    fn rs_get_seacells;
}
