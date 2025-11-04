use extendr_api::*;
use faer::Mat;
use rustc_hash::FxHashMap;
use std::time::Instant;

use crate::core::data::sparse_io::ParallelSparseReader;
use crate::core::data::sparse_structures::CompressedSparseData;
use crate::single_cell::cell_aggregations::*;
use crate::single_cell::dge_aucs::*;
use crate::single_cell::methods::seacells::*;
use crate::single_cell::sc_knn_snn::*;
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

/// Generate meta cells (hdWGCNA method)
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
/// @param cells_to_keep Optional indices of the cells to keep, i.e., the
/// cells used for the generation of the embedding.
/// @param cells_to_use Optional indices of cells to use for meta cell
/// generation. Useful if you wish to generate meta cells in specific cell
/// types.
/// @param meta_cell_params A list containing the meta cell parameters.
/// @param target_size Numeric. Target library size for re-normalisation of
/// the meta cells. Typically `1e4`.
/// @param seed Integer. For reproducibility purposes.
/// @param verbose Boolean. Controls verbosity of the function.
///
/// @returns A list with the following elements:
/// \itemize{
///  \item assignments - A list containing assignment information with elements:
///    assignments (vector), metacells (list), unassigned (vector), n_metacells,
///    n_cells, n_unassigned
///  \item aggregated - A list with indptr, indices, raw_counts, norm_counts,
///    nrow, ncol in sparse format.
/// }
///
/// @export
#[extendr]
#[allow(clippy::too_many_arguments)]
fn rs_get_metacells(
    f_path: String,
    knn_mat: Option<RMatrix<i32>>,
    embd: Option<RMatrix<f64>>,
    cells_to_keep: Option<Vec<i32>>,
    cells_to_use: Option<Vec<i32>>,
    meta_cell_params: List,
    target_size: f64,
    seed: usize,
    verbose: bool,
) -> extendr_api::Result<List> {
    let meta_cell_params = MetaCellParams::from_r_list(meta_cell_params);

    // If subsetting, we need both the QC cells and the subset to use
    if cells_to_use.is_some() && (cells_to_keep.is_none() || embd.is_none()) {
        return Err(
            "When using 'cells_to_use', both 'cells_to_keep' and 'embd' must be provided".into(),
        );
    }

    let (subset_to_orig, n_total_cells, knn_graph) = match cells_to_use {
        Some(ref use_cells) => {
            let cells_to_keep = cells_to_keep.unwrap();
            let embd = embd.unwrap();

            let qc_cells: Vec<usize> = cells_to_keep.iter().map(|&x| x as usize).collect();
            let use_cells: Vec<usize> = use_cells.iter().map(|&x| x as usize).collect();

            // Create mapping: original index -> PCA row index
            let orig_to_pca: FxHashMap<usize, usize> = qc_cells
                .iter()
                .enumerate()
                .map(|(pca_row, &orig_idx)| (orig_idx, pca_row))
                .collect();

            // Find which PCA rows correspond to cells_to_use
            let mut pca_rows_to_use = Vec::new();
            let mut subset_to_orig = Vec::new();

            for &orig_idx in &use_cells {
                if let Some(&pca_row) = orig_to_pca.get(&orig_idx) {
                    pca_rows_to_use.push(pca_row);
                    subset_to_orig.push(orig_idx);
                }
            }

            let n_total = use_cells.iter().max().map(|&x| x + 1).unwrap_or(0);

            if verbose {
                println!(
                    "Subsetting to {} cells (from {} QC-passing cells) and regenerating kNN graph",
                    pca_rows_to_use.len(),
                    qc_cells.len()
                );
            }

            // Subset the embedding using PCA row indices
            let ncol = embd.ncols();
            let nrow = embd.nrows();
            let data = embd.data();

            let subset_data: Vec<f64> = pca_rows_to_use
                .iter()
                .flat_map(|&row| (0..ncol).map(move |col| data[row + col * nrow]))
                .collect();

            let embd_subset = Mat::from_fn(pca_rows_to_use.len(), ncol, |i, j| {
                subset_data[i * ncol + j] as f32
            });

            // Generate kNN on subset
            let knn_method = parse_knn_method(&meta_cell_params.knn_method).ok_or_else(|| {
                format!("Invalid KNN search method: {}", meta_cell_params.knn_method)
            })?;

            let knn = match knn_method {
                KnnSearch::Hnsw => generate_knn_hnsw(
                    embd_subset.as_ref(),
                    &meta_cell_params.ann_dist,
                    meta_cell_params.k,
                    seed,
                    verbose,
                ),
                KnnSearch::Annoy => generate_knn_annoy(
                    embd_subset.as_ref(),
                    &meta_cell_params.ann_dist,
                    meta_cell_params.k,
                    meta_cell_params.n_trees,
                    meta_cell_params.search_budget,
                    seed,
                    verbose,
                ),
                KnnSearch::NNDescent => generate_knn_nndescent(
                    embd_subset.as_ref(),
                    &meta_cell_params.ann_dist,
                    meta_cell_params.k,
                    meta_cell_params.nn_max_iter,
                    meta_cell_params.delta,
                    meta_cell_params.rho,
                    seed,
                    verbose,
                ),
            };

            (subset_to_orig, n_total, knn)
        }
        None => {
            // No subsetting - original logic
            let n_total = match (&knn_mat, &embd) {
                (Some(mat), _) => mat.nrows(),
                (_, Some(em)) => em.nrows(),
                _ => return Err("Must provide either 'knn_mat' or 'embd' parameter".into()),
            };

            let knn = match (knn_mat, embd) {
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
                (None, Some(embd)) => {
                    if verbose {
                        println!("Calculating the kNN matrix from the provided data.");
                    }

                    let embd = r_matrix_to_faer_fp32(&embd);
                    let knn_method =
                        parse_knn_method(&meta_cell_params.knn_method).ok_or_else(|| {
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
                (None, None) => {
                    return Err("Must provide either 'knn_mat' or 'embd' parameter".into());
                }
            };

            ((0..n_total).collect(), n_total, knn)
        }
    };

    // Rest of the function remains the same...
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

    let n_subset_cells = nn_map.len();
    let assignments_subset = metacells_to_assignments(&meta_cell_indices, n_subset_cells);

    let is_subset = cells_to_use.is_some();

    let assignments_full = if is_subset {
        remap_assignments_to_original(&assignments_subset, &subset_to_orig, n_total_cells)
    } else {
        assignments_subset
    };

    let assignment_list = assignments_to_r_list(&assignments_full, n_total_cells);

    if verbose {
        println!("Aggregating meta cells.");
    }

    let reader = ParallelSparseReader::new(&f_path).unwrap();
    let n_genes = reader.get_header().total_genes;

    let metacells_original: Vec<Vec<usize>> = if is_subset {
        remap_metacells_to_original(&meta_cell_indices, &subset_to_orig)
    } else {
        meta_cell_indices
            .iter()
            .map(|&slice| slice.to_vec())
            .collect()
    };

    let metacells_refs: Vec<&[usize]> = metacells_original.iter().map(|v| v.as_slice()).collect();

    let aggregated: CompressedSparseData<u32, f32> =
        aggregate_meta_cells(&reader, &metacells_refs, target_size as f32, n_genes);

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
/// @param cells_to_keep Optional indices of the cells to keep, i.e., the
/// cells used for the generation of the embedding.
/// @param cells_to_use Optional indices of cells to use for meta cell
/// generation. Useful if you wish to generate meta cells in specific cell
/// types.
/// @param seacells_params A list containing the SEACells parameters.
/// @param target_size Numeric. Target library size for re-normalisation of
/// the meta cells. Typically `1e4`.
/// @param seed Integer. For reproducibility purposes.
/// @param verbose Boolean. Controls verbosity of the function.
///
/// @returns A list with the following elements:
/// \itemize{
///  \item assignments - A list containing assignment information with elements:
///    assignments (vector), metacells (list), unassigned (vector), n_metacells,
///    n_cells, n_unassigned
///  \item aggregated - A list with indptr, indices, raw_counts, norm_counts,
///    nrow, ncol in sparse format.
///  \item rss - Vector of RSS values from each iteration.
///  \item archetypes - Vector of cell indices selected as archetypes.
/// }
///
/// @export
///
/// @references Persad, et al., Nat. Biotechnol., 2023.
#[extendr]
#[allow(clippy::too_many_arguments)]
fn rs_get_seacells(
    f_path: String,
    embd: RMatrix<f64>,
    cells_to_keep: Option<Vec<i32>>,
    cells_to_use: Option<Vec<i32>>,
    seacells_params: List,
    target_size: f64,
    seed: usize,
    verbose: bool,
) -> extendr_api::Result<List> {
    let start_seacell = Instant::now();

    let seacells_params = SEACellsParams::from_r_list(seacells_params);

    let (subset_to_orig, n_total_cells, embd_mat) = match cells_to_use {
        Some(ref use_cells) => {
            if cells_to_keep.is_none() {
                return Err("When using 'cells_to_use', 'cells_to_keep' must be provided".into());
            }

            let cells_to_keep = cells_to_keep.unwrap();
            let qc_cells: Vec<usize> = cells_to_keep.iter().map(|&x| x as usize).collect();
            let use_cells: Vec<usize> = use_cells.iter().map(|&x| x as usize).collect();

            // Create mapping: original index -> PCA row index
            let orig_to_pca: FxHashMap<usize, usize> = qc_cells
                .iter()
                .enumerate()
                .map(|(pca_row, &orig_idx)| (orig_idx, pca_row))
                .collect();

            // Find which PCA rows correspond to cells_to_use
            let mut pca_rows_to_use = Vec::new();
            let mut subset_to_orig = Vec::new();

            for &orig_idx in &use_cells {
                if let Some(&pca_row) = orig_to_pca.get(&orig_idx) {
                    pca_rows_to_use.push(pca_row);
                    subset_to_orig.push(orig_idx);
                }
            }

            let n_total = use_cells.iter().max().map(|&x| x + 1).unwrap_or(0);

            if verbose {
                println!(
                    "Subsetting to {} cells (from {} QC-passing cells)",
                    pca_rows_to_use.len(),
                    qc_cells.len()
                );
            }

            // Subset the embedding using PCA row indices
            let ncol = embd.ncols();
            let nrow = embd.nrows();
            let data = embd.data();

            let subset_data: Vec<f64> = pca_rows_to_use
                .iter()
                .flat_map(|&row| (0..ncol).map(move |col| data[row + col * nrow]))
                .collect();

            let embd_subset = Mat::from_fn(pca_rows_to_use.len(), ncol, |i, j| {
                subset_data[i * ncol + j] as f32
            });

            (subset_to_orig, n_total, embd_subset)
        }
        None => {
            let n = embd.nrows();
            let embd_full = r_matrix_to_faer_fp32(&embd);
            ((0..n).collect(), n, embd_full)
        }
    };

    let is_subset = cells_to_use.is_some();

    let knn_method = parse_knn_method(&seacells_params.knn_method)
        .ok_or_else(|| format!("Invalid KNN search method: {}", seacells_params.knn_method))
        .unwrap();

    let start_knn = Instant::now();

    let (knn_indices, knn_dist) = match knn_method {
        KnnSearch::Annoy => {
            let annoy_index = build_annoy_index(embd_mat.as_ref(), seacells_params.n_trees, seed);
            query_annoy_index(
                embd_mat.as_ref(),
                &annoy_index,
                &seacells_params.ann_dist,
                seacells_params.k,
                seacells_params.search_budget,
                true,
                verbose,
            )
        }
        KnnSearch::Hnsw => {
            let hnsw_index = build_hnsw_index(embd_mat.as_ref(), &seacells_params.ann_dist, seed);
            query_hnsw_index(
                embd_mat.as_ref(),
                &hnsw_index,
                &seacells_params.ann_dist,
                seacells_params.k,
                true,
                verbose,
            )
        }
        KnnSearch::NNDescent => generate_knn_nndescent_with_dist(
            embd_mat.as_ref(),
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

    let end_knn = start_knn.elapsed();

    if verbose {
        println!(
            "kNN generation done in : {:.2?} with {}",
            end_knn, seacells_params.knn_method
        );
    }

    let mut seacell = SEACells::new(embd_mat.nrows(), &seacells_params);

    let knn_dist = knn_dist.unwrap();

    seacell.construct_kernel_mat(embd_mat.as_ref(), &knn_indices, &knn_dist, verbose);
    seacell.initialise_archetypes(&knn_indices, &knn_dist, verbose, seed as u64);

    seacell.fit(seed, verbose);

    let assignments_subset = seacell.get_hard_assignments();
    let archetypes_subset = seacell.get_archetypes();
    let k = seacells_params.n_sea_cells;

    let rss = seacell.get_rss_history();

    let assignments_subset_opt: Vec<Option<usize>> =
        assignments_subset.iter().map(|&x| Some(x)).collect();

    let assignments_full = if is_subset {
        remap_assignments_to_original(&assignments_subset_opt, &subset_to_orig, n_total_cells)
    } else {
        assignments_subset_opt
    };

    let archetypes_original: Vec<usize> = if is_subset {
        archetypes_subset
            .iter()
            .map(|&idx| subset_to_orig[idx])
            .collect()
    } else {
        archetypes_subset
    };

    let assignment_list = assignments_to_r_list(&assignments_full, n_total_cells);

    if verbose {
        println!("Aggregating meta cells.");
    }

    let meta_cell_indices = assignments_to_metacells(&assignments_subset, k);

    let metacells_original: Vec<Vec<usize>> = if is_subset {
        remap_metacells_to_original(
            &meta_cell_indices
                .iter()
                .map(|v| v.as_slice())
                .collect::<Vec<_>>(),
            &subset_to_orig,
        )
    } else {
        meta_cell_indices
    };

    let metacells_refs: Vec<&[usize]> = metacells_original.iter().map(|v| v.as_slice()).collect();

    let reader = ParallelSparseReader::new(&f_path).unwrap();
    let n_genes = reader.get_header().total_genes;

    let aggregated: CompressedSparseData<u32, f32> =
        aggregate_meta_cells(&reader, &metacells_refs, target_size as f32, n_genes);

    let end_seacell = start_seacell.elapsed();

    if verbose {
        println!("SEACells found in : {:.2?}", end_seacell);
    }

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
        rss = rss.to_vec(),
        archetypes = archetypes_original.r_int_convert()
    ))
}

extendr_module! {
    mod r_sc_analysis;
    fn rs_calculate_dge_mann_whitney;
    fn rs_aucell;
    fn rs_get_metacells;
    fn rs_get_seacells;
}
