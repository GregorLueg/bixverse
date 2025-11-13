use extendr_api::*;
use faer::Mat;

use crate::core::base::stats::calc_fdr;
use crate::single_cell::dge_pathway_scores::*;
use crate::single_cell::methods::vision_hotspot::*;
use crate::single_cell::sc_knn_snn::*;
use crate::utils::r_rust_interface::*;
use crate::utils::traits::*;

extendr_module! {
    mod r_sc_analysis;
    fn rs_aucell;
    fn rs_calculate_dge_mann_whitney;
    fn rs_hotspot_autocor;
    fn rs_hotspot_gene_cor;
    fn rs_vision;
    fn rs_vision_with_autocorrelation;
}

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
/// @return A matrix of cells x gene sets with the values representing the
/// AUC.
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

    let cells_to_keep = cells_to_keep.r_int_convert();

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

    let auc_mat = Mat::from_fn(res[0].len(), res.len(), |i, j| res[j][i] as f64);
    Ok(faer_to_r_matrix(auc_mat.as_ref()))
}

/// Calculate VISION pathway scores in Rust
///
/// @description
/// The function will take in a list of gene sets that contains lists of `"pos"`
/// and `"neg"` gene indices (0-indexed). You don't have to provide the `"neg"`,
/// but it can be useful to classify the delta of two stats (EMT, Th1; Th2) etc.
///
/// @param f_path String. Path to the `counts_cells.bin` file.
/// @param gs_list Nested list. Each sublist contains the (0-indexed!) positive
/// and negative gene indices of that specific gene set.
/// @param cells_to_keep Integer. Vector of indices of the cells to keep.
/// @param streaming Boolean. Shall the data be streamed.
/// @param verbose Boolean. Controls verbosity of the function.
///
/// @return A matrix of cells x vision scores per gene set.
///
/// @export
#[extendr]
fn rs_vision(
    f_path: String,
    gs_list: List,
    cells_to_keep: Vec<i32>,
    streaming: bool,
    verbose: bool,
) -> extendr_api::Result<extendr_api::RArray<f64, [usize; 2]>> {
    let cells_to_keep = cells_to_keep.r_int_convert();
    let gene_signatures = r_list_to_sig_genes(gs_list)?;

    let res = if streaming {
        calculate_vision_streaming(&f_path, &gene_signatures, &cells_to_keep, verbose)
    } else {
        calculate_vision(&f_path, &gene_signatures, &cells_to_keep, verbose)
    };

    let vision_mat = Mat::from_fn(res.len(), res[0].len(), |i, j| res[i][j] as f64);

    Ok(faer_to_r_matrix(vision_mat.as_ref()))
}

/// Calculate VISION pathway scores in Rust with auto-correlation
///
/// @description
/// The function will take in a list of gene sets that contains lists of `"pos"`
/// and `"neg"` gene indices (0-indexed). You don't have to provide the `"neg"`,
/// but it can be useful to classify the delta of two stats (EMT, Th1; Th2) etc.
/// Additionally, it will take a random gene list and calculate an
/// auto-correlation score based on Gaery's C to identify pathways that show
/// significant patterns on the kNN graph generate on the provided embedding.
///
/// @param f_path String. Path to the `counts_cells.bin` file.
/// @param embd Numerical matrix. The embedding matrix to use to generate the
/// kNN graph.
/// @param gs_list Nested list. Each sublist contains the (0-indexed!) positive
/// and negative gene indices of that specific gene set.
/// @param random_gs_list Double-nested list. The outer list represents the
/// clusters of clusters and the inner list represents the permutations within
/// that cluster.
/// @param vision_params List. Contains various parameters to use in terms
/// of the kNN generation.
/// @param cells_to_keep Integer. Vector of indices of the cells to keep.
/// @param cluster_membership Integer. Vector that indicates to which of the
/// permuted gene set clusters the given gene set belongs.
/// @param streaming Boolean. Shall the data be streamed.
/// @param verbose Boolean. Controls verbosity of the function.
/// @param seed Integer. Random seed for reproducibility.
///
/// @return A list with the following items:
/// \itemize{
///   \item autocor_res - Auto-correlation results, i.e., 1 - C, p-value and
///   FDR.
///   \item vision_mat - A matrix of cells x vision scores per gene set.
/// }
///
/// @export
#[extendr]
#[allow(clippy::too_many_arguments)]
fn rs_vision_with_autocorrelation(
    f_path: String,
    embd: RMatrix<f64>,
    gs_list: List,
    random_gs_list: List,
    vision_params: List,
    cells_to_keep: Vec<i32>,
    cluster_membership: Vec<i32>,
    streaming: bool,
    verbose: bool,
    seed: usize,
) -> extendr_api::Result<List> {
    if verbose {
        println!("Calculating the VISION scores of the actual gene sets.")
    }
    let cells_to_keep = cells_to_keep.r_int_convert();
    let gene_signatures = r_list_to_sig_genes(gs_list).unwrap();

    let res = if streaming {
        calculate_vision_streaming(&f_path, &gene_signatures, &cells_to_keep, verbose)
    } else {
        calculate_vision(&f_path, &gene_signatures, &cells_to_keep, verbose)
    };

    if verbose {
        println!("Calculating the VISION scores of the permuted gene sets.")
    }

    let mut random_scores_by_cluster: Vec<Vec<Vec<f32>>> = Vec::with_capacity(random_gs_list.len());

    for cluster_idx in 0..random_gs_list.len() {
        let cluster_sigs = random_gs_list
            .elt(cluster_idx)?
            .as_list()
            .ok_or("Cluster element not a list")?;

        let cluster_signatures = r_list_to_sig_genes(cluster_sigs).unwrap();

        let cluster_random_scores = if streaming {
            calculate_vision_streaming(&f_path, &cluster_signatures, &cells_to_keep, verbose)
        } else {
            calculate_vision(&f_path, &cluster_signatures, &cells_to_keep, verbose)
        };

        random_scores_by_cluster.push(cluster_random_scores);

        if verbose {
            println!(
                "Completed random signatures for cluster {} / {}",
                cluster_idx + 1,
                random_gs_list.len()
            );
        }
    }

    let embd = r_matrix_to_faer_fp32(&embd);
    let knn_params = KnnParams::from_r_list(vision_params);

    let (knn_indices, knn_dist) =
        generate_knn_with_dist(embd.as_ref(), &knn_params, true, seed, verbose);

    let cluster_membership = cluster_membership.r_int_convert_shift();

    let auto_cor_res = calc_autocorr_with_clusters(
        &res,
        &random_scores_by_cluster,
        &cluster_membership,
        knn_indices,
        knn_dist.unwrap(),
        verbose,
    );

    let gaery_c = auto_cor_res.0;
    let p_val = auto_cor_res.1;
    let fdr = calc_fdr(&p_val);

    let vision_mat = Mat::from_fn(res.len(), res[0].len(), |i, j| res[i][j] as f64);

    Ok(list!(
        autocor_res = list!(auto_cor = gaery_c, p_val = p_val, fdr = fdr),
        vision_mat = faer_to_r_matrix(vision_mat.as_ref())
    ))
}

/// @export
#[extendr]
#[allow(clippy::too_many_arguments)]
fn rs_hotspot_autocor(
    f_path_genes: String,
    f_path_cells: String,
    embd: RMatrix<f64>,
    hotspot_params: List,
    cells_to_keep: Vec<i32>,
    genes_to_use: Vec<i32>,
    streaming: bool,
    verbose: bool,
    seed: usize,
) -> extendr_api::Result<List> {
    let embd = r_matrix_to_faer_fp32(&embd);
    let cells_to_keep = cells_to_keep.r_int_convert();
    let genes_to_use = genes_to_use.r_int_convert();

    let hotspot_params = HotSpotParams::from_r_list(hotspot_params);
    let knn_params = KnnParams::from_hotspot_params(&hotspot_params);

    if verbose {
        println!("Generating kNN graph...")
    }

    let (knn_indices, knn_dist) =
        generate_knn_with_dist(embd.as_ref(), &knn_params, true, seed, verbose);

    let mut knn_dist = knn_dist.unwrap();

    let mut hotspot = Hotspot::new(
        f_path_genes,
        f_path_cells,
        &cells_to_keep,
        &knn_indices,
        &mut knn_dist,
    );

    let res: HotSpotGeneRes = if streaming {
        hotspot.compute_all_genes_streaming(
            &genes_to_use,
            &hotspot_params.model,
            hotspot_params.normalise,
            verbose,
        )
    } else {
        hotspot.compute_all_genes(
            &genes_to_use,
            &hotspot_params.model,
            hotspot_params.normalise,
            verbose,
        )
    }?;

    Ok(list!(
        gene_idx = res.gene_idx,
        gaerys_c = res.c,
        z_score = res.z,
        pval = res.pval,
        fdr = res.fdr
    ))
}

/// @export
#[extendr]
#[allow(clippy::too_many_arguments)]
fn rs_hotspot_gene_cor(
    f_path_genes: String,
    f_path_cells: String,
    embd: RMatrix<f64>,
    hotspot_params: List,
    cells_to_keep: Vec<i32>,
    genes_to_use: Vec<i32>,
    streaming: bool,
    verbose: bool,
    seed: usize,
) -> extendr_api::Result<List> {
    let embd = r_matrix_to_faer_fp32(&embd);
    let cells_to_keep = cells_to_keep.r_int_convert();
    let genes_to_use = genes_to_use.r_int_convert();

    let hotspot_params = HotSpotParams::from_r_list(hotspot_params);
    let knn_params = KnnParams::from_hotspot_params(&hotspot_params);

    if verbose {
        println!("Generating kNN graph...")
    }

    let (knn_indices, knn_dist) =
        generate_knn_with_dist(embd.as_ref(), &knn_params, true, seed, verbose);

    let mut knn_dist = knn_dist.unwrap();

    let mut hotspot = Hotspot::new(
        f_path_genes,
        f_path_cells,
        &cells_to_keep,
        &knn_indices,
        &mut knn_dist,
    );

    let res: HotSpotPairRes = if streaming {
        hotspot.compute_gene_cor_streaming(&genes_to_use, &hotspot_params.model, verbose)
    } else {
        hotspot.compute_gene_cor(&genes_to_use, &hotspot_params.model, verbose)
    }
    .unwrap();

    Ok(list!(
        cor = faer_to_r_matrix(res.cor.as_ref()),
        z = faer_to_r_matrix(res.z_scores.as_ref())
    ))
}
