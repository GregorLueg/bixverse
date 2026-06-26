use bixverse_rs::core::math::stats::calc_fdr;
use bixverse_rs::methods::nmf_hals::HalsOpts;
use bixverse_rs::prelude::*;
use bixverse_rs::single_cell::sc_analysis::{
    dge_pathway_scores::*,
    hotspot::*,
    meld::*,
    milo_r::*,
    module_scoring::*,
    nichenet::activity_scoring::*,
    nichenet::ligand_regulatory_potential::*,
    nmf_sc::{nmf_multiple_run_sc, nmf_single_run_sc},
    scenic::*,
    vision::*,
};
use extendr_api::*;
use faer::Mat;
use rand::prelude::*;
use rayon::prelude::*;
use std::cmp::Ordering;

use crate::single_cell::utils::{knn_data_to_rust, panel_size_from_mem, prep_nichenet_network};

////////////////////
// extendr Module //
////////////////////

extendr_module! {
    mod r_sc_analysis;
    // dge
    fn rs_calculate_dge_mann_whitney;
    // aucell
    fn rs_aucell;
    // hotspot
    fn rs_hotspot_autocor;
    fn rs_hotspot_cluster_genes;
    fn rs_hotspot_gene_cor;
    // module scoring
    fn rs_module_scoring;
    // miloR
    fn rs_make_milor_nhoods;
    // MELD
    fn rs_meld_sc;
    // vision
    fn rs_vision;
    fn rs_vision_with_autocorrelation;
    // scenic
    fn rs_scenic_gene_filter;
    fn rs_scenic_grn;
    fn rs_scenic_grn_streaming;
    fn rs_top_k_targets;
    fn rs_importance_threshold;
    // nmf
    fn rs_nmf_single_sc;
    fn rs_nmf_multi_sc;
    // nichenet
    fn rs_generate_ligand_target_influence;
    fn rs_ligand_activity_scores;
}

//////////
// DGEs //
//////////

/// Calculate DGEs between cells based on Mann Whitney stats
///
/// @description
/// `r lifecycle::badge("experimental")`
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
/// @param verbose Integer. `0L` - quiet; `1L` - normal verbosity; `2L` -
/// detailed verbosity.
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
///
/// @keywords internal
#[extendr]
fn rs_calculate_dge_mann_whitney(
    f_path: String,
    cell_indices_1: &[i32],
    cell_indices_2: &[i32],
    min_prop: f64,
    alternative: String,
    verbose: usize,
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
    )
    .to_extendr()?;

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

/// Calculate module activity scores in Rust
///
/// @description
/// `r lifecycle::badge("experimental")`
/// Calculates module activity scores following Seurat's `AddModuleScore`.
/// For each module (gene set), computes the average expression of genes in the
/// set minus the average expression of randomly selected control genes from the
/// same expression bins. Genes are binned based on their average expression
/// across cells to ensure controls are expression-matched.
///
/// @param f_path_cells String. Path to the cell-based binary file.
/// @param f_path_genes String. Path to the gene-based binary file.
/// @param gs_list List. List of integer vectors, where each vector contains
/// gene indices (0-based) for a module/gene set.
/// @param cells_to_keep Integer. Vector of indices of the cells to keep.
/// @param nbin Integer. Number of bins for gene stratification.
/// @param ctrl Integer. Number of control genes to sample per gene in each
/// module.
/// @param seed Integer. Random seed for reproducible control gene sampling.
/// @param streaming Logical. If TRUE, processes cells and genes are read in in
/// chunks to reduce memory usage.
/// @param verbose Integer. `0L` - quiet; `1L` - normal verbosity; `2L` -
/// detailed verbosity.
///
/// @return Matrix of module scores (modules x cells). Each row corresponds to a
/// module from gs_list, each column to a cell from cells_to_keep.
///
/// @references
/// Tirosh et al, Science (2016)
///
/// @export
///
/// @keywords internal
#[extendr]
#[allow(clippy::too_many_arguments)]
fn rs_module_scoring(
    f_path_cells: String,
    f_path_genes: String,
    gs_list: List,
    cells_to_keep: Vec<i32>,
    nbin: usize,
    ctrl: usize,
    seed: usize,
    streaming: bool,
    verbose: usize,
) -> Result<RMatrix<f64>> {
    let cells_to_keep = cells_to_keep.r_int_convert();

    let mut gs_indices: Vec<Vec<usize>> = Vec::with_capacity(gs_list.len());
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

    let module_scores = calculate_module_scores_main(
        &f_path_genes,
        &f_path_cells,
        &gs_indices,
        &cells_to_keep,
        nbin,
        ctrl,
        streaming,
        seed,
        verbose,
    )
    .to_extendr()?;

    let nrows = module_scores[0].len(); // cells
    let ncols = module_scores.len(); // gene sets

    Ok(RMatrix::new_matrix(nrows, ncols, |r, c| {
        module_scores[c][r] as f64
    }))
}

/// Calculate AUCell in Rust
///
/// @description
/// `r lifecycle::badge("experimental")`
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
/// @param verbose Integer. `0L` - quiet; `1L` - normal verbosity; `2L` -
/// detailed verbosity.
///
/// @return A matrix of cells x gene sets with the values representing the
/// AUC.
///
/// @export
///
/// @keywords internal
#[extendr]
fn rs_aucell(
    f_path: String,
    gs_list: List,
    cells_to_keep: Vec<i32>,
    auc_type: &str,
    streaming: bool,
    verbose: usize,
) -> Result<RArray<f64, 2>> {
    let cells_to_keep = cells_to_keep.r_int_convert();

    let mut gs_indices: Vec<Vec<usize>> = Vec::with_capacity(gs_list.len());
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
        calculate_aucell_streaming(&f_path, &gs_indices, &cells_to_keep, auc_type, verbose)
            .to_extendr()?
    } else {
        calculate_aucell(&f_path, &gs_indices, &cells_to_keep, auc_type, verbose).to_extendr()?
    };

    let auc_mat = Mat::from_fn(res[0].len(), res.len(), |i, j| res[j][i] as f64);
    Ok(faer_to_r_matrix(auc_mat.as_ref()))
}

/// Calculate VISION pathway scores in Rust
///
/// @description
/// `r lifecycle::badge("experimental")`
/// The function will take in a list of gene sets that contains lists of `"pos"`
/// and `"neg"` gene indices (0-indexed). You don't have to provide the `"neg"`,
/// but it can be useful to classify the delta of two stats (EMT, Th1; Th2) etc.
///
/// @param f_path String. Path to the `counts_cells.bin` file.
/// @param gs_list Nested list. Each sublist contains the (0-indexed!) positive
/// and negative gene indices of that specific gene set.
/// @param cells_to_keep Integer. Vector of indices of the cells to keep.
/// @param streaming Boolean. Shall the data be streamed.
/// @param verbose Integer. `0L` - quiet; `1L` - normal verbosity; `2L` -
/// detailed verbosity.
///
/// @return A matrix of cells x vision scores per gene set.
///
/// @export
///
/// @keywords internal
#[extendr]
fn rs_vision(
    f_path: String,
    gs_list: List,
    cells_to_keep: Vec<i32>,
    streaming: bool,
    verbose: usize,
) -> Result<RArray<f64, 2>> {
    let cells_to_keep = cells_to_keep.r_int_convert();
    let gene_signatures = r_list_to_sig_genes(gs_list)?;

    let res = if streaming {
        calculate_vision_streaming(&f_path, &gene_signatures, &cells_to_keep, verbose)
            .to_extendr()?
    } else {
        calculate_vision(&f_path, &gene_signatures, &cells_to_keep, verbose).to_extendr()?
    };

    let vision_mat = Mat::from_fn(res.len(), res[0].len(), |i, j| res[i][j] as f64);

    Ok(faer_to_r_matrix(vision_mat.as_ref()))
}

/// Calculate VISION pathway scores in Rust with auto-correlation
///
/// @description
/// `r lifecycle::badge("experimental")`
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
/// @param knn_data Optional list. This contains pre-computed kNN data
/// (including distances). The user has to ensure consistency! If provided,
/// this will be used.
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
/// @param verbose Integer. `0L` - quiet; `1L` - normal verbosity; `2L` -
/// detailed verbosity.
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
///
/// @keywords internal
#[extendr]
#[allow(clippy::too_many_arguments)]
fn rs_vision_with_autocorrelation(
    f_path: String,
    embd: RMatrix<f64>,
    knn_data: Nullable<List>,
    gs_list: List,
    random_gs_list: List,
    vision_params: List,
    cells_to_keep: Vec<i32>,
    cluster_membership: Vec<i32>,
    streaming: bool,
    verbose: usize,
    seed: usize,
) -> extendr_api::Result<List> {
    assert!(
        embd.nrows() == cells_to_keep.len(),
        "The embedding matrix need to have the same nrow as the the cells to use"
    );

    let verbosity = parse_verbosity_level(verbose);
    let knn_provided = knn_data != extendr_api::Nullable::Null;

    if verbosity.normal_verbosity() {
        println!("Calculating the VISION scores of the actual gene sets.")
    }
    let cells_to_keep = cells_to_keep.r_int_convert();
    let gene_signatures = r_list_to_sig_genes(gs_list).unwrap();

    let res = if streaming {
        calculate_vision_streaming(&f_path, &gene_signatures, &cells_to_keep, verbose)
            .to_extendr()?
    } else {
        calculate_vision(&f_path, &gene_signatures, &cells_to_keep, verbose).to_extendr()?
    };

    if verbosity.normal_verbosity() {
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
                .to_extendr()?
        } else {
            calculate_vision(&f_path, &cluster_signatures, &cells_to_keep, verbose).to_extendr()?
        };

        random_scores_by_cluster.push(cluster_random_scores);

        if verbosity.normal_verbosity() {
            println!(
                "Completed random signatures for cluster {} / {}",
                cluster_idx + 1,
                random_gs_list.len()
            );
        }
    }

    let embd = r_matrix_to_faer_fp32(&embd);
    let knn_params = KnnParams::from_r_list(vision_params)?;

    let (knn_indices, knn_dist) = if knn_provided {
        if verbosity.normal_verbosity() {
            println!("Using generated kNN graph.")
        }

        let knn_data = knn_data
            .into_robj()
            .as_list()
            .ok_or_else(|| Error::Other("'knn_data' is not a list".into()))?;
        let (knn_indices, knn_dist, _, _) = knn_data_to_rust(knn_data)?;

        (knn_indices, knn_dist)
    } else {
        if verbosity.normal_verbosity() {
            println!("Using generated kNN graph.")
        }

        let (knn_indices, knn_dist) = generate_knn_with_dist(
            embd.as_ref(),
            &knn_params,
            true,
            false,
            seed,
            verbosity.detailed_verbosity(),
        )
        .to_extendr()?;

        (knn_indices, knn_dist.unwrap())
    };

    let cluster_membership = cluster_membership.r_int_convert_shift();

    let auto_cor_res = calc_autocorr_with_clusters(
        &res,
        &random_scores_by_cluster,
        &cluster_membership,
        knn_indices,
        knn_dist,
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

///////////////////////
// Gene module stuff //
///////////////////////

/// Calculate gene spatial auto-correlations
///
/// @description
/// `r lifecycle::badge("experimental")`
/// This function implements the HotSpot auto-correlation functionality and
/// will return to what extent a given gene shows auto-correlation in the
/// generated kNN-graph from the embeddings. For details see DeTomaso, et al.
///
/// @param f_path_genes Path to the `counts_genes.bin` file.
/// @param f_path_cells Path to the `counts_cells.bin` file.
/// @param embd Numerical matrix. The embedding matrix from which to generate
/// the kNN graph.
/// @param hotspot_params List. The HotSpot parameter list.
/// @param cells_to_keep Integer vector. 0-index vector indicating which cells
/// to include in the analysis. Ensure that this is of same order/length
/// as the embedding matrix.
/// @param knn_data Optional list. This contains pre-computed kNN data
/// (including distances). The user has to ensure consistency! If provided,
/// this will be used.
/// @param genes_to_use Integer vector. 0-index vector indicating which genes
/// to include.
/// @param streaming Boolean. Shall the data be streamed in chunks. Useful
/// for large data sets.
/// @param verbose Integer. `0L` - quiet; `1L` - normal verbosity; `2L` -
/// detailed verbosity.
/// @param seed Integer. Random seed for reproducibility.
///
/// @returns A list with the following elements.
/// \itemize{
///   \item gene_idx - 0-based integer indicating the gene index.
///   \item gaerys_c - Gaery's C calculation for the autocorrelation
///   coefficient.
///   \item z_score - Z-score of the auto-correlation.
///   \item pval - P-value derived from the Z-score.
///   \item fdr - False discovery rate based on the p-value.
/// }
///
/// @export
///
/// @references DeTomaso, et al., Cell Systems, 2021
///
/// @keywords internal
#[extendr]
#[allow(clippy::too_many_arguments)]
fn rs_hotspot_autocor(
    f_path_genes: String,
    f_path_cells: String,
    embd: RMatrix<f64>,
    knn_data: Nullable<List>,
    hotspot_params: List,
    cells_to_keep: Vec<i32>,
    genes_to_use: Vec<i32>,
    streaming: bool,
    verbose: usize,
    seed: usize,
) -> extendr_api::Result<List> {
    assert!(
        embd.nrows() == cells_to_keep.len(),
        "The embedding matrix need to have the same nrow as the the cells to use."
    );

    let verbosity = parse_verbosity_level(verbose);
    let knn_provided = knn_data != extendr_api::Nullable::Null;

    let hotspot_params = HotSpotParams::from_r_list(hotspot_params)?;

    let embd = r_matrix_to_faer_fp32(&embd);
    let cells_to_keep = cells_to_keep.r_int_convert();
    let genes_to_use = genes_to_use.r_int_convert();

    let (knn_indices, mut knn_dist) = if knn_provided {
        if verbosity.normal_verbosity() {
            println!("Using provided kNN graph...")
        }
        let knn_data = knn_data
            .into_robj()
            .as_list()
            .ok_or_else(|| Error::Other("'knn_data' is not a list".into()))?;
        let (knn_indices, knn_dist, _, _) = knn_data_to_rust(knn_data)?;

        (knn_indices, knn_dist)
    } else {
        if verbosity.normal_verbosity() {
            println!("Generating a kNN graph from scratch")
        }

        let (knn_indices, knn_dist) = generate_knn_with_dist(
            embd.as_ref(),
            &hotspot_params.knn_params,
            true,
            false,
            seed,
            verbosity.detailed_verbosity(),
        )
        .to_extendr()?;

        (knn_indices, knn_dist.unwrap())
    };

    let mut hotspot = Hotspot::new(
        f_path_genes,
        f_path_cells,
        &cells_to_keep,
        &knn_indices,
        &mut knn_dist,
    )
    .to_extendr()?;

    let res: HotSpotGeneRes = if streaming {
        hotspot
            .compute_all_genes_streaming(
                &genes_to_use,
                &hotspot_params.model,
                hotspot_params.normalise,
                verbose,
            )
            .to_extendr()?
    } else {
        hotspot
            .compute_all_genes(
                &genes_to_use,
                &hotspot_params.model,
                hotspot_params.normalise,
                verbose,
            )
            .to_extendr()?
    };

    Ok(list!(
        gene_idx = res.gene_idx,
        gaerys_c = res.c,
        z_score = res.z,
        pval = res.pval,
        fdr = res.fdr
    ))
}

/// Calculate gene to gene spatial correlations
///
/// @description
/// `r lifecycle::badge("experimental")`
/// This function implements the HotSpot gene <> gene local correlation
/// functionality from HotSpot, see DeTomaso, et al.
///
/// @param f_path_genes Path to the `counts_genes.bin` file.
/// @param f_path_cells Path to the `counts_cells.bin` file.
/// @param embd Numerical matrix. The embedding matrix from which to generate
/// the kNN graph.
/// @param knn_data Optional list. This contains pre-computed kNN data
/// (including distances). The user has to ensure consistency! If provided,
/// this will be used.
/// @param hotspot_params List. The HotSpot parameter list.
/// @param cells_to_keep Integer vector. 0-index vector indicating which cells
/// to include in the analysis. Ensure that this is of same order/length
/// as the embedding matrix.
/// @param genes_to_use Integer vector. 0-index vector indicating which genes
/// to include.
/// @param working_mem_gb Numeric. Approximate working memory (GB) the streaming
/// pair path may use for resident gene panels. Ignored when `streaming` is
/// `FALSE`. Larger values mean fewer disk re-reads. Note this excludes the two
/// dense N_genes x N_genes output matrices, which scale with `genes_to_use`.
/// @param streaming Boolean. Shall the data be streamed in chunks. Useful
/// for large data sets.
/// @param verbose Integer. `0L` - quiet; `1L` - normal verbosity; `2L` -
/// detailed verbosity.
/// @param seed Integer. Random seed for reproducibility.
///
/// @returns A list with the following elements.
/// \itemize{
///   \item cor - A matrix of the N x N genes_to_use length with the auto-
///   correlation coefficients.
///   \item z - A matrix of N x N genes_to_use length with the Z-scores of the
///   local correlations between two genes.
/// }
///
/// @export
///
/// @references DeTomaso, et al., Cell Systems, 2021
///
/// @keywords internal
#[extendr]
#[allow(clippy::too_many_arguments)]
fn rs_hotspot_gene_cor(
    f_path_genes: String,
    f_path_cells: String,
    embd: RMatrix<f64>,
    knn_data: Nullable<List>,
    hotspot_params: List,
    cells_to_keep: Vec<i32>,
    genes_to_use: Vec<i32>,
    working_mem_gb: f64,
    streaming: bool,
    verbose: usize,
    seed: usize,
) -> extendr_api::Result<List> {
    let embd = r_matrix_to_faer_fp32(&embd);
    let cells_to_keep = cells_to_keep.r_int_convert();
    let genes_to_use = genes_to_use.r_int_convert();

    let hotspot_params = HotSpotParams::from_r_list(hotspot_params)?;

    let verbosity = parse_verbosity_level(verbose);
    let knn_provided = knn_data != extendr_api::Nullable::Null;

    let (knn_indices, mut knn_dist) = if knn_provided {
        if verbosity.normal_verbosity() {
            println!("Using provided kNN graph...")
        }
        let knn_data = knn_data
            .into_robj()
            .as_list()
            .ok_or_else(|| Error::Other("'knn_data' is not a list".into()))?;
        let (knn_indices, knn_dist, _, _) = knn_data_to_rust(knn_data)?;

        (knn_indices, knn_dist)
    } else {
        if verbosity.normal_verbosity() {
            println!("Generating a kNN graph from scratch")
        }

        let (knn_indices, knn_dist) = generate_knn_with_dist(
            embd.as_ref(),
            &hotspot_params.knn_params,
            true,
            false,
            seed,
            verbosity.detailed_verbosity(),
        )
        .to_extendr()?;

        (knn_indices, knn_dist.unwrap())
    };

    let mut hotspot = Hotspot::new(
        f_path_genes,
        f_path_cells,
        &cells_to_keep,
        &knn_indices,
        &mut knn_dist,
    )
    .to_extendr()?;

    let res: HotSpotPairRes = if streaming {
        const STREAM_BATCH_SIZE: usize = 1000;
        let n_cells = knn_indices.len();
        let panel_size = panel_size_from_mem(working_mem_gb, n_cells, genes_to_use.len());
        hotspot
            .compute_gene_cor_streaming(
                &genes_to_use,
                &hotspot_params.model,
                STREAM_BATCH_SIZE,
                panel_size,
                verbose,
            )
            .to_extendr()?
    } else {
        hotspot
            .compute_gene_cor(&genes_to_use, &hotspot_params.model, verbose)
            .to_extendr()?
    };

    Ok(list!(
        cor = faer_to_r_matrix(res.cor.as_ref()),
        z = faer_to_r_matrix(res.z_scores.as_ref())
    ))
}

/// Cluster the genes by Z-score together
///
/// @description
/// `r lifecycle::badge("experimental")`
///
/// @param z_matrix Numerical matrix representing the Z-scores.
/// @param fdr_threshold Float. The FDR thresholds in terms of the Z-scores.
/// @param min_size Integer. Minimum cluster size.
///
/// @returns An assignment vector. NA indicates that the gene did not pass the
/// thresholds and has not been assigned.
///
/// @export
///
/// @keywords internal
#[extendr]
fn rs_hotspot_cluster_genes(
    z_matrix: RMatrix<f64>,
    fdr_threshold: f64,
    min_size: usize,
) -> Vec<f64> {
    let z_matrix = r_matrix_to_faer(&z_matrix);

    hotspot_gene_clusters(z_matrix, fdr_threshold, min_size)
}

////////////////////////////
// Differential abundance //
////////////////////////////

/// Generate the neighbourhoods akin to the miloR approach
///
/// @description
/// `r lifecycle::badge("experimental")`
///
/// @param embd Numeric matrix. Represents the matrix used to generate the kNN
/// graph and will be used to refine the neighbourhoods.
/// @param knn_indices Integer matrix. Each row represents a given cell and
/// the columns the neighbours. (0-indexed!)
/// @param milor_params Named list. Contains the parameters for running the
/// miloR approach.
/// @param seed Integer. Seed for reproducibility.
/// @param verbose Integer. `0L` - quiet; `1L` - normal verbosity; `2L` -
/// detailed verbosity.
///
/// @returns A list with the following elements:
/// \itemize{
///  \item index_cell - Integer. 0-indexed positions of the cells defining the
///  neighbourhood.
///  \item nhoods_i - Integer. 0-indexed positions of the cells in the
///  neighbourhood.
///  \item nhoods_j - Integer. To which neighbourhood the cell belongs.
///  \item nhoods_x - Numeric. The x-value of the COO type matrix, i.e.,
///  defaults to `1.0`.
///  \item nrows - Integer. Number of cells in the matrix
///  \item ncols - Integer. Number of refined neighbourhoods.
///  \item kth_distances - The k-th distances for spatial FDR calculations.
/// }
///
/// @export
///
/// @keywords internal
#[extendr]
fn rs_make_milor_nhoods(
    embd: RMatrix<f64>,
    knn_indices: RMatrix<i32>,
    milor_params: List,
    seed: usize,
    verbose: usize,
) -> Result<List> {
    let milor_params = MiloRParams::from_r_list(milor_params)?;
    let embd = r_matrix_to_faer_fp32(&embd);
    let k_original = knn_indices.ncols();
    let n_cells = embd.nrows();
    let knn_data = knn_indices.data();

    let verbosity = parse_verbosity_level(verbose);

    let knn_indices: Vec<Vec<usize>> = (0..n_cells)
        .map(|i| {
            (0..k_original)
                .map(|j| knn_data[i + j * n_cells] as usize)
                .collect()
        })
        .collect();

    let refinement_strategy = parse_refinement_strategy(&milor_params.refinement_strategy)
        .unwrap_or(RefinementStrategy::Approximate);

    let knn_idx = if matches!(refinement_strategy, RefinementStrategy::IndexBased) {
        let index_type =
            parse_index_type(&milor_params.index_type).unwrap_or(KnnIndexType::AnnoyIndex);

        Some(
            KnnIndex::new(
                embd.as_ref(),
                index_type,
                &milor_params.knn_params,
                seed,
                verbosity.detailed_verbosity(),
            )
            .to_extendr()?,
        )
    } else {
        None
    };

    // random sampling
    let n_cells = embd.nrows();
    let n_sample = ((n_cells as f64) * milor_params.prop).floor() as usize;
    let mut rng = StdRng::seed_from_u64(seed as u64);

    let mut random_indices: Vec<usize> = (0..n_cells).collect();
    random_indices.shuffle(&mut rng);
    random_indices.truncate(n_sample);

    if verbosity.normal_verbosity() {
        println!("Sampled {} vertices from {} cells", n_sample, n_cells);
    }

    // refinement
    let indices_refined = refine_sampling_with_strategy(
        embd.as_ref(),
        &knn_indices,
        &random_indices,
        milor_params.k_refine,
        &milor_params.knn_params,
        &refinement_strategy,
        knn_idx.as_ref(),
        verbose,
    )
    .to_extendr()?;

    // deduplicate
    let mut unique_indices = indices_refined;
    unique_indices.sort_unstable();
    unique_indices.dedup();

    let len_unique_indices = unique_indices.len();

    if verbosity.normal_verbosity() {
        println!(
            "Refined to {} unique neighbourhood indices",
            unique_indices.len()
        );
    }

    // build neighbourhood matrix
    let nhoods_triplets = build_nhood_matrix(&knn_indices, &unique_indices);

    // compute k-th distance
    let kth_distances = compute_kth_distances_from_matrix(
        embd.as_ref(),
        &knn_indices,
        &unique_indices,
        knn_indices[0].len() - 1,
    );

    Ok(list!(
        index_cell = unique_indices.r_int_convert(),
        nhoods_i = nhoods_triplets.0,
        nhoods_j = nhoods_triplets.1,
        nhoods_x = nhoods_triplets.2,
        nrows = n_cells as i32,
        ncols = len_unique_indices,
        kth_distances = kth_distances
    ))
}

//////////
// GRNs //
//////////

/// Identifies genes to include into a SCENIC analysis
///
/// @description
/// `r lifecycle::badge("experimental")`
///
/// @param f_path_genes Path to the `counts_genes.bin` file.
/// @param cell_indices Integer vector. 0-indexed(!) positions of cells to
/// include in the analysis
/// @param scenic_params Named list. Contains all of the parameters need for
/// SCENIC.
/// @param verbose Integer. `0L` - quiet; `1L` - normal verbosity; `2L` -
/// detailed verbosity.
///
/// @returns The 0-indexed positions of the genes to include in the scenic
/// analysis.
///
/// @export
///
/// @keywords internal
#[extendr]
fn rs_scenic_gene_filter(
    f_path_genes: String,
    cell_indices: Vec<i32>,
    scenic_params: List,
    verbose: usize,
) -> Result<Vec<i32>> {
    let cell_indices = cell_indices.r_int_convert();

    let scenic_params = ScenicParams::from_r_list(scenic_params)?;

    let genes_to_use =
        scenic_gene_filter(&f_path_genes, &cell_indices, &scenic_params, verbose).to_extendr()?;

    Ok(genes_to_use.r_int_convert())
}

/// SCENIC: Generating gene-regulatory networks
///
/// @description
/// `r lifecycle::badge("experimental")`
///
/// @param f_path_genes Path to the `counts_genes.bin` file.
/// @param cell_indices Integer vector. 0-indexed(!) positions of cells to
/// include in the analysis
/// @param gene_indices Integer vector. 0-indexed(!) positions of the genes
/// to include.
/// @param tf_indices Integer vector. 0-indexed(!) positions of the TF
/// predictor variables to use in the generation of the regression learners.
/// @param scenic_params Named list. Contains all of the parameters need for
/// SCENIC.
/// @param seed Integer. Controls reproducibility of the function.
/// @param verbose Integer. `0L` - quiet; `1L` - normal verbosity; `2L` -
/// detailed verbosity.
///
/// @returns A gene x TF importance matrix
///
/// @export
///
/// @keywords internal
#[extendr]
fn rs_scenic_grn(
    f_path_genes: String,
    cell_indices: Vec<i32>,
    gene_indices: Vec<i32>,
    tf_indices: Vec<i32>,
    scenic_params: List,
    seed: usize,
    verbose: usize,
) -> Result<RArray<f64, 2>> {
    let cell_indices = cell_indices.r_int_convert();
    let gene_indices = gene_indices.r_int_convert();
    let tf_indices = tf_indices.r_int_convert();

    let scenic_params = ScenicParams::from_r_list(scenic_params)?;

    let grn_matrix = run_scenic_grn(
        &f_path_genes,
        &cell_indices,
        &gene_indices,
        &tf_indices,
        &scenic_params,
        seed,
        verbose,
    )
    .to_extendr()?;

    Ok(faer_to_r_matrix(grn_matrix.as_ref()))
}

/// SCENIC: Generating gene-regulatory networks (streaming version)
///
/// @description
/// `r lifecycle::badge("experimental")`
/// Loads the genes in as chunks to avoid high memory pressure.
///
/// @param f_path_genes Path to the `counts_genes.bin` file.
/// @param cell_indices Integer vector. 0-indexed(!) positions of cells to
/// include in the analysis
/// @param gene_indices Integer vector. 0-indexed(!) positions of the genes
/// to include.
/// @param tf_indices Integer vector. 0-indexed(!) positions of the TF
/// predictor variables to use in the generation of the regression learners.
/// @param scenic_params Named list. Contains all of the parameters need for
/// SCENIC.
/// @param seed Integer. Controls reproducibility of the function.
/// @param verbose Integer. `0L` - quiet; `1L` - normal verbosity; `2L` -
/// detailed verbosity.
///
/// @returns A gene x TF importance matrix
///
/// @export
///
/// @keywords internal
#[extendr]
fn rs_scenic_grn_streaming(
    f_path_genes: String,
    cell_indices: Vec<i32>,
    gene_indices: Vec<i32>,
    tf_indices: Vec<i32>,
    scenic_params: List,
    seed: usize,
    verbose: usize,
) -> Result<RArray<f64, 2>> {
    let cell_indices = cell_indices.r_int_convert();
    let gene_indices = gene_indices.r_int_convert();
    let tf_indices = tf_indices.r_int_convert();

    let scenic_params = ScenicParams::from_r_list(scenic_params)?;

    let grn_matrix = run_scenic_grn_streaming(
        &f_path_genes,
        &cell_indices,
        &gene_indices,
        &tf_indices,
        &scenic_params,
        seed,
        verbose,
    )
    .to_extendr()?;

    Ok(faer_to_r_matrix(grn_matrix.as_ref()))
}

/// SCENIC: Select the Top TF <> Gene pairs
///
/// @description
/// `r lifecycle::badge("experimental")`
///
/// @param matrix Numeric matrix with genes x TF importance values
/// @param k Integer. Number of top genes / TFs to extract.
/// @param margin If set to 1, the top k TFs per gene are used. If set to 2, the
/// top k genes per TF are used. Both versions were used in the original paper.
/// @param min_value Float. An
///
/// @returns A list with three vectors: tf, gene, importance
///
/// @export
///
/// @keywords internal
#[extendr]
fn rs_top_k_targets(matrix: RMatrix<f64>, k: i32, margin: i32, min_value: Option<f64>) -> List {
    let nrow = matrix.nrows();
    let ncol = matrix.ncols();
    let data = matrix.data();

    let dimnames = matrix.get_attrib("dimnames").unwrap();
    let dimnames_list: List = dimnames.try_into().unwrap();
    let rownames: Strings = dimnames_list.elt(0).unwrap().try_into().unwrap();
    let colnames: Strings = dimnames_list.elt(1).unwrap().try_into().unwrap();

    let (outer_len, inner_len) = match margin {
        1 => (nrow, ncol),
        2 => (ncol, nrow),
        _ => panic!("margin must be 1 (rows) or 2 (cols)"),
    };
    let k = (k as usize).min(inner_len);

    let mut tfs: Vec<String> = Vec::new();
    let mut genes: Vec<String> = Vec::new();
    let mut importances: Vec<f64> = Vec::new();

    for i in 0..outer_len {
        let mut idx: Vec<usize> = (0..inner_len)
            .filter(|&j| {
                let v = match margin {
                    1 => data[j * nrow + i],
                    _ => data[i * nrow + j],
                };
                min_value.is_none_or(|min| v >= min)
            })
            .collect();

        let take = k.min(idx.len());
        if take == 0 {
            continue;
        }

        idx.select_nth_unstable_by(take - 1, |&a, &b| {
            let va = match margin {
                1 => data[a * nrow + i],
                _ => data[i * nrow + a],
            };
            let vb = match margin {
                1 => data[b * nrow + i],
                _ => data[i * nrow + b],
            };
            vb.partial_cmp(&va).unwrap_or(Ordering::Equal)
        });

        for &j in &idx[..take] {
            let (tf_name, gene_name, val) = match margin {
                1 => (
                    colnames[j].as_ref(),
                    rownames[i].as_ref(),
                    data[j * nrow + i],
                ),
                _ => (
                    colnames[i].as_ref(),
                    rownames[j].as_ref(),
                    data[i * nrow + j],
                ),
            };
            tfs.push(tf_name.to_string());
            genes.push(gene_name.to_string());
            importances.push(val);
        }
    }

    list!(tf = tfs, gene = genes, importance = importances)
}

/// SCENIC: Select TF-gene pairs by per-gene importance threshold
///
/// @description
/// `r lifecycle::badge("experimental")`
/// For each gene (row), computes mean + n_sd * SD of the importance scores
/// across all TFs and retains only pairs exceeding that threshold.
///
/// @param matrix Numeric matrix with genes (rows) x TFs (columns) importance
/// values.
/// @param n_sd Float. Number of standard deviations above the mean to use as
/// the per-gene threshold.
/// @param min_value Optional float. Absolute minimum importance score. Pairs
/// below this are excluded even if they pass the per-gene threshold.
///
/// @returns A list with three vectors: tf, gene, importance
///
/// @export
///
/// @keywords internal
#[extendr]
fn rs_importance_threshold(matrix: RMatrix<f64>, n_sd: f64, min_value: Option<f64>) -> List {
    let nrow = matrix.nrows();
    let ncol = matrix.ncols();
    let data = matrix.data();

    let dimnames = matrix.get_attrib("dimnames").unwrap();
    let dimnames_list: List = dimnames.try_into().unwrap();
    let rownames: Strings = dimnames_list.elt(0).unwrap().try_into().unwrap();
    let colnames: Strings = dimnames_list.elt(1).unwrap().try_into().unwrap();

    let mut tfs: Vec<String> = Vec::new();
    let mut genes: Vec<String> = Vec::new();
    let mut importances: Vec<f64> = Vec::new();

    for i in 0..nrow {
        let mut sum = 0.0;
        let mut sum_sq = 0.0;
        for j in 0..ncol {
            let v = data[j * nrow + i];
            sum += v;
            sum_sq += v * v;
        }
        let mean = sum / ncol as f64;
        let var = sum_sq / ncol as f64 - mean * mean;
        let sd = var.max(0.0).sqrt();
        let threshold = mean + n_sd * sd;

        for j in 0..ncol {
            let v = data[j * nrow + i];
            if v >= threshold && min_value.is_none_or(|min| v >= min) {
                tfs.push(colnames[j].as_ref().to_string());
                genes.push(rownames[i].as_ref().to_string());
                importances.push(v);
            }
        }
    }

    list!(tf = tfs, gene = genes, importance = importances)
}

//////////
// MELD //
//////////

/// Run MELD
///
/// @description
/// `r lifecycle::badge("experimental")`
/// This implements a Rust-based version of the MELD algorithm, see Burkhardt,
/// et al. Nat. Biotechnol., 2021.
///
/// @param embd Numeric matrix. The original embedding that was used to generate
/// the kNN graph.
/// @param knn_data Optional named list. This contains pre-computed kNN data
/// (including distances). The user has to ensure consistency! If provided, this
/// will be used.
/// @param meld_params Named list. Contains the parameters to use for MELD.
/// @param landmark Boolean. Shall a landmark method be used for accelerated
/// MELD.
/// @param n_landmarks Integer. If `landmark = TRUE`, how many landmarks to use.
/// @param labels Integer. The labels of the different groups. (1-indexed!)
/// @param seed Integer. For reproducibility.
/// @param verbose Integer. `0L` - quiet; `1L` - normal verbosity; `2L` -
/// detailed verbosity.
///
/// @returns A list with the following items
/// \itemize{
///   \item raw_scores - The raw MELD scores
///   \item norm_scores - Negative values were clamped to 0 and the rows L1
///   normalised. This yields probability-like values.
/// }
///
/// @export
///
/// @keywords internal
#[extendr]
#[allow(clippy::too_many_arguments)]
fn rs_meld_sc(
    embd: RMatrix<f64>,
    knn_data: Nullable<List>,
    meld_params: List,
    labels: &[i32],
    n_labels: usize,
    seed: usize,
    verbose: usize,
) -> Result<List> {
    let embd = r_matrix_to_faer_fp32(&embd);
    let meld_params = MeldParams::from_r_list(meld_params)?;
    let labels = labels.r_int_convert_shift();
    let verbosity = parse_verbosity_level(verbose);

    // deal with kNN
    let knn_provided = knn_data != extendr_api::Nullable::Null;

    let (knn_indices, knn_dist, dist) = if knn_provided {
        if verbosity.normal_verbosity() {
            println!("Using provided kNN graph...")
        }
        let knn_data = knn_data
            .into_robj()
            .as_list()
            .ok_or_else(|| Error::Other("'knn_data' is not a list".into()))?;
        let (knn_indices, knn_dist, _, dist) = knn_data_to_rust(knn_data)?;

        (knn_indices, knn_dist, dist)
    } else {
        if verbosity.normal_verbosity() {
            println!("Generating a kNN graph from scratch")
        }

        let (knn_indices, knn_dist) = generate_knn_with_dist(
            embd.as_ref(),
            &meld_params.knn_params,
            true,
            false,
            seed,
            verbosity.detailed_verbosity(),
        )
        .to_extendr()?;

        (
            knn_indices,
            knn_dist.unwrap(),
            meld_params.knn_params.ann_dist.clone(),
        )
    };

    let is_squared_distance = dist == "euclidean";

    let (meld_raw, meld_norm) = meld(
        &knn_indices,
        &knn_dist,
        &labels,
        n_labels,
        is_squared_distance,
        &meld_params,
        seed as u64,
        verbose,
    )
    .to_extendr()?;

    Ok(list!(
        raw_scores = faer_to_r_matrix(meld_raw.as_ref()),
        norm_scores = faer_to_r_matrix(meld_norm.as_ref())
    ))
}

/////////
// NMF //
/////////

/// Run NMF (HALS) over a set of single cells and genes
///
/// @description
/// `r lifecycle::badge("experimental")`
/// Runs a single NMF (HALS) run with the specified initialisation.
///
/// @param f_path_gene Path to the `counts_genes.bin` file.
/// @param gene_indices Integer vector. 0-indexed(!) positions of the genes
/// to include.
/// @param cell_indices Integer vector. 0-indexed(!) positions of cells to
/// include in the analysis.
/// @param k Integer. Number of latent factors to return.
/// @param preprocessing String. One of `c("none", "sd", "sqrt_sd")`. Takes the
/// data as is, or scales by standard deviation or squared standard deviation
/// per feature.
/// @param use_second_layer Boolean. If `TRUE`, runs NMF on the normalised
/// counts; if `FALSE`, on the raw counts.
/// @param nmf_hals_params Named list. Contains the NMF parameters.
/// @param seed Integer. Random seed for initialisation.
/// @param verbose Integer. `0L` - quiet; `1L` - normal verbosity; `2L` -
/// detailed verbosity.
///
/// @returns A list with the following items
/// \itemize{
///   \item w - The left factor matrix (n_features x k)
///   \item h - The right factor matrix (k x n_samples)
///   \item final_loss - Loss at the final iteration
///   \item n_iter - Number of iterations the algorithm run for
///   \item converged - Did the NMF algorithm converge
/// }
///
/// @export
///
/// @keywords internal
#[extendr]
#[allow(clippy::too_many_arguments)]
fn rs_nmf_single_sc(
    f_path_gene: &str,
    gene_indices: &[i32],
    cell_indices: &[i32],
    k: usize,
    preprocessing: &str,
    use_second_layer: bool,
    nmf_hals_params: List,
    seed: usize,
    verbose: usize,
) -> Result<List> {
    let cell_indices = cell_indices.r_int_convert();
    let gene_indices = gene_indices.r_int_convert();
    let nmf_hals_opt: HalsOpts<f32> = HalsOpts::from_r_list(nmf_hals_params, seed).to_extendr()?;
    let nmf_res = nmf_single_run_sc(
        f_path_gene,
        &gene_indices,
        &cell_indices,
        k,
        preprocessing,
        use_second_layer,
        Some(nmf_hals_opt),
        verbose,
    )
    .to_extendr()?;
    Ok(list!(
        w = faer_to_r_matrix(nmf_res.w.as_ref()),
        h = faer_to_r_matrix(nmf_res.h.as_ref()),
        final_loss = nmf_res.final_loss as f64,
        n_iter = nmf_res.n_iter,
        converged = nmf_res.converged
    ))
}

/// Run multiple NMF (HALS) restarts over a set of single cells and genes
///
/// @description
/// `r lifecycle::badge("experimental")`
/// Runs `n_runs` HALS NMF with random initialisations seeded by `seed + i`.
/// The `nmf_init` field in `nmf_hals_params` is ignored; random init is always
/// used. The returned `w_all` is the column-binding of all run W matrices.
///
/// @param f_path_gene Path to the `counts_genes.bin` file.
/// @param gene_indices Integer vector. 0-indexed(!) positions of the genes
/// to include.
/// @param cell_indices Integer vector. 0-indexed(!) positions of cells to
/// include in the analysis.
/// @param k Integer. Number of latent factors per run.
/// @param preprocessing String. One of `c("none", "sd", "sqrt_sd")`.
/// @param use_second_layer Boolean. If `TRUE`, runs NMF on the normalised
/// counts; if `FALSE`, on the raw counts.
/// @param nmf_hals_params Named list. Contains the NMF parameters.
/// @param n_runs Integer. Number of random restarts.
/// @param seed Integer. Base random seed. Run `i` uses `seed + i`.
/// @param verbose Integer. `0L` - quiet; `1L` - normal verbosity; `2L` -
/// detailed verbosity.
///
/// @returns A list with the following items
/// \itemize{
///   \item w_all - Column-bound W matrices across all runs,
///   shape `n_features x (k * n_runs)`. Columns `i*k+1..(i+1)*k` are run `i`'s
///   components (1-indexed).
///   \item h_per_run - List of H matrices, each `k x n_cells`.
///   \item losses - Numeric vector. Final reconstruction loss per run.
///   \item converged - Logical vector. Convergence flag per run.
///   \item best_idx - Integer. 1-indexed position of the run with the lowest
///   final loss.
/// }
///
/// @export
///
/// @keywords internal
#[extendr]
#[allow(clippy::too_many_arguments)]
fn rs_nmf_multi_sc(
    f_path_gene: &str,
    gene_indices: &[i32],
    cell_indices: &[i32],
    k: usize,
    preprocessing: &str,
    use_second_layer: bool,
    nmf_hals_params: List,
    n_runs: usize,
    seed: usize,
    verbose: usize,
) -> Result<List> {
    let gene_indices = gene_indices.r_int_convert();
    let cell_indices = cell_indices.r_int_convert();
    let nmf_hals_opt: HalsOpts<f32> = HalsOpts::from_r_list(nmf_hals_params, seed).to_extendr()?;
    let nmf_res = nmf_multiple_run_sc(
        f_path_gene,
        &gene_indices,
        &cell_indices,
        k,
        preprocessing,
        use_second_layer,
        Some(nmf_hals_opt),
        n_runs,
        seed,
        verbose,
    )
    .to_extendr()?;
    let h_per_run: List = nmf_res
        .h_per_run
        .iter()
        .map(|h| faer_to_r_matrix(h.as_ref()))
        .collect();
    Ok(list!(
        w_all = faer_to_r_matrix(nmf_res.w_all.as_ref()),
        h_per_run = h_per_run,
        losses = nmf_res.losses.r_float_convert(),
        converged = nmf_res.converged,
        best_idx = (nmf_res.best_idx + 1) as i32
    ))
}

//////////////
// NicheNet //
//////////////

/// Generate the ligand to target influence matrices
///
/// @description
/// `r lifecycle::badge("experimental")`
/// Helper function to generate the ligand to target influence matrix for the
/// NicheNet like approach.
///
/// @param ligand_seeds List. Contains the indices of the seeds, i.e., ligands.
/// @param ppi_network Named list. Contains the PPI network with the ligand
/// to receptor to signalling to TFs. Must contain from (indices), to
/// (indices), and edge weights.
/// @param grn_network Named list. Contains the gene regulatory network with the
/// TF to target gene network. Must contain from (indices), to (indices), and
/// edge weights.
/// @param n_nodes Integer. Number of total nodes.
/// @param params Named list.
///
/// @returns A dense matrix of ligands x genes that contains the influence
/// scores of each
///
/// @export
#[extendr]
fn rs_generate_ligand_target_influence(
    ligand_seeds: List,
    ppi_network: List,
    grn_network: List,
    n_nodes: usize,
    params: List,
) -> Result<RMatrix<f64>> {
    let (sig_from, sig_to, sig_weight) = prep_nichenet_network(ppi_network)?;
    let (grn_from, grn_to, grn_weight) = prep_nichenet_network(grn_network)?;
    let mut seed_vec = Vec::with_capacity(ligand_seeds.len());

    for i in 0..ligand_seeds.len() {
        let elem_i = ligand_seeds.elt(i)?;
        let vec_i = elem_i
            .as_integer_vector()
            .ok_or_else(|| {
                Error::Other("One of the ligand seeds could not be transformed to integers.".into())
            })?
            .iter()
            .map(|x| *x as u32)
            .collect();
        seed_vec.push(vec_i);
    }

    let params: LigandTargetParams<f64> = LigandTargetParams::from_r_list(params)?;

    let ligand_influence_matrix = construct_ligand_target_mat(
        n_nodes,
        &sig_from,
        &sig_to,
        &sig_weight,
        &grn_from,
        &grn_to,
        &grn_weight,
        &seed_vec,
        &params,
    )
    .to_extendr()?;

    Ok(faer_to_r_matrix(ligand_influence_matrix.as_ref()))
}

/// [LigandActivityScores] to R list Helper
///
/// ### Params
///
/// * `scores` - The [LigandActivityScores] structure to transform
///
/// ### Returns
///
/// The list with the results
fn activity_to_list(scores: &LigandActivityScores<f64>) -> List {
    list!(
        auroc = scores.auroc.clone(),
        aupr = scores.aupr.clone(),
        aupr_corrected = scores.aupr_corrected.clone(),
        pearson = scores.pearson.clone(),
        spearman = scores.spearman.clone()
    )
}

/// Calculate the NicheNet ligand activity scores
///
/// @description
/// `r lifecycle::badge("experimental")`
///
/// @param ligand_influence A ligand x background genes matrix that measures the
/// ligand to target gene influence.
/// @param in_gene_sets A list of logicals with the genes of interest being set
/// to `TRUE` and the background genes set to `FALSE`.
///
/// @returns A list with internal lists with:
/// \itemize{
///   \item `auroc` - The Area Under the Receiver Operating Characteristic for
///   that ligand
///   \item `aupr` - The Area Under the Precision-Recall curve for that ligand.
///   \item `aupr_corrected` - The corrected AUPR
///   \item `pearson` - The Pearson correlations
///   \item `spearman` - The Spearman correlations
/// }
///
/// @export
#[extendr]
fn rs_ligand_activity_scores(ligand_influence: RMatrix<f64>, in_gene_sets: List) -> Result<List> {
    let ligand_influence = r_matrix_to_faer(&ligand_influence);

    let mut in_gene_sets_vec = Vec::with_capacity(in_gene_sets.len());

    for i in 0..in_gene_sets.len() {
        let elem_i = in_gene_sets.elt(i)?;
        let vec_i: Vec<bool> = elem_i
            .as_logical_vector()
            .ok_or_else(|| {
                Error::Other("One of the ligand seeds could not be transformed to integers.".into())
            })?
            .iter()
            .map(|x| x.to_bool())
            .collect();
        in_gene_sets_vec.push(vec_i);
    }

    let scores: Vec<LigandActivityScores<f64>> = in_gene_sets_vec
        .par_iter()
        .map(|x| ligand_activity_scores(&ligand_influence.as_ref(), x))
        .collect();

    let mut res = List::new(scores.len());

    for i in 0..res.len() {
        let res_i = activity_to_list(&scores[i]);
        res.set_elt(i, res_i.into())?;
    }

    Ok(res)
}
