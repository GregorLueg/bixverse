use bixverse_rs::methods::lda::{LdaParams, lda_fit, lda_k_sweep};
use bixverse_rs::methods::methods_r_wrapper::{lda_result_to_r_list, lda_sweep_to_r_list};
use bixverse_rs::prelude::*;
use extendr_api::prelude::*;

/////////////
// extendR //
/////////////

extendr_module! {
    mod r_lda;
    fn rs_lda;
    fn rs_lda_k_sweep;
}

///////////////
// Functions //
///////////////

/// Fit a latent Dirichlet allocation model to a document-term matrix
///
/// @description
/// `r lifecycle::badge("experimental")`
/// Variational Bayes, following Hoffman, et al. The matrix is documents x
/// terms; for a cisTopic-style run that is a binarised cells x regions (or
/// cells x regulons) matrix, but any count matrix works.
///
/// @param sparse_data A named list that needs to have `data`, `indptr`,
/// `indices`, `nrow`, `ncol` and `cs_type`. Documents x terms.
/// @param k Integer. Number of topics.
/// @param lda_params List. The LDA parameters, see [bixverse::params_lda()].
/// @param seed Integer. Seed for the variational initialisation.
/// @param verbose Integer. `0L` - quiet; `1L` - normal verbosity; `2L` -
/// detailed verbosity.
///
/// @returns A list with the following items
/// \itemize{
///   \item cell_topic - Topic proportions per document, `k x n_documents`.
///   \item topic_region - Term probabilities per topic, `n_terms x k`.
///   \item bound - Final variational bound (ELBO). Higher is better.
///   \item perplexity - Per-token perplexity. Lower is better.
///   \item n_iter - Outer iterations run.
///   \item converged - Whether the relative bound change fell below `tol`.
/// }
///
/// @references Hoffman, Blei and Bach, NIPS, 2010
///
/// @export
///
/// @keywords internal
#[extendr]
fn rs_lda(
    sparse_data: List,
    k: usize,
    lda_params: List,
    seed: usize,
    verbose: usize,
) -> Result<List, extendr_api::Error> {
    let sparse: CompressedSparseData2<f64, f64> =
        list_to_sparse_matrix(sparse_data, true).to_extendr()?;
    let lda_opt: LdaParams<f64> = LdaParams::from_r_list(lda_params, seed)?;
    let model = lda_fit(&sparse, k, Some(lda_opt), verbose).to_extendr()?;

    lda_result_to_r_list(&model)
}

/// Fit LDA across a range of topic counts and score each fit
///
/// @description
/// `r lifecycle::badge("experimental")`
/// Fits one model per requested topic count and evaluates each with the Arun,
/// Cao Juan and Mimno metrics. The corpus is built once and shared, and each
/// fit is already parallel over documents, so the topic counts run
/// sequentially rather than nesting Rayon pools.
///
/// @param sparse_data A named list that needs to have `data`, `indptr`,
/// `indices`, `nrow`, `ncol` and `cs_type`. Documents x terms.
/// @param k_range Integer vector. Topic counts to evaluate.
/// @param lda_params List. The LDA parameters, see [bixverse::params_lda()].
/// @param top_topics_coh Optional integer. Number of top-scoring topics
/// averaged into the reported coherence. If `NULL`, defaults to `5L`.
/// @param seed Integer. Seed for the variational initialisation.
/// @param verbose Integer. `0L` - quiet; `1L` - normal verbosity; `2L` -
/// detailed verbosity.
///
/// @returns A list with the following items
/// \itemize{
///   \item k - The topic counts that were tried.
///   \item models - One fitted model per topic count, see [rs_lda()].
///   \item metrics - One metric list per topic count, with `arun_2010`,
///   `cao_juan_2009`, `mimno_2011`, `coherence_per_topic`, `bound` and
///   `perplexity`.
///   \item combined_score - Rescaled combined score per entry. `NaN` where the
///   entry was excluded from selection by the coherence topic-count floor.
///   \item best_k - Topic count with the highest combined score.
/// }
///
/// @references Arun, et al., PAKDD, 2010; Cao, et al., Neurocomputing, 2009;
/// Mimno, et al., EMNLP, 2011
///
/// @export
///
/// @keywords internal
#[extendr]
#[allow(clippy::too_many_arguments)]
fn rs_lda_k_sweep(
    sparse_data: List,
    k_range: &[i32],
    lda_params: List,
    top_topics_coh: Option<usize>,
    seed: usize,
    verbose: usize,
) -> Result<List, extendr_api::Error> {
    let k_range = k_range.r_int_convert();
    let sparse: CompressedSparseData2<f64, f64> =
        list_to_sparse_matrix(sparse_data, true).to_extendr()?;
    let lda_opt: LdaParams<f64> = LdaParams::from_r_list(lda_params, seed)?;
    let sweep = lda_k_sweep(&sparse, &k_range, Some(lda_opt), top_topics_coh, verbose)
        .to_extendr()?;

    lda_sweep_to_r_list(&sweep)
}
