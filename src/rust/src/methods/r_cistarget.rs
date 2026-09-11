use bixverse_rs::methods::cis_target::*;
use bixverse_rs::methods::methods_r_wrapper::motif_enrichments_to_r_list;
use bixverse_rs::prelude::*;
use extendr_api::prelude::*;

use rayon::prelude::*;
use std::sync::atomic::{AtomicUsize, Ordering};

/////////////
// extendR //
/////////////

extendr_module! {
    mod r_cistarget;
    fn rs_cistarget;
}

///////////////
// Functions //
///////////////

/// Run CisTarget motif enrichment analysis
///
/// @description
/// `r lifecycle::badge("experimental")`
/// Core Rust function for motif enrichment analysis using recovery curves.
/// Gene sets are processed in parallel; only motifs with an NES at or above
/// `nes_threshold` are returned.
///
/// @param rankings Integer matrix with motif rankings for genes (genes in rows,
/// motifs in columns). Lower ranks indicate higher regulatory potential.
/// @param gs_list List of integer vectors. Each element contains 1-based
/// indices of genes in the gene set (matching row indices in rankings).
/// @param auc_threshold Integer. Absolute number of top-ranked genes to use
/// for the AUC calculation (e.g., for 5% of 10000 genes, use 500).
/// @param nes_threshold Numeric. Normalised Enrichment Score threshold for
/// filtering significant motifs.
/// @param max_rank Integer. Maximum rank to consider for the recovery curves
/// (at most `nrow(rankings)`).
/// @param method String. Recovery curve calculation method, one of
/// `c("approx", "icistarget")`. Anything else falls back to `"approx"`.
/// @param n_mean Integer. Window size for the smoothing in the approximate
/// method.
/// @param verbose Boolean. Report progress per decile of gene sets.
///
/// @returns List of lists, one per gene set, each containing
/// \itemize{
///   \item motif_idx - 1-based column index of the motif in `rankings`.
///   \item nes - Normalised enrichment score.
///   \item auc - Area under the recovery curve.
///   \item rank_at_max - Rank at which the leading edge is reached.
///   \item n_enriched - Number of genes in the leading edge.
///   \item leading_edge - List of 1-based row indices of the leading edge
///   genes, one element per motif.
/// }
///
/// @export
#[extendr]
#[allow(clippy::too_many_arguments)]
fn rs_cistarget(
    rankings: RMatrix<i32>,
    gs_list: List,
    auc_threshold: i32,
    nes_threshold: f64,
    max_rank: i32,
    method: String,
    n_mean: usize,
    verbose: bool,
) -> List {
    let rankings = r_matrix_to_faer(&rankings);
    let mut gs_indices: Vec<Vec<usize>> = Vec::with_capacity(gs_list.len());
    let rcc_method = parse_rcc_type(&method).unwrap_or(RccType::Approx);
    for i in 0..gs_list.len() {
        let list_elem = gs_list.elt(i).unwrap();
        let elem = list_elem
            .as_integer_vector()
            .unwrap()
            .iter()
            .map(|x| (*x - 1) as usize)
            .collect();
        gs_indices.push(elem);
    }

    let total = gs_indices.len();
    let done = AtomicUsize::new(0);
    let last_decile = AtomicUsize::new(0);

    let results: Vec<Vec<MotifEnrichment<f64>>> = gs_indices
        .par_iter()
        .map(|gs_idx| {
            let res = process_gene_set(
                rankings,
                gs_idx,
                auc_threshold,
                nes_threshold,
                max_rank,
                &rcc_method,
                n_mean,
            );
            if verbose {
                let n = done.fetch_add(1, Ordering::Relaxed) + 1;
                let decile = (n * 10) / total;
                if decile > last_decile.swap(decile, Ordering::Relaxed) {
                    println!(" cistarget: {}% ({}/{})", decile * 10, n, total);
                }
            }
            res
        })
        .collect();

    let mut r_results = List::new(results.len());
    for (i, enrichments) in results.iter().enumerate() {
        r_results
            .set_elt(i, Robj::from(motif_enrichments_to_r_list(enrichments)))
            .unwrap();
    }
    r_results
}
