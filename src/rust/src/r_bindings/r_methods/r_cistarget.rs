use extendr_api::prelude::*;

use crate::core::methods::cistarget::*;
use crate::utils::r_rust_interface::*;

/// Run CisTarget motif enrichment analysis
///
/// @description Core Rust function for motif enrichment analysis using recovery
/// curves.
///
/// @param rankings Integer matrix with motif rankings for genes (genes in rows,
/// motifs in columns). Lower ranks indicate higher regulatory potential.
/// @param gs_list List of integer vectors. Each element contains 1-based
/// indices of genes in the gene set (matching row indices in rankings).
/// @param auc_threshold Absolute number of top-ranked genes to use for AUC
/// calculation (e.g., for 5% of 10000 genes, use 500).
/// @param nes_threshold Normalised Enrichment Score threshold for filtering
/// significant motifs
/// @param max_rank Maximum rank to consider (typically nrow(rankings)).
/// @param method Recovery curve calculation method: "approx" or "icistarget".
/// @param n_mean Number of points for averaging in approximate method.
///
/// @return List of lists, one per gene set, each containing motif_idx, nes, auc,
/// rank_at_max, n_enriched, and leading_edge.
///
/// @export
#[extendr]
fn rs_cistarget(
    rankings: RMatrix<i32>,
    gs_list: List,
    auc_threshold: i32,
    nes_threshold: f64,
    max_rank: i32,
    method: String,
    n_mean: usize,
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

    let results: Vec<Vec<MotifEnrichment>> = gs_indices
        .iter()
        .map(|gs_idx| {
            process_gene_set(
                rankings,
                gs_idx,
                auc_threshold,
                nes_threshold,
                max_rank,
                &rcc_method,
                n_mean,
            )
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

extendr_module! {
    mod r_cistarget;
    fn rs_cistarget;
}
