use extendr_api::prelude::*;

use crate::core::methods::cistarget::*;
use crate::utils::r_rust_interface::*;

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
