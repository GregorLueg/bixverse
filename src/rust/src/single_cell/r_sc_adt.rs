use bixverse_rs::prelude::*;
use extendr_api::*;
use rayon::prelude::*;

////////////////////
// extendr Module //
////////////////////

extendr_module! {
    mod r_sc_adt;
    // processing
    fn rs_adt_clr;
}

////////////////
// Processing //
////////////////

/// Applies CLR normalisation on ADT counts (Seurat-style, per cell)
///
/// @param counts R matrix of shape cells x features.
///
/// @returns CLR-transformed matrix.
///
/// @export
#[extendr]
fn rs_adt_clr(counts: RMatrix<f64>) -> RMatrix<f64> {
    let counts = r_matrix_to_faer(&counts);
    let nrow = counts.nrows();
    let ncol = counts.ncols();

    // per-row geometric-mean factor
    let g: Vec<f64> = (0..nrow)
        .into_par_iter()
        .map(|i| {
            let mut s = 0.0_f64;
            for j in 0..ncol {
                let v = counts[(i, j)];
                if v > 0.0 {
                    s += v.ln_1p();
                }
            }
            (s / ncol as f64).exp()
        })
        .collect();

    RMatrix::new_matrix(nrow, ncol, |i, j| (counts[(i, j)] / g[i]).ln_1p())
}
