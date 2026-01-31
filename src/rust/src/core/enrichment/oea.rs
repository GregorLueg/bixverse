////////////////////
// Result structs //
////////////////////

///////////////
// Functions //
///////////////

/// Calculate odds ratios
///
/// ### Params
///
/// * `a1_b1` - In both gene set and target set
/// * `a0_b1` - In gene set, but not in target set
/// * `a1_b0` - In target set, but not in gene set
/// * `a0_b0` - Not in either
///
/// ### Return
///
/// The odds ratio. Pending values, can become infinity.
#[inline]
pub fn hypergeom_odds_ratio(a1_b1: usize, a0_b1: usize, a1_b0: usize, a0_b0: usize) -> f64 {
    (a1_b1 as f64 / a0_b1 as f64) / (a1_b0 as f64 / a0_b0 as f64)
}
