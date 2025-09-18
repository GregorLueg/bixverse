use extendr_api::*;

use crate::core::base::loess::*;

#[extendr]
fn rs_2d_loess(x: &[f64], y: &[f64], span: f64, degree: usize) -> List {
    let loess = LoessRegression::new(span, degree);
    let loess_res: LoessRes = loess.fit(x, y);

    list!(
        predicted = loess_res.fitted_vals,
        residuals = loess_res.residuals,
        valid_idx = loess_res.valid_indices,
    )
}

extendr_module! {
  mod r_loess;
  fn rs_2d_loess;
}
