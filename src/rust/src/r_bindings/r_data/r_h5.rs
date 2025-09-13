use extendr_api::prelude::*;

use crate::core::data::sparse_io_h5::*;

/// @export
#[extendr]
fn rs_h5ad_data(
    f_path: String,
    cs_type: String,
    nrows: usize,
    ncols: usize,
) -> extendr_api::Result<List> {
    let res = read_h5ad_x_data(&f_path, &cs_type, (nrows, ncols))
        .map_err(|e| extendr_api::Error::from(format!("HDF5 error: {}", e)))?;

    Ok(list!(
        data = res.data.iter().map(|x| *x as i32).collect::<Vec<i32>>(),
        indices = res.indices.iter().map(|x| *x as i32).collect::<Vec<i32>>(),
        indptr = res.indptr.iter().map(|x| *x as i32).collect::<Vec<i32>>()
    ))
}

extendr_module! {
    mod r_h5;
    fn rs_h5ad_data;
}
