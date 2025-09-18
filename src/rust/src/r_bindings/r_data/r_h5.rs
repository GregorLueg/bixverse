use extendr_api::prelude::*;

use crate::core::data::sparse_io::*;
use crate::core::data::sparse_io_h5::*;
use crate::single_cell::processing::*;

/// @export
#[extendr]
fn rs_h5ad_data(
    f_path: String,
    cs_type: String,
    nrows: usize,
    ncols: usize,
    cell_quality: List,
) -> extendr_api::Result<List> {
    let cell_quality = MinCellQuality::from_r_list(cell_quality);

    let file_quality: CellOnFileQuality =
        parse_h5_csr_quality(&f_path, (nrows, ncols), &cell_quality)
            .map_err(|e| extendr_api::Error::from(format!("HDF5 error: {}", e)))?;

    println!("No cells passing: {:?}", file_quality.cells_to_keep.len());
    println!("No genes passing: {:?}", file_quality.genes_to_keep.len());

    let res = read_h5ad_x_data(&f_path, &cs_type, &file_quality)
        .map_err(|e| extendr_api::Error::from(format!("HDF5 error: {}", e)))?;

    Ok(list!(
        data = res.data.iter().map(|x| *x as i32).collect::<Vec<i32>>(),
        indices = res.indices.iter().map(|x| *x as i32).collect::<Vec<i32>>(),
        indptr = res.indptr.iter().map(|x| *x as i32).collect::<Vec<i32>>(),
        no_genes = file_quality.genes_to_keep.len(),
        no_cells = file_quality.cells_to_keep.len()
    ))
}

extendr_module! {
    mod r_h5;
    fn rs_h5ad_data;
}
