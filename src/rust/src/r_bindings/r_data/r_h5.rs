use extendr_api::prelude::*;

use crate::core::data::sparse_io_h5::*;
use crate::core::data::sparse_structures::*;
use crate::single_cell::processing::*;

/// @export
#[extendr]
fn rs_h5ad_data(
    f_path: String,
    cs_type: String,
    nrows: usize,
    ncols: usize,
    cell_quality: List,
    verbose: bool,
) -> extendr_api::Result<List> {
    let cell_quality = MinCellQuality::from_r_list(cell_quality);

    let file_format = parse_cs_format(&cs_type).unwrap();

    let file_quality = match file_format {
        CompressedSparseFormat::Csc => {
            parse_h5_csc_quality(&f_path, (nrows, ncols), &cell_quality, verbose).unwrap()
        }
        CompressedSparseFormat::Csr => {
            parse_h5_csr_quality(&f_path, (nrows, ncols), &cell_quality, verbose).unwrap()
        }
    };

    println!("No cells passing: {:?}", file_quality.cells_to_keep.len());
    println!("No genes passing: {:?}", file_quality.genes_to_keep.len());

    let file_data: CompressedSparseData<u16> = match file_format {
        CompressedSparseFormat::Csc => {
            read_h5ad_x_data_csc(&f_path, &file_quality, verbose).unwrap()
        }
        CompressedSparseFormat::Csr => {
            let data = read_h5ad_x_data_csr(&f_path, &file_quality, verbose).unwrap();
            data.transpose_and_convert()
        }
    };

    Ok(list!(
        data = file_data
            .data
            .iter()
            .map(|x| *x as i32)
            .collect::<Vec<i32>>(),
        indices = file_data
            .indices
            .iter()
            .map(|x| *x as i32)
            .collect::<Vec<i32>>(),
        indptr = file_data
            .indptr
            .iter()
            .map(|x| *x as i32)
            .collect::<Vec<i32>>(),
        no_genes = file_quality.genes_to_keep.len(),
        no_cells = file_quality.cells_to_keep.len()
    ))
}

extendr_module! {
    mod r_h5;
    fn rs_h5ad_data;
}
