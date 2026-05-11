use extendr_api::prelude::*;

use bixverse_rs::prelude::*;
use bixverse_rs::single_cell::sc_data::h5ad_io::*;

extendr_module! {
    mod r_h5;
    fn rs_h5ad_data;
}

/// Load in h5ad data via Rust
///
/// @description
/// Loads in h5ad data within Rust and automatically converts the data into
/// CSR with cells x genes.
///
/// @param f_path File path. The path to the h5ad file.
/// @param cs_type String. Is the data stored in CSC or CSR.
/// @param nrows Integer. Number of rows in the file.
/// @param ncols Integer. Number of columns in the file.
/// @param cell_quality List. Specifiying the cell quality. Please refer
/// to [bixverse::params_sc_min_quality()].
/// @param verbose Boolean. Controls verbosity of the function
///
/// @returns A list with:
/// \itemize{
///   \item data - The data of the sparse matrix stored on the h5ad file.
///   \item indices - The indices of the sparse matrix stored in the h5ad file.
///   \item indptr - The indptr of the sparse matrix stored in the h5ad file.
///   \item no_genes - No of genes in the sparse matrix (i.e., ncol).
///   \item no_cells - No of cells in the sparse matrix (i.e., nrow).
/// }
///
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
    let cell_quality = MinCellQuality::from_r_list(cell_quality)?;

    let file_format = parse_compressed_sparse_format(&cs_type)
        .ok_or_else(|| BixverseErrors::UnknownSparseFormat(cs_type.to_string()))
        .to_extendr()?;

    let file_quality = match file_format {
        CompressedSparseFormat::Csr => {
            parse_h5_csr_quality(&f_path, (nrows, ncols), &cell_quality, verbose).unwrap()
        }
        CompressedSparseFormat::Csc => {
            parse_h5_csc_quality(&f_path, (nrows, ncols), &cell_quality, verbose).unwrap()
        }
    };

    let file_data: CompressedSparseData2<u32> = match file_format {
        CompressedSparseFormat::Csr => {
            read_h5ad_x_data_csr(&f_path, &file_quality, verbose).to_extendr()?
        }
        CompressedSparseFormat::Csc => {
            let data = read_h5ad_x_data_csc(&f_path, &file_quality, verbose).to_extendr()?;
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
