use extendr_api::prelude::*;

use bixverse_rs::prelude::*;
use bixverse_rs::single_cell::sc_data::h5ad_io::*;

/////////////
// extendR //
/////////////

extendr_module! {
    mod r_h5;
    fn rs_h5ad_data;
}

///////////////
// Fucntions //
///////////////

/// Load in h5ad data via Rust
///
/// @description
/// `r lifecycle::badge("experimental")`
/// Loads in h5ad data within Rust and automatically converts the data into
/// CSR with cells x genes.
///
/// @param f_path File path. The path to the h5ad file.
/// @param cs_type String. One of `c("csr", "csc")`. How the data is stored in
/// the file. Other values raise an error.
/// @param nrows Integer. Number of rows in the file.
/// @param ncols Integer. Number of columns in the file.
/// @param cell_quality List. Specifying the cell quality. Please refer
/// to [bixverse::params_sc_min_quality()].
/// @param slot String. In which slot the raw data can be found. One of
/// `c("X", "raw", "layers.counts")`. Unknown strings default to `"X"`.
/// @param verbose Boolean. Controls verbosity of the function
///
/// @returns A list with the CSR data (cells x genes) of the cells and genes
/// passing `cell_quality`:
/// \itemize{
///   \item data - The counts of the sparse matrix.
///   \item indices - The 0-based gene indices of the sparse matrix.
///   \item indptr - The index pointers of the sparse matrix.
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
    slot: String,
    verbose: bool,
) -> extendr_api::Result<List> {
    let cell_quality = MinCellQuality::from_r_list(cell_quality)?;

    let raw_slot = parse_raw_slot(&slot).unwrap_or_else(|| {
        println!(
            "The provided string ({:?}) could not be matched. Defaulting to X",
            slot
        );
        RawDataSlot::default()
    });

    let file_format = parse_compressed_sparse_format(&cs_type)
        .ok_or_else(|| BixverseErrors::UnknownSparseFormat(cs_type.to_string()))
        .to_extendr()?;

    let file_quality = match file_format {
        CompressedSparseFormat::Csr => {
            parse_h5_csr_quality(&f_path, (nrows, ncols), &raw_slot, &cell_quality, verbose)
                .to_extendr()?
        }
        CompressedSparseFormat::Csc => {
            parse_h5_csc_quality(&f_path, (nrows, ncols), &cell_quality, &raw_slot, verbose)
                .to_extendr()?
        }
    };

    let file_data: CompressedSparseData2<u32> = match file_format {
        CompressedSparseFormat::Csr => {
            read_h5ad_x_data_csr(&f_path, &file_quality, &raw_slot, verbose).to_extendr()?
        }
        CompressedSparseFormat::Csc => {
            let data =
                read_h5ad_x_data_csc(&f_path, &file_quality, &raw_slot, verbose).to_extendr()?;
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
