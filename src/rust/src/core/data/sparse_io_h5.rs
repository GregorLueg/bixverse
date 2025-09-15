use hdf5::{File, Result};
use std::path::Path;

use crate::core::data::sparse_io::*;
use crate::core::data::sparse_structures::*;
use crate::single_cell::processing::*;

/// Helper function to parse compressed sparse format
///
/// ### Params
///
/// * `s` - String specifying the type
///
/// ### Return
///
/// Returns an `Option<CompressedSparseFormat>`
pub fn parse_cs_format(s: &str) -> Option<CompressedSparseFormat> {
    let res = match s.to_lowercase().as_str() {
        "csr" => Some(CompressedSparseFormat::Csr),
        "csc" => Some(CompressedSparseFormat::Csc),
        _ => None,
    };

    res
}

/////////////
// Writers //
/////////////

/// Writes h5ad data to disk in the binarised format
///
/// ### Params
///
/// * `h5_path` - Path to the h5 object.
/// * `bin_path` - Path to the binarised object on disk to write to
/// * `cs_type` - Was the h5ad data stored in "csc" or "csr". Important! h5ad
///   stores data in genes x cells; bixverse stores in cells x genes!
/// * `no_cells` - Total number of obversations in the data.
/// * `no_genes` - Total number of vars in the data.
/// * `min_genes` - Minimum number of genes to be detected in the cell.
/// * `size_factor` - To which library size to normalise the data
///
/// ### Returns
///
/// A boolean vector indicating which cells have sufficient genes.
pub fn write_h5_counts<P: AsRef<Path>>(
    h5_path: P,
    bin_path: P,
    cs_type: &str,
    no_cells: usize,
    no_genes: usize,
    cell_quality: MinCellQuality,
) -> CellQuality {
    let file_data: CompressedSparseData<u16> =
        read_h5ad_x_data(h5_path, cs_type, (no_genes, no_cells)).unwrap();

    let file_data = file_data.transpose_and_convert();

    let mut writer = CellGeneSparseWriter::new(bin_path, true, no_cells, no_genes).unwrap();

    let (cell_chunk_vec, cell_qc) =
        CsrCellChunk::generate_chunks_sparse_data(file_data, cell_quality);

    let mut cells_writen = 0_usize;

    for cell_chunk in cell_chunk_vec {
        if cell_chunk.to_keep {
            writer.write_cell_chunk(cell_chunk).unwrap();
            cells_writen += 1;
        }
    }

    writer.update_header_no_cells(cells_writen);
    writer.finalise().unwrap();

    cell_qc
}

/////////////
// Helpers //
/////////////

/// Helper function that reads in full CSR data from an h5 file
///
/// ### Params
///
/// * `file_path` - Path to the h5ad file
///
/// ### Returns
///
/// The `CsEitherData` with the counts stored as u16
pub fn read_h5ad_x_data<P: AsRef<Path>>(
    file_path: P,
    cs_type: &str,
    shape: (usize, usize),
) -> Result<CompressedSparseData<u16>> {
    let cs_type = parse_cs_format(cs_type)
        .ok_or_else(|| format!("Invalid Compressed Sparse type: {}", cs_type))?;

    let file = File::open(file_path)?;
    let data_ds = file.dataset("X/data")?;
    let indices_ds = file.dataset("X/indices")?;
    let indptr_ds = file.dataset("X/indptr")?;

    // Read data as f32 and convert to u16
    let data_raw: Vec<f32> = data_ds.read_1d()?.to_vec();
    let data: Vec<u16> = data_raw.iter().map(|&x| x as u16).collect();

    // Read indices as u32 and convert to usize
    let indices_raw: Vec<u32> = indices_ds.read_1d()?.to_vec();
    let indices: Vec<usize> = indices_raw.iter().map(|&x| x as usize).collect();

    // Read indptr as u32 and convert to usize
    let indptr_raw: Vec<u32> = indptr_ds.read_1d()?.to_vec();
    let indptr: Vec<usize> = indptr_raw.iter().map(|&x| x as usize).collect();

    Ok(CompressedSparseData {
        data,
        indices,
        indptr,
        cs_type,
        data_2: None::<Vec<u16>>,
        shape,
    })
}

#[cfg(test)]
mod tests {
    use hdf5::File;

    #[test]
    fn test_read_h5ad_file() {
        let file_path = "/Users/gregorlueg/Downloads/ERX11148735.h5ad";

        // Try to open the file
        let file = match File::open(file_path) {
            Ok(f) => f,
            Err(e) => {
                println!("Failed to open file: {}", e);
                return;
            }
        };

        println!("File opened successfully");

        // Try to access datasets
        let data_ds = match file.dataset("X/data") {
            Ok(ds) => ds,
            Err(e) => {
                println!("Failed to get data dataset: {}", e);
                return;
            }
        };

        let indices_ds = match file.dataset("X/indices") {
            Ok(ds) => ds,
            Err(e) => {
                println!("Failed to get indices dataset: {}", e);
                return;
            }
        };

        let indptr_ds = match file.dataset("X/indptr") {
            Ok(ds) => ds,
            Err(e) => {
                println!("Failed to get indptr dataset: {}", e);
                return;
            }
        };

        println!("All datasets accessed successfully");

        // Print dataset info
        println!("Data shape: {:?}", data_ds.shape());
        println!("Data size: {:?}", data_ds.size());
        println!("Data dtype: {:?}", data_ds.dtype().unwrap().to_descriptor());

        println!("Indices shape: {:?}", indices_ds.shape());
        println!("Indices size: {:?}", indices_ds.size());
        println!(
            "Indices dtype: {:?}",
            indices_ds.dtype().unwrap().to_descriptor()
        );

        println!("Indptr shape: {:?}", indptr_ds.shape());
        println!("Indptr size: {:?}", indptr_ds.size());
        println!(
            "Indptr dtype: {:?}",
            indptr_ds.dtype().unwrap().to_descriptor()
        );

        // Try different reading methods
        println!("Trying to read data...");

        // Method 1: read_dyn
        match data_ds.read_dyn::<f32>() {
            Ok(arr) => println!("read_dyn() succeeded, shape: {:?}", arr.shape()),
            Err(e) => println!("read_dyn() failed: {}", e),
        }

        // Method 2: read_1d
        match data_ds.read_1d::<f32>() {
            Ok(arr) => println!("read_1d::<f32>() succeeded, len: {}", arr.len()),
            Err(e) => println!("read_1d::<f32>() failed: {}", e),
        }

        // Method 3: as_reader().read_1d()
        match data_ds.as_reader().read_1d::<f32>() {
            Ok(arr) => println!("as_reader().read_1d::<f32>() succeeded, len: {}", arr.len()),
            Err(e) => println!("as_reader().read_1d::<f32>() failed: {}", e),
        }

        // Method 4: read_raw
        match data_ds.read_raw::<f32>() {
            Ok(vec) => println!("read_raw::<f32>() succeeded, len: {}", vec.len()),
            Err(e) => println!("read_raw::<f32>() failed: {}", e),
        }
    }
}
