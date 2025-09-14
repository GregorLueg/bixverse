use hdf5::{File, Result};

use crate::core::data::sparse_io::*;
use crate::core::data::sparse_structures::*;

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

pub fn write_h5_counts(
    h5_path: &str,
    bin_path: &str,
    cs_type: &str,
    no_cells: usize,
    no_genes: usize,
    min_genes: usize,
    size_factor: f32,
) -> Vec<bool> {
    let file_data: CompressedSparseData<u16> =
        read_h5ad_x_data(h5_path, cs_type, (no_genes, no_cells)).unwrap();

    println!("Data length: {:?}", file_data.data.len());
    println!("Indptr length: {:?}", file_data.indptr.len());

    let file_data = if file_data.cs_type.is_csc() {
        file_data.transform()
    } else {
        file_data
    };

    println!("Data length: {:?}", file_data.data.len());
    println!("Indptr length: {:?}", file_data.indptr.len());

    let mut writer = CellGeneSparseWriter::new(bin_path, true, no_cells, no_genes).unwrap();

    let (cell_chunk_vec, to_keep) =
        CsrCellChunk::generate_chunks_sparse_data(file_data, min_genes, size_factor);

    for cell_chunk in cell_chunk_vec {
        writer.write_cell_chunk(cell_chunk).unwrap();
    }

    writer.finalise().unwrap();

    to_keep
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
pub fn read_h5ad_x_data(
    file_path: &str,
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
