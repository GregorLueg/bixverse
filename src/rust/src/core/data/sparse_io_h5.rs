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
/// A tuple with `(no_cells, no_genes, cell quality metrics)`
pub fn write_h5_counts<P: AsRef<Path>>(
    h5_path: P,
    bin_path: P,
    cs_type: &str,
    no_cells: usize,
    no_genes: usize,
    cell_quality: MinCellQuality,
    verbose: bool,
) -> (usize, usize, CellQuality) {
    let file_quality = parse_h5_csr_quality(&h5_path, (no_genes, no_cells), &cell_quality).unwrap();

    if verbose {
        println!(
            "A total of {:?} genes pass the threshold",
            file_quality.genes_to_keep.len()
        );
        println!(
            "A total of {:?} cells pass the threshold",
            file_quality.cells_to_keep.len()
        );
    }

    let file_data: CompressedSparseData<u16> =
        read_h5ad_x_data(&h5_path, cs_type, &file_quality).unwrap();

    let file_data = file_data.transpose_and_convert();

    let mut writer = CellGeneSparseWriter::new(bin_path, true, no_cells, no_genes).unwrap();

    let (cell_chunk_vec, cell_qc) =
        CsrCellChunk::generate_chunks_sparse_data(file_data, cell_quality);

    for cell_chunk in cell_chunk_vec {
        writer.write_cell_chunk(cell_chunk).unwrap();
    }

    // update the header with the actual files written
    writer.update_header_no_cells(file_quality.cells_to_keep.len());
    writer.update_header_no_genes(file_quality.genes_to_keep.len());
    writer.finalise().unwrap();

    (
        file_quality.cells_to_keep.len(),
        file_quality.genes_to_keep.len(),
        cell_qc,
    )
}

/////////////
// Helpers //
/////////////

/// Helper function that reads in full CSR data from an h5 file
///
/// The function assumes that the data is stored as genes x cells.
///
/// ### Params
///
/// * `file_path` - Path to the h5ad file.
/// * `cs_type` - Which type of sparse compression format the file is stored in.
/// * `shape` - The final dimension of the matrix.
///
/// ### Returns
///
/// The `CsEitherData` with the counts stored as u16
pub fn read_h5ad_x_data<P: AsRef<Path>>(
    file_path: P,
    cs_type: &str,
    quality: &CellOnFileQuality,
) -> Result<CompressedSparseData<u16>> {
    let cs_type = parse_cs_format(cs_type)
        .ok_or_else(|| format!("Invalid Compressed Sparse type: {}", cs_type))?;
    let file = File::open(file_path)?;

    let data_ds = file.dataset("X/data")?;
    let indices_ds = file.dataset("X/indices")?;
    let indptr_ds = file.dataset("X/indptr")?;

    let data_raw: Vec<f32> = data_ds.read_1d()?.to_vec();
    let indices_raw: Vec<u32> = indices_ds.read_1d()?.to_vec();
    let indptr_raw: Vec<u32> = indptr_ds.read_1d()?.to_vec();

    let mut new_data: Vec<u16> = Vec::new();
    let mut new_indices: Vec<usize> = Vec::new();
    let mut new_indptr: Vec<usize> = Vec::with_capacity(quality.genes_to_keep.len() + 1);
    new_indptr.push(0);

    for &gene_idx in &quality.genes_to_keep {
        let start = indptr_raw[gene_idx] as usize;
        let end = indptr_raw[gene_idx + 1] as usize;

        for idx in start..end {
            let old_cell_idx = indices_raw[idx] as usize;
            if let Some(&new_cell_idx) = quality.cell_old_to_new.get(&old_cell_idx) {
                new_data.push(data_raw[idx] as u16);
                new_indices.push(new_cell_idx);
            }
        }
        new_indptr.push(new_data.len());
    }

    let shape = (quality.genes_to_keep.len(), quality.cells_to_keep.len());

    Ok(CompressedSparseData {
        data: new_data,
        indices: new_indices,
        indptr: new_indptr,
        cs_type,
        data_2: None::<Vec<u16>>,
        shape,
    })
}

/// Get the cell quality data from a CSR file
///
/// This file assumes that the rows are representing the genes and the columns
/// the cells and the data was stored in CSR type format.
///
/// ### Params
///
/// * `file_path` - Path to the h5ad file
/// * `cell_quality` -
pub fn parse_h5_csr_quality<P: AsRef<Path>>(
    file_path: P,
    shape: (usize, usize),
    cell_quality: &MinCellQuality,
) -> Result<CellOnFileQuality> {
    let file = File::open(file_path)?;
    let data_ds = file.dataset("X/data")?;
    let indices_ds = file.dataset("X/indices")?;
    let indptr_ds = file.dataset("X/indptr")?;

    // First parse over the indptr to understand in how many cells the gene is
    // expressed
    let mut no_cells_exp_gene: Vec<usize> = Vec::with_capacity(shape.0);

    let data: Vec<f32> = data_ds.read_1d()?.to_vec();
    let indices: Vec<u32> = indices_ds.read_1d()?.to_vec();
    let indptr: Vec<u32> = indptr_ds.read_1d()?.to_vec();

    for i in 0..shape.0 {
        let start_ptr_i = indptr[i] as usize;
        let end_ptr_i = indptr[i + 1] as usize;

        no_cells_exp_gene.push(end_ptr_i - start_ptr_i);
    }

    let genes_to_keep_bool: Vec<bool> = no_cells_exp_gene
        .iter()
        .map(|x| *x >= cell_quality.min_cells)
        .collect();

    // Calculate the cell metrics using the kept genes
    // Bit inefficient here because we data is in CSR, but oh well...
    let mut cell_unique_genes = vec![0usize; shape.1];
    let mut cell_lib_size = vec![0.0f32; shape.1];

    for gene_idx in 0..shape.0 {
        if !genes_to_keep_bool[gene_idx] {
            continue;
        }

        let start = indptr[gene_idx] as usize;
        let end = indptr[gene_idx + 1] as usize;

        for idx in start..end {
            let cell_idx = indices[idx] as usize;
            cell_unique_genes[cell_idx] += 1;
            cell_lib_size[cell_idx] += data[idx];
        }
    }

    let cells_to_keep: Vec<usize> = (0..shape.1)
        .filter(|&i| {
            cell_unique_genes[i] >= cell_quality.min_unique_genes
                && cell_lib_size[i] >= cell_quality.min_lib_size as f32
        })
        .collect();

    let genes_to_keep: Vec<usize> = (0..shape.0).filter(|&i| genes_to_keep_bool[i]).collect();

    let mut file_quality_data = CellOnFileQuality::new(cells_to_keep, genes_to_keep);
    file_quality_data.generate_maps_sets();

    Ok(file_quality_data)
}

///////////
// Tests //
///////////

#[cfg(test)]
mod tests {
    use hdf5::File;

    #[test]
    fn test_read_h5ad_file() {
        // TODO - need to generate testing data and add to the repo itself!

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
