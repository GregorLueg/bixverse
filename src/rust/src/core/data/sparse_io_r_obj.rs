use std::ops::{Add, Mul};
use std::path::Path;
use thousands::Separable;

use crate::core::data::sparse_io::*;
use crate::core::data::sparse_structures::*;
use crate::single_cell::processing::*;

//////////
// Main //
//////////

/// Write R counts to binarised file
///
/// ### Params
///
/// * `bin_path` - Path to the h5 object.
/// * `compressed_data` - The original R data stored as a CompressedSparseData
///   structure.
/// * `cell_quality` - Structure containing information on the desired minimum
///   cell and gene quality + target size for library normalisation.
/// * `verbose` - Controls verbosity of the function
///
/// ### Returns
///
/// A tuple with `(no_cells, no_genes, cell quality metrics)`
pub fn write_r_counts<P: AsRef<Path>, T>(
    bin_path: P,
    compressed_data: CompressedSparseData<T>,
    cell_quality: MinCellQuality,
    verbose: bool,
) -> (usize, usize, CellQuality)
where
    T: Clone + Default + Into<u32> + Sync + Into<f64> + Add + PartialEq + Mul,
{
    let (no_cells, no_genes) = compressed_data.shape();

    if verbose {
        println!(
            "Processing R sparse matrix (shape: {} x {})...",
            no_cells.separate_with_underscores(),
            no_genes.separate_with_underscores()
        );
    }

    match compressed_data.cs_type {
        CompressedSparseFormat::Csr => {
            write_r_counts_csr(bin_path, compressed_data, cell_quality, verbose)
        }
        CompressedSparseFormat::Csc => {
            if verbose {
                println!("Converting CSC to CSR...");
            }
            let csr_data = csc_to_csr(&compressed_data);
            write_r_counts_csr(bin_path, csr_data, cell_quality, verbose)
        }
    }
}

/////////////
// Helpers //
/////////////

/// Write R counts to binarised file helper
///
/// ### Params
///
/// * `bin_path` - Path to the h5 object.
/// * `compressed_data` - The original R data stored as a CompressedSparseData
///   structure.
/// * `cell_quality` - Structure containing information on the desired minimum
///   cell and gene quality + target size for library normalisation.
/// * `verbose` - Controls verbosity of the function
///
/// ### Returns
///
/// A tuple with `(no_cells, no_genes, cell quality metrics)`
pub fn write_r_counts_csr<P: AsRef<Path>, T>(
    bin_path: P,
    compressed_data: CompressedSparseData<T>,
    cell_quality: MinCellQuality,
    verbose: bool,
) -> (usize, usize, CellQuality)
where
    T: Clone + Default + Into<u32> + Sync + Into<f64> + Add + PartialEq + Mul,
{
    let (no_cells, no_genes) = compressed_data.shape();

    if verbose {
        println!("Generating cell chunks...");
    }

    let (cell_chunk_vec, cell_qc): (Vec<CsrCellChunk>, CellQuality) =
        CsrCellChunk::generate_chunks_sparse_data(compressed_data, cell_quality);

    let cells_passing = cell_qc.cell_indices.len();
    let genes_passing = cell_qc.gene_indices.len();

    if verbose {
        println!(
            "Cells passing QC: {} / {}",
            cells_passing.separate_with_underscores(),
            no_cells.separate_with_underscores()
        );
        println!(
            "Genes passing QC: {} / {}",
            genes_passing.separate_with_underscores(),
            no_genes.separate_with_underscores()
        );
        println!("Writing to binary format...");
    }

    let mut writer =
        CellGeneSparseWriter::new(bin_path, true, cells_passing, genes_passing).unwrap();

    // Filter to passing cells and remap indices to sequential
    let mut passing_chunks: Vec<_> = cell_chunk_vec
        .into_iter()
        .filter(|chunk| chunk.to_keep)
        .collect();

    // CRITICAL: Remap original_index to sequential 0, 1, 2, ... for correct transpose
    for (new_idx, chunk) in passing_chunks.iter_mut().enumerate() {
        chunk.original_index = new_idx;
    }

    for (i, cell_chunk) in passing_chunks.into_iter().enumerate() {
        writer.write_cell_chunk(cell_chunk).unwrap();

        if verbose && (i + 1) % 100000 == 0 {
            println!(
                "  Written {} / {} cells to disk.",
                (i + 1).separate_with_underscores(),
                cells_passing.separate_with_underscores()
            );
        }
    }

    if verbose {
        println!(
            "  Written {} / {} cells (complete).",
            cells_passing.separate_with_underscores(),
            cells_passing.separate_with_underscores()
        );
        println!("Finalising file...");
    }

    writer.finalise().unwrap();

    // Return cell_qc as-is - it already contains correct original indices
    // from generate_chunks_sparse_data
    (cells_passing, genes_passing, cell_qc)
}
