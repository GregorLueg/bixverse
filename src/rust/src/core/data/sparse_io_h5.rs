use hdf5::{File, Result};
use std::io::Result as IoResult;
use std::path::Path;
use thousands::Separable;

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
/// Function to take in h5 data and write it to disk (first the cells) in a
/// binarised format for fast retrieval of cells. Pending on how the file was
/// stored (CSC or CSR) different paths will be used to process the data.
///
/// ### Params
///
/// * `h5_path` - Path to the h5 object.
/// * `bin_path` - Path to the binarised object on disk to write to
/// * `cs_type` - Was the h5ad data stored in "csc" or "csr". Important! h5ad
///   stores data in genes x cells; bixverse stores in cells x genes!
/// * `no_cells` - Total number of obversations in the data.
/// * `no_genes` - Total number of vars in the data.
/// * `cell_quality` - Structure containing information on the desired minimum
///   cell and gene quality + target size for library normalisation.
/// * `verbose` - Controls verbosity of the function.
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
    if verbose {
        println!("Step 1/4: Analysing data structure, calculating QC metrics and identifying cells/genes to take...");
    }

    let file_format = parse_cs_format(cs_type).unwrap();

    let file_quality = match file_format {
        CompressedSparseFormat::Csr => {
            parse_h5_csr_quality(&h5_path, (no_cells, no_genes), &cell_quality, verbose).unwrap()
        }
        CompressedSparseFormat::Csc => {
            parse_h5_csc_quality(&h5_path, (no_cells, no_genes), &cell_quality, verbose).unwrap()
        }
    };

    if verbose {
        println!("Step 2/4: QC Results:");
        println!(
            "  Genes passing QC (i.e., getting included): {} / {}",
            file_quality.genes_to_keep.len().separate_with_underscores(),
            no_genes.separate_with_underscores()
        );
        println!(
            "  Cells passing QC (i.e., getting included): {} / {}",
            file_quality.cells_to_keep.len().separate_with_underscores(),
            no_cells.separate_with_underscores()
        );
        println!("Step 3/4: Loading filtered data from h5...");
    }

    let file_data: CompressedSparseData<u16> = match file_format {
        CompressedSparseFormat::Csr => {
            read_h5ad_x_data_csr(&h5_path, &file_quality, verbose).unwrap()
        }
        CompressedSparseFormat::Csc => {
            let data = read_h5ad_x_data_csc(&h5_path, &file_quality, verbose).unwrap();
            data.transpose_and_convert()
        }
    };

    if verbose {
        println!("Step 4/4: Writing to binary format...");
    }

    let mut writer = CellGeneSparseWriter::new(bin_path, true, no_cells, no_genes).unwrap();

    let (cell_chunk_vec, mut cell_qc): (Vec<CsrCellChunk>, CellQuality) =
        CsrCellChunk::generate_chunks_sparse_data(file_data, cell_quality);

    let total_chunks = cell_chunk_vec.len();
    for (i, cell_chunk) in cell_chunk_vec.into_iter().enumerate() {
        writer.write_cell_chunk(cell_chunk).unwrap();

        if verbose && (i + 1) % 100000 == 0 {
            println!(
                "  Written {} / {} cells to disk.",
                (i + 1).separate_with_underscores(),
                total_chunks.separate_with_underscores()
            );
        }
    }

    if verbose {
        println!(
            "  Written {} / {} cells (complete).",
            total_chunks.separate_with_underscores(),
            total_chunks.separate_with_underscores()
        );
        println!("Finalising file...");
    }

    writer.update_header_no_cells(file_quality.cells_to_keep.len());
    writer.update_header_no_genes(file_quality.genes_to_keep.len());
    writer.finalise().unwrap();

    cell_qc.set_cell_indices(&file_quality.cells_to_keep);
    cell_qc.set_gene_indices(&file_quality.genes_to_keep);

    (
        file_quality.cells_to_keep.len(),
        file_quality.genes_to_keep.len(),
        cell_qc,
    )
}

/// Function that streams the h5 counts to disk
///
/// This function will avoid (as far as possible) loading in the data into
/// memory and leverage direct streaming to disk where possible. Works best
/// with data that is already stored in CSR on disk (assuming genes * cells).
///
/// ### Params
///
/// * `h5_path` - Path to the h5 object.
/// * `bin_path` - Path to the binarised object on disk to write to
/// * `cs_type` - Was the h5ad data stored in "csc" or "csr". Important! h5ad
///   stores data in genes x cells; bixverse stores in cells x genes!
/// * `no_cells` - Total number of obversations in the data.
/// * `no_genes` - Total number of vars in the data.
/// * `cell_quality` - Structure containing information on the desired minimum
///   cell and gene quality + target size for library normalisation.
/// * `verbose` - Controls verbosity of the function.
///
/// ### Returns
///
/// A tuple with `(no_cells, no_genes, cell quality metrics)`
pub fn stream_h5_counts<P: AsRef<Path>>(
    h5_path: P,
    bin_path: P,
    cs_type: &str,
    no_cells: usize,
    no_genes: usize,
    cell_quality: MinCellQuality,
    verbose: bool,
) -> (usize, usize, CellQuality) {
    if verbose {
        println!("Step 1/3: Analysing data structure and calculating QC metrics...");
    }

    let file_format = parse_cs_format(cs_type).unwrap();

    let file_quality = match file_format {
        CompressedSparseFormat::Csr => {
            parse_h5_csr_quality(&h5_path, (no_cells, no_genes), &cell_quality, verbose).unwrap()
        }
        CompressedSparseFormat::Csc => {
            parse_h5_csc_quality(&h5_path, (no_cells, no_genes), &cell_quality, verbose).unwrap()
        }
    };

    if verbose {
        println!("Step 2/3: QC Results:");
        println!(
            "  Genes passing QC: {} / {}",
            file_quality.genes_to_keep.len().separate_with_underscores(),
            no_genes.separate_with_underscores()
        );
        println!(
            "  Cells passing QC: {} / {}",
            file_quality.cells_to_keep.len().separate_with_underscores(),
            no_cells.separate_with_underscores()
        );
        println!("Step 3/3: Writing cells to CSR format...");
    }

    let mut cell_qc = match file_format {
        CompressedSparseFormat::Csr => {
            write_h5_csr_streaming(&h5_path, &bin_path, &file_quality, cell_quality, verbose)
                .unwrap()
        }
        CompressedSparseFormat::Csc => {
            if verbose {
                println!("  Pass 1/2: Scanning library sizes...");
            }
            let cell_lib_sizes = scan_h5_csc_library_sizes(&h5_path, &file_quality).unwrap();

            if verbose {
                println!("  Pass 2/2: Writing cells with normalization...");
            }
            write_h5_csc_to_csr_streaming(
                &h5_path,
                &bin_path,
                &file_quality,
                &cell_lib_sizes,
                cell_quality.target_size,
                verbose,
            )
            .unwrap()
        }
    };

    cell_qc.set_cell_indices(&file_quality.cells_to_keep);
    cell_qc.set_gene_indices(&file_quality.genes_to_keep);

    (
        file_quality.cells_to_keep.len(),
        file_quality.genes_to_keep.len(),
        cell_qc,
    )
}

//////////////
// CSC data //
//////////////

/// Get the cell quality data from a CSC file
///
/// This file assumes that the rows are representing cells and the columns
/// represent genes (typical for h5ad).
///
/// ### Params
///
/// * `file_path` - Path to the h5ad file
/// * `shape` - Tuple with `(no_cells, no_genes)`.
/// * `cell_quality` - Structure defining the minimum quality values that are
///   expected here.
/// * `verbose` - Controls verbosity of the function.
///
/// ### Returns
///
/// `CellOnFileQuality` structure that contains all of the information about
/// which cells and genes to include.
pub fn parse_h5_csc_quality<P: AsRef<Path>>(
    file_path: P,
    shape: (usize, usize),
    cell_quality: &MinCellQuality,
    verbose: bool,
) -> Result<CellOnFileQuality> {
    let file = File::open(file_path)?;
    let data_ds = file.dataset("X/data")?;
    let indices_ds = file.dataset("X/indices")?;
    let indptr_ds = file.dataset("X/indptr")?;

    if verbose {
        println!(
            "  Reading CSC matrix structure (shape: {} x {})...",
            shape.0.separate_with_underscores(),
            shape.1.separate_with_underscores()
        );
        println!(
            "  Data size: {}, Indices size: {}, Indptr size: {}",
            data_ds.size().separate_with_underscores(),
            indices_ds.size().separate_with_underscores(),
            indptr_ds.size().separate_with_underscores()
        );
    }

    let indptr: Vec<u32> = indptr_ds.read_1d()?.to_vec();

    if (indptr.len() != shape.1 + 1) && verbose {
        println!(
            "  WARNING: indptr length {} doesn't match expected {} (shape.1 + 1)",
            indptr.len().separate_with_underscores(),
            (shape.1 + 1).separate_with_underscores()
        );
    }

    // Count cells per gene
    let mut no_cells_exp_gene: Vec<usize> = Vec::with_capacity(shape.1);
    for i in 0..shape.1 {
        let start_ptr = indptr[i] as usize;
        let end_ptr = indptr[i + 1] as usize;
        no_cells_exp_gene.push(end_ptr - start_ptr);
    }

    if verbose {
        let max_expr = no_cells_exp_gene.iter().max().unwrap_or(&0);
        let min_expr = no_cells_exp_gene.iter().min().unwrap_or(&0);
        let avg_expr = if shape.1 > 0 {
            no_cells_exp_gene.iter().sum::<usize>() / shape.1
        } else {
            0
        };
        println!(
            "  Gene expression stats: min = {} | max = {} | avg = {} cells per gene",
            min_expr.separate_with_underscores(),
            max_expr.separate_with_underscores(),
            avg_expr.separate_with_underscores()
        );
    }

    // Filter genes first
    let genes_to_keep: Vec<usize> = (0..shape.1)
        .filter(|&i| no_cells_exp_gene[i] >= cell_quality.min_cells)
        .collect();

    if verbose {
        println!(
            "  Genes passing filter: {} / {}",
            genes_to_keep.len().separate_with_underscores(),
            shape.1.separate_with_underscores()
        );
    }

    // Calculate cell metrics using only kept genes
    let mut cell_unique_genes = vec![0usize; shape.0];
    let mut cell_lib_size = vec![0.0f32; shape.0];

    if verbose {
        println!("  Calculating cell metrics in chunks...");
    }

    const GENE_CHUNK_SIZE: usize = 10000;

    for (chunk_idx, gene_chunk) in genes_to_keep.chunks(GENE_CHUNK_SIZE).enumerate() {
        if verbose && chunk_idx % 10 == 0 {
            let processed = chunk_idx * GENE_CHUNK_SIZE;
            println!(
                "   Processing genes {} / {}",
                processed
                    .min(genes_to_keep.len())
                    .separate_with_underscores(),
                genes_to_keep.len().separate_with_underscores()
            );
        }

        let chunk_start_gene = gene_chunk[0];
        let chunk_end_gene = gene_chunk[gene_chunk.len() - 1];

        let data_start = indptr[chunk_start_gene] as usize;
        let data_end = indptr[chunk_end_gene + 1] as usize;

        if data_start >= data_end {
            continue;
        }

        let chunk_data: Vec<f32> = data_ds.read_slice_1d(data_start..data_end)?.to_vec();
        let chunk_indices: Vec<u32> = indices_ds.read_slice_1d(data_start..data_end)?.to_vec();

        for &gene_idx in gene_chunk {
            let gene_data_start = indptr[gene_idx] as usize;
            let gene_data_end = indptr[gene_idx + 1] as usize;

            for idx in gene_data_start..gene_data_end {
                let local_idx = idx - data_start;
                let cell_idx = chunk_indices[local_idx] as usize;

                if cell_idx < shape.0 {
                    cell_unique_genes[cell_idx] += 1;
                    cell_lib_size[cell_idx] += chunk_data[local_idx];
                }
            }
        }
    }

    if verbose {
        println!(
            "   Processing genes {} / {} (complete)",
            genes_to_keep.len().separate_with_underscores(),
            genes_to_keep.len().separate_with_underscores()
        );

        let max_genes = cell_unique_genes.iter().max().unwrap_or(&0);
        let min_genes = cell_unique_genes.iter().min().unwrap_or(&0);
        let max_lib = cell_lib_size.iter().fold(0.0f32, |a, &b| a.max(b));
        let min_lib = cell_lib_size.iter().fold(f32::INFINITY, |a, &b| a.min(b));

        println!(
            "  Cell stats: genes per cell: min = {} | max={}",
            min_genes.separate_with_underscores(),
            max_genes.separate_with_underscores()
        );
        println!(
            "  Cell stats: library size: min = {:.1} | max={:.1}",
            min_lib.separate_with_underscores(),
            max_lib.separate_with_underscores()
        );
    }

    let cells_to_keep: Vec<usize> = (0..shape.0)
        .filter(|&i| {
            cell_unique_genes[i] >= cell_quality.min_unique_genes
                && cell_lib_size[i] >= cell_quality.min_lib_size as f32
        })
        .collect();

    let mut file_quality_data = CellOnFileQuality::new(cells_to_keep, genes_to_keep);
    file_quality_data.generate_maps_sets();

    Ok(file_quality_data)
}

/// Helper function to scan a CSC file and get the lib sizes
///
/// ### Params
///
/// * `file_path` - Path to the h5ad file.
/// * `quality` - Information on the which cells and genes to include after a
///   first pass of the file.
///
/// ### Returns
///
/// A vector of library sizes in the cells.
pub fn scan_h5_csc_library_sizes<P: AsRef<Path>>(
    file_path: P,
    quality: &CellOnFileQuality,
) -> Result<Vec<u32>> {
    let file = File::open(file_path)?;
    let data_ds = file.dataset("X/data")?;
    let indices_ds = file.dataset("X/indices")?;
    let indptr_ds = file.dataset("X/indptr")?;
    let indptr_raw: Vec<u32> = indptr_ds.read_1d()?.to_vec();

    // Only track library sizes - very light in terms of memory
    let mut cell_lib_sizes = vec![0u32; quality.cells_to_keep.len()];

    const GENE_CHUNK_SIZE: usize = 5000;

    for gene_chunk in quality.genes_to_keep.chunks(GENE_CHUNK_SIZE) {
        let start_pos = gene_chunk
            .iter()
            .map(|&g| indptr_raw[g] as usize)
            .min()
            .unwrap_or(0);
        let end_pos = gene_chunk
            .iter()
            .map(|&g| indptr_raw[g + 1] as usize)
            .max()
            .unwrap_or(0);

        if start_pos >= end_pos {
            continue;
        }

        // only read data values, not indices - saves memory
        let chunk_data: Vec<f32> = data_ds.read_slice_1d(start_pos..end_pos)?.to_vec();
        let chunk_indices: Vec<u32> = indices_ds.read_slice_1d(start_pos..end_pos)?.to_vec();

        for &gene_idx in gene_chunk {
            let gene_start = indptr_raw[gene_idx] as usize;
            let gene_end = indptr_raw[gene_idx + 1] as usize;

            for idx in gene_start..gene_end {
                let local_idx = idx - start_pos;
                let old_cell_idx = chunk_indices[local_idx] as usize;

                if let Some(&new_cell_idx) = quality.cell_old_to_new.get(&old_cell_idx) {
                    cell_lib_sizes[new_cell_idx] += chunk_data[local_idx] as u32;
                }
            }
        }
    }

    Ok(cell_lib_sizes)
}

/// Helper function that reads in full CSC data from an h5 file
///
/// The function assumes that the data is stored as cells x genes.
///
/// ### Params
///
/// * `file_path` - Path to the h5ad file.
/// * `quality` - Information on the which cells and genes to include after a
///   first pass of the file.
/// * `verbose` - Controls verbosity of the function.
///
/// ### Returns
///
/// The `CompressedSparseData` in CSR format with the counts stored as u16.
pub fn read_h5ad_x_data_csc<P: AsRef<Path>>(
    file_path: P,
    quality: &CellOnFileQuality,
    verbose: bool,
) -> Result<CompressedSparseData<u16>> {
    let file = File::open(file_path)?;

    let data_ds = file.dataset("X/data")?;
    let indices_ds = file.dataset("X/indices")?;
    let indptr_ds = file.dataset("X/indptr")?;

    // Read indptr first (small array)
    let indptr_raw: Vec<u32> = indptr_ds.read_1d()?.to_vec();

    let mut new_data: Vec<u16> = Vec::new();
    let mut new_indices: Vec<usize> = Vec::new();
    let mut new_indptr: Vec<usize> = Vec::with_capacity(quality.genes_to_keep.len() + 1);
    new_indptr.push(0);

    let total_genes = quality.genes_to_keep.len();

    if verbose {
        println!(
            "  Processing {} genes in chunks...",
            total_genes.separate_with_underscores()
        );
    }

    // Process genes in chunks to reduce memory usage
    // should think about exposing this as a function parameter... ?
    const CHUNK_SIZE: usize = 1000;

    for (chunk_idx, gene_chunk) in quality.genes_to_keep.chunks(CHUNK_SIZE).enumerate() {
        if verbose && chunk_idx % 10 == 0 {
            let processed = chunk_idx * CHUNK_SIZE;
            println!(
                "   Processed {} / {} genes",
                processed.min(total_genes).separate_with_underscores(),
                total_genes.separate_with_underscores()
            );
        }

        // calculate range for this chunk
        let start_pos = gene_chunk
            .iter()
            .map(|&g| indptr_raw[g] as usize)
            .min()
            .unwrap_or(0);
        let end_pos = gene_chunk
            .iter()
            .map(|&g| indptr_raw[g + 1] as usize)
            .max()
            .unwrap_or(0);

        if start_pos >= end_pos {
            continue;
        }

        // read only the data range needed for this chunk
        let chunk_data: Vec<f32> = data_ds.read_slice_1d(start_pos..end_pos)?.to_vec();
        let chunk_indices: Vec<u32> = indices_ds.read_slice_1d(start_pos..end_pos)?.to_vec();

        // process each gene in the chunk
        for &gene_idx in gene_chunk {
            let gene_start = indptr_raw[gene_idx] as usize;
            let gene_end = indptr_raw[gene_idx + 1] as usize;

            for idx in gene_start..gene_end {
                let local_idx = idx - start_pos;
                let old_cell_idx = chunk_indices[local_idx] as usize;

                if let Some(&new_cell_idx) = quality.cell_old_to_new.get(&old_cell_idx) {
                    new_data.push(chunk_data[local_idx] as u16);
                    new_indices.push(new_cell_idx);
                }
            }
            new_indptr.push(new_data.len());
        }
    }

    if verbose {
        println!(
            "   Processed {} / {} genes (complete)",
            total_genes.separate_with_underscores(),
            total_genes.separate_with_underscores()
        );
    }

    let shape = (quality.genes_to_keep.len(), quality.cells_to_keep.len());

    Ok(CompressedSparseData {
        data: new_data,
        indices: new_indices,
        indptr: new_indptr,
        cs_type: CompressedSparseFormat::Csr,
        data_2: None::<Vec<u16>>,
        shape,
    })
}

/// Helper function that streams CSC data to CSR format on disk
///
/// ### Params
///
/// * `file_path` - Path to the h5ad file.
/// * `bin_path` - Path to the binary file on disk to write to.
/// * `quality` - Information on the which cells and genes to include after a
///   first pass of the file.
/// * `cell_lib_sizes` - Vector with the pre-calculated library sizes, see
///   scan_h5_csc_library_sizes().
/// * `target_size` - Target size for the library normalisation.
///
/// ### Returns
///
/// After writing, information on the CellQuality with cell and gene indices,
/// NNZ and lib size per gene.
pub fn write_h5_csc_to_csr_streaming<P: AsRef<Path>>(
    file_path: P,
    bin_path: P,
    quality: &CellOnFileQuality,
    cell_lib_sizes: &[u32],
    target_size: f32,
    verbose: bool,
) -> IoResult<CellQuality> {
    let file = File::open(file_path)?;
    let data_ds = file.dataset("X/data")?;
    let indices_ds = file.dataset("X/indices")?;
    let indptr_ds = file.dataset("X/indptr")?;
    let indptr_raw: Vec<u32> = indptr_ds.read_1d()?.to_vec();

    // accumulate cells in memory (necessary for CSR)
    let mut cell_data: Vec<Vec<(u16, u16)>> = vec![Vec::new(); quality.cells_to_keep.len()];

    const GENE_CHUNK_SIZE: usize = 5000;
    let total_genes = quality.genes_to_keep.len();

    if verbose {
        println!(
            "  Processing {} genes in chunks...",
            total_genes.separate_with_underscores()
        );
    }

    for (chunk_idx, gene_chunk) in quality.genes_to_keep.chunks(GENE_CHUNK_SIZE).enumerate() {
        if verbose && chunk_idx % 10 == 0 {
            let processed = chunk_idx * GENE_CHUNK_SIZE;
            println!(
                "   Processed {} / {} genes",
                processed.min(total_genes).separate_with_underscores(),
                total_genes.separate_with_underscores()
            );
        }

        let start_pos = gene_chunk
            .iter()
            .map(|&g| indptr_raw[g] as usize)
            .min()
            .unwrap_or(0);
        let end_pos = gene_chunk
            .iter()
            .map(|&g| indptr_raw[g + 1] as usize)
            .max()
            .unwrap_or(0);

        if start_pos >= end_pos {
            continue;
        }

        let chunk_data: Vec<f32> = data_ds.read_slice_1d(start_pos..end_pos)?.to_vec();
        let chunk_indices: Vec<u32> = indices_ds.read_slice_1d(start_pos..end_pos)?.to_vec();

        for &gene_idx in gene_chunk {
            let new_gene_idx = quality.gene_old_to_new[&gene_idx] as u16;
            let gene_start = indptr_raw[gene_idx] as usize;
            let gene_end = indptr_raw[gene_idx + 1] as usize;

            for idx in gene_start..gene_end {
                let local_idx = idx - start_pos;
                let old_cell_idx = chunk_indices[local_idx] as usize;

                if let Some(&new_cell_idx) = quality.cell_old_to_new.get(&old_cell_idx) {
                    let raw_count = chunk_data[local_idx] as u16;
                    cell_data[new_cell_idx].push((new_gene_idx, raw_count));
                }
            }
        }
    }

    if verbose {
        println!(
            "   Processed {} / {} genes (complete)",
            total_genes.separate_with_underscores(),
            total_genes.separate_with_underscores()
        );
        println!(
            "  Writing {} cells to disk...",
            cell_lib_sizes.len().separate_with_underscores()
        );
    }

    // now write cells with correct normalisation
    let mut writer = CellGeneSparseWriter::new(
        bin_path,
        true,
        cell_lib_sizes.len(),
        quality.genes_to_keep.len(),
    )?;
    let mut lib_size = Vec::with_capacity(cell_lib_sizes.len());
    let mut nnz = Vec::with_capacity(cell_lib_sizes.len());

    let total_cells = cell_data.len();

    for (cell_idx, mut data) in cell_data.into_iter().enumerate() {
        if verbose && (cell_idx + 1) % 100000 == 0 {
            println!(
                "   Written {} / {} cells to disk.",
                (cell_idx + 1).separate_with_underscores(),
                total_cells.separate_with_underscores()
            );
        }

        data.sort_by_key(|(gene_idx, _)| *gene_idx);

        let gene_indices: Vec<u16> = data.iter().map(|(g, _)| *g).collect();
        let gene_counts: Vec<u16> = data.iter().map(|(_, c)| *c).collect();

        let cell_chunk =
            CsrCellChunk::from_data(&gene_counts, &gene_indices, cell_idx, target_size, true);

        let (nnz_i, lib_size_i) = cell_chunk.get_qc_info();
        nnz.push(nnz_i);
        lib_size.push(lib_size_i);

        writer.write_cell_chunk(cell_chunk)?;
    }

    writer.finalise()?;

    if verbose {
        println!(
            "   Written {} / {} cells (complete).",
            total_cells.separate_with_underscores(),
            total_cells.separate_with_underscores()
        );
    }

    Ok(CellQuality {
        cell_indices: Vec::new(),
        gene_indices: Vec::new(),
        lib_size,
        no_genes: nnz,
    })
}

/////////
// CSR //
/////////

/// Get the cell quality data from a CSR file
///
/// This file assumes that the rows are representing the cells and the columns
/// the genes and the data was stored in CSR type format.
///
/// ### Params
///
/// * `file_path` - Path to the h5ad file
/// * `shape` - Tuple with `(no_cells, no_genes)`.
/// * `cell_quality` - Structure defining the minimum quality values that are
///   expected here.
/// * `verbose` - Controls verbosity of the function.
///
/// ### Returns
///
/// `CellOnFileQuality` structure that contains all of the information about
/// which cells and genes to include.
pub fn parse_h5_csr_quality<P: AsRef<Path>>(
    file_path: P,
    shape: (usize, usize),
    cell_quality: &MinCellQuality,
    verbose: bool,
) -> Result<CellOnFileQuality> {
    let file = File::open(file_path)?;
    let data_ds = file.dataset("X/data")?;
    let indices_ds = file.dataset("X/indices")?;
    let indptr_ds = file.dataset("X/indptr")?;

    if verbose {
        println!(
            "  Reading CSR matrix structure (shape: {} x {} )...",
            shape.0.separate_with_underscores(),
            shape.1.separate_with_underscores()
        );
        println!(
            "  Data size: {}, Indices size: {}, Indptr size: {}",
            data_ds.size().separate_with_underscores(),
            indices_ds.size().separate_with_underscores(),
            indptr_ds.size().separate_with_underscores()
        );
    }

    let indptr: Vec<u32> = indptr_ds.read_1d()?.to_vec();

    if (indptr.len() != shape.0 + 1) && verbose {
        println!(
            "  WARNING: indptr length {} doesn't match expected {} (no_cells + 1)",
            indptr.len().separate_with_underscores(),
            (shape.0 + 1).separate_with_underscores()
        );
    }

    // first pass - count how many cells express each gene
    let mut no_cells_exp_gene = vec![0usize; shape.1];

    if verbose {
        println!("  Pass 1: Calculating gene expression in chunks...");
    }

    const CELL_CHUNK_SIZE: usize = 10000;

    for chunk_start_cell in (0..shape.0).step_by(CELL_CHUNK_SIZE) {
        let chunk_end_cell = (chunk_start_cell + CELL_CHUNK_SIZE).min(shape.0) - 1;

        if verbose && (chunk_start_cell / CELL_CHUNK_SIZE) % 10 == 0 {
            println!(
                "   Processing cells {} / {}",
                chunk_start_cell.separate_with_underscores(),
                shape.0.separate_with_underscores()
            );
        }

        let data_start = indptr[chunk_start_cell] as usize;
        let data_end = indptr[chunk_end_cell + 1] as usize;

        if data_start >= data_end {
            continue;
        }

        let chunk_indices: Vec<u32> = indices_ds.read_slice_1d(data_start..data_end)?.to_vec();

        for cell_idx in chunk_start_cell..=chunk_end_cell {
            let cell_data_start = indptr[cell_idx] as usize;
            let cell_data_end = indptr[cell_idx + 1] as usize;

            let local_start = cell_data_start - data_start;
            let local_end = cell_data_end - data_start;

            #[allow(clippy::needless_range_loop)]
            for local_idx in local_start..local_end {
                let gene_idx = chunk_indices[local_idx] as usize;
                if gene_idx < shape.1 {
                    no_cells_exp_gene[gene_idx] += 1;
                }
            }
        }
    }

    if verbose {
        let max_expr = no_cells_exp_gene.iter().max().unwrap_or(&0);
        let min_expr = no_cells_exp_gene.iter().min().unwrap_or(&0);
        let avg_expr = if shape.1 > 0 {
            no_cells_exp_gene.iter().sum::<usize>() / shape.1
        } else {
            0
        };

        println!(
            "  Gene expression stats: min = {} | max = {} | avg = {} cells per gene",
            min_expr.separate_with_underscores(),
            max_expr.separate_with_underscores(),
            avg_expr.separate_with_underscores()
        );
    }

    // Filter genes first
    let genes_to_keep: Vec<usize> = (0..shape.1)
        .filter(|&i| no_cells_exp_gene[i] >= cell_quality.min_cells)
        .collect();

    if verbose {
        println!(
            "  Genes passing filter: {} / {}",
            genes_to_keep.len().separate_with_underscores(),
            shape.1.separate_with_underscores()
        );
    }

    // Create a boolean lookup vector for faster gene filtering
    let mut genes_to_keep_lookup = vec![false; shape.1];
    for &gene_idx in &genes_to_keep {
        genes_to_keep_lookup[gene_idx] = true;
    }

    // second pass - calculate cell metrics using only kept genes
    let mut cell_unique_genes = vec![0usize; shape.0];
    let mut cell_lib_size = vec![0.0f32; shape.0];

    if verbose {
        println!("  Pass 2: Calculating cell metrics from kept genes in chunks...");
    }

    for chunk_start_cell in (0..shape.0).step_by(CELL_CHUNK_SIZE) {
        let chunk_end_cell = (chunk_start_cell + CELL_CHUNK_SIZE).min(shape.0) - 1;

        if verbose && (chunk_start_cell / CELL_CHUNK_SIZE) % 10 == 0 {
            println!(
                "   Processing cells {} / {}",
                chunk_start_cell.separate_with_underscores(),
                shape.0.separate_with_underscores()
            );
        }

        let data_start = indptr[chunk_start_cell] as usize;
        let data_end = indptr[chunk_end_cell + 1] as usize;

        if data_start >= data_end {
            continue;
        }

        let chunk_data: Vec<f32> = data_ds.read_slice_1d(data_start..data_end)?.to_vec();
        let chunk_indices: Vec<u32> = indices_ds.read_slice_1d(data_start..data_end)?.to_vec();

        for cell_idx in chunk_start_cell..=chunk_end_cell {
            let cell_data_start = indptr[cell_idx] as usize;
            let cell_data_end = indptr[cell_idx + 1] as usize;

            let local_start = cell_data_start - data_start;
            let local_end = cell_data_end - data_start;

            for local_idx in local_start..local_end {
                let gene_idx = chunk_indices[local_idx] as usize;

                if genes_to_keep_lookup[gene_idx] {
                    cell_unique_genes[cell_idx] += 1;
                    cell_lib_size[cell_idx] += chunk_data[local_idx];
                }
            }
        }
    }

    if verbose {
        let max_genes = cell_unique_genes.iter().max().unwrap_or(&0);
        let min_genes = cell_unique_genes.iter().min().unwrap_or(&0);
        let max_lib = cell_lib_size.iter().fold(0.0f32, |a, &b| a.max(b));
        let min_lib = cell_lib_size.iter().fold(f32::INFINITY, |a, &b| a.min(b));

        println!(
            "  Cell stats: genes per cell: min = {} | max={}",
            min_genes.separate_with_underscores(),
            max_genes.separate_with_underscores()
        );
        println!(
            "  Cell stats: library size: min = {:.1} | max = {:.1}",
            min_lib.separate_with_underscores(),
            max_lib.separate_with_underscores()
        );
    }

    // Filter cells based on kept genes
    let cells_to_keep: Vec<usize> = (0..shape.0)
        .filter(|&i| {
            cell_unique_genes[i] >= cell_quality.min_unique_genes
                && cell_lib_size[i] >= cell_quality.min_lib_size as f32
        })
        .collect();

    if verbose {
        println!(
            "  Cells passing filter: {} / {}",
            cells_to_keep.len().separate_with_underscores(),
            shape.0.separate_with_underscores()
        );
    }

    let mut file_quality_data = CellOnFileQuality::new(cells_to_keep, genes_to_keep);
    file_quality_data.generate_maps_sets();

    Ok(file_quality_data)
}

/// Helper function that reads in full CSR data from an h5 file
///
/// The function assumes that the data is stored as cells x genes.
///
/// ### Params
///
/// * `file_path` - Path to the h5ad file.
/// * `quality` - Information on the which cells and genes to include after a
///   first pass of the file.
/// * `shape` - The final dimension of the matrix.
///
/// ### Returns
///
/// The `CompressedSparseData` in CSR format with the counts stored as u16.
pub fn read_h5ad_x_data_csr<P: AsRef<Path>>(
    file_path: P,
    quality: &CellOnFileQuality,
    verbose: bool,
) -> Result<CompressedSparseData<u16>> {
    let file = File::open(file_path)?;

    let data_ds = file.dataset("X/data")?;
    let indices_ds = file.dataset("X/indices")?;
    let indptr_ds = file.dataset("X/indptr")?;

    let indptr_raw: Vec<u32> = indptr_ds.read_1d()?.to_vec();

    // Build CSC format directly: indptr = cells, indices = genes
    let mut new_data: Vec<u16> = Vec::new();
    let mut new_indices: Vec<usize> = Vec::new();
    let mut new_indptr: Vec<usize> = Vec::with_capacity(quality.cells_to_keep.len() + 1);
    new_indptr.push(0);

    let total_cells = quality.cells_to_keep.len();

    if verbose {
        println!(
            "  Processing {} cells in chunks (CSC format)...",
            total_cells.separate_with_underscores()
        );
    }

    // should think about moving this into a function parameter ... ?
    const CHUNK_SIZE: usize = 1000;

    for (chunk_idx, cell_chunk) in quality.cells_to_keep.chunks(CHUNK_SIZE).enumerate() {
        if verbose && chunk_idx % 100 == 0 {
            let processed = chunk_idx * CHUNK_SIZE;
            println!(
                "   Processed {} / {} cells",
                processed.min(total_cells).separate_with_underscores(),
                total_cells.separate_with_underscores()
            );
        }

        let start_pos = cell_chunk
            .iter()
            .map(|&c| indptr_raw[c] as usize)
            .min()
            .unwrap_or(0);
        let end_pos = cell_chunk
            .iter()
            .map(|&c| indptr_raw[c + 1] as usize)
            .max()
            .unwrap_or(0);

        if start_pos >= end_pos {
            // Add empty cells to maintain indptr structure
            for _ in cell_chunk {
                new_indptr.push(new_data.len());
            }
            continue;
        }

        let chunk_data: Vec<f32> = data_ds.read_slice_1d(start_pos..end_pos)?.to_vec();
        let chunk_indices: Vec<u32> = indices_ds.read_slice_1d(start_pos..end_pos)?.to_vec();

        for &old_cell_idx in cell_chunk {
            let cell_start = indptr_raw[old_cell_idx] as usize;
            let cell_end = indptr_raw[old_cell_idx + 1] as usize;

            // Collect this cell's data with gene filtering
            let mut cell_data: Vec<(usize, u16)> = Vec::new();

            for idx in cell_start..cell_end {
                let local_idx = idx - start_pos;
                let old_gene_idx = chunk_indices[local_idx] as usize;

                if let Some(&new_gene_idx) = quality.gene_old_to_new.get(&old_gene_idx) {
                    cell_data.push((new_gene_idx, chunk_data[local_idx] as u16));
                }
            }

            // Sort by gene index to maintain consistent ordering
            // IMPORTANT!
            cell_data.sort_by_key(|&(gene_idx, _)| gene_idx);

            // Add to CSC arrays
            for (gene_idx, value) in cell_data {
                new_data.push(value);
                new_indices.push(gene_idx);
            }

            new_indptr.push(new_data.len());
        }
    }

    if verbose {
        println!(
            "   Processed {} / {} cells (complete)",
            total_cells.separate_with_underscores(),
            total_cells.separate_with_underscores()
        );
    }

    let shape = (quality.genes_to_keep.len(), quality.cells_to_keep.len());

    Ok(CompressedSparseData {
        data: new_data,
        indices: new_indices,
        indptr: new_indptr,
        cs_type: CompressedSparseFormat::Csc, // Return CSC format
        data_2: None::<Vec<u16>>,
        shape,
    })
}

/// Stream h5 CSR data directly to disk with batched reading
///
/// This is memory-efficient because we process cells in batches and write
/// immediately. CSR format is ideal since we can calculate library size and
/// normalisation per-cell.
///
/// ### Params
///
/// * `file_path` - Path to the h5ad file
/// * `bin_path` - Path to the to-be-written binary file for the cell data
/// * `quality` - Information on the which cells and genes to include after a
///   first pass of the file.
/// * `cell_qc` - Structure containing the information on which minimum criteria
///   cells and genes need to pass.
/// * `verbose` - Controls verbosity of the function
///
/// ### Returns
///
/// After writing, information on the CellQuality with cell and gene indices,
/// NNZ and lib size per gene.
pub fn write_h5_csr_streaming<P: AsRef<Path>>(
    file_path: P,
    bin_path: P,
    quality: &CellOnFileQuality,
    cell_qc: MinCellQuality,
    verbose: bool,
) -> IoResult<CellQuality> {
    let file = File::open(file_path)?;
    let data_ds = file.dataset("X/data")?;
    let indices_ds = file.dataset("X/indices")?;
    let indptr_ds = file.dataset("X/indptr")?;
    let indptr_raw: Vec<u32> = indptr_ds.read_1d()?.to_vec();

    let mut writer = CellGeneSparseWriter::new(
        bin_path,
        true,
        quality.cells_to_keep.len(),
        quality.genes_to_keep.len(),
    )?;

    let mut lib_size = Vec::with_capacity(quality.cells_to_keep.len());
    let mut nnz = Vec::with_capacity(quality.cells_to_keep.len());

    const CELL_BATCH_SIZE: usize = 1000;
    let total_cells = quality.cells_to_keep.len();

    if verbose {
        println!(
            "  Processing {} cells in batches of {}...",
            total_cells.separate_with_underscores(),
            CELL_BATCH_SIZE.separate_with_underscores()
        );
    }

    for (batch_idx, cell_batch) in quality.cells_to_keep.chunks(CELL_BATCH_SIZE).enumerate() {
        if verbose && batch_idx % 100 == 0 {
            let processed = batch_idx * CELL_BATCH_SIZE;
            println!(
                "   Processed {} / {} cells",
                processed.min(total_cells).separate_with_underscores(),
                total_cells.separate_with_underscores()
            );
        }

        // Find data range for this batch
        let start_pos = cell_batch
            .iter()
            .map(|&c| indptr_raw[c] as usize)
            .min()
            .unwrap_or(0);
        let end_pos = cell_batch
            .iter()
            .map(|&c| indptr_raw[c + 1] as usize)
            .max()
            .unwrap_or(0);

        if start_pos >= end_pos {
            // Empty cells in this batch
            for _ in cell_batch {
                lib_size.push(0);
                nnz.push(0);
                let empty_chunk = CsrCellChunk::from_data(
                    &[] as &[u16],
                    &[] as &[u16],
                    0,
                    cell_qc.target_size,
                    true,
                );
                writer.write_cell_chunk(empty_chunk)?;
            }
            continue;
        }

        // load batch data
        let chunk_data: Vec<f32> = data_ds.read_slice_1d(start_pos..end_pos)?.to_vec();
        let chunk_indices: Vec<u32> = indices_ds.read_slice_1d(start_pos..end_pos)?.to_vec();

        // process each cell in batch
        for &old_cell_idx in cell_batch {
            let cell_start = indptr_raw[old_cell_idx] as usize;
            let cell_end = indptr_raw[old_cell_idx + 1] as usize;

            let mut cell_data: Vec<(usize, u16)> = Vec::new();

            for idx in cell_start..cell_end {
                let local_idx = idx - start_pos;
                let old_gene_idx = chunk_indices[local_idx] as usize;

                if let Some(&new_gene_idx) = quality.gene_old_to_new.get(&old_gene_idx) {
                    let raw_val = chunk_data[local_idx] as u16;
                    cell_data.push((new_gene_idx, raw_val));
                }
            }

            // Sort by gene index for consistent ordering
            cell_data.sort_by_key(|&(gene_idx, _)| gene_idx);

            let gene_indices: Vec<u16> = cell_data.iter().map(|(g, _)| *g as u16).collect();
            let gene_counts: Vec<u16> = cell_data.iter().map(|(_, c)| *c).collect();

            let new_cell_idx = quality.cell_old_to_new[&old_cell_idx];
            let cell_chunk = CsrCellChunk::from_data(
                &gene_counts,
                &gene_indices,
                new_cell_idx,
                cell_qc.target_size,
                true,
            );

            let (nnz_i, lib_size_i) = cell_chunk.get_qc_info();
            nnz.push(nnz_i);
            lib_size.push(lib_size_i);

            // Write immediately - memory is freed after this
            writer.write_cell_chunk(cell_chunk)?;
        }
        // Batch data dropped here - memory freed
    }

    writer.finalise()?;

    if verbose {
        println!(
            "   Processed {} / {} cells (complete)",
            total_cells.separate_with_underscores(),
            total_cells.separate_with_underscores()
        );
    }

    Ok(CellQuality {
        cell_indices: Vec::new(),
        gene_indices: Vec::new(),
        lib_size,
        no_genes: nnz,
    })
}
