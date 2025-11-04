use hdf5::{File, Result};
use rayon::prelude::*;
use std::io::Result as IoResult;
use std::path::Path;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;
use std::time::Instant;
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

    let n_cells = file_data.indptr.len() - 1;
    let mut writer = CellGeneSparseWriter::new(bin_path, true, no_cells, no_genes).unwrap();

    let mut lib_size = Vec::with_capacity(n_cells);
    let mut nnz = Vec::with_capacity(n_cells);

    for i in 0..n_cells {
        let start_i = file_data.indptr[i];
        let end_i = file_data.indptr[i + 1];

        let cell_data = &file_data.data[start_i..end_i];
        let cell_indices = &file_data.indices[start_i..end_i];

        let cell_chunk =
            CsrCellChunk::from_data(cell_data, cell_indices, i, cell_quality.target_size, true);

        let (nnz_i, lib_size_i) = cell_chunk.get_qc_info();
        nnz.push(nnz_i);
        lib_size.push(lib_size_i);

        writer.write_cell_chunk(cell_chunk).unwrap();

        if verbose && (i + 1) % 100000 == 0 {
            println!(
                "  Written {} / {} cells to disk.",
                (i + 1).separate_with_underscores(),
                n_cells.separate_with_underscores()
            );
        }
    }

    if verbose {
        println!(
            "  Written {} / {} cells (complete).",
            n_cells.separate_with_underscores(),
            n_cells.separate_with_underscores()
        );
        println!("Finalising file...");
    }

    writer.update_header_no_cells(file_quality.cells_to_keep.len());
    writer.update_header_no_genes(file_quality.genes_to_keep.len());
    writer.finalise().unwrap();

    let cell_qc = CellQuality {
        cell_indices: file_quality.cells_to_keep.clone(),
        gene_indices: file_quality.genes_to_keep.clone(),
        lib_size,
        no_genes: nnz,
    };

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
                println!("  Pass 2/2: Writing cells with normalisation...");
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
    let file_path = file_path.as_ref();

    if verbose {
        println!(
            "  Reading CSC matrix structure (shape: {} x {})...",
            shape.0.separate_with_underscores(),
            shape.1.separate_with_underscores()
        );
    }

    let file = File::open(file_path)?;
    let indptr: Vec<u32> = file.dataset("X/indptr")?.read_1d()?.to_vec();

    let mut no_cells_exp_gene: Vec<usize> = Vec::with_capacity(shape.1);
    for i in 0..shape.1 {
        no_cells_exp_gene.push((indptr[i + 1] - indptr[i]) as usize);
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

    let genes_to_keep: Vec<usize> = (0..shape.1)
        .filter(|&i| no_cells_exp_gene[i] >= cell_quality.min_cells)
        .collect();

    if verbose {
        println!(
            "  Genes passing filter: {} / {}",
            genes_to_keep.len().separate_with_underscores(),
            shape.1.separate_with_underscores()
        );
        println!("Calculating cell metrics in parallel...");
    }

    const GENE_CHUNK_SIZE: usize = 10000;
    let gene_chunks: Vec<&[usize]> = genes_to_keep.chunks(GENE_CHUNK_SIZE).collect();
    let num_chunks = gene_chunks.len();

    let calc_time = Instant::now();
    let completed_chunks = Arc::new(AtomicUsize::new(0));
    let report_interval = (num_chunks / 10).max(1);

    let cell_stats: Vec<(Vec<usize>, Vec<f32>)> = gene_chunks
        .par_iter()
        .map(|gene_chunk| {
            let mut local_unique = vec![0usize; shape.0];
            let mut local_lib_size = vec![0.0f32; shape.0];

            if let Ok(file) = File::open(file_path) {
                if let (Ok(data_ds), Ok(indices_ds)) =
                    (file.dataset("X/data"), file.dataset("X/indices"))
                {
                    let chunk_start_gene = gene_chunk[0];
                    let chunk_end_gene = gene_chunk[gene_chunk.len() - 1];
                    let data_start = indptr[chunk_start_gene] as usize;
                    let data_end = indptr[chunk_end_gene + 1] as usize;

                    if data_start < data_end {
                        if let (Ok(chunk_data), Ok(chunk_indices)) = (
                            data_ds.read_slice_1d(data_start..data_end),
                            indices_ds.read_slice_1d(data_start..data_end),
                        ) {
                            let chunk_data: Vec<f32> = chunk_data.to_vec();
                            let chunk_indices: Vec<u32> = chunk_indices.to_vec();

                            for &gene_idx in gene_chunk.iter() {
                                let gene_data_start = indptr[gene_idx] as usize - data_start;
                                let gene_data_end = indptr[gene_idx + 1] as usize - data_start;

                                for idx in gene_data_start..gene_data_end {
                                    let cell_idx = chunk_indices[idx] as usize;
                                    if cell_idx < shape.0 {
                                        local_unique[cell_idx] += 1;
                                        local_lib_size[cell_idx] += chunk_data[idx];
                                    }
                                }
                            }
                        }
                    }
                }
            }

            if verbose {
                let completed = completed_chunks.fetch_add(1, Ordering::Relaxed) + 1;
                if completed.is_multiple_of(report_interval) || completed == num_chunks {
                    let progress =
                        ((completed as f64 / num_chunks as f64 * 10.0).round() as usize) * 10;
                    println!(
                        "  Processed {}% of chunks ({}/{})",
                        progress, completed, num_chunks
                    );
                }
            }

            (local_unique, local_lib_size)
        })
        .collect();

    let mut cell_unique_genes = vec![0usize; shape.0];
    let mut cell_lib_size = vec![0.0f32; shape.0];

    for (local_unique, local_lib) in cell_stats {
        for i in 0..shape.0 {
            cell_unique_genes[i] += local_unique[i];
            cell_lib_size[i] += local_lib[i];
        }
    }

    let calc_elapsed = calc_time.elapsed();

    if verbose {
        println!("Cell metrics calculation done: {:.2?}", calc_elapsed);

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
    let file_path = file_path.as_ref();

    if verbose {
        println!(
            "  Reading CSR matrix structure (shape: {} x {} )...",
            shape.0.separate_with_underscores(),
            shape.1.separate_with_underscores()
        );
    }

    let file = File::open(file_path)?;
    let indptr: Vec<u32> = file.dataset("X/indptr")?.read_1d()?.to_vec();

    const CELL_CHUNK_SIZE: usize = 10000;
    let chunks: Vec<usize> = (0..shape.0).step_by(CELL_CHUNK_SIZE).collect();
    let num_chunks = chunks.len();

    if verbose {
        println!("First pass - gene expression statistics:");
    }

    let first_pass_time = Instant::now();
    let completed_chunks = Arc::new(AtomicUsize::new(0));
    let report_interval = (num_chunks / 10).max(1);

    let gene_counts: Vec<Vec<usize>> = chunks
        .par_iter()
        .map(|&chunk_start_cell| {
            let mut local_counts = vec![0usize; shape.1];

            if let Ok(file) = File::open(file_path) {
                if let Ok(indices_ds) = file.dataset("X/indices") {
                    let chunk_end_cell = (chunk_start_cell + CELL_CHUNK_SIZE).min(shape.0) - 1;
                    let data_start = indptr[chunk_start_cell] as usize;
                    let data_end = indptr[chunk_end_cell + 1] as usize;

                    if data_start < data_end {
                        if let Ok(chunk_indices) = indices_ds.read_slice_1d(data_start..data_end) {
                            let chunk_indices: Vec<u32> = chunk_indices.to_vec();

                            for cell_idx in chunk_start_cell..=chunk_end_cell {
                                let cell_data_start = indptr[cell_idx] as usize - data_start;
                                let cell_data_end = indptr[cell_idx + 1] as usize - data_start;

                                for local_idx in cell_data_start..cell_data_end {
                                    let gene_idx = chunk_indices[local_idx] as usize;
                                    if gene_idx < shape.1 {
                                        local_counts[gene_idx] += 1;
                                    }
                                }
                            }
                        }
                    }
                }
            }

            if verbose {
                let completed = completed_chunks.fetch_add(1, Ordering::Relaxed) + 1;
                if completed.is_multiple_of(report_interval) || completed == num_chunks {
                    let progress =
                        ((completed as f64 / num_chunks as f64 * 10.0).round() as usize) * 10;
                    println!(
                        "  Processed {}% of chunks ({}/{})",
                        progress, completed, num_chunks
                    );
                }
            }

            local_counts
        })
        .collect();

    let mut no_cells_exp_gene = vec![0usize; shape.1];
    for local_counts in gene_counts {
        for (i, count) in local_counts.into_iter().enumerate() {
            no_cells_exp_gene[i] += count;
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

    let genes_to_keep: Vec<usize> = (0..shape.1)
        .filter(|&i| no_cells_exp_gene[i] >= cell_quality.min_cells)
        .collect();

    let first_pass_elapsed = first_pass_time.elapsed();

    if verbose {
        println!("First pass done: {:.2?}", first_pass_elapsed);
        println!(
            "  Genes passing filter: {} / {}",
            genes_to_keep.len().separate_with_underscores(),
            shape.1.separate_with_underscores()
        );
        println!("Second pass - cell statistics:");
    }

    let mut genes_to_keep_lookup = vec![false; shape.1];
    for &gene_idx in &genes_to_keep {
        genes_to_keep_lookup[gene_idx] = true;
    }

    let second_pass_time = Instant::now();
    let completed_chunks = Arc::new(AtomicUsize::new(0));

    let cell_stats: Vec<(Vec<usize>, Vec<f32>)> = chunks
        .par_iter()
        .map(|&chunk_start_cell| {
            let chunk_end_cell = (chunk_start_cell + CELL_CHUNK_SIZE).min(shape.0) - 1;
            let mut local_unique = vec![0usize; chunk_end_cell - chunk_start_cell + 1];
            let mut local_lib_size = vec![0.0f32; chunk_end_cell - chunk_start_cell + 1];

            if let Ok(file) = File::open(file_path) {
                if let (Ok(data_ds), Ok(indices_ds)) =
                    (file.dataset("X/data"), file.dataset("X/indices"))
                {
                    let data_start = indptr[chunk_start_cell] as usize;
                    let data_end = indptr[chunk_end_cell + 1] as usize;

                    if data_start < data_end {
                        if let (Ok(chunk_data), Ok(chunk_indices)) = (
                            data_ds.read_slice_1d(data_start..data_end),
                            indices_ds.read_slice_1d(data_start..data_end),
                        ) {
                            let chunk_data: Vec<f32> = chunk_data.to_vec();
                            let chunk_indices: Vec<u32> = chunk_indices.to_vec();

                            for cell_idx in chunk_start_cell..=chunk_end_cell {
                                let cell_data_start = indptr[cell_idx] as usize - data_start;
                                let cell_data_end = indptr[cell_idx + 1] as usize - data_start;
                                let local_cell_idx = cell_idx - chunk_start_cell;

                                for local_idx in cell_data_start..cell_data_end {
                                    let gene_idx = chunk_indices[local_idx] as usize;
                                    if genes_to_keep_lookup[gene_idx] {
                                        local_unique[local_cell_idx] += 1;
                                        local_lib_size[local_cell_idx] += chunk_data[local_idx];
                                    }
                                }
                            }
                        }
                    }
                }
            }

            if verbose {
                let completed = completed_chunks.fetch_add(1, Ordering::Relaxed) + 1;
                if completed.is_multiple_of(report_interval) || completed == num_chunks {
                    let progress =
                        ((completed as f64 / num_chunks as f64 * 10.0).round() as usize) * 10;
                    println!(
                        "  Processed {}% of chunks ({}/{})",
                        progress, completed, num_chunks
                    );
                }
            }

            (local_unique, local_lib_size)
        })
        .collect();

    let mut cell_unique_genes = vec![0usize; shape.0];
    let mut cell_lib_size = vec![0.0f32; shape.0];

    for (chunk_idx, (local_unique, local_lib)) in cell_stats.into_iter().enumerate() {
        let chunk_start = chunks[chunk_idx];
        for (i, (unique, lib)) in local_unique.into_iter().zip(local_lib).enumerate() {
            cell_unique_genes[chunk_start + i] = unique;
            cell_lib_size[chunk_start + i] = lib;
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
            min_lib, max_lib
        );
    }

    let cells_to_keep: Vec<usize> = (0..shape.0)
        .filter(|&i| {
            cell_unique_genes[i] >= cell_quality.min_unique_genes
                && cell_lib_size[i] >= cell_quality.min_lib_size as f32
        })
        .collect();

    let second_pass_elapsed = second_pass_time.elapsed();

    if verbose {
        println!("Second pass done: {:.2?}", second_pass_elapsed);
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

    let start_write = Instant::now();

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

    let end_write = start_write.elapsed();

    if verbose {
        println!(
            "   Processed {} / {} cells (complete) in {:.2?}.",
            total_cells.separate_with_underscores(),
            total_cells.separate_with_underscores(),
            end_write
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
    let file = File::open(&file_path)?;
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
    let num_batches = total_cells.div_ceil(CELL_BATCH_SIZE);

    if verbose {
        println!(
            "  Processing {} cells in batches of {}...",
            total_cells.separate_with_underscores(),
            CELL_BATCH_SIZE.separate_with_underscores()
        );
    }

    let start_write = Instant::now();

    // Reusable buffers to avoid allocations
    let mut cell_data: Vec<(usize, u16)> = Vec::with_capacity(10000);
    let mut gene_indices: Vec<u16> = Vec::with_capacity(10000);
    let mut gene_counts: Vec<u16> = Vec::with_capacity(10000);

    for (batch_idx, cell_batch) in quality.cells_to_keep.chunks(CELL_BATCH_SIZE).enumerate() {
        if verbose && (batch_idx % ((num_batches / 10).max(1)) == 0 || batch_idx == num_batches - 1)
        {
            let progress = ((batch_idx as f64 / num_batches as f64 * 10.0).round() as usize) * 10;
            let processed = (batch_idx + 1) * CELL_BATCH_SIZE;
            println!(
                "  Processed {}% ({} / {} cells)",
                progress,
                processed.min(total_cells).separate_with_underscores(),
                total_cells.separate_with_underscores()
            );
        }

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
            for &old_cell_idx in cell_batch {
                lib_size.push(0);
                nnz.push(0);
                let new_cell_idx = quality.cell_old_to_new[&old_cell_idx];
                let empty_chunk = CsrCellChunk::from_data(
                    &[] as &[u16],
                    &[] as &[u16],
                    new_cell_idx,
                    cell_qc.target_size,
                    true,
                );
                writer.write_cell_chunk(empty_chunk)?;
            }
            continue;
        }

        let chunk_data: Vec<f32> = data_ds.read_slice_1d(start_pos..end_pos)?.to_vec();
        let chunk_indices: Vec<u32> = indices_ds.read_slice_1d(start_pos..end_pos)?.to_vec();

        for &old_cell_idx in cell_batch {
            let cell_start = indptr_raw[old_cell_idx] as usize;
            let cell_end = indptr_raw[old_cell_idx + 1] as usize;

            cell_data.clear();
            gene_indices.clear();
            gene_counts.clear();

            for idx in cell_start..cell_end {
                let local_idx = idx - start_pos;
                let old_gene_idx = chunk_indices[local_idx] as usize;

                if let Some(&new_gene_idx) = quality.gene_old_to_new.get(&old_gene_idx) {
                    let raw_val = chunk_data[local_idx] as u16;
                    cell_data.push((new_gene_idx, raw_val));
                }
            }

            if !cell_data.is_empty() {
                // Check if already sorted (common in CSR)
                let needs_sort = cell_data.windows(2).any(|w| w[0].0 > w[1].0);
                if needs_sort {
                    cell_data.sort_unstable_by_key(|&(gene_idx, _)| gene_idx);
                }

                gene_indices.extend(cell_data.iter().map(|(g, _)| *g as u16));
                gene_counts.extend(cell_data.iter().map(|(_, c)| *c));
            }

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

            writer.write_cell_chunk(cell_chunk)?;
        }
    }

    writer.finalise()?;

    let end_write = start_write.elapsed();

    if verbose {
        println!("  Writing complete in {:.2?}", end_write);
    }

    Ok(CellQuality {
        cell_indices: Vec::new(),
        gene_indices: Vec::new(),
        lib_size,
        no_genes: nnz,
    })
}
