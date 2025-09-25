use rustc_hash::FxHashSet;
use std::fs::File;
use std::io::{BufRead, BufReader, Result as IoResult, Seek};
use std::path::Path;
use std::time::Instant;

use crate::core::data::sparse_io::*;
use crate::single_cell::processing::*;

/////////
// MTX //
/////////

/// MTX file metadata
///
/// ### Fields
///
/// * `total_cells` - Number of cells identified in the .mtx header.
/// * `total_genes` - Number of genes identified in the .mtx header.
/// * `total_entries` - Number of entries identified in the .mtx header.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct MtxHeader {
    pub total_cells: usize,
    pub total_genes: usize,
    pub total_entries: usize,
}

/// MTX final data
///
/// Structure to store final results after reading in the .mtx file
///
/// ### Fields
///
/// * `to_keep` - Vector of cells that passed the threshold as booleans.
/// * `no_genes` - Number of genes that were read in.
/// * `no_cells` - Number of cells that were read in.
#[derive(Debug, Clone)]
pub struct MtxFinalData {
    pub cell_qc: CellQuality,
    pub no_genes: usize,
    pub no_cells: usize,
}

/// MTX Reader for bixverse
///
/// ### Fields
///
/// * `reader` - Buffered reader of the mtx file
/// * `header` - The header of the mtx file
/// * `min_genes` - Minimum number of unique genes identified in the cell to
///   be included in the binarised file.
/// * `min_lib_size` - Minimum library size in the cell to be included in the
///   binarised file.
/// * `target_size` - Target size for the normalised counts.
pub struct MtxReader {
    reader: BufReader<File>,
    header: MtxHeader,
    qc_params: MinCellQuality,
}

impl MtxReader {
    /// Generate a new instance of the reader
    ///
    /// ### Params
    ///
    /// * `path` - Path to the mtx file.
    /// * `min_genes` - Minimum number of genes to that have to be found in the
    ///   cell to be included.
    /// * `min_lib_size` - Minimum library size in the cell to be included.
    /// * `target_size` - Target size after normalisation.
    pub fn new<P: AsRef<Path>>(path: P, qc_params: MinCellQuality) -> IoResult<Self> {
        let file = File::open(path)?;
        let mut reader = BufReader::new(file);

        let header = Self::parse_header(&mut reader)?;

        Ok(Self {
            reader,
            header,
            qc_params,
        })
    }

    /// Parse the header of the mtx file
    ///
    /// ### Returns
    ///
    /// The `MtxHeader`
    fn parse_header(reader: &mut BufReader<File>) -> IoResult<MtxHeader> {
        let mut line = String::new();

        loop {
            line.clear();
            reader.read_line(&mut line)?;
            if !line.starts_with('%') {
                break;
            }
        }

        let parts: Vec<&str> = line.split_whitespace().collect();

        if parts.len() != 3 {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "Invalid MTX header format",
            ));
        }

        let total_cells = parts[0].parse().map_err(|_| {
            std::io::Error::new(std::io::ErrorKind::InvalidData, "Invalid cell count")
        })?;

        let total_genes = parts[1].parse().map_err(|_| {
            std::io::Error::new(std::io::ErrorKind::InvalidData, "Invalid gene count")
        })?;

        let total_entries = parts[2].parse().map_err(|_| {
            std::io::Error::new(std::io::ErrorKind::InvalidData, "Invalid entry count")
        })?;

        Ok(MtxHeader {
            total_cells,
            total_genes,
            total_entries,
        })
    }

    /// Helper to parse the file to understand which cells to keep
    ///
    /// ### Returns
    ///
    /// The CellOnFileQuality file containing the indices and mappings
    /// for the cells/genes to keep.
    pub fn parse_mtx_quality(&mut self, verbose: bool) -> IoResult<CellOnFileQuality> {
        let mut gene_cell_count = vec![0u32; self.header.total_genes];
        let mut line_buffer = Vec::with_capacity(128);

        self.reader.rewind()?;
        Self::skip_header(&mut self.reader)?;

        if verbose {
            println!("First file pass - getting gene statistics:");
        }

        let first_scan_time = Instant::now();

        // Pass 1: Gene statistics only
        while {
            line_buffer.clear();
            self.reader.read_until(b'\n', &mut line_buffer)? > 0
        } {
            if line_buffer.len() < 3 {
                continue;
            }
            if line_buffer[line_buffer.len() - 1] == b'\n' {
                line_buffer.pop();
            }
            if !line_buffer.is_empty() && line_buffer[line_buffer.len() - 1] == b'\r' {
                line_buffer.pop();
            }

            if let Some((_, col, _)) = parse_mtx_line(&line_buffer) {
                let gene_idx = (col - 1) as usize;
                if gene_idx < self.header.total_genes {
                    gene_cell_count[gene_idx] += 1;
                }
            }
        }

        // Filter genes
        let genes_to_keep_set: FxHashSet<usize> = (0..self.header.total_genes)
            .filter(|&i| gene_cell_count[i] as usize >= self.qc_params.min_cells)
            .collect();

        let first_scan_end = first_scan_time.elapsed();

        if verbose {
            println!("First file scan done: {:.2?}", first_scan_end);
            println!("Second file pass - getting cell statistics");
        }

        let second_scan_time = Instant::now();

        // Pass 2: Cell statistics with filtered genes
        let mut cell_gene_count = vec![0u32; self.header.total_cells];
        let mut cell_lib_size = vec![0u32; self.header.total_cells];

        self.reader.rewind()?;
        Self::skip_header(&mut self.reader)?;

        while {
            line_buffer.clear();
            self.reader.read_until(b'\n', &mut line_buffer)? > 0
        } {
            if line_buffer.len() < 3 {
                continue;
            }
            if line_buffer[line_buffer.len() - 1] == b'\n' {
                line_buffer.pop();
            }
            if !line_buffer.is_empty() && line_buffer[line_buffer.len() - 1] == b'\r' {
                line_buffer.pop();
            }

            if let Some((row, col, value)) = parse_mtx_line(&line_buffer) {
                let gene_idx = (col - 1) as usize;
                let cell_idx = (row - 1) as usize;

                if genes_to_keep_set.contains(&gene_idx) && cell_idx < self.header.total_cells {
                    cell_gene_count[cell_idx] += 1;
                    cell_lib_size[cell_idx] += value as u32;
                }
            }
        }

        // Filter cells
        let cells_to_keep: Vec<usize> = (0..self.header.total_cells)
            .filter(|&i| {
                cell_gene_count[i] as usize >= self.qc_params.min_unique_genes
                    && cell_lib_size[i] as f32 >= self.qc_params.min_lib_size as f32
            })
            .collect();

        let genes_to_keep: Vec<usize> = genes_to_keep_set.into_iter().collect();
        let mut quality = CellOnFileQuality::new(cells_to_keep, genes_to_keep);
        quality.generate_maps_sets();

        let second_scan_end = second_scan_time.elapsed();

        if verbose {
            println!("Second file scan done: {:.2?}", second_scan_end);
            println!(
                "Genes passing QC: {}/{}",
                quality.genes_to_keep.len(),
                self.header.total_genes
            );
            println!(
                "Cells passing QC: {}/{}",
                quality.cells_to_keep.len(),
                self.header.total_cells
            );
        }

        Ok(quality)
    }

    /// Process the mtx file and write to binarised Rust file
    ///
    /// ### Params
    ///
    /// * `bin_path` - Where to save the binarised file
    ///
    /// ### Returns
    ///
    /// The `MtxFinalData` with information how many cells were written to
    /// file, how many genes were included in which cells did not parse
    /// the thresholds.
    pub fn process_mtx_and_write_bin(
        mut self,
        bin_path: &str,
        quality: &CellOnFileQuality,
        verbose: bool,
    ) -> IoResult<MtxFinalData> {
        let mut writer = CellGeneSparseWriter::new(
            bin_path,
            true,
            quality.cells_to_keep.len(),
            quality.genes_to_keep.len(),
        )?;

        // Collect all data first
        let mut cell_data: Vec<Vec<(u16, u16)>> = vec![Vec::new(); quality.cells_to_keep.len()];
        let mut line_buffer = Vec::with_capacity(64);

        self.reader.rewind()?;
        Self::skip_header(&mut self.reader)?;

        if verbose {
            println!(
                "Starting to write high quality cells, genes in a cell-friendly format to disk."
            )
        }

        while {
            line_buffer.clear();
            self.reader.read_until(b'\n', &mut line_buffer)? > 0
        } {
            if line_buffer.last() == Some(&b'\n') {
                line_buffer.pop();
            }
            if line_buffer.last() == Some(&b'\r') {
                line_buffer.pop();
            }
            if line_buffer.is_empty() {
                continue;
            }

            let (row, col, value) = match parse_mtx_line(&line_buffer) {
                Some(parsed) => parsed,
                None => continue,
            };

            let old_gene_idx = (col - 1) as usize;
            let old_cell_idx = (row - 1) as usize;

            if !quality.genes_to_keep_set.contains(&old_gene_idx)
                || !quality.cells_to_keep_set.contains(&old_cell_idx)
            {
                continue;
            }

            let new_cell_idx = quality.cell_old_to_new[&old_cell_idx];
            let new_gene_idx = quality.gene_old_to_new[&old_gene_idx] as u16;

            cell_data[new_cell_idx].push((new_gene_idx, value));
        }

        // Write cells in order
        let mut lib_size = Vec::with_capacity(quality.cells_to_keep.len());
        let mut nnz = Vec::with_capacity(quality.cells_to_keep.len());

        for (cell_idx, data) in cell_data.iter_mut().enumerate() {
            if data.is_empty() {
                continue;
            }

            data.sort_by_key(|(gene_idx, _)| *gene_idx);

            let gene_indices: Vec<u16> = data.iter().map(|(g, _)| *g).collect();
            let gene_counts: Vec<u16> = data.iter().map(|(_, c)| *c).collect();

            let total_umi: u32 = gene_counts.iter().map(|&x| x as u32).sum();
            let n_genes = gene_counts.len();

            lib_size.push(total_umi as usize);
            nnz.push(n_genes);

            let cell_chunk = CsrCellChunk::from_data(
                &gene_counts,
                &gene_indices,
                cell_idx,
                self.qc_params.target_size,
                true,
            );
            writer.write_cell_chunk(cell_chunk)?;
        }

        let mut to_include = vec![false; self.header.total_cells];
        for &cell_idx in &quality.cells_to_keep {
            to_include[cell_idx] = true;
        }

        writer.finalise()?;

        let cell_quality = CellQuality {
            cell_indices: quality.cells_to_keep.to_vec(),
            gene_indices: quality.genes_to_keep.to_vec(),
            lib_size,
            no_genes: nnz,
        };

        Ok(MtxFinalData {
            cell_qc: cell_quality,
            no_genes: quality.genes_to_keep.len(),
            no_cells: quality.cells_to_keep.len(),
        })
    }

    fn skip_header(reader: &mut BufReader<File>) -> IoResult<()> {
        let mut line = String::new();
        loop {
            line.clear();
            reader.read_line(&mut line)?;
            if !line.starts_with('%') {
                break;
            }
        }
        Ok(())
    }
}

/////////////
// Helpers //
/////////////

/// Parse mtx line from bytes
///
/// ### Params
///
/// * `line` - The file line as bytes
///
/// ### Return
///
/// Returns an Option of a tuple representing
/// `<cell_index, gene_index, raw_count>`
#[inline]
fn parse_mtx_line(line: &[u8]) -> Option<(u32, u32, u16)> {
    let mut i = 0;
    let len = line.len();

    // Parse first number (row)
    let mut row = 0u32;
    while i < len && line[i].is_ascii_digit() {
        row = row * 10 + (line[i] - b'0') as u32;
        i += 1;
    }
    if i == 0 {
        return None;
    }

    // Skip whitespace
    while i < len && (line[i] == b' ' || line[i] == b'\t') {
        i += 1;
    }
    if i >= len {
        return None;
    }

    // Parse second number (col)
    let mut col = 0u32;
    while i < len && line[i].is_ascii_digit() {
        col = col * 10 + (line[i] - b'0') as u32;
        i += 1;
    }

    // Skip whitespace
    while i < len && (line[i] == b' ' || line[i] == b'\t') {
        i += 1;
    }
    if i >= len {
        return None;
    }

    // Parse third number (value)
    let mut val = 0u16;
    while i < len && line[i].is_ascii_digit() {
        val = val * 10 + (line[i] - b'0') as u16;
        i += 1;
    }

    Some((row, col, val))
}
