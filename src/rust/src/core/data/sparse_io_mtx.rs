use rayon::prelude::*;
use rustc_hash::FxHashSet;
use std::fs::File;
use std::io::{BufRead, BufReader, Read, Result as IoResult, Seek};
use std::path::{Path, PathBuf};
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
/// * `cell_qc` - Structure containing the information on which cells/genes to
///   keep and library size and NNZ for cells.
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
/// * `qc_params` - The min quality parameters that genes and cells have to
///   reach and the target library size
/// * `cells_as_rows` - Boolean. Are the cells the rows (= true) or columns in
///   the mtx file.
pub struct MtxReader {
    path: PathBuf,
    reader: BufReader<File>,
    header: MtxHeader,
    qc_params: MinCellQuality,
    cells_as_rows: bool,
}

impl MtxReader {
    /// Generate a new instance of the reader
    ///
    /// ### Params
    ///
    /// * `path` - Path to the mtx file.
    /// * `qc_params` - The min quality parameters that genes and cells have to
    ///   reach and the target library size
    /// * `cells_as_rows` - Boolean. Are the cells the rows (= true) or columns in
    ///   the mtx file.
    pub fn new<P: AsRef<Path>>(
        path: P,
        qc_params: MinCellQuality,
        cells_as_rows: bool,
    ) -> IoResult<Self> {
        let path = path.as_ref().to_path_buf();
        let file = File::open(&path)?;
        let mut reader = BufReader::with_capacity(1024 * 1024, file);

        let header = Self::parse_header(&mut reader, cells_as_rows)?;

        Ok(Self {
            path,
            reader,
            header,
            qc_params,
            cells_as_rows,
        })
    }

    /// Parse the header of the mtx file
    ///
    /// ### Returns
    ///
    /// The `MtxHeader`
    fn parse_header(reader: &mut BufReader<File>, cells_as_rows: bool) -> IoResult<MtxHeader> {
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

        let (total_cells, total_genes) = if cells_as_rows {
            // Header format: cells genes entries
            (
                parts[0].parse().map_err(|_| {
                    std::io::Error::new(std::io::ErrorKind::InvalidData, "Invalid cell count")
                })?,
                parts[1].parse().map_err(|_| {
                    std::io::Error::new(std::io::ErrorKind::InvalidData, "Invalid gene count")
                })?,
            )
        } else {
            // Header format: genes cells entries
            (
                parts[1].parse().map_err(|_| {
                    std::io::Error::new(std::io::ErrorKind::InvalidData, "Invalid cell count")
                })?,
                parts[0].parse().map_err(|_| {
                    std::io::Error::new(std::io::ErrorKind::InvalidData, "Invalid gene count")
                })?,
            )
        };

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
    /// ### Params
    ///
    /// * `verbose` - Controls verbosity of the function.
    ///
    /// ### Returns
    ///
    /// The CellOnFileQuality file containing the indices and mappings
    /// for the cells/genes to keep.
    pub fn parse_mtx_quality(&mut self, verbose: bool) -> IoResult<CellOnFileQuality> {
        const CHUNK_SIZE: u64 = 64 * 1024 * 1024;
        let file_size = self.reader.get_ref().metadata()?.len();
        let num_chunks = ((file_size / CHUNK_SIZE) as usize).max(1);

        if verbose {
            println!("First file pass - getting gene statistics:");
        }

        let first_scan_time = Instant::now();

        let boundaries = self.find_chunk_boundaries(num_chunks)?;

        let results: Vec<_> = boundaries
            .par_iter()
            .map(|&(start, end)| {
                let mut local_gene_cells = vec![FxHashSet::default(); self.header.total_genes];

                if let Ok(file) = File::open(&self.path) {
                    let mut reader = BufReader::with_capacity(256 * 1024, file);
                    if reader.seek(std::io::SeekFrom::Start(start)).is_ok() {
                        let mut line_buffer = Vec::with_capacity(64);
                        let mut bytes_read = 0u64;

                        while bytes_read < (end - start) {
                            line_buffer.clear();
                            if let Ok(n) = reader.read_until(b'\n', &mut line_buffer) {
                                if n == 0 {
                                    break;
                                }
                                bytes_read += n as u64;

                                let len = line_buffer.len();
                                if len < 3 {
                                    continue;
                                }
                                let trim_end = if line_buffer[len - 1] == b'\n' {
                                    if len > 1 && line_buffer[len - 2] == b'\r' {
                                        len - 2
                                    } else {
                                        len - 1
                                    }
                                } else {
                                    len
                                };

                                if let Some((row, col, _)) =
                                    parse_mtx_line(&line_buffer[..trim_end])
                                {
                                    let (cell_idx, gene_idx) = if self.cells_as_rows {
                                        ((row - 1) as usize, (col - 1) as usize)
                                    } else {
                                        ((col - 1) as usize, (row - 1) as usize)
                                    };

                                    if gene_idx < self.header.total_genes {
                                        local_gene_cells[gene_idx].insert(cell_idx);
                                    }
                                }
                            }
                        }
                    }
                }

                local_gene_cells
            })
            .collect();

        let mut gene_cells = vec![FxHashSet::default(); self.header.total_genes];
        for local_genes in results {
            for (i, local_set) in local_genes.into_iter().enumerate() {
                gene_cells[i].extend(local_set);
            }
        }

        // Filter genes
        let genes_to_keep_set: FxHashSet<usize> = (0..self.header.total_genes)
            .filter(|&i| gene_cells[i].len() >= self.qc_params.min_cells)
            .collect();

        let first_scan_end = first_scan_time.elapsed();

        if verbose {
            println!("First pass done: {:.2?}", first_scan_end);
            println!("Second pass - cell statistics:");
        }

        let second_scan_time = Instant::now();

        // Parallel second pass - cell stats with filtered genes
        let results: Vec<_> = boundaries
            .par_iter()
            .map(|&(start, end)| {
                let mut local_cell_stats = vec![(0u32, 0u32); self.header.total_cells];

                if let Ok(file) = File::open(&self.path) {
                    let mut reader = BufReader::with_capacity(256 * 1024, file);
                    if reader.seek(std::io::SeekFrom::Start(start)).is_ok() {
                        let mut line_buffer = Vec::with_capacity(64);
                        let mut bytes_read = 0u64;

                        while bytes_read < (end - start) {
                            line_buffer.clear();
                            if let Ok(n) = reader.read_until(b'\n', &mut line_buffer) {
                                if n == 0 {
                                    break;
                                }
                                bytes_read += n as u64;

                                let len = line_buffer.len();
                                if len < 3 {
                                    continue;
                                }
                                let trim_end = if line_buffer[len - 1] == b'\n' {
                                    if len > 1 && line_buffer[len - 2] == b'\r' {
                                        len - 2
                                    } else {
                                        len - 1
                                    }
                                } else {
                                    len
                                };

                                if let Some((row, col, value)) =
                                    parse_mtx_line(&line_buffer[..trim_end])
                                {
                                    let (cell_idx, gene_idx) = if self.cells_as_rows {
                                        ((row - 1) as usize, (col - 1) as usize)
                                    } else {
                                        ((col - 1) as usize, (row - 1) as usize)
                                    };

                                    if genes_to_keep_set.contains(&gene_idx)
                                        && cell_idx < self.header.total_cells
                                    {
                                        local_cell_stats[cell_idx].0 += 1;
                                        local_cell_stats[cell_idx].1 += value as u32;
                                    }
                                }
                            }
                        }
                    }
                }

                local_cell_stats
            })
            .collect();

        // Merge cell results
        let mut cell_gene_count = vec![0u32; self.header.total_cells];
        let mut cell_lib_size = vec![0u32; self.header.total_cells];

        for local_cells in results {
            for (i, (count, size)) in local_cells.into_iter().enumerate() {
                cell_gene_count[i] += count;
                cell_lib_size[i] += size;
            }
        }

        // Filter cells
        let cells_to_keep: Vec<usize> = (0..self.header.total_cells)
            .filter(|&i| {
                cell_gene_count[i] as usize >= self.qc_params.min_unique_genes
                    && cell_lib_size[i] as f32 >= self.qc_params.min_lib_size as f32
            })
            .collect();

        let mut genes_to_keep: Vec<usize> = genes_to_keep_set.into_iter().collect();
        genes_to_keep.sort_unstable();

        let mut quality = CellOnFileQuality::new(cells_to_keep, genes_to_keep);
        quality.generate_maps_sets();

        let second_scan_end = second_scan_time.elapsed();

        if verbose {
            println!("Second pass done: {:.2?}", second_scan_end);
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
    /// * `bin_path` - Where to save the binarised file.
    /// * `quality` - Structure indicating which cells and genes to keep from
    ///   the mtx file.
    /// * `verbose` - Controls verbosity of the function.
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

        let mut cell_data: Vec<Vec<(u16, u16)>> = vec![Vec::new(); quality.cells_to_keep.len()];
        let mut line_buffer = Vec::with_capacity(64);

        self.reader.rewind()?;
        Self::skip_header(&mut self.reader)?;

        if verbose {
            println!(
                "Starting to write high quality cells, genes in a cell-friendly format to disk."
            )
        }

        let mut lines_read = 0usize;
        let report_interval = (self.header.total_entries / 10).max(1);

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

            let (old_cell_idx, old_gene_idx) = if self.cells_as_rows {
                ((row - 1) as usize, (col - 1) as usize)
            } else {
                ((col - 1) as usize, (row - 1) as usize)
            };

            if !quality.genes_to_keep_set.contains(&old_gene_idx)
                || !quality.cells_to_keep_set.contains(&old_cell_idx)
            {
                continue;
            }

            let new_cell_idx = quality.cell_old_to_new[&old_cell_idx];
            let new_gene_idx = quality.gene_old_to_new[&old_gene_idx] as u16;

            cell_data[new_cell_idx].push((new_gene_idx, value));

            lines_read += 1;
            if verbose && lines_read % report_interval == 0 {
                let progress =
                    (lines_read as f64 / self.header.total_entries as f64 * 100.0) as usize;
                println!("  Processed {}% of entries", progress);
            }
        }

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

    /// Helper function to skip the header
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

    /// Generate chunk boundaries for parallel processing
    ///
    /// ### Params
    ///
    /// * `num_chunks` - Number of desired chunks
    ///
    /// ### Results
    ///
    /// A vector of tuples indicating the chunk boundaries.
    fn find_chunk_boundaries(&mut self, num_chunks: usize) -> IoResult<Vec<(u64, u64)>> {
        let file_size = self.reader.get_ref().metadata()?.len();

        self.reader.rewind()?;
        let mut line = String::new();
        loop {
            line.clear();
            self.reader.read_line(&mut line)?;
            if !line.starts_with('%') {
                break;
            }
        }

        let data_start = self.reader.stream_position()?;

        let data_size = file_size - data_start;
        let chunk_size = data_size / num_chunks as u64;

        let mut boundaries = vec![(data_start, data_start)];

        // idea is to chunk the data into a vector of boundaries that can
        // be dealt with in parallel via Rayon
        for i in 1..num_chunks {
            let target_pos = data_start + (chunk_size * i as u64);
            if target_pos >= file_size {
                break;
            }

            self.reader.seek(std::io::SeekFrom::Start(target_pos))?;

            let mut byte = [0u8; 1];
            while self.reader.read(&mut byte)? > 0 {
                if byte[0] == b'\n' {
                    break;
                }
            }

            let boundary = self.reader.stream_position()?;
            boundaries.push((boundary, boundary));
        }

        boundaries.push((file_size, file_size));

        for i in 0..boundaries.len() - 1 {
            boundaries[i].1 = boundaries[i + 1].0;
        }
        boundaries.pop();

        Ok(boundaries)
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
