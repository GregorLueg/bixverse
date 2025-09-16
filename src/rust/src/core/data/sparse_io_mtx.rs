use std::fs::File;
use std::io::{BufRead, BufReader, Result as IoResult};
use std::path::Path;

use crate::core::data::sparse_io::{CellGeneSparseWriter, CsrCellChunk};
use crate::single_cell::processing::*;

////////////////////////////////
// Cell accumulator structure //
////////////////////////////////

/// Structure to keep Cell information during MTX parsing
///
/// ### Fields
///
/// * `current_idx` - Current index position.
/// * `gene_indices` - The gene indices for this cell.
/// * `genen_counts` - The raw counts for this cell.
/// * `total_umi` - Total library size of the cell.
#[derive(Debug)]
struct CellAccumulator {
    current_idx: usize,
    gene_indices: Vec<u16>,
    gene_counts: Vec<u16>,
    total_umi: u32,
}

impl CellAccumulator {
    /// Generate a new instace of the accumulator
    ///
    /// ### Returns
    ///
    /// Returns initialised `CellAccumulator`.
    fn new() -> Self {
        Self {
            current_idx: 0,
            gene_indices: Vec::new(),
            gene_counts: Vec::new(),
            total_umi: 0,
        }
    }

    /// Reset the accumulator
    ///
    /// ### Params
    ///
    /// * `new_index` - To which new index to set the accumulator
    #[inline]
    fn reset(&mut self, new_index: usize) {
        self.current_idx = new_index;
        self.gene_indices.clear();
        self.gene_counts.clear();
        self.total_umi = 0;
    }

    /// Add a gene
    ///
    /// * `gene_idx` - Index position of the gene
    /// * `count` - Raw counts for this gene
    #[inline]
    fn add_gene(&mut self, gene_idx: u16, count: u16) {
        self.gene_indices.push(gene_idx);
        self.gene_counts.push(count);
        self.total_umi += count as u32;
    }

    /// Check if the cell has any associated genes
    ///
    /// ### Returns
    ///
    /// Boolean if cell has no associated genes
    #[inline]
    fn is_empty(&self) -> bool {
        self.gene_indices.is_empty()
    }

    /// Does the cell pass thresholds
    ///
    /// Should the cell not pass thresholds based on minimum genes per cell
    /// and library size
    ///
    /// ### Params
    ///
    /// * `min_genes` - Minimum genes that the cell needs to have to pass the
    ///   threshold.
    /// * `min_umi` - Minimum library size that the cell needs to have to pass
    ///   the threshold.
    ///
    /// ### Return
    ///
    /// Boolean indicating if the cell passes the thresholds.
    #[inline]
    fn to_include(&self, min_genes: usize, min_umi: u32) -> bool {
        self.gene_indices.len() >= min_genes && self.total_umi >= min_umi
    }

    /// Return the cell data
    ///
    /// ### Returns
    ///
    /// A tuple of `(library_size, no_genes_expressed)`
    fn get_cell_data(&self) -> (u32, usize) {
        (self.total_umi, self.gene_indices.len())
    }
}

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
    pub fn process_mtx_and_write_bin(mut self, bin_path: &str) -> IoResult<MtxFinalData> {
        let mut writer = CellGeneSparseWriter::new(
            bin_path,
            true, // cell_based
            self.header.total_cells,
            self.header.total_genes,
        )?;

        let mut to_include = vec![false; self.header.total_cells];
        let mut lib_size = Vec::new();
        let mut nnz = Vec::new();
        let mut accumulator = CellAccumulator::new();
        let mut line = String::new();
        let mut cells_written = 0_usize;
        let mut last_processed_cell = 0_usize;

        let mut line_buffer = Vec::with_capacity(64);

        while {
            line_buffer.clear();
            self.reader.read_until(b'\n', &mut line_buffer)? > 0
        } {
            // Remove trailing newline
            if line_buffer.last() == Some(&b'\n') {
                line_buffer.pop();
            }
            if line_buffer.last() == Some(&b'\r') {
                line_buffer.pop();
            }

            if line_buffer.is_empty() {
                continue;
            }

            // Fast parsing without string allocation
            let (row, col, value) = match parse_mtx_line(&line_buffer) {
                Some(parsed) => parsed,
                None => continue, // Skip malformed lines
            };

            let cell_idx = (row - 1) as usize;
            let gene_idx = (col - 1) as u16;

            if !accumulator.is_empty() && cell_idx != accumulator.current_idx {
                let (lib_size_i, nnz_i) = accumulator.get_cell_data();
                lib_size.push(lib_size_i);
                nnz.push(nnz_i);
                cells_written += self.process_cell(
                    &mut accumulator,
                    &mut writer,
                    &mut to_include,
                    &cells_written,
                )?;
                last_processed_cell = accumulator.current_idx;
            }

            if accumulator.is_empty() || cell_idx != accumulator.current_idx {
                accumulator.reset(cell_idx);
            }

            accumulator.add_gene(gene_idx, value);
            line.clear();
        }

        // Process the last accumulated cell
        if !accumulator.is_empty() {
            let (lib_size_i, nnz_i) = accumulator.get_cell_data();
            lib_size.push(lib_size_i);
            nnz.push(nnz_i);
            cells_written += self.process_cell(
                &mut accumulator,
                &mut writer,
                &mut to_include,
                &cells_written,
            )?;
            last_processed_cell = accumulator.current_idx;
        }

        #[allow(clippy::needless_range_loop)]
        for cell_idx in (last_processed_cell + 1)..self.header.total_cells {
            // These cells have no entries, so they automatically fail QC
            to_include[cell_idx] = false;
        }

        writer.update_header_no_cells(cells_written);

        writer.finalise()?;

        Ok(MtxFinalData {
            cell_qc: CellQuality {
                to_keep: to_include,
                lib_size: Some(lib_size.iter().map(|x| *x as usize).collect()),
                no_genes: Some(nnz),
            },
            no_genes: self.header.total_genes,
            no_cells: cells_written,
        })
    }

    /// Helper function to process a given cell
    ///
    /// ### Params
    ///
    /// * `accumulator` - The `CellAccumulator` structure with the cell data
    ///   to write to disk.
    /// * `writer` - The `CellGeneSparseWriter` structure to write the cell
    ///   data to disk.
    /// * `to_include` - Reference to the mutable vector indicating if the cell
    ///   will be included (or not).
    ///
    /// ### Returns
    ///
    /// A usize indicating that the cell was included (`1`) or not (`0`).
    fn process_cell(
        &self,
        accumulator: &mut CellAccumulator,
        writer: &mut CellGeneSparseWriter,
        to_include: &mut [bool],
        cell_written: &usize,
    ) -> IoResult<usize> {
        let cell_idx = accumulator.current_idx;
        let passes_qc = accumulator.to_include(
            self.qc_params.min_unique_genes,
            self.qc_params.min_lib_size as u32,
        );

        to_include[cell_idx] = passes_qc;

        if passes_qc {
            let cell_chunk = CsrCellChunk::from_data(
                &accumulator.gene_counts,
                &accumulator.gene_indices,
                *cell_written,
                self.qc_params.target_size,
                true, // to_keep = true since it passed QC
            );

            writer.write_cell_chunk(cell_chunk)?;
            Ok(1)
        } else {
            Ok(0)
        }
    }
}

/////////////
// Helpers //
/////////////

/// Fast integer parsing from byte slice
///
/// ### Params
///
/// * `bytes` - Slice of the bytes
///
/// ### Return
///
/// Optional u32 if the bytes were an ASCII digit
#[inline]
fn fast_parse_u32(bytes: &[u8]) -> Option<u32> {
    let mut result = 0u32;
    for &byte in bytes {
        if byte.is_ascii_digit() {
            result = result.wrapping_mul(10).wrapping_add((byte - b'0') as u32);
        } else {
            return None;
        }
    }
    Some(result)
}

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
    let mut parts = line.split(|&b| b == b' ' || b == b'\t');

    let row = fast_parse_u32(parts.next()?)?;
    let col = fast_parse_u32(parts.next()?)?;
    let val = fast_parse_u32(parts.next()?)? as u16;

    Some((row, col, val))
}
