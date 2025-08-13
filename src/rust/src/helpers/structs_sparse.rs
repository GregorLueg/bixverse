use bincode::{config, decode_from_slice, serde::encode_to_vec, Decode, Encode};
use faer::traits::ComplexField;
use faer::{Mat, MatRef};
use half::f16;
use rustc_hash::FxHashMap;
use serde::{Deserialize, Serialize};
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use std::iter::Sum;

////////////////
// Structures //
////////////////

/// Structure for SparseColumnMatrices
///
/// ### Fields
///
/// * `data` - Vector with the data.
/// * `row_indices` - The row indices of the data.
/// * `col_ptrs` - The column pointers of the data.
/// * `ncol` - Original number of columns.
/// * `nrow` - Original number of rows.
#[derive(Debug, Clone)]
pub struct SparseColumnMatrix<T> {
    pub data: Vec<T>,
    pub row_indices: Vec<usize>,
    pub col_ptrs: Vec<usize>,
    pub ncol: usize,
    pub nrow: usize,
}

#[allow(dead_code)]
impl<T> SparseColumnMatrix<T>
where
    T: Clone + Default + PartialEq + ComplexField + From<f64>,
{
    /// Generate a new sparse column matrix from values pre-computed data
    ///
    /// ### Params
    ///
    /// * `data` - Slice of the data.
    /// * `row_indices` - Slice of the row indices of the data.
    /// * `col_ptrs` - Slice of the column pointers of the data.
    /// * `ncol` - Original number of columns.
    /// * `nrow` - Original number of rows.
    pub fn new(
        data: &[T],
        row_indices: &[usize],
        col_ptrs: &[usize],
        ncol: usize,
        nrow: usize,
    ) -> Self {
        Self {
            data: data.to_vec(),
            row_indices: row_indices.to_vec(),
            col_ptrs: col_ptrs.to_vec(),
            ncol,
            nrow,
        }
    }

    /// Convert a faer dense matrix to sparse column format
    ///
    /// ### Params
    ///
    /// * `dense` - The original dense matrix.
    pub fn from_dense_matrix(dense: MatRef<T>) -> Self {
        let ncol = dense.ncols();
        let nrow = dense.nrows();

        let mut values = Vec::new();
        let mut row_indices = Vec::new();
        let mut col_ptrs = Vec::with_capacity(ncol + 1);

        col_ptrs.push(0_usize);

        for col in 0..ncol {
            for row in 0..nrow {
                let value = dense.get(row, col).clone();
                if value != T::default() {
                    values.push(value);
                    row_indices.push(row);
                }
            }
            col_ptrs.push(values.len());
        }

        Self {
            data: values,
            row_indices,
            col_ptrs,
            ncol,
            nrow,
        }
    }

    /// To a dense faer matrix
    ///
    /// ### Returns
    ///
    /// Returns a dense faer matrix.
    pub fn to_dense_matrix(&self) -> Mat<T> {
        let mut dense = Mat::zeros(self.nrow, self.ncol);

        for col in 0..self.ncol {
            let start = self.col_ptrs[col];
            let end = self.col_ptrs[col + 1];

            for idx in start..end {
                let row = self.row_indices[idx];
                let value = &self.data[idx];
                *dense.get_mut(row, col) = value.clone();
            }
        }

        dense
    }

    /// Return the number of non-zero values
    ///
    /// ### Returns
    ///
    /// Return the total number of NNZ values in the data
    pub fn nnz(&self) -> usize {
        self.data.len()
    }

    /// Create a sparse matrix from a row major upper triangle stored value
    ///
    /// This is a helper function to transform potentially sparse symmetric matrices
    /// stored as upper-triangles into a sparse matrix format in Rust.
    ///
    /// ### Params
    ///
    /// * `upper_triangle` - Represents the values of the upper triangle in
    ///   row major formant
    /// * `n` - Original nrows and ncols.
    /// * `include_diagonal` - Are the diagonal values included.
    pub fn from_upper_triangle_sym(upper_triangle: &[T], n: usize, include_diagonal: bool) -> Self {
        let mut values = Vec::new();
        let mut row_indices = Vec::new();
        let mut col_ptrs = Vec::with_capacity(n + 1);

        col_ptrs.push(0);

        for col in 0..n {
            for row in 0..n {
                let value = if row == col && !include_diagonal {
                    T::from(1.0)
                } else if row < col {
                    // Upper triangle (row < col)
                    let offset = if include_diagonal {
                        row * n - row * (row + 1) / 2 + col
                    } else {
                        row * (n - 1) - row * (row + 1) / 2 + col - 1
                    };
                    upper_triangle[offset].clone()
                } else if row == col {
                    // Diagonal (when included)
                    let offset = row * n - row * (row + 1) / 2 + col;
                    upper_triangle[offset].clone()
                } else {
                    // Lower triangle (col < row) - symmetric
                    let offset = if include_diagonal {
                        col * n - col * (col + 1) / 2 + row
                    } else {
                        col * (n - 1) - col * (col + 1) / 2 + row - 1
                    };
                    upper_triangle[offset].clone()
                };

                if value != T::default() {
                    values.push(value);
                    row_indices.push(row);
                }
            }
            col_ptrs.push(values.len());
        }

        Self {
            data: values,
            row_indices,
            col_ptrs,
            ncol: n,
            nrow: n,
        }
    }
}

////////////////
// Conversion //
////////////////

/// A type alias representing effect size results
///
/// ### Fields
///
/// * `0` - The data in the CSR format
/// * `1` - The row pointers
/// * `2` - The column indices
/// * `3` - An optional second data layer
#[allow(dead_code)]
pub type CsrData<T, U = T> = (Vec<T>, Vec<usize>, Vec<usize>, Option<Vec<U>>);

/// A type alias representing effect size results
///
/// ### Fields
///
/// * `0` - The data in the CSR format
/// * `1` - The column pointers
/// * `2` - The row indices
/// * `3` - An optional second data layer
pub type CscData<T, U = T> = (Vec<T>, Vec<usize>, Vec<usize>, Option<Vec<U>>);

/// Transform CSC stored data into CSR stored data
///
/// This version does a full memory copy of the data
///
/// ### Params
///
/// * `data` - The data stored in CSC format.
/// * `row_ind` - The row indices from the CSC format.
/// * `col_ptr` - The column pointers from the CSC format.
/// * `nrows` - The number of rows in the data.
/// * `data2` - An optional second data layer.
///
/// ### Returns
///
/// The data in CSR format, i.e., `CsrData`
#[allow(dead_code)]
pub fn csc_to_csr<T: Clone + Default, U: Clone + Default>(
    data: &[T],
    row_ind: &[usize],
    col_ptr: &[usize],
    nrows: usize,
    data2: Option<&[U]>,
) -> CsrData<T, U> {
    let nnz = data.len();
    let mut row_ptr = vec![0; nrows + 1];

    for &r in row_ind {
        row_ptr[r + 1] += 1;
    }

    for i in 0..nrows {
        row_ptr[i + 1] += row_ptr[i];
    }

    let mut csr_data = vec![T::default(); nnz];
    let mut csr_data2 = data2.map(|_| vec![U::default(); nnz]);
    let mut csr_col_ind = vec![0; nnz];
    let mut next = row_ptr[..nrows].to_vec();

    for col in 0..(col_ptr.len() - 1) {
        for idx in col_ptr[col]..col_ptr[col + 1] {
            let row = row_ind[idx];
            let pos = next[row];
            csr_data[pos] = data[idx].clone();
            csr_col_ind[pos] = col;

            // Add the second layer
            if let (Some(d2), Some(ref mut csr_d2)) = (data2, &mut csr_data2) {
                csr_d2[pos] = d2[idx].clone();
            }

            next[row] += 1;
        }
    }

    (csr_data, row_ptr, csr_col_ind, csr_data2)
}

/// Transform CSR stored data into CSC stored data
///
/// This version does a full memory copy of the data
///
/// ### Params
///
/// * `data` - The data stored in CSR format.
/// * `col_ind` - The column indices from the CSR format.
/// * `row_ptr` - The row pointers from the CSR format.
/// * `ncols` - The number of columns in the data.
/// * `data2` - An optional second data layer.
///
/// ### Returns
///
/// The data in CSC format, i.e., `CscData`
pub fn csr_to_csc<T: Clone + Default, U: Clone + Default>(
    data: &[T],
    col_ind: &[usize],
    row_ptr: &[usize],
    ncols: usize,
    data2: Option<&[U]>,
) -> CscData<T, U> {
    let nnz = data.len();
    let mut col_ptr = vec![0; ncols + 1];

    // Count occurrences per column
    for &c in col_ind {
        col_ptr[c + 1] += 1;
    }

    // Cumulative sum to get column pointers
    for i in 0..ncols {
        col_ptr[i + 1] += col_ptr[i];
    }

    let mut csc_data = vec![T::default(); nnz];
    let mut csc_data2 = data2.map(|_| vec![U::default(); nnz]);
    let mut csc_row_ind = vec![0; nnz];
    let mut next = col_ptr[..ncols].to_vec();

    // Iterate through rows and place data in CSC format
    for row in 0..(row_ptr.len() - 1) {
        for idx in row_ptr[row]..row_ptr[row + 1] {
            let col = col_ind[idx];
            let pos = next[col];
            csc_data[pos] = data[idx].clone();
            csc_row_ind[pos] = row;

            // Add second layer if there
            if let (Some(d2), Some(ref mut csc_d2)) = (data2, &mut csc_data2) {
                csc_d2[pos] = d2[idx].clone();
            }

            next[col] += 1;
        }
    }

    (csc_data, col_ptr, csc_row_ind, csc_data2)
}

///////////////////////////
// Sparse data streaming //
///////////////////////////

////////////
// Traits //
////////////

#[derive(Encode, Decode, Serialize, Deserialize, Debug, Clone, Copy, Default)]
pub struct F16(u16);

impl From<half::f16> for F16 {
    fn from(f: half::f16) -> Self {
        F16(f.to_bits())
    }
}

impl From<F16> for half::f16 {
    fn from(f: F16) -> Self {
        half::f16::from_bits(f.0)
    }
}

impl Sum for F16 {
    fn sum<I: Iterator<Item = F16>>(iter: I) -> Self {
        let sum: half::f16 = iter.map(half::f16::from).sum();
        F16::from(sum)
    }
}

impl<'a> Sum<&'a F16> for F16 {
    fn sum<I: Iterator<Item = &'a F16>>(iter: I) -> Self {
        let sum: half::f16 = iter.map(|&f| half::f16::from(f)).sum();
        F16::from(sum)
    }
}

/////////////////////
// Data structures //
/////////////////////

/// CsrCellChunk
///
/// This structure is designed to store the data of a single cell in a
/// CSR-like format optimised for rapid access on disk.
///
/// ### Fields
///
/// * `data_raw` - Array of the raw counts of this cell.
/// * `data_norm` - Array of the normalised counts of this cell.
/// * `library_size` - Total library size of the cell.
/// * `col_indices` - The col indices of the genes.
/// * `original_index` - Original (row) index of the cell.
/// * `to_keep` - Flat if the cell should be included in certain analysis.
#[derive(Encode, Decode, Serialize, Deserialize, Debug)]
pub struct CsrCellChunk {
    pub data_raw: Vec<u16>,
    pub data_norm: Vec<F16>,
    pub library_size: usize,
    pub col_indices: Vec<u16>,
    pub original_index: usize,
    pub to_keep: bool,
}

impl CsrCellChunk {
    /// Function to generate the chunk from R data
    pub fn from_r_data(
        data: &[i32],
        col_idx: &[i32],
        original_index: usize,
        target_size: f32,
    ) -> Self {
        let data_f32 = data.iter().map(|x| *x as f32).collect::<Vec<f32>>();
        let sum = data_f32.iter().sum::<f32>();
        let data_norm: Vec<F16> = data_f32
            .into_iter()
            .map(|x| {
                let norm = (x / sum * target_size).ln_1p();
                F16::from(f16::from_f32(norm))
            })
            .collect();

        Self {
            data_raw: data.iter().map(|x| *x as u16).collect::<Vec<u16>>(),
            data_norm,
            library_size: sum as usize,
            col_indices: col_idx.iter().map(|x| *x as u16).collect::<Vec<u16>>(),
            original_index,
            to_keep: true,
        }
    }
}

/// CscGeneChunk
///
/// This structure is designed to store the data of a single gene in a
/// CSC-like format optimised for rapid access on disk.
///
/// ### Fields
///
/// * `data_raw` - Vector with the raw data.
/// * `data_norm` - Vector with the normalised data (library size adjusted and
///   log-normalised).
/// * `avg_exp` - Vector with average expression.
/// * `nnz` - Number non-zero values.
/// * `row_indices` - The column indices of the data.
/// * `original_index` - Original index of the gene.
/// * `to_keep` - Boolean if the gene should be included into anything.
#[derive(Encode, Decode, Serialize, Deserialize, Debug)]
pub struct CscGeneChunk {
    pub data_raw: Vec<u16>,
    pub data_norm: Vec<F16>,
    pub avg_exp: F16,
    pub nnz: usize,
    pub row_indices: Vec<u32>,
    pub original_index: usize,
    pub to_keep: bool,
}

#[allow(dead_code)]
impl CscGeneChunk {
    pub fn from_conversion(
        data: &[u16],
        data_2: &[F16],
        col_idx: &[usize],
        original_index: usize,
    ) -> Self {
        let avg_exp = data_2.iter().sum::<F16>();
        let nnz = data.len();

        Self {
            data_raw: data.to_vec(),
            data_norm: data_2.to_vec(),
            avg_exp,
            nnz,
            row_indices: col_idx.iter().map(|x| *x as u32).collect::<Vec<u32>>(),
            original_index,
            to_keep: true,
        }
    }
}

/// SparseDataHeader
///
/// Stores the information in terms of total cells, total genes, number of
/// chunks in terms of cells and genes and the offset vectors
///
/// ### Params
///
/// * `total_cells` - Total number of cells in the experiment.
/// * `total_genes` - Total number of genes in the experiemnt.
/// * `cell_based` - Boolean. If `true` the data stores cells; if `false` the
///   data stores genes.
/// * `no_chunks` - No of chunks that store either cell or gene data.
/// * `chunk_offsets` - Vector containing the offsets for the cell or gene
///   chunks.
/// * `cell_index_map` - FxHashMap with the original index -> chunk info
#[derive(Encode, Decode, Serialize, Deserialize)]
pub struct SparseDataHeader {
    pub total_cells: usize,
    pub total_genes: usize,
    pub cell_based: bool,
    pub no_chunks: usize,
    pub chunk_offsets: Vec<u64>,
    pub index_map: FxHashMap<usize, usize>,
}

/// Fixed-size file header that points to the main header location
///
/// ### Params
///
/// * `magic` - Magic string as bytes to recognise the file
/// * `version` - Version of the file
/// * `main_header_offset` - Offset of the main header, i.e., 64 bytes
/// * `_reserved_1` - 32 additional reserved bytes for the future
/// * `_reserved_2` - 4 additional reserved bytes for the future
#[repr(C)]
#[derive(Encode, Decode, Serialize, Deserialize)]
struct FileHeader {
    magic: [u8; 8],
    version: u32,
    main_header_offset: u64,
    cell_based: bool,
    // Needs to be split into two arrays to get to 64 bytes
    _reserved_1: [u8; 32],
    _reserved_2: [u8; 3],
}

impl FileHeader {
    /// Generate a new header
    ///
    /// ### Params
    ///
    /// * `cell_based` - Is the data stored for fast cell retrieval.
    ///
    /// ### Returns
    ///
    /// A new object of `FileHeader`
    fn new(cell_based: bool) -> Self {
        Self {
            magic: *b"SCRNASEQ",
            version: 1,
            main_header_offset: 0,
            cell_based,
            _reserved_1: [0; 32],
            _reserved_2: [0; 3],
        }
    }
}

//////////////////////
// Streaming writer //
//////////////////////

/// CellGeneSparseWriter
///
/// Implementation of a structure for writing in a streamed manner two different
/// types of sparse stored data.
///
/// ### Params
///
/// * `header` - The header of the file.
/// * `writer` - BufWriter to the file.
/// * `chunks_start_pos` - The current position of the chunks.
/// * `cell_based` - Boolean indicating if the writer is designed to write in
///   an efficient manner for cells.
pub struct CellGeneSparseWriter {
    header: SparseDataHeader,
    writer: BufWriter<File>,
    chunks_start_pos: u64,
    cell_based: bool,
}

#[allow(dead_code)]
impl CellGeneSparseWriter {
    /// Create a new sparse writer instance
    ///
    /// This writer assumes that rows represent genes and columns represent
    /// cells.
    ///
    /// ### Params
    ///
    /// * `path_f` - Path to the .bin file to which to write to.
    /// * `cell_based` - Shall the writer be set up for writing cell-based
    ///   (`true`) or gene-based chunks.
    /// * `total_cells` - Total cells in the data.
    /// * `total_genes` - Total genes in the data.
    pub fn new(
        path_f: &str,
        cell_based: bool,
        total_cells: usize,
        total_genes: usize,
    ) -> std::io::Result<Self> {
        let file = File::create(path_f)?;
        let mut writer = BufWriter::new(file);

        let file_header = FileHeader::new(cell_based);
        let file_header_enc = encode_to_vec(&file_header, config::standard()).unwrap();
        if file_header_enc.len() < 64 {
            writer.write_all(&file_header_enc)?;
            writer.write_all(&vec![0u8; 64 - file_header_enc.len()])?;
        } else {
            writer.write_all(&file_header_enc[..64])?;
        }
        writer.flush()?;

        let chunks_start_pos = 64;

        let header = SparseDataHeader {
            total_cells,
            total_genes,
            cell_based,
            no_chunks: 0,
            chunk_offsets: Vec::new(),
            index_map: FxHashMap::default(),
        };

        Ok(Self {
            header,
            writer,
            chunks_start_pos,
            cell_based,
        })
    }

    /// Write a Cell (Chunk) to the file
    ///
    /// This function will panic if the file was not set to cell-based! The
    /// data is represented in a CSR-type format.
    ///
    /// ### Params
    ///
    /// * `cell_chunk` - The data representing that specific cell.
    pub fn write_cell_chunk(&mut self, cell_chunk: CsrCellChunk) -> std::io::Result<()> {
        assert!(
            self.cell_based,
            "The writer is not set to write in a cell-based manner!"
        );

        // get current position in file
        let current_pos = self.writer.stream_position()?;

        // store offset relative to chunks start
        let chunk_offset = current_pos - self.chunks_start_pos;
        self.header.chunk_offsets.push(chunk_offset);

        self.header
            .index_map
            .insert(cell_chunk.original_index, self.header.no_chunks);

        let encoded = encode_to_vec(&cell_chunk, config::standard()).unwrap();
        let chunk_size = encoded.len() as u64;

        self.writer.write_all(&chunk_size.to_le_bytes())?;
        self.writer.write_all(&encoded)?;
        self.writer.flush()?;

        self.header.no_chunks += 1;

        Ok(())
    }

    /// Write a Gene to the file
    ///
    /// This function will panic if the file was set to cell-based!
    ///
    /// ### Params
    ///
    /// * `gene_chunk` - The data representing that specific gene.
    pub fn write_gene_chunk(&mut self, gene_chunk: CscGeneChunk) -> std::io::Result<()> {
        assert!(
            !self.cell_based,
            "The writer is not set to write in a gene-based manner!"
        );

        // get current position in file
        let current_pos = self.writer.stream_position()?;

        // store offset relative to chunks start
        let chunk_offset = current_pos - self.chunks_start_pos;
        self.header.chunk_offsets.push(chunk_offset);

        self.header
            .index_map
            .insert(gene_chunk.original_index, self.header.no_chunks);

        let encoded = encode_to_vec(&gene_chunk, config::standard()).unwrap();
        let chunk_size = encoded.len() as u64;

        self.writer.write_all(&chunk_size.to_le_bytes())?;
        self.writer.write_all(&encoded)?;
        self.writer.flush()?;

        self.header.no_chunks += 1;
        Ok(())
    }

    /// Finalise the file
    pub fn finalise(mut self) -> std::io::Result<()> {
        // write header size and header
        let main_header_offset = self.writer.stream_position()?;

        let header_data = encode_to_vec(&self.header, config::standard()).unwrap();
        let header_size = header_data.len() as u64;

        self.writer.write_all(&header_size.to_le_bytes())?;
        self.writer.write_all(&header_data)?;

        self.writer.seek(SeekFrom::Start(0))?;
        let mut file_header = FileHeader::new(self.cell_based);
        file_header.main_header_offset = main_header_offset;
        let file_header_enc = encode_to_vec(&file_header, config::standard()).unwrap();

        // ensure it's exactly 64 bytes
        if file_header_enc.len() < 64 {
            self.writer.write_all(&file_header_enc)?;
            self.writer
                .write_all(&vec![0u8; 64 - file_header_enc.len()])?;
        } else {
            self.writer.write_all(&file_header_enc[..64])?;
        }

        self.writer.flush()?;

        Ok(())
    }
}

//////////////////////
// Streaming reader //
//////////////////////

/// Iterator for reading cell chunks sequentially
///
/// ### Fields
///
/// * `reader` - A borrowed reference to the `StreamingSparseReader`.
/// * `current_chunk` - The current chunk
pub struct CellChunkIterator<'a> {
    reader: &'a mut StreamingSparseReader,
    current_chunk: usize,
}

impl<'a> Iterator for CellChunkIterator<'a> {
    type Item = std::io::Result<CsrCellChunk>;

    /// Get the next element of the iterator
    fn next(&mut self) -> Option<Self::Item> {
        if self.current_chunk >= self.reader.header.no_chunks {
            return None;
        }

        let result = self.reader.read_cell_chunk_at(self.current_chunk);
        self.current_chunk += 1;

        match result {
            Ok(Some(chunk)) => Some(Ok(chunk)),
            Ok(None) => None,
            Err(e) => Some(Err(e)),
        }
    }
}

/// Iterator for reading gene chunks sequentially
///
/// ### Fields
///
/// * `reader` - A borrowed reference to the `StreamingSparseReader`.
/// * `current_chunk` - The current chunk
pub struct GeneChunkIterator<'a> {
    reader: &'a mut StreamingSparseReader,
    current_chunk: usize,
}

impl<'a> Iterator for GeneChunkIterator<'a> {
    type Item = std::io::Result<CscGeneChunk>;

    /// Get the next element of the iterator
    fn next(&mut self) -> Option<Self::Item> {
        if self.current_chunk >= self.reader.header.no_chunks {
            return None;
        }

        let result = self.reader.read_gene_chunk_at(self.current_chunk);
        self.current_chunk += 1;

        match result {
            Ok(Some(chunk)) => Some(Ok(chunk)),
            Ok(None) => None,
            Err(e) => Some(Err(e)),
        }
    }
}

/// StreamingSparseReader
///
/// Implementation of a structure for writing in a streamed manner
///
/// ### Params
///
/// * `header` - The header of the file.
/// * `reader` - BufReader of the file.
/// * `chunks_start` - Position of where chunks start.
#[allow(dead_code)]
pub struct StreamingSparseReader {
    header: SparseDataHeader,
    reader: BufReader<File>,
    chunks_start: u64,
}

#[allow(dead_code)]
impl StreamingSparseReader {
    /// Generate a new StreamingReader
    ///
    /// ### Params
    ///
    /// * `new`
    pub fn new(f_path: &str) -> std::io::Result<Self> {
        let file = File::open(f_path)?;
        let mut reader = BufReader::new(file);

        // read file header (64 bytes)
        let mut file_header_buf = [0u8; 64];
        reader.read_exact(&mut file_header_buf)?;
        let (file_header, _) = decode_from_slice::<FileHeader, _>(
            &file_header_buf,
            config::standard(),
        )
        .map_err(|_| {
            std::io::Error::new(std::io::ErrorKind::InvalidData, "File header decode failed")
        })?;

        // verify magic number <- nice touch by Claude!
        if &file_header.magic != b"SCRNASEQ" {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "Invalid file format",
            ));
        }

        // seek to main header
        reader.seek(SeekFrom::Start(file_header.main_header_offset))?;

        // read main header size
        let mut size_buf = [0u8; 8];
        reader.read_exact(&mut size_buf)?;
        let header_size = u64::from_le_bytes(size_buf) as usize;

        // read main header
        let mut header_buf = vec![0u8; header_size];
        reader.read_exact(&mut header_buf)?;

        let (header, _) = decode_from_slice::<SparseDataHeader, _>(&header_buf, config::standard())
            .map_err(|_| {
                std::io::Error::new(std::io::ErrorKind::InvalidData, "Header decode failed")
            })?;

        let chunks_start = 64; // after file header

        Ok(Self {
            header,
            reader,
            chunks_start,
        })
    }

    ///////////
    // Cells //
    ///////////

    /// Read a specific cell by its original index
    ///
    /// This function will panic if the file was not designed for cell-based
    /// read/write!
    ///
    /// ### Param
    ///
    /// * `original_index` - Original index position of the cell.
    ///
    /// ### Returns
    ///
    /// `Result<Option<CsrCellChunk>>`
    pub fn read_cell_by_index(
        &mut self,
        original_index: usize,
    ) -> std::io::Result<Option<CsrCellChunk>> {
        assert!(
            self.header.cell_based,
            "The file was not designed for cell-based read/write!"
        );

        // find chunk index for this cell
        let chunk_index = match self.header.index_map.get(&original_index) {
            Some(idx) => *idx,
            None => return Ok(None),
        };

        self.read_cell_chunk_at(chunk_index)
    }

    /// Read multiple cells by their original indices
    ///
    /// This function will panic if the file was not designed for cell-based
    /// read/write!
    ///
    /// ### Params
    ///
    /// * `indices` - A slice of positions to retrieve
    ///
    /// ### Returns
    ///
    /// A vector of `CsrCellChunk`
    pub fn read_cells_by_indices(
        &mut self,
        indices: &[usize],
    ) -> std::io::Result<Vec<CsrCellChunk>> {
        assert!(
            self.header.cell_based,
            "The file was not designed for cell-based read/write!"
        );
        let mut cells = Vec::with_capacity(indices.len());
        for &idx in indices {
            match self.read_cell_by_index(idx)? {
                Some(cell) => cells.push(cell),
                None => {
                    return Err(std::io::Error::new(
                        std::io::ErrorKind::NotFound,
                        format!("Cell index {} not found", idx),
                    ))
                }
            }
        }
        Ok(cells)
    }

    /// Read a cell chunk at a specific chunk index
    ///
    /// This function will panic if the file was not designed for cell-based
    /// read/write! If the index is larger than the what is written in the file,
    /// it will return `None`.
    ///
    /// ### Params
    ///
    /// * `chunk_index` - Index postion of the chunk.
    ///
    /// ### Returns
    ///
    /// `Result<Option<CsrCellChunk>>`
    pub fn read_cell_chunk_at(
        &mut self,
        chunk_index: usize,
    ) -> std::io::Result<Option<CsrCellChunk>> {
        assert!(
            self.header.cell_based,
            "The file was not designed for cell-based read/write!"
        );

        if chunk_index >= self.header.no_chunks {
            return Ok(None);
        }

        // Calculate absolute position
        let chunk_offset = self.chunks_start + self.header.chunk_offsets[chunk_index];
        self.reader.seek(SeekFrom::Start(chunk_offset))?;

        // Read chunk size
        let mut size_buf = [0u8; 8];
        self.reader.read_exact(&mut size_buf)?;
        let chunk_size = u64::from_le_bytes(size_buf) as usize;

        // Read chunk data
        let mut chunk_buf = vec![0u8; chunk_size];
        self.reader.read_exact(&mut chunk_buf)?;

        let (chunk, _) = decode_from_slice(&chunk_buf, config::standard()).unwrap();

        Ok(Some(chunk))
    }

    /// Iterator over all cell chunks
    ///
    /// ### Returns
    ///
    /// An iterator (`CellChunkIterator`) that can be subsequently consumed
    pub fn iter_cells(&mut self) -> CellChunkIterator<'_> {
        assert!(
            self.header.cell_based,
            "The file was not designed for cell-based read/write!"
        );

        CellChunkIterator {
            reader: self,
            current_chunk: 0,
        }
    }

    ///////////
    // Genes //
    ///////////

    /// Read a specific gene by its original index
    ///
    /// This function will panic if the file was not designed for gene-based
    /// read/write!
    ///
    /// ### Param
    ///
    /// * `original_index` - Original index position of the cell.
    ///
    /// ### Returns
    ///
    /// `Result<Option<CscGeneChunk>>`
    pub fn read_gene_by_index(
        &mut self,
        original_index: usize,
    ) -> std::io::Result<Option<CscGeneChunk>> {
        assert!(
            !self.header.cell_based,
            "The file was not designed for gene-based read/write!"
        );

        // find chunk index for this cell
        let chunk_index = match self.header.index_map.get(&original_index) {
            Some(idx) => *idx,
            None => return Ok(None),
        };

        self.read_gene_chunk_at(chunk_index)
    }

    /// Read multiple cells by their original indices
    ///
    /// This function will panic if the file was not designed for cell-based
    /// read/write!
    ///
    /// ### Params
    ///
    /// * `indices` - A slice of positions to retrieve
    ///
    /// ### Returns
    ///
    /// A vector of `CscGeneChunk`
    pub fn read_genes_by_indices(
        &mut self,
        indices: &[usize],
    ) -> std::io::Result<Vec<CscGeneChunk>> {
        assert!(
            !self.header.cell_based,
            "The file was not designed for cell-based read/write!"
        );

        let mut genes = Vec::new();

        for &idx in indices {
            if let Some(cell) = self.read_gene_chunk_at(idx)? {
                genes.push(cell);
            }
        }

        Ok(genes)
    }

    /// Read a gene chunk at a specific chunk index
    ///
    /// This function will panic if the file was not designed for gene-based
    /// read/write! If the index is larger than the what is written in the file,
    /// it will return `None`.
    ///
    /// ### Params
    ///
    /// * `chunk_index` - Index postion of the chunk.
    ///
    /// ### Returns
    ///
    /// `Result<Option<CscGeneChunk>>`
    pub fn read_gene_chunk_at(
        &mut self,
        chunk_index: usize,
    ) -> std::io::Result<Option<CscGeneChunk>> {
        assert!(
            !self.header.cell_based,
            "The file was not designed for gene-based read/write!"
        );

        if chunk_index >= self.header.no_chunks {
            return Ok(None);
        }

        // Calculate absolute position
        let chunk_offset = self.chunks_start + self.header.chunk_offsets[chunk_index];
        self.reader.seek(SeekFrom::Start(chunk_offset))?;

        // Read chunk size
        let mut size_buf = [0u8; 8];
        self.reader.read_exact(&mut size_buf)?;
        let chunk_size = u64::from_le_bytes(size_buf) as usize;

        // Read chunk data
        let mut chunk_buf = vec![0u8; chunk_size];
        self.reader.read_exact(&mut chunk_buf)?;

        let (chunk, _) = decode_from_slice(&chunk_buf, config::standard()).unwrap();

        Ok(Some(chunk))
    }

    /// Iterator over all gene chunks
    ///
    /// ### Returns
    ///
    /// An iterator (`GeneChunkIterator`) that can be subsequently consumed
    pub fn iter_genes(&mut self) -> GeneChunkIterator<'_> {
        assert!(
            !self.header.cell_based,
            "The file was not designed for gene-based read/write!"
        );

        GeneChunkIterator {
            reader: self,
            current_chunk: 0,
        }
    }

    /////////////
    // General //
    /////////////

    /// Get the header of the file
    ///
    /// ### Returns
    ///
    /// The borrowed `SparseDataHeader`
    pub fn get_header(&self) -> &SparseDataHeader {
        &self.header
    }

    /// Get all available cell or gene indices
    ///
    /// ### Returns
    ///
    /// A vector of usizes with all indices.
    pub fn get_available_indices(&self) -> Vec<usize> {
        let mut indices: Vec<usize> = self.header.index_map.keys().copied().collect();
        indices.sort();
        indices
    }
}

///////////
// Tests //
///////////

#[cfg(test)]
mod tests {
    use super::*;
    use faer::mat;

    #[test]
    fn test_dense_to_sparse_conversion() {
        let dense = mat![[1.0, 0.0, 3.0], [0.0, 2.0, 0.0], [4.0, 0.0, 5.0]];

        let sparse_obj = SparseColumnMatrix::from_dense_matrix(dense.as_ref());

        assert_eq!(sparse_obj.nrow, 3);
        assert_eq!(sparse_obj.ncol, 3);
        assert_eq!(sparse_obj.nnz(), 5);

        let expected_values = vec![1.0, 4.0, 2.0, 3.0, 5.0];

        assert_eq!(sparse_obj.data, expected_values);
    }

    #[test]
    fn test_dense_to_sparse_to_dense_conversion() {
        let dense = mat![[1.0, 0.0, 3.0], [0.0, 2.0, 0.0], [4.0, 0.0, 5.0]];

        let sparse_obj = SparseColumnMatrix::from_dense_matrix(dense.as_ref());

        let redense = sparse_obj.to_dense_matrix();

        assert_eq!(dense, redense);
    }

    #[test]
    fn test_raw_to_dense_conversion() {
        let data = vec![1.0, 4.0, 2.0, 3.0, 5.0];
        let row_indices: Vec<usize> = vec![0, 2, 1, 0, 2];
        let col_ptr: Vec<usize> = vec![0, 2, 3, 5];

        let sparse_obj = SparseColumnMatrix::new(&data, &row_indices, &col_ptr, 3, 3);

        let dense = mat![[1.0, 0.0, 3.0], [0.0, 2.0, 0.0], [4.0, 0.0, 5.0]];

        let redense = sparse_obj.to_dense_matrix();

        assert_eq!(dense, redense);
    }

    #[test]
    fn test_from_upper_triangle_vec_with_diagonal() {
        // Upper triangle with diagonal: [1, 0.8, 0.6, 1, 0.3, 1]
        // Represents correlation matrix:
        // [1.0, 0.8, 0.6]
        // [0.8, 1.0, 0.3]
        // [0.6, 0.3, 1.0]
        let upper_tri = vec![1.0, 0.8, 0.6, 1.0, 0.3, 1.0];
        let sparse = SparseColumnMatrix::from_upper_triangle_sym(&upper_tri, 3, true);

        let dense = sparse.to_dense_matrix();
        let expected = mat![[1.0, 0.8, 0.6], [0.8, 1.0, 0.3], [0.6, 0.3, 1.0]];

        assert_eq!(dense, expected);
    }

    #[test]
    fn test_from_upper_triangle_vec_without_diagonal() {
        // Upper triangle without diagonal: [0.8, 0.6, 0.3]
        // With implied diagonal of 1's, represents:
        // [1.0, 0.8, 0.6]
        // [0.8, 1.0, 0.3]
        // [0.6, 0.3, 1.0]
        let upper_tri = vec![0.8, 0.6, 0.3];
        let sparse = SparseColumnMatrix::from_upper_triangle_sym(&upper_tri, 3, false);

        let dense = sparse.to_dense_matrix();
        let expected = mat![[1.0, 0.8, 0.6], [0.8, 1.0, 0.3], [0.6, 0.3, 1.0]];

        assert_eq!(dense, expected);
    }

    #[test]
    fn test_csc_to_csr_conversion_single_layer() {
        // Test matrix:
        // [1 0 2]
        // [0 3 0]
        // [4 0 5]
        // CSC format
        let data = vec![1, 4, 3, 2, 5];
        let row_ind = vec![0, 2, 1, 0, 2];
        let col_ptr = vec![0, 2, 3, 5];
        let (csr_data, csr_row_ptr, csr_col_ind, csr_data2) =
            csc_to_csr(&data, &row_ind, &col_ptr, 3, None::<&[i32]>);

        // Expected CSR format
        assert_eq!(csr_data, vec![1, 2, 3, 4, 5]);
        assert_eq!(csr_col_ind, vec![0, 2, 1, 0, 2]);
        assert_eq!(csr_row_ptr, vec![0, 2, 3, 5]);
        assert_eq!(csr_data2, None);
    }

    #[test]
    fn test_csc_to_csr_conversion_dual_layer() {
        // Test matrix:
        // [1 0 2]    [10 0  20]
        // [0 3 0]    [0  30 0 ]
        // [4 0 5]    [40 0  50]
        // CSC format for both layers
        let data = vec![1, 4, 3, 2, 5];
        let data2 = vec![1.0, 4.0, 3.0, 2.0, 5.0];
        let row_ind = vec![0, 2, 1, 0, 2];
        let col_ptr = vec![0, 2, 3, 5];
        let (csr_data, csr_row_ptr, csr_col_ind, csr_data2) =
            csc_to_csr(&data, &row_ind, &col_ptr, 3, Some(&data2));

        // Expected CSR format
        assert_eq!(csr_data, vec![1, 2, 3, 4, 5]);
        assert_eq!(csr_col_ind, vec![0, 2, 1, 0, 2]);
        assert_eq!(csr_row_ptr, vec![0, 2, 3, 5]);
        assert_eq!(csr_data2, Some(vec![1.0, 2.0, 3.0, 4.0, 5.0]));
    }

    #[test]
    fn test_csr_to_csc_conversion_single_layer() {
        // Test matrix:
        // [1 0 2]
        // [0 3 0]
        // [4 0 5]
        // CSR format
        let data = vec![1, 2, 3, 4, 5];
        let col_ind = vec![0, 2, 1, 0, 2];
        let row_ptr = vec![0, 2, 3, 5];
        let (csc_data, csc_col_ptr, csc_row_ind, csc_data2) =
            csr_to_csc(&data, &col_ind, &row_ptr, 3, None::<&[i32]>);

        // Expected CSC format
        assert_eq!(csc_data, vec![1, 4, 3, 2, 5]);
        assert_eq!(csc_row_ind, vec![0, 2, 1, 0, 2]);
        assert_eq!(csc_col_ptr, vec![0, 2, 3, 5]);
        assert_eq!(csc_data2, None);
    }

    #[test]
    fn test_csr_to_csc_conversion_dual_layer() {
        // Test matrix:
        // [1 0 2]    [10 0  20]
        // [0 3 0]    [0  30 0 ]
        // [4 0 5]    [40 0  50]
        // CSR format for both layers
        let data = vec![1, 2, 3, 4, 5];
        let data2 = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let col_ind = vec![0, 2, 1, 0, 2];
        let row_ptr = vec![0, 2, 3, 5];
        let (csc_data, csc_col_ptr, csc_row_ind, csc_data2) =
            csr_to_csc(&data, &col_ind, &row_ptr, 3, Some(&data2));

        // Expected CSC format
        assert_eq!(csc_data, vec![1, 4, 3, 2, 5]);
        assert_eq!(csc_row_ind, vec![0, 2, 1, 0, 2]);
        assert_eq!(csc_col_ptr, vec![0, 2, 3, 5]);
        assert_eq!(csc_data2, Some(vec![1.0, 4.0, 3.0, 2.0, 5.0]));
    }

    #[test]
    fn test_streaming_write_read() {
        // Write some cells
        let mut writer = CellGeneSparseWriter::new("test.bin", true, 1000, 2000).unwrap();

        // Write cells one by one (not keeping them in memory)
        for i in 0..100 {
            let cell = CsrCellChunk::from_r_data(&[1, 2, 3], &[0, 5, 10], i, 1e5);
            writer.write_cell_chunk(cell).unwrap();
        }

        writer.finalise().unwrap();

        // Read specific cells
        let mut reader = StreamingSparseReader::new("test.bin").unwrap();

        // Read cells 5, 10, and 50
        let cells = reader.read_cells_by_indices(&[5, 10, 50]).unwrap();
        assert_eq!(cells.len(), 3);
        assert_eq!(cells[0].original_index, 5);
        assert_eq!(cells[1].original_index, 10);
        assert_eq!(cells[2].original_index, 50);

        // Iterate over all cells
        let mut reader = StreamingSparseReader::new("test.bin").unwrap();
        let all_cells: Vec<_> = reader.iter_cells().map(|r| r.unwrap()).collect();
        assert_eq!(all_cells.len(), 100);
    }
}
