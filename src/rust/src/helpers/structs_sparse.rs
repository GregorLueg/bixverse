use bincode::{config, decode_from_slice, serde::encode_to_vec, Decode, Encode};
use faer::traits::ComplexField;
use faer::{Mat, MatRef};
use serde::{Deserialize, Serialize};
#[allow(unused_imports)]
use std::fs::{File, OpenOptions};
#[allow(unused_imports)]
use std::io::{BufReader, BufWriter, Read, Seek, SeekFrom, Write};

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
    ///                      row major formant
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

///////////////////////////
// Sparse data streaming //
///////////////////////////

////////////
// Traits //
////////////

#[derive(Encode, Decode, Serialize, Deserialize, Debug)]
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

/////////////////////
// Data structures //
/////////////////////

/// CscCellChunk
///
/// This structure is designed to store the data for fast access in terms of
/// cells for single cell experiments.
///
/// ### Fields
///
/// * `raw_data` - Vector with the raw data (likely `u16`).
/// * `norm_data` - Vector with the normalised data (likely `f16`).
/// * `library_size` - Vector with the library size (likely `u32`).
/// * `to_keep` - Boolean if the cell should be included into anything.
/// * `row_indices` - The row indices of the data.
/// * `col_ptrs` - The column pointers of the data.
/// * `chunk_id` - Id of the chunk.
#[derive(Encode, Decode, Serialize, Deserialize, Debug)]
pub struct CscCellChunk {
    pub raw_data: Vec<u16>,
    pub norm_data: Vec<F16>,
    pub library_size: Vec<u32>,
    pub to_keep: Vec<bool>,
    pub row_indices: Vec<usize>,
    pub col_ptrs: Vec<usize>,
    pub chunk_id: usize,
}

/// CsrGeneChunk
///
/// This structure is designed to store the data for fast access in terms of
/// genes for single cell experiments.
///
/// ### Fields
///
/// * `raw_data` - Vector with the raw data (likely `u16`).
/// * `norm_data` - Vector with the normalised data (likely `f16`).
/// * `average_exp` - Vector with average expression (likely `f16`).
/// * `no_expressed` - Vector with the number of cells expressing this gene.
/// * `to_keep` - Boolean if the gene should be included into anything.
/// * `row_ptrs` - The row pointers of the data.
/// * `col_indices` - The column indices of the data.
/// * `chunk_id` - Id of the chunk.
#[derive(Encode, Decode, Serialize, Deserialize, Debug)]
pub struct CsrGeneChunk {
    pub raw_data: Vec<u16>,
    pub norm_data: Vec<F16>,
    pub average_exp: Vec<F16>,
    pub no_expressed: Vec<usize>,
    pub to_keep: Vec<bool>,
    pub row_ptrs: Vec<usize>,
    pub col_indices: Vec<usize>,
    pub chunk_id: usize,
}

/// CompressedDataHeader
///
/// Stores the information in terms of total cells, total genes, number of
/// chunks in terms of cells and genes and the offset vectors
///
/// ### Params
///
/// * `total_cells` - Total number of cells in the experiment.
/// * `total_genes` - Total number of genes in the experiemnt.
/// * `cell_based` - Boolean indicating that if `true` the data is stored
///                  in a cell-favourable format. If `false`, data is stored
///                  in a gene-favourable format.
/// * `no_chunks` - No of chunks that store either cell or gene data.
/// * `chunk_offsets` - Vector containing the offsets for the cell or gene
///                     chunks.
#[derive(Encode, Decode, Serialize, Deserialize)]
pub struct CompressedDataHeader {
    pub total_cells: usize,
    pub total_genes: usize,
    pub cell_based: bool,
    pub no_chunks: usize,
    pub chunk_offsets: Vec<usize>,
}

//////////////////////
// Streaming writer //
//////////////////////

/// StreamingWriter
///
/// Implementation of a structure for writing in a streamed manner two different
/// types of sparse stored data.
///
/// ### Params
///
/// * `header` - The header of the file.
/// * `writer` - BufWriter to the file.
/// * `current_pos` - The current position of the chunks.
#[allow(dead_code)]
pub struct StreamingSparseWriter {
    header: CompressedDataHeader,
    writer: BufWriter<File>,
    current_pos: usize,
}

#[allow(dead_code)]
impl StreamingSparseWriter {
    /// Generate a new StreamingSparseWriter
    ///
    /// ### Params
    ///
    /// * `path_f` - Path to the file storing the compressed sparse column
    ///                  data for rapid access for cells.
    /// * `cell_based` - Boolean indicating if you are about to store data in
    ///                  cell favourable, i.e., CSC format.
    /// * `total_cells` - Total number of cells.
    /// * `total_genes` - Total number of genes.
    ///
    /// ### Returns
    ///
    /// Returns the ready `StreamingSparseWriter`.
    pub fn new(
        path_f: &str,
        cell_based: bool,
        total_cells: usize,
        total_genes: usize,
    ) -> std::io::Result<Self> {
        // Cells
        let file = File::create(path_f)?;
        let mut writer = BufWriter::new(file);
        let placeholder_cell_header = CompressedDataHeader {
            total_cells,
            total_genes,
            cell_based,
            no_chunks: 0_usize,
            chunk_offsets: Vec::new(),
        };
        let header_size = encode_to_vec(&placeholder_cell_header, config::standard())
            .unwrap()
            .len();
        writer.write_all(&vec![0u8; header_size as usize])?;

        Ok(Self {
            header: placeholder_cell_header,
            writer,
            current_pos: header_size,
        })
    }

    /// Write a Cell Chunk to disk
    ///
    /// * `csc_chunk` - The chunk with the cell data, i.e., `CscCellChunk`
    pub fn write_csc_chunk(&mut self, csc_chunk: CscCellChunk) -> std::io::Result<()> {
        assert!(
            self.header.cell_based,
            "This is not set for a cell-based (CSC format) writing!"
        );

        self.header.chunk_offsets.push(self.current_pos);

        let encoded = encode_to_vec(&csc_chunk, config::standard()).unwrap();
        let chunk_size = encoded.len() as u64;

        self.writer.write_all(&chunk_size.to_le_bytes())?;
        self.writer.write_all(&encoded)?;

        self.current_pos += 8 + chunk_size as usize;
        self.header.no_chunks += 1;

        Ok(())
    }

    /// Write a Gene Chunk to disk
    ///
    /// * `csr_chunk` - The chunk with the cell data, i.e., `CsrGeneChunk`
    pub fn write_csr_chunk(&mut self, csr_chunk: CsrGeneChunk) -> std::io::Result<()> {
        assert!(
            !self.header.cell_based,
            "This is not set for a gene-based (CSR format) writing!"
        );
        self.header.chunk_offsets.push(self.current_pos);

        let encoded = encode_to_vec(&csr_chunk, config::standard()).unwrap();
        let chunk_size = encoded.len();

        self.writer.write_all(&chunk_size.to_le_bytes())?;
        self.writer.write_all(&encoded)?;

        self.current_pos += 8 + chunk_size as usize;
        self.header.no_chunks += 1;

        Ok(())
    }

    /// Finalise writing
    pub fn finalise(mut self) -> std::io::Result<()> {
        self.writer.seek(SeekFrom::Start(0))?;
        let header_data = encode_to_vec(&self.header, config::standard()).unwrap();
        self.writer.write_all(&header_data)?;
        self.writer.flush()?;

        Ok(())
    }
}

//////////////////////
// Streaming reader //
//////////////////////

/// StreamingSparseReader
///
/// Implementation of a structure for writing in a streamed manner
///
/// ### Params
///
/// * `header` - The header of the file.
/// * `writer` - BufReader of the file.
/// * `current_pos` - The current position of the chunks.
#[allow(dead_code)]
pub struct StreamingSparseReader {
    header: CompressedDataHeader,
    reader: BufReader<File>,
    current_chunk: usize,
}

#[allow(dead_code)]
impl StreamingSparseReader {
    /// Initialise the `StreamingSparseReader`
    ///
    /// ### Params
    ///
    /// * `f_path` - Path to the file.
    pub fn new(f_path: &str) -> std::io::Result<Self> {
        let file = File::open(f_path)?;
        let mut reader = BufReader::new(file);

        // Read in the header
        let mut header_size_buf = [0u8; 8];
        reader.read_exact(&mut header_size_buf)?;

        reader.seek(SeekFrom::Start(0))?;

        let mut header_buf = vec![
            0u8;
            encode_to_vec(
                &CompressedDataHeader {
                    total_cells: 0,
                    total_genes: 0,
                    cell_based: true,
                    no_chunks: 0,
                    chunk_offsets: Vec::new(),
                },
                config::standard()
            )
            .unwrap()
            .len() as usize
        ];

        reader.read_exact(&mut header_buf)?;

        let (header, _): (CompressedDataHeader, usize) =
            decode_from_slice(&header_buf, config::standard()).unwrap();

        Ok(Self {
            reader,
            header,
            current_chunk: 0,
        })
    }

    /// Read in the cell chunk
    ///
    /// The function will panic if the underlying file does not seem to contain
    /// `CscCellChunk`s.
    ///
    /// ### Returns
    ///
    /// `Result<Option<CscCellChunk>>`
    pub fn read_cell_chunk(&mut self) -> std::io::Result<Option<CscCellChunk>> {
        assert!(
            self.header.cell_based,
            "The data you are trying to read in is not set up for CscCellChunk."
        );

        if self.current_chunk >= self.header.no_chunks {
            return Ok(None);
        }

        // get to the chunk position
        let chunk_offset = self.header.chunk_offsets[self.current_chunk] as u64;
        self.reader.seek(SeekFrom::Start(chunk_offset))?;

        // read chunk size
        let mut size_buf = [0u8; 8];
        self.reader.read_exact(&mut size_buf)?;
        let chunk_size = u64::from_le_bytes(size_buf) as usize;

        // read chunk data
        let mut chunk_buf = vec![0u8; chunk_size];
        self.reader.read_exact(&mut chunk_buf)?;

        let (chunk, _): (CscCellChunk, usize) =
            decode_from_slice(&chunk_buf, config::standard()).unwrap();

        self.current_chunk += 1;

        Ok(Some(chunk))
    }

    /// Read in the gene chunk
    ///
    /// The function will panic if the underlying file does not seem to contain
    /// `CsrGeneChunk`s.
    ///
    /// ### Returns
    ///
    /// `Result<Option<CsrGeneChunk>>`
    pub fn read_gene_chunk(&mut self) -> std::io::Result<Option<CsrGeneChunk>> {
        assert!(
            !self.header.cell_based,
            "The data you are trying to read in is not set up for CscCellChunk."
        );

        if self.current_chunk >= self.header.no_chunks {
            return Ok(None);
        }

        // get to the chunk position
        let chunk_offset = self.header.chunk_offsets[self.current_chunk] as u64;
        self.reader.seek(SeekFrom::Start(chunk_offset))?;

        // read chunk size
        let mut size_buf = [0u8; 8];
        self.reader.read_exact(&mut size_buf)?;
        let chunk_size = u64::from_le_bytes(size_buf) as usize;

        // read chunk data
        let mut chunk_buf = vec![0u8; chunk_size];
        self.reader.read_exact(&mut chunk_buf)?;

        let (chunk, _): (CsrGeneChunk, usize) =
            decode_from_slice(&chunk_buf, config::standard()).unwrap();

        self.current_chunk += 1;

        Ok(Some(chunk))
    }

    /// Get the header
    ///
    /// ### Returns
    ///
    /// The `CompressedDataHeader` with information about the file.
    pub fn get_header(&self) -> &CompressedDataHeader {
        &self.header
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
}
