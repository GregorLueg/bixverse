use bincode::{config, decode_from_slice, serde::encode_to_vec, Decode, Encode};
use half::f16;
use rustc_hash::FxHashMap;
use serde::{Deserialize, Serialize};
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Seek, SeekFrom, Write};

use crate::utils::traits::F16;

///////////////////////////
// Sparse data streaming //
///////////////////////////

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
/// * `data_norm` - Array of the normalised counts of this cell. This will do a
///   CPM type transformation and then calculate the ln_1p.
/// * `library_size` - Total library size/UMI counts of the cell.
/// * `col_indices` - The col indices of the genes.
/// * `original_index` - Original (row) index of the cell.
/// * `to_keep` - Flat if the cell should be included in certain analysis.
///   Future feature.
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
    ///
    /// Assumes columns = genes, rows = cells.
    ///
    /// ### Params
    ///
    /// * `data` - The raw counts present in this cell
    /// * `col_idx` - The column indices where the gene is expressed.
    /// * `original_index` - Original row index in the matrix
    /// * `size_factor` - To which size to normalise to. 1e6 -> CPM normalisation.
    ///
    /// ### Returns
    ///
    /// The `CsrCellChunk` for this cell.
    pub fn from_r_data(
        data: &[i32],
        col_idx: &[i32],
        original_index: usize,
        size_factor: f32,
    ) -> Self {
        let data_f32 = data.iter().map(|x| *x as f32).collect::<Vec<f32>>();
        let sum = data_f32.iter().sum::<f32>();
        let data_norm: Vec<F16> = data_f32
            .into_iter()
            .map(|x| {
                let norm = (x / sum * size_factor).ln_1p();
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
///   Future feature.
#[derive(Encode, Decode, Serialize, Deserialize, Debug)]
pub struct CscGeneChunk {
    pub data_raw: Vec<u16>,
    pub data_norm: Vec<F16>,
    pub avg_exp: F16,
    pub nnz: usize,
    // u32 as there might be clearly more than 65_535 cells in the data
    // 4_294_967_295 should be enough however...
    pub row_indices: Vec<u32>,
    pub original_index: usize,
    pub to_keep: bool,
}

impl CscGeneChunk {
    /// Helper function to generate the CscGeneChunk from converted data
    ///
    /// ### Params
    ///
    /// * `data_raw` - The raw counts for this gene
    /// * `data_norm` - The normalised counts for this gene.
    /// * `col_idx` - The column indices for which cells this gene is expressed.
    /// * `original_index` - Original row index
    ///
    /// ### Returns
    ///
    /// The `CscGeneChunk` for this gene.
    pub fn from_conversion(
        data_raw: &[u16],
        data_norm: &[F16],
        col_idx: &[usize],
        original_index: usize,
    ) -> Self {
        let avg_exp = data_norm.iter().sum::<F16>();
        let nnz = data_raw.len();

        Self {
            data_raw: data_raw.to_vec(),
            data_norm: data_norm.to_vec(),
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
/// * `index_map` - FxHashMap with the original index -> chunk info
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
    ///
    /// ### Returns
    ///
    /// The `CellGeneSparseWriter`.
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

impl Iterator for CellChunkIterator<'_> {
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

impl Iterator for GeneChunkIterator<'_> {
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
