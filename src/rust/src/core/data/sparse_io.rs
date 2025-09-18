use bincode::{config, decode_from_slice, serde::encode_to_vec, Decode, Encode};
use half::f16;
use memmap2::MmapOptions;
use rayon::iter::*;
use rustc_hash::{FxHashMap, FxHashSet};
use serde::{Deserialize, Serialize};
use std::fs::File;
use std::io::{BufWriter, Seek, SeekFrom, Write};
use std::marker::Sync;
use std::path::Path;
use std::sync::Arc;

use crate::core::data::sparse_structures::*;
use crate::single_cell::processing::*;
use crate::utils::traits::*;

///////////////////////////
// Sparse data streaming //
///////////////////////////

/////////////////////
// Data structures //
/////////////////////

/// CellOnFileQuality
///
/// This structure is being generate after a first scan of the file on disk and
/// defining which cells and genes to actually read in.
///
/// ### Fields
///
/// * `cells_to_keep` - Vector of indices of the cells to keep.
/// * `genes_to_keep` - Vector of indices of the genes to keep.
/// * `cells_to_keep_set` - HashSet of the indices to keep.
/// * `genes_to_keep_set` - HashSet of the genes to keep.
/// * `cell_old_to_new` - Mapping of the old indices to the new indices for the
///   cells.
/// * `gene_old_new` - Mapping of old indices to new indices for the genes.
#[derive(Debug, Clone)]
pub struct CellOnFileQuality {
    pub cells_to_keep: Vec<usize>,
    pub genes_to_keep: Vec<usize>,
    pub cells_to_keep_set: FxHashSet<usize>,
    pub genes_to_keep_set: FxHashSet<usize>,
    pub cell_old_to_new: FxHashMap<usize, usize>,
    pub gene_old_to_new: FxHashMap<usize, usize>,
}

impl CellOnFileQuality {
    /// Create a new CellOnFileQuality structure
    ///
    /// ### Params
    ///
    /// * cells_to_keep - Index positions of the cells to keep
    /// * genes_to_keep - Index positions of the genes to keep
    ///
    /// ### Returns
    ///
    /// Initiliased self
    pub fn new(cells_to_keep: Vec<usize>, genes_to_keep: Vec<usize>) -> Self {
        Self {
            cells_to_keep,
            genes_to_keep,
            cells_to_keep_set: FxHashSet::default(),
            genes_to_keep_set: FxHashSet::default(),
            cell_old_to_new: FxHashMap::default(),
            gene_old_to_new: FxHashMap::default(),
        }
    }

    /// Generate internally the sets and maps
    pub fn generate_maps_sets(&mut self) {
        let cells_to_keep_set: FxHashSet<usize> = self.cells_to_keep.iter().cloned().collect();
        let genes_to_keep_set: FxHashSet<usize> = self.genes_to_keep.iter().cloned().collect();
        let cell_old_to_new: FxHashMap<usize, usize> = self
            .cells_to_keep
            .iter()
            .enumerate()
            .map(|(new_idx, &old_idx)| (old_idx, new_idx))
            .collect();
        let gene_old_to_new: FxHashMap<usize, usize> = self
            .genes_to_keep
            .iter()
            .enumerate()
            .map(|(new_idx, &old_idx)| (old_idx, new_idx))
            .collect();

        self.cells_to_keep_set = cells_to_keep_set;
        self.genes_to_keep_set = genes_to_keep_set;
        self.cell_old_to_new = cell_old_to_new;
        self.gene_old_to_new = gene_old_to_new
    }
}

/// CsrCellChunk
///
/// This structure is designed to store the data of a single cell in a
/// CSR-like format optimised for rapid access on disk.
///
/// ### Fields
///
/// * `data_raw` - Array of the raw counts of this cell.
/// * `data_norm` - Array of the normalised counts of this cell. This will do a
///   CPM-type transformation and then calculate the ln_1p.
/// * `library_size` - Total library size/UMI counts of the cell.
/// * `col_indices` - The col indices of the genes.
/// * `original_index` - Original (row) index of the cell.
/// * `to_keep` - Flat if the cell should be included in certain analysis.
///   Future feature.
#[derive(Debug)]
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
    pub fn from_data<T, U>(
        data: &[T],
        col_idx: &[U],
        original_index: usize,
        size_factor: f32,
        to_keep: bool,
    ) -> Self
    where
        T: ToF32AndU16,
        U: ToF32AndU16,
    {
        let data_f32 = data.iter().map(|&x| x.to_f32()).collect::<Vec<f32>>();
        let sum = data_f32.iter().sum::<f32>();
        let data_norm: Vec<F16> = data_f32
            .into_iter()
            .map(|x| {
                let norm = (x / sum * size_factor).ln_1p();
                F16::from(f16::from_f32(norm))
            })
            .collect();
        Self {
            data_raw: data.iter().map(|&x| x.to_u16()).collect::<Vec<u16>>(),
            data_norm,
            library_size: sum as usize,
            col_indices: col_idx.iter().map(|&x| x.to_u16()).collect::<Vec<u16>>(),
            original_index,
            to_keep,
        }
    }

    /// Write directly to bytes on disk
    ///
    /// ### Params
    ///
    /// * `writer` - Something that can write, i.e., has the implementation
    ///   `Write`
    pub fn write_to_bytes(&self, writer: &mut impl Write) -> std::io::Result<()> {
        writer.write_all(&(self.data_raw.len() as u32).to_le_bytes())?;
        writer.write_all(&(self.data_norm.len() as u32).to_le_bytes())?;
        writer.write_all(&(self.col_indices.len() as u32).to_le_bytes())?;
        writer.write_all(&(self.library_size as u64).to_le_bytes())?;
        writer.write_all(&(self.original_index as u64).to_le_bytes())?;
        // include some padding
        writer.write_all(&[self.to_keep as u8, 0, 0, 0])?;

        // unsafe fun with direct write to disk
        let data_raw_bytes = unsafe {
            std::slice::from_raw_parts(self.data_raw.as_ptr() as *const u8, self.data_raw.len() * 2)
        };
        writer.write_all(data_raw_bytes)?;

        let data_norm_bytes = unsafe {
            std::slice::from_raw_parts(
                self.data_norm.as_ptr() as *const u8,
                self.data_norm.len() * 2,
            )
        };
        writer.write_all(data_norm_bytes)?;

        let col_indices_bytes = unsafe {
            std::slice::from_raw_parts(
                self.col_indices.as_ptr() as *const u8,
                self.col_indices.len() * 2,
            )
        };
        writer.write_all(col_indices_bytes)?;

        Ok(())
    }

    /// Read data from buffer
    ///
    /// ### Params
    ///
    /// * `buffer` - A slice of u8's representing the buffer
    ///
    /// ### Return
    ///
    /// The `CsrCellChunk`
    pub fn read_from_buffer(buffer: &[u8]) -> std::io::Result<Self> {
        if buffer.len() < 32 {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "Buffer too small for header",
            ));
        }

        // parse header
        let header = &buffer[0..32];

        let data_raw_len =
            u32::from_le_bytes([header[0], header[1], header[2], header[3]]) as usize;
        let data_norm_len =
            u32::from_le_bytes([header[4], header[5], header[6], header[7]]) as usize;
        let col_indices_len =
            u32::from_le_bytes([header[8], header[9], header[10], header[11]]) as usize;
        let library_size = u64::from_le_bytes([
            header[12], header[13], header[14], header[15], header[16], header[17], header[18],
            header[19],
        ]) as usize;
        let original_index = u64::from_le_bytes([
            header[20], header[21], header[22], header[23], header[24], header[25], header[26],
            header[27],
        ]) as usize;
        let to_keep = header[28] != 0;

        let data_start = 32;
        let data_end = data_start + data_raw_len * 2;
        let norm_end = data_end + data_norm_len * 2;

        // direct transmutation for raw counts - no intermediate allocation
        let data_raw = unsafe {
            let ptr = buffer.as_ptr().add(data_start) as *const u16;
            std::slice::from_raw_parts(ptr, data_raw_len).to_vec()
        };

        // For F16, I need to convert
        let data_norm = unsafe {
            let ptr = buffer.as_ptr().add(data_end) as *const u16;
            let slice = std::slice::from_raw_parts(ptr, data_norm_len);
            let mut norm = Vec::with_capacity(data_norm_len);
            norm.extend(slice.iter().map(|&bits| F16::from_bits(bits)));
            norm
        };

        let col_indices = unsafe {
            let ptr = buffer.as_ptr().add(norm_end) as *const u16;
            std::slice::from_raw_parts(ptr, col_indices_len).to_vec()
        };

        Ok(Self {
            data_raw,
            data_norm,
            library_size,
            col_indices,
            original_index,
            to_keep,
        })
    }

    /// Generate a vector of Chunks from CompressedSparseData
    ///
    /// ### Params
    ///
    /// * `sparse_data` - The `CompressedSparseData` (in CSR format!)
    /// * `min_genes` - Number of genes per cell to be included
    /// * `size_factor` - Size factor for normalisation. 1e6 -> CPM
    ///
    /// ### Returns
    ///
    /// A tuple of the `Vec<CsrCellChunk>` and if the cell should be kept.
    pub fn generate_chunks_sparse_data<T, U>(
        sparse_data: CompressedSparseData<T, U>,
        cell_qc: MinCellQuality,
    ) -> (Vec<CsrCellChunk>, CellQuality)
    where
        T: Clone + Default + Into<u32> + Sync,
        U: Clone + Default,
    {
        let n_cells = sparse_data.indptr.len() - 1;

        let (nnz, to_keep) = filter_by_nnz(&sparse_data.indptr, cell_qc.min_unique_genes);

        let (res, lib_size): (Vec<CsrCellChunk>, Vec<usize>) = (0..n_cells)
            .into_par_iter()
            .map(|i| {
                let start_i = sparse_data.indptr[i];
                let end_i = sparse_data.indptr[i + 1];
                let indices_i = &sparse_data.indices[start_i..end_i];
                let data_i: &Vec<u32> = &sparse_data.data[start_i..end_i]
                    .iter()
                    .map(|i| i.clone().into())
                    .collect();
                let sum_data_i = data_i.iter().sum::<u32>() as usize;
                let to_keep_i = to_keep[i] & (sum_data_i >= cell_qc.min_lib_size);

                let chunk_i =
                    CsrCellChunk::from_data(data_i, indices_i, i, cell_qc.target_size, to_keep_i);

                (chunk_i, sum_data_i)
            })
            .unzip();

        let qc_data = CellQuality {
            to_keep,
            lib_size: Some(lib_size),
            no_genes: Some(nnz),
        };

        (res, qc_data)
    }

    /// Helper function to update the column index
    ///
    /// ### Params
    ///
    /// * `new_index` - Which is the new column index to set to
    pub fn update_index(&mut self, new_index: &usize) {
        self.original_index = *new_index;
    }

    /// Helper function to check if chunk passes QC
    ///
    /// ### Params
    ///
    /// * `cell_qc` - Structure containiner the required library size and min
    ///   number of genes
    pub fn passes_qc(&self, cell_qc: &MinCellQuality) -> bool {
        self.col_indices.len() >= cell_qc.min_unique_genes
            && self.library_size >= cell_qc.min_lib_size
    }

    /// Helper function to get QC parameters for this cell
    ///
    /// ### Reutrns
    ///
    /// A tuple of `(no_genes, library_size)`
    pub fn get_qc_info(&self) -> (usize, usize) {
        (self.col_indices.len(), self.library_size)
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
    /// * `to_keep` - Shall the gene be included in later analysis.
    ///
    /// ### Returns
    ///
    /// The `CscGeneChunk` for this gene.
    pub fn from_conversion(
        data_raw: &[u16],
        data_norm: &[F16],
        col_idx: &[usize],
        original_index: usize,
        to_keep: bool,
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
            to_keep,
        }
    }

    /// Write directly to bytes on disk
    ///
    /// ### Params
    ///
    /// * `writer` - Something that can write, i.e., has the implementation
    ///   `Write`
    pub fn write_to_bytes(&self, writer: &mut impl Write) -> std::io::Result<()> {
        writer.write_all(&(self.data_raw.len() as u32).to_le_bytes())?;
        writer.write_all(&(self.data_norm.len() as u32).to_le_bytes())?;
        writer.write_all(&(self.row_indices.len() as u32).to_le_bytes())?;
        writer.write_all(&self.avg_exp.to_le_bytes())?;
        // first padding
        writer.write_all(&[0, 0])?;
        writer.write_all(&(self.nnz as u64).to_le_bytes())?;
        writer.write_all(&(self.original_index as u64).to_le_bytes())?;
        // bit padding
        writer.write_all(&[self.to_keep as u8, 0, 0, 0])?;

        let data_raw_bytes = unsafe {
            std::slice::from_raw_parts(self.data_raw.as_ptr() as *const u8, self.data_raw.len() * 2)
        };
        writer.write_all(data_raw_bytes)?;

        let data_norm_bytes = unsafe {
            std::slice::from_raw_parts(
                self.data_norm.as_ptr() as *const u8,
                self.data_norm.len() * 2,
            )
        };
        writer.write_all(data_norm_bytes)?;

        let row_indices_bytes = unsafe {
            std::slice::from_raw_parts(
                self.row_indices.as_ptr() as *const u8,
                self.row_indices.len() * 4,
            )
        };
        writer.write_all(row_indices_bytes)?;

        Ok(())
    }

    /// Read data from buffer
    ///
    /// ### Params
    ///
    /// * `buffer` - A slice of u8's representing the buffer
    ///
    /// ### Return
    ///
    /// The `CscGeneChunk`
    pub fn read_from_buffer(buffer: &[u8]) -> std::io::Result<Self> {
        if buffer.len() < 36 {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "Buffer too small for header",
            ));
        }

        // parse header
        let header = &buffer[0..36];
        let data_raw_len =
            u32::from_le_bytes([header[0], header[1], header[2], header[3]]) as usize;
        let data_norm_len =
            u32::from_le_bytes([header[4], header[5], header[6], header[7]]) as usize;
        let row_indices_len =
            u32::from_le_bytes([header[8], header[9], header[10], header[11]]) as usize;
        let avg_exp = F16::from_le_bytes([header[12], header[13]]);
        // skip 2 bytes padding
        let nnz = u64::from_le_bytes([
            header[16], header[17], header[18], header[19], header[20], header[21], header[22],
            header[23],
        ]) as usize;
        let original_index = u64::from_le_bytes([
            header[24], header[25], header[26], header[27], header[28], header[29], header[30],
            header[31],
        ]) as usize;
        let to_keep = header[32] != 0;

        let data_start = 36;
        let data_end = data_start + data_raw_len * 2;
        let norm_end = data_end + data_norm_len * 2;

        // direct transmutation for raw counts - no intermediate allocation
        let data_raw = unsafe {
            let ptr = buffer.as_ptr().add(data_start) as *const u16;
            std::slice::from_raw_parts(ptr, data_raw_len).to_vec()
        };

        // For F16, I need to convert
        let data_norm = unsafe {
            let ptr = buffer.as_ptr().add(data_end) as *const u16;
            let slice = std::slice::from_raw_parts(ptr, data_norm_len);
            let mut norm = Vec::with_capacity(data_norm_len);
            norm.extend(slice.iter().map(|&bits| F16::from_bits(bits)));
            norm
        };

        let row_indices = unsafe {
            let ptr = buffer.as_ptr().add(norm_end) as *const u32;
            std::slice::from_raw_parts(ptr, row_indices_len).to_vec()
        };

        Ok(Self {
            data_raw,
            data_norm,
            avg_exp,
            nnz,
            row_indices,
            original_index,
            to_keep,
        })
    }

    /// Helper function that allows to filter out cells
    ///
    /// ### Params
    ///
    /// * `cells_to_keep` - HashSet with cell index positions to keep.
    pub fn filter_selected_cells(&mut self, cells_to_keep: &FxHashSet<u32>) {
        let mut indices_to_keep = Vec::new();

        for (i, &cell_index) in self.row_indices.iter().enumerate() {
            if cells_to_keep.contains(&cell_index) {
                indices_to_keep.push(i);
            }
        }

        // Update internal values
        self.row_indices = indices_to_keep
            .iter()
            .map(|&i| self.row_indices[i])
            .collect();

        self.data_raw = indices_to_keep.iter().map(|&i| self.data_raw[i]).collect();

        self.data_norm = indices_to_keep.iter().map(|&i| self.data_norm[i]).collect();

        self.nnz = self.row_indices.len();
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
#[derive(Encode, Decode, Serialize, Deserialize, Clone)]
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
    chunks_since_flush: usize,
    flush_frequency: usize,
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
    pub fn new<P: AsRef<Path>>(
        path_f: P,
        cell_based: bool,
        total_cells: usize,
        total_genes: usize,
    ) -> std::io::Result<Self> {
        let file = File::create(path_f)?;
        let mut writer = BufWriter::with_capacity(10 * 1024 * 1024, file);

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

        let flush_frequency = if cell_based { 100000_usize } else { 1000_usize };

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
            chunks_since_flush: 0_usize,
            flush_frequency,
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

        let current_pos = self.writer.stream_position()?;
        let chunk_offset = current_pos - self.chunks_start_pos;
        self.header.chunk_offsets.push(chunk_offset);

        self.header
            .index_map
            .insert(cell_chunk.original_index, self.header.no_chunks);

        // Calculate size first (for the size prefix)
        let size = 32 + // header size
               (cell_chunk.data_raw.len() * 2) +
               (cell_chunk.data_norm.len() * 2) +
               (cell_chunk.col_indices.len() * 2);

        self.writer.write_all(&(size as u64).to_le_bytes())?;

        // Use our custom serialization
        cell_chunk.write_to_bytes(&mut self.writer)?;

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

        let current_pos = self.writer.stream_position()?;
        let chunk_offset = current_pos - self.chunks_start_pos;
        self.header.chunk_offsets.push(chunk_offset);
        self.header
            .index_map
            .insert(gene_chunk.original_index, self.header.no_chunks);

        // Calculate size first (for the size prefix)
        let size = 36 + // header size
        (gene_chunk.data_raw.len() * 2) +
        (gene_chunk.data_norm.len() * 2) +
        (gene_chunk.row_indices.len() * 4);

        self.writer.write_all(&(size as u64).to_le_bytes())?;

        // Use our custom serialization
        gene_chunk.write_to_bytes(&mut self.writer)?;

        self.header.no_chunks += 1;
        self.chunks_since_flush += 1;

        if self.chunks_since_flush >= self.flush_frequency {
            self.writer.flush()?;
            self.chunks_since_flush = 0;
        }

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

    /// Update the number of cells in the header
    ///
    /// Helper function to update the number of cells to the number that was
    /// finally written on disk.
    ///
    /// ### Params
    ///
    /// * `no_cells` - New number of cells
    pub fn update_header_no_cells(&mut self, no_cells: usize) {
        self.header.total_cells = no_cells;
    }

    /// Update the number of cells in the header
    ///
    /// Helper function to update the number of genes to the number that was
    /// finally written on disk.
    ///
    /// ### Params
    ///
    /// * `no_genes` - New number of genes
    pub fn update_header_no_genes(&mut self, no_genes: usize) {
        self.header.total_genes = no_genes;
    }
}

//////////////////////
// Streaming reader //
//////////////////////

pub struct ParallelSparseReader {
    header: SparseDataHeader,
    mmap: Arc<memmap2::Mmap>,
    chunks_start: u64,
}

#[allow(dead_code)]
impl ParallelSparseReader {
    /// Generate a new parallelised streaming reader
    ///
    /// ### Params
    ///
    /// * `f_path` - Path to the file
    ///
    /// ### Returns
    ///
    /// Initialised `ParallelSparseReader`.
    pub fn new(f_path: &str) -> std::io::Result<Self> {
        let file = File::open(f_path)?;
        let mmap = unsafe { MmapOptions::new().map(&file)? };

        // Parse headers from mmap
        let file_header_bytes = &mmap[0..64];
        let (file_header, _) = decode_from_slice::<FileHeader, _>(
            file_header_bytes,
            config::standard(),
        )
        .map_err(|_| {
            std::io::Error::new(std::io::ErrorKind::InvalidData, "File header decode failed")
        })?;

        // Read main header
        let main_header_offset = file_header.main_header_offset as usize;
        let header_size = u64::from_le_bytes(
            mmap[main_header_offset..main_header_offset + 8]
                .try_into()
                .unwrap(),
        ) as usize;

        let header_bytes = &mmap[main_header_offset + 8..main_header_offset + 8 + header_size];
        let (header, _) =
            decode_from_slice::<SparseDataHeader, _>(header_bytes, config::standard()).map_err(
                |_| std::io::Error::new(std::io::ErrorKind::InvalidData, "Header decode failed"),
            )?;

        Ok(Self {
            header,
            mmap: Arc::new(mmap),
            chunks_start: 64,
        })
    }

    /// Read in cells by indices in a multithreaded manner
    ///
    /// ### Params
    ///
    /// * `indices` - Slice of index positions of the cells to retrieve
    ///
    /// ### Returns
    ///
    /// Returns an array of `CsrCellChunk`.
    pub fn read_cells_parallel(&self, indices: &[usize]) -> Vec<CsrCellChunk> {
        assert!(
            self.header.cell_based,
            "The file is not set up for CellChunks."
        );

        indices
            .par_iter()
            .map(|&original_index| {
                let chunk_index = *self.header.index_map.get(&original_index).unwrap();
                let chunk_offset =
                    (self.chunks_start + self.header.chunk_offsets[chunk_index]) as usize;

                // Read size
                let size = u64::from_le_bytes(
                    self.mmap[chunk_offset..chunk_offset + 8]
                        .try_into()
                        .unwrap(),
                ) as usize;

                // Parse chunk
                let buffer = &self.mmap[chunk_offset + 8..chunk_offset + 8 + size];
                CsrCellChunk::read_from_buffer(buffer).unwrap()
            })
            .collect()
    }

    /// Read in genes by indices in a multithreaded manner
    ///
    /// ### Params
    ///
    /// * `indices` - Slice of index positions of the genes to retrieve
    ///
    /// ### Returns
    ///
    /// Returns an array of `CscGeneChunk`.
    pub fn read_gene_parallel(&self, indices: &[usize]) -> Vec<CscGeneChunk> {
        assert!(
            !self.header.cell_based,
            "The file is not set up for CellChunks."
        );

        indices
            .par_iter()
            .map(|&original_index| {
                let chunk_index = *self.header.index_map.get(&original_index).unwrap();
                let chunk_offset =
                    (self.chunks_start + self.header.chunk_offsets[chunk_index]) as usize;

                // Read size
                let size = u64::from_le_bytes(
                    self.mmap[chunk_offset..chunk_offset + 8]
                        .try_into()
                        .unwrap(),
                ) as usize;

                // Parse chunk
                let buffer = &self.mmap[chunk_offset + 8..chunk_offset + 8 + size];
                CscGeneChunk::read_from_buffer(buffer).unwrap()
            })
            .collect()
    }

    /// Return all cells
    ///
    /// ### Returns
    ///
    /// Returns an array of `CsrCellChunk` containing all cells on disk.
    pub fn get_all_cells(&self) -> Vec<CsrCellChunk> {
        let iter: Vec<usize> = (0..self.header.total_cells).collect();

        self.read_cells_parallel(&iter)
    }

    /// Return all genes
    ///
    /// ### Returns
    ///
    /// Returns an array of `CscGeneChunk` containing all genes on disk.
    pub fn get_all_genes(&self) -> Vec<CscGeneChunk> {
        let iter: Vec<usize> = (0..self.header.total_genes).collect();

        self.read_gene_parallel(&iter)
    }

    /// Helper to return the header
    ///
    /// ### Returns
    ///
    /// Returns the header file
    pub fn get_header(&self) -> SparseDataHeader {
        self.header.clone()
    }
}
