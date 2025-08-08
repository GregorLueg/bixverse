use extendr_api::prelude::*;

use crate::{helpers::structs_sparse::*, utils::general::flatten_vector};

/// @export
#[extendr]
pub fn rs_csc_to_binary_f(
    f_path: &str,
    no_cells: usize,
    no_genes: usize,
    data: &[i32],
    col_ptr: &[i32],
    row_idx: &[i32],
) {
    let ptr_range = 1..col_ptr.len();

    let mut writer = CellGeneSparseWriter::new(f_path, true, no_cells, no_genes).unwrap();

    for i in ptr_range {
        // get the index position
        let end_i = col_ptr[i] as usize;
        let start_i = col_ptr[i - 1] as usize;
        // create chunk from raw data
        let chunk_i =
            CscCellChunk::from_r_data(&data[start_i..end_i], &row_idx[start_i..end_i], end_i);

        writer.write_cell_chunk(chunk_i).unwrap();
    }

    writer.finalise().unwrap();
}

pub fn rs_binary_f_to_csc(f_path: &str) -> List {
    let mut reader = StreamingSparseReader::new(f_path).unwrap();

    let total_chunks = reader.get_header().no_chunks;
    let no_cells = reader.get_header().total_cells;
    let no_genes = reader.get_header().total_genes;

    let mut data: Vec<Vec<u16>> = Vec::new();
    let mut row_idx: Vec<Vec<u16>> = Vec::new();
    let mut col_ptr: Vec<usize> = Vec::new();

    // Add zero
    col_ptr.push(0);

    for i in 0..total_chunks {
        let chunk_i = reader.read_cell_chunk().unwrap().unwrap();
        data.push(chunk_i.data_raw);
        row_idx.push(chunk_i.row_indices);
        col_ptr.push(chunk_i.original_index);
    }

    let data = flatten_vector(data);
    let row_idx = flatten_vector(row_idx);

    let data_final = data.iter().map(|x| *x as i32).collect::<Vec<i32>>();
    let row_idx_final = row_idx.iter().map(|x| *x as i32).collect::<Vec<i32>>();

    list!(
        col_ptr = col_ptr,
        row_idx = row_idx_final,
        data = data_final,
        no_cells = no_cells,
        no_genes = no_genes
    )
}

extendr_module! {
    mod r_single_cell_obj;
    fn rs_csc_to_binary_f;
}
