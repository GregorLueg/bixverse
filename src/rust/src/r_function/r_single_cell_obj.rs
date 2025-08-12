use extendr_api::prelude::*;

use crate::{helpers::structs_sparse::*, utils::general::flatten_vector};

/// A class for handling single cell count data
/// @export  
#[extendr]
struct SingeCellCountData {
    pub f_path: String,
}

#[extendr]
impl SingeCellCountData {
    pub fn new(f_path: String) -> Self {
        Self { f_path }
    }

    pub fn r_csc_mat_to_file(
        &self,
        no_cells: usize,
        no_genes: usize,
        data: &[i32],
        col_ptr: &[i32],
        row_idx: &[i32],
        target_size: f64,
    ) {
        let ptr_range = 1..col_ptr.len();

        let mut writer = CellGeneSparseWriter::new(&self.f_path, true, no_cells, no_genes).unwrap();

        for i in ptr_range {
            // get the index position
            let end_i = col_ptr[i] as usize;
            let start_i = col_ptr[i - 1] as usize;
            // create chunk from raw data
            let chunk_i = CscCellChunk::from_r_data(
                &data[start_i..end_i],
                &row_idx[start_i..end_i],
                i - 1,
                target_size as f32,
            );

            writer.write_cell_chunk(chunk_i).unwrap();
        }

        writer.finalise().unwrap();
    }

    pub fn file_to_r_csc_mat(&self) -> List {
        let mut reader = StreamingSparseReader::new(&self.f_path).unwrap();

        let no_cells = reader.get_header().total_cells;
        let no_genes = reader.get_header().total_genes;

        let mut data: Vec<Vec<u16>> = Vec::new();
        let mut row_idx: Vec<Vec<u16>> = Vec::new();
        let mut col_ptr: Vec<usize> = Vec::new();

        let all_cells: Vec<_> = reader.iter_cells().map(|r| r.unwrap()).collect();

        let mut current_pol_ptr = 0_usize;
        col_ptr.push(current_pol_ptr);

        for cell in all_cells {
            let data_i = cell.data_raw;
            let len_data_i = data_i.len();
            current_pol_ptr += len_data_i;
            // Add data
            data.push(data_i);
            row_idx.push(cell.row_indices);
            col_ptr.push(current_pol_ptr);
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

    pub fn get_cells_by_indices(&self, indices: &[i32]) -> List {
        // TODO better error handling here
        let mut reader = StreamingSparseReader::new(&self.f_path).unwrap();

        let no_genes = reader.get_header().total_genes;

        let indices = indices
            .iter()
            .map(|x| (*x - 1_i32) as usize)
            .collect::<Vec<usize>>();
        let cells = reader.read_cells_by_indices(&indices).unwrap();

        let mut data: Vec<Vec<u16>> = Vec::new();
        let mut row_idx: Vec<Vec<u16>> = Vec::new();
        let mut col_ptr: Vec<usize> = Vec::new();

        let mut current_pol_ptr = 0_usize;
        col_ptr.push(current_pol_ptr);

        for cell in &cells {
            let data_i = cell.data_raw.clone();
            let len_data_i = data_i.len();
            current_pol_ptr += len_data_i;
            // Add data
            data.push(data_i);
            row_idx.push(cell.row_indices.clone());
            col_ptr.push(current_pol_ptr);
        }

        let data = flatten_vector(data);
        let row_idx = flatten_vector(row_idx);

        let data_final = data.iter().map(|x| *x as i32).collect::<Vec<i32>>();
        let row_idx_final = row_idx.iter().map(|x| *x as i32).collect::<Vec<i32>>();

        list!(
            col_ptr = col_ptr,
            row_idx = row_idx_final,
            data = data_final,
            no_cells = cells.len(),
            no_genes = no_genes
        )
    }
}

extendr_module! {
    mod r_single_cell_obj;
    impl SingeCellCountData;
}
