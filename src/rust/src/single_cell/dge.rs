use crate::core::data::sparse_io::*;

pub fn calculate_dge_grps(f_path: &str, grp_1_indices: &[usize], grp_2_indices: &[usize]) {
    let reader = ParallelSparseReader::new(f_path).unwrap();

    let cell_chunks_1 = reader.read_cells_parallel(grp_1_indices);
    let cell_chunks_2 = reader.read_cells_parallel(grp_2_indices);
}
