use extendr_api::{list, List, Robj};
use faer::Mat;
use rand::prelude::IndexedRandom;
use rand::{Rng, SeedableRng};
use rayon::prelude::*;
use rustc_hash::{FxHashMap, FxHashSet};
use std::time::Instant;

use crate::core::data::sparse_io::*;
use crate::core::data::sparse_structures::*;
use crate::core::graph::community_detection::*;
use crate::core::graph::knn::*;
use crate::single_cell::sc_knn_snn::KnnParams;

/////////////
// Helpers //
/////////////

/// Helper function to aggregate the meta cells
///
/// The function will generate the metacells based on the provided indices.
/// Per meta-cell it will aggregate the raw counts and recalculate the norm
/// counts based on the aggregated counts.
///
/// ### Params
///
/// * `reader` - The reader structure to get the cells from disk
/// * `meta_cells` - The indices of the meta cells
/// * `target_size` - Float defining the target size for the normalisation
///   procedure. Usually defaults to `1e4` in single cell.
/// * `n_genes` - Total number of genes in the data
///
/// ### Return
///
/// `CompressedSparseData` in CSR format with aggregated raw counts and re-
/// normalised counts per meta cell.
pub fn aggregate_meta_cells(
    reader: &ParallelSparseReader,
    metacells: &[&[usize]],
    target_size: f32,
    n_genes: usize,
) -> CompressedSparseData<u32, f32> {
    let n_metacells = metacells.len();
    let mut all_data: Vec<u32> = Vec::new();
    let mut all_data_norm: Vec<f32> = Vec::new();
    let mut all_indices: Vec<usize> = Vec::new();
    let mut all_indptr: Vec<usize> = vec![0];

    const CHUNK_SIZE: usize = 1000;

    for chunk_start in (0..n_metacells).step_by(CHUNK_SIZE) {
        let chunk_end = (chunk_start + CHUNK_SIZE).min(n_metacells);
        let chunk = &metacells[chunk_start..chunk_end];

        let results: Vec<(Vec<usize>, Vec<u32>, Vec<f32>)> = chunk
            .par_iter()
            .map(|cell_idx| {
                let cells = reader.read_cells_parallel(cell_idx);
                let mut gene_counts: FxHashMap<usize, u32> = FxHashMap::default();
                let mut library_size: u32 = 0;

                for cell in &cells {
                    for (idx, &count) in cell.indices.iter().zip(cell.data_raw.iter()) {
                        *gene_counts.entry(*idx as usize).or_insert(0) += count as u32;
                        library_size += count as u32;
                    }
                }

                // sort for CSR format
                let mut entries: Vec<(usize, u32)> = gene_counts.into_iter().collect();
                entries.sort_by_key(|(idx, _)| *idx);

                let indices: Vec<usize> = entries.iter().map(|(idx, _)| *idx).collect();
                let raw_counts: Vec<u32> = entries.iter().map(|(_, count)| *count).collect();
                let norm_counts: Vec<f32> = entries
                    .iter()
                    .map(|(_, count)| {
                        let norm = (*count as f32 / library_size as f32) * target_size;
                        (norm + 1.0).ln()
                    })
                    .collect();

                (indices, raw_counts, norm_counts)
            })
            .collect();

        for (indices, raw_counts, norm_counts) in results {
            all_indices.extend(indices);
            all_data.extend(raw_counts);
            all_data_norm.extend(norm_counts);
            all_indptr.push(all_indices.len());
        }
    }

    CompressedSparseData::new_csr(
        &all_data,
        &all_indices,
        &all_indptr,
        Some(&all_data_norm),
        (n_metacells, n_genes),
    )
}

/// Convert metacell groups to flat assignments, handling unassigned cells
///
/// ### Params
///
/// * `metacells` - Vector of cell groups (metacell → [cells])
/// * `n_cells` - Total number of cells
///
/// ### Returns
///
/// Flat assignment vector where assignments[cell_id] = Some(metacell_id)
/// or None if cell is unassigned
pub fn metacells_to_assignments(metacells: &[&[usize]], n_cells: usize) -> Vec<Option<usize>> {
    let mut assignments = vec![None; n_cells];

    for (metacell_id, &cells) in metacells.iter().enumerate() {
        for &cell_id in cells {
            if cell_id < n_cells {
                assignments[cell_id] = Some(metacell_id);
            }
        }
    }

    assignments
}

// /// Transform the community memberships to metacells
// ///
// /// ### Params
// ///
// /// * `membership` - Community membership vector
// /// * `k` - Number of member cells
// ///
// /// ### Returns
// ///
// /// The community to metacell assignements
// pub fn communities_to_metacells(membership: &[usize], k: usize) -> Vec<Vec<usize>> {
//     let mut metacells = vec![Vec::new(); k];
//     for (cell_id, &metacell_id) in membership.iter().enumerate() {
//         if metacell_id < k {
//             metacells[metacell_id].push(cell_id);
//         }
//     }
//     metacells
// }

/// Convert assignments to R-friendly list format with unassigned cell handling
///
/// Returns -1 for unassigned cells (R convention for missing/unassigned).
///
/// ### Params
///
/// * `assignments` - Vector where assignments[cell_id] = Some(metacell_id) or None
/// * `n_cells` - Total number of cells
/// * `k` - Number of metacells
///
/// ### Returns
///
/// R List with -1 indicating unassigned cells
pub fn assignments_to_r_list(assignments: &[Option<usize>], n_cells: usize) -> List {
    let r_assignments: Vec<i32> = assignments
        .iter()
        .map(|&x| match x {
            Some(id) => (id + 1) as i32,
            None => -1,
        })
        .collect();

    let n_unassigned = assignments.iter().filter(|x| x.is_none()).count();

    let actual_k = assignments
        .iter()
        .filter_map(|&x| x)
        .max()
        .map(|x| x + 1)
        .unwrap_or(0);

    let mut metacells = vec![Vec::new(); actual_k];

    for (cell_id, &metacell_id) in assignments.iter().enumerate() {
        if let Some(id) = metacell_id {
            metacells[id].push((cell_id + 1) as i32);
        }
    }

    let unassigned: Vec<i32> = assignments
        .iter()
        .enumerate()
        .filter_map(|(cell_id, &x)| {
            if x.is_none() {
                Some((cell_id + 1) as i32)
            } else {
                None
            }
        })
        .collect();

    let metacells_list: List = metacells.into_iter().map(Robj::from).collect();

    list!(
        assignments = r_assignments,
        metacells = metacells_list,
        unassigned = unassigned,
        n_metacells = actual_k,
        n_cells = n_cells,
        n_unassigned = n_unassigned
    )
}

/// Remap subset assignments back to original index space
///
/// Takes metacell assignments computed on a subset of cells and maps them
/// back to the full original index space. Cells not in the subset will have
/// `None` assignments.
///
/// ### Params
///
/// * `subset_assignments` - Vector of metacell assignments in subset index space
/// * `subset_to_orig` - Mapping from subset indices to original indices
/// * `n_total` - Total number of cells in original space
///
/// ### Return
///
/// Vector of assignments in original index space with `None` for cells not
/// in the subset
pub fn remap_assignments_to_original(
    subset_assignments: &[Option<usize>],
    subset_to_orig: &[usize],
    n_total: usize,
) -> Vec<Option<usize>> {
    let mut full_assignments = vec![None; n_total];
    for (subset_idx, &metacell_id) in subset_assignments.iter().enumerate() {
        if let Some(orig_idx) = subset_to_orig.get(subset_idx) {
            full_assignments[*orig_idx] = metacell_id;
        }
    }
    full_assignments
}

/// Remap metacell indices from subset space to original space
///
/// Takes metacell groups where cell indices are in subset space and transforms
/// all indices back to original space. This is used when metacells are computed
/// on a subset of cells but need to be aggregated from the full dataset.
///
/// ### Params
///
/// * `metacells` - Vector of metacell groups with cell indices in subset space
/// * `subset_to_orig` - Mapping from subset indices to original indices
///
/// ### Return
///
/// Vector of metacell groups with cell indices in original space
pub fn remap_metacells_to_original(
    metacells: &[&[usize]],
    subset_to_orig: &[usize],
) -> Vec<Vec<usize>> {
    metacells
        .iter()
        .map(|&cells| cells.iter().map(|&idx| subset_to_orig[idx]).collect())
        .collect()
}

////////////////////////
// hdWGCNA meta cells //
////////////////////////

/// Structure for the MetaCell parameters
///
/// ### Fields
///
/// ** Meta cell params**
///
/// * `max_shared` - Maximum number of shared cells for the meta cell
///   aggregation
/// * `target_no_metacells` - Number of target meta cells.
/// * `max_iter` - Maximum iterations for the algorithm.
///
/// **General kNN params**
///
/// * `k` - Number of neighbours for the kNN algorithm.
/// * `knn_method` - Which method to use for the generation of the kNN graph.
///   One of `"hnsw"`, `"annoy"` or `"nndescent"`
/// * `ann_dist` - The distance metric for the approximate nearest neighbour
///   search. One of `"cosine"` or `"euclidean"`.
///
/// **Annoy**
///
/// * `n_tree` - Number of trees for the generation of the index
/// * `search_budget` - Search budget during querying
///
/// **NN Descent**
///
/// * `max_iter` - Maximum iterations for the algorithm
/// * `rho` - Sampling rate for the algorithm
/// * `delta` - Early termination criterium
#[derive(Clone, Debug)]
pub struct MetaCellParams {
    // meta cell params
    pub max_shared: usize,
    pub target_no_metacells: usize,
    pub max_iter: usize,
    // general knn params
    pub knn_params: KnnParams,
}

impl MetaCellParams {
    /// Generate the MetaCellParams from an R list
    ///
    /// ### Params
    ///
    /// * `r_list` - The R list with the parameters.
    ///
    /// ### Return
    ///
    /// The `MetaCellParams` structure.
    pub fn from_r_list(r_list: List) -> Self {
        let knn_params = KnnParams::from_r_list(r_list.clone());
        let meta_cell_params = r_list.into_hashmap();

        // meta cell
        let max_shared = meta_cell_params
            .get("max_shared")
            .and_then(|v| v.as_integer())
            .unwrap_or(15) as usize;
        let target_no_metacells = meta_cell_params
            .get("target_no_metacells")
            .and_then(|v| v.as_integer())
            .unwrap_or(1000) as usize;
        let max_iter = meta_cell_params
            .get("max_iter")
            .and_then(|v| v.as_integer())
            .unwrap_or(5000) as usize;

        Self {
            max_shared,
            target_no_metacells,
            max_iter,
            knn_params,
        }
    }
}

/// Select meta cells
///
/// ### Params
///
/// * `nn_map` - Nearest neighbours with self
/// * `max_shared` - Maximum number of shared neighbours to allow
/// * `target_no` - Target number of meta cells
/// * `max_iter` - Maximum iterations for the algorithm
/// * `seed` - seed for reproducibility purposes
///
/// ### Returns
///
/// Borrowed slice of the selected meta cell indices
pub fn identify_meta_cells(
    nn_map: &[Vec<usize>],
    max_shared: usize,
    target_no: usize,
    max_iter: usize,
    seed: usize,
    verbose: bool,
) -> Vec<&[usize]> {
    let mut rng = rand::rngs::StdRng::seed_from_u64(seed as u64);
    let k = nn_map[0].len();
    let k2 = k * 2;

    let mut good_choices: Vec<usize> = (0..nn_map.len()).collect();
    let mut chosen: Vec<usize> = Vec::new();

    if let Some(&first) = good_choices.choose(&mut rng) {
        chosen.push(first);
        good_choices.retain(|&x| x != first);
    }

    let mut it = 0;

    // cache the HashSets to avoid regeneration during loops
    let mut set_cache: FxHashMap<usize, FxHashSet<usize>> = FxHashMap::default();

    // bootstrap meta cells
    while !good_choices.is_empty() && chosen.len() < target_no && it < max_iter {
        it += 1;

        // sample remaining cells
        let choice_idx = rng.random_range(0..good_choices.len());
        let candidate = good_choices[choice_idx];
        good_choices.remove(choice_idx);

        // check overlap with existing meta cells
        set_cache
            .entry(candidate)
            .or_insert_with(|| nn_map[candidate].iter().copied().collect());

        let mut max_overlap = 0;
        for &existing in &chosen {
            set_cache
                .entry(existing)
                .or_insert_with(|| nn_map[existing].iter().copied().collect());

            let candidate_set = &set_cache[&candidate];
            let existing_set = &set_cache[&existing];

            let shared = k2 - candidate_set.union(existing_set).count();
            max_overlap = max_overlap.max(shared);
        }

        if verbose && it % 10000 == 0 {
            println!(
                "Meta cell neighbour search - iter {} out of {} max iters",
                it, max_iter
            );
        }

        if max_overlap <= max_shared {
            chosen.push(candidate);
        }
    }

    chosen
        .iter()
        .map(|&center| nn_map[center].as_slice())
        .collect()
}

///////////////
// SuperCell //
///////////////

/// Structure for the SuperCell parameters
///
/// ### Fields
///
/// **SuperCell params**
///
/// * `walk_length` - Walk length for the Walktrap algorithm
/// * `graining_factor` - Graining level of data (proportion of number of single
///   cells in the initial dataset to the number of metacells in the final
///   dataset)
/// * `linkage_dist` - Which type of distance metric to use for the linkage.
///
/// **General kNN params**
///
/// * `knn_params` - All of the kNN parameters
#[derive(Clone, Debug)]
pub struct SuperCellParams {
    // supercell params
    pub walk_length: usize,
    pub graining_factor: f64,
    pub linkage_dist: String,
    // knn
    pub knn_params: KnnParams,
}

impl SuperCellParams {
    /// Generate the SuperCellParams from an R list
    ///
    /// ### Params
    ///
    /// * `r_list` - The R list with the parameters
    ///
    /// ### Return
    ///
    /// The `SuperCellParams` structure
    pub fn from_r_list(r_list: List) -> Self {
        let knn_params = KnnParams::from_r_list(r_list.clone());

        let params = r_list.into_hashmap();

        // supercell
        let walk_length = params
            .get("walk_length")
            .and_then(|v| v.as_integer())
            .unwrap_or(3) as usize;

        let graining_factor = params
            .get("graining_factor")
            .and_then(|v| v.as_real())
            .unwrap_or(50.0);

        let linkage_dist = params
            .get("linkage_dist")
            .and_then(|v| v.as_str())
            .unwrap_or("average")
            .to_string();

        Self {
            walk_length,
            graining_factor,
            linkage_dist,
            knn_params,
        }
    }
}

/// SuperCell algorithm
///
/// ### Params
///
/// * `knn_mat` - The kNN matrix
/// * `walk_length` - Walk length for the Walktrap algorithm
/// * `no_meta_cells` - Number of communities, i.e., metacells to identify
/// * `linkage_dist` - The distance metric to use for the linkage.
/// * `verbose` - Controls the verbosity of the function
///
/// ### Returns
///
/// Membership of included cells to MetaCells.
pub fn supercell(
    knn_mat: &[Vec<usize>],
    walk_length: usize,
    no_meta_cells: usize,
    linkage_dist: &str,
    verbose: bool,
) -> Vec<usize> {
    let knn_graph = knn_to_sparse_graph(knn_mat);
    walktrap_sparse_graph(
        &knn_graph,
        walk_length,
        no_meta_cells,
        linkage_dist,
        verbose,
    )
}

////////////////////
// Pseudo-bulking //
////////////////////

/// Pseudo-bulk data across cells based on cell indices
///
/// ### Params
///
/// * `f_path` - File path to the cell-based binary file.
/// * `cell_indices` - Slice of indices to pseudo-bulk.
/// * `verbose` - Controls verbosity of the function.
///
/// ### Return
///
/// Matrix of samples x genes pseudo-bulked.
#[allow(dead_code)]
pub fn get_pseudo_bulked_counts(
    f_path: &str,
    cell_indices: &[Vec<usize>],
    verbose: bool,
) -> Mat<f64> {
    let reader = ParallelSparseReader::new(f_path).unwrap();

    let n_genes = reader.get_header().total_genes;
    let n_groups = cell_indices.len();

    // Initialize matrix: groups (rows) × genes (columns)
    let mut result = Mat::zeros(n_groups, n_genes);

    for (group_idx, indices) in cell_indices.iter().enumerate() {
        let start_group = Instant::now();
        let chunks = reader.read_cells_parallel(indices);

        for chunk in chunks {
            for (value, &gene_idx) in chunk.data_raw.iter().zip(chunk.indices.iter()) {
                result[(group_idx, gene_idx as usize)] += *value as f64;
            }
        }

        if verbose && (group_idx + 1) % 10 == 0 {
            let elapsed = start_group.elapsed();
            let pct_complete = ((group_idx + 1) as f32 / n_groups as f32) * 100.0;
            println!(
                "Processed group {} out of {} (took {:.2?}, completed {:.1}%)",
                group_idx + 1,
                n_groups,
                elapsed,
                pct_complete
            );
        }
    }

    result
}
