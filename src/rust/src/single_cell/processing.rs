use extendr_api::*;

/// Structure to store QC information on cells
///
/// ### Fields
///
/// * `to_keep` - Boolean vector indicating if the cells passes thresholds
/// * `lib_size` - Optional library size of the cells
/// * `no_genes` - Optional number of genes of the cells
/// * `mt_perc` - Optional number of mitochondrial reads per cell
#[derive(Clone, Debug)]
#[allow(dead_code)]
pub struct CellQuality {
    pub to_keep: Vec<bool>,
    pub lib_size: Option<Vec<usize>>,
    pub no_genes: Option<Vec<usize>>,
}

/// Structure that stores minimum QC thresholds/info for single cell
///
/// ### Fields
///
/// * `min_unique_genes` - Minimum number of unique genes per cell/spot.
/// * `min_lib_size` - Minimum library size per cell/spot.
/// * `min_cells` - Minimum cells per gene.
/// * `target_size` - Target size for normalisation.
#[derive(Clone, Debug)]
pub struct MinCellQuality {
    pub min_unique_genes: usize,
    pub min_lib_size: usize,
    pub min_cells: usize,
    pub target_size: f32,
}

impl MinCellQuality {
    /// Generate the MinCellQuality params from an R list
    ///
    /// or default to sensible defaults
    ///
    /// ### Params
    ///
    /// * `r_list` - R list with the parameters
    ///
    /// ### Returns
    ///
    /// Self with the specified parameters.
    pub fn from_r_list(r_list: List) -> Self {
        let min_qc = r_list.into_hashmap();

        let min_unique_genes = min_qc
            .get("min_unique_genes")
            .and_then(|v| v.as_integer())
            .unwrap_or(100) as usize;

        let min_lib_size = min_qc
            .get("min_lib_size")
            .and_then(|v| v.as_integer())
            .unwrap_or(250) as usize;

        let min_cells = min_qc
            .get("min_cells")
            .and_then(|v| v.as_integer())
            .unwrap_or(10) as usize;

        let target_size = min_qc
            .get("target_size")
            .and_then(|v| v.as_real())
            .unwrap_or(1e5) as f32;

        MinCellQuality {
            min_unique_genes,
            min_lib_size,
            min_cells,
            target_size,
        }
    }
}
