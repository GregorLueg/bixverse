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
    pub mt_perc: Option<Vec<f32>>,
}
