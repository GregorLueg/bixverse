use extendr_api::prelude::*;

use bixverse_rs::single_cell::sc_data::depracated_conversion::migrate_v2_to_v3_pair;

////////////////////
// extendr Module //
////////////////////

extendr_module! {
    mod r_sc_data;
}

/////////////
// Helpers //
/////////////

/// Convert old v2 files to v3 files
///
/// @param cell_input Path to the `counts_cells.bin` file.
/// @param cell_output Path to the `counts_cells_new.bin` file.
/// @param gene_input Path to the `counts_genes.bin` file.
/// @param gene_output Path to the `counts_genes_new.bin` file.
/// @param verbose Boolean. Controls verbosity
///
/// @returns A potential error if it occurs
///
/// @export
#[extendr]
fn rs_data_v2_3_conversion(
    cell_input: &str,
    cell_output: &str,
    gene_input: &str,
    gene_output: &str,
    verbose: bool,
) -> extendr_api::Result<()> {
    migrate_v2_to_v3_pair(cell_input, cell_output, gene_input, gene_output, verbose)
        .map_err(|e| e.to_string())?;
    Ok(())
}
