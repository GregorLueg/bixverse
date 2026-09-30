//! Rust to R interface for Bonsai trees over single cell and metacell counts.
//!
//! Every cell, gene and node index that crosses this boundary, in and out, is
//! 0-indexed. The root's parent is `-1`.

use bixverse_rs::prelude::*;
use bixverse_rs::single_cell::mc_analysis::bonsai_mc::sanity_bonsai_mc;
use bixverse_rs::single_cell::sc_analysis::bonsai::{
    bonsai_layout, parse_bonsai_layout, sanity_bonsai_sc, BonsaiScParams,
};
use bixverse_rs::single_cell::sc_r_wrappers::{bonsai_sc_to_r_list, parents_from_r, parents_to_r};
use extendr_api::*;

use crate::meta_cell::utils::mc_list_to_sparse_u32;

////////////////////
// extendr Module //
////////////////////

extendr_module! {
    mod r_sc_bonsai;
    fn rs_sc_bonsai;
    fn rs_mc_bonsai;
    fn rs_bonsai_layout;
}

////////////
// Bonsai //
////////////

/// Build a Bonsai tree from single cell counts
///
/// @description
/// `r lifecycle::badge("experimental")`
/// Reads the raw counts of the given genes and cells from the binary stores,
/// runs Sanity on them for posterior log fold changes with error bars, builds
/// a Bonsai tree over the cells and lays it out in 2D. Sanity runs on the CPU.
///
/// @param f_path_gene String. Path to the `counts_genes.bin` file.
/// @param f_path_cell String. Path to the `counts_cells.bin` file. Supplies the
/// library sizes.
/// @param cell_indices Integer. The cell indices to use. (0-indexed!) Sets the
/// leaf order.
/// @param gene_indices Integer. The gene indices to use. (0-indexed!)
/// @param bonsai_params List. Parameter list, see
/// [bixverse::params_sc_bonsai()].
/// @param verbose Integer. `0L` - quiet; `1L` - normal verbosity; `2L` -
/// detailed verbosity.
///
/// @returns A list with
/// \itemize{
///  \item parent - Integer. Parent of each node (0-indexed!), `-1` for the
///  root. Nodes `0` to `n_leaves - 1` are the cells in `cell_indices` order.
///  \item branch - Numeric. Branch length above each node.
///  \item x - Numeric. Horizontal coordinate of each node.
///  \item y - Numeric. Vertical coordinate of each node.
///  \item n_leaves - Integer. Number of leaves.
///  \item loglik - Numeric. Final tree loglikelihood.
///  \item steps - List with `step` and `loglik` after each search step.
///  \item genes_used - Integer. Genes the tree was built on. (0-indexed!)
/// }
///
/// @export
///
/// @references de Groot, et al., Nat Biotechnol, 2026; Breda, et al., Nat
/// Biotechnol, 2021.
///
/// @keywords internal
#[extendr]
fn rs_sc_bonsai(
    f_path_gene: &str,
    f_path_cell: &str,
    cell_indices: Vec<i32>,
    gene_indices: Vec<i32>,
    bonsai_params: List,
    verbose: usize,
) -> extendr_api::Result<List> {
    let verbosity = parse_verbosity_level(verbose);
    let cell_indices = cell_indices.r_int_convert();
    let gene_indices = gene_indices.r_int_convert();
    let params = BonsaiScParams::from_r_list(bonsai_params)?;

    let gene_reader = ParallelSparseReader::new(f_path_gene).to_extendr()?;
    let cell_reader = ParallelSparseReader::new(f_path_cell).to_extendr()?;

    let res = sanity_bonsai_sc(
        &gene_reader,
        &cell_reader,
        &cell_indices,
        &gene_indices,
        &params,
        verbosity,
    )
    .to_extendr()?;

    Ok(bonsai_sc_to_r_list(res))
}

/// Build a Bonsai tree from metacell counts
///
/// @description
/// `r lifecycle::badge("experimental")`
/// As [bixverse::rs_sc_bonsai()], with the metacells' aggregated raw counts
/// from memory in place of the binary files. Each metacell is a leaf, and its
/// total counts over all genes are its library size.
///
/// @param sparse_data List. The raw metacell counts, see
/// [bixverse::mc_counts_to_list()] with `assay = "raw"`.
/// @param gene_indices Integer. The candidate genes. (0-indexed!)
/// @param bonsai_params List. Parameter list, see
/// [bixverse::params_sc_bonsai()].
/// @param verbose Integer. `0L` - quiet; `1L` - normal verbosity; `2L` -
/// detailed verbosity.
///
/// @returns The same list as [bixverse::rs_sc_bonsai()], with the metacells
/// as the leaves in their row order.
///
/// @export
///
/// @references de Groot, et al., Nat Biotechnol, 2026; Breda, et al., Nat
/// Biotechnol, 2021.
///
/// @keywords internal
#[extendr]
fn rs_mc_bonsai(
    sparse_data: List,
    gene_indices: Vec<i32>,
    bonsai_params: List,
    verbose: usize,
) -> extendr_api::Result<List> {
    let verbosity = parse_verbosity_level(verbose);
    let gene_indices = gene_indices.r_int_convert();
    let params = BonsaiScParams::from_r_list(bonsai_params)?;
    let counts = mc_list_to_sparse_u32(sparse_data)?;

    let res = sanity_bonsai_mc(&counts, &gene_indices, &params, verbosity).to_extendr()?;

    Ok(bonsai_sc_to_r_list(res))
}

/// Lay out an existing Bonsai tree
///
/// @description
/// `r lifecycle::badge("experimental")`
/// Computes a new 2D layout for a finished tree without searching again. The
/// tree is renumbered internally, so inferred ancestors can come back with
/// different indices than they went in with; leaves keep theirs. Replace the
/// whole node table with the output.
///
/// @param parent Integer. Parent of each node (0-indexed!), negative for the
/// root.
/// @param branch Numeric. Branch length above each node.
/// @param n_leaves Integer. Number of leaves, which occupy the first
/// `n_leaves` nodes.
/// @param layout String. One of `c("equal_angle", "equal_daylight",
/// "dendrogram")`.
/// @param hyperbolic Boolean. Project onto the hyperbolic disk.
///
/// @returns A list with `parent` (0-indexed!, `-1` for the root), `branch`, `x`
/// and `y`, all indexed by the tree's own node numbering.
///
/// @export
///
/// @keywords internal
#[extendr]
fn rs_bonsai_layout(
    parent: Vec<i32>,
    branch: Vec<f64>,
    n_leaves: usize,
    layout: &str,
    hyperbolic: bool,
) -> extendr_api::Result<List> {
    let layout = parse_bonsai_layout(layout)
        .ok_or_else(|| Error::Other(format!("Invalid Bonsai layout: {layout}")))?;

    let (parent, branch, coords) = bonsai_layout(
        parents_from_r(&parent),
        branch,
        n_leaves,
        layout,
        hyperbolic,
    )
    .to_extendr()?;

    Ok(list!(
        parent = parents_to_r(&parent),
        branch = branch,
        x = coords.x,
        y = coords.y
    ))
}
