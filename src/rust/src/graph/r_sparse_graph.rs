//! The neighbour graph as an opaque Rust handle.
//!
//! R only ever sees `ExternalPtr<SparseGraph<f32>>`. The graph stays in Rust
//! between the kNN symmetrisation and its consumers (clustering, label
//! smoothing, export), so it is never flattened into R vectors on the hot path.
//! The pointer does not survive `saveRDS()`; persistence goes through
//! `rs_graph_to_bytes()` / `rs_graph_from_bytes()`.

use bixverse_r::graph::graph_structures::{edge_pairs_to_sparse_graph, SparseGraph};
use bixverse_r::graph::leiden::{leiden_sparse_graph, LeidenParams};
use bixverse_r::prelude::*;
use extendr_api::prelude::*;

/////////////
// extendR //
/////////////

extendr_module! {
    mod r_sparse_graph;
    fn rs_graph_info;
    fn rs_graph_to_bytes;
    fn rs_graph_from_bytes;
    fn rs_graph_to_edges;
    fn rs_graph_from_edges;
    fn rs_graph_to_csr;
    fn rs_graph_leiden;
}

/// The graph type behind every pointer this module hands out.
pub type GraphPtr = ExternalPtr<SparseGraph<f32>>;

///////////////
// Functions //
///////////////

/// Summary of a Rust neighbour graph
///
/// @description
/// `r lifecycle::badge("experimental")`
/// Reads the shape of the graph behind an external pointer.
///
/// @param ptr External pointer to the Rust graph.
///
/// @returns A list with
/// \itemize{
///   \item n_nodes - Integer. Number of nodes.
///   \item n_edges - Numeric. Number of (undirected) edges. Numeric, since a
///   very large graph can hold more than `.Machine$integer.max` edges.
///   \item directed - Logical. Whether the graph is directed.
/// }
///
/// @export
///
/// @keywords internal
#[extendr]
fn rs_graph_info(ptr: GraphPtr) -> List {
    list!(
        n_nodes = ptr.get_node_number() as i32,
        n_edges = ptr.get_edge_number() as f64,
        directed = ptr.is_directed()
    )
}

/// Serialise a Rust neighbour graph to bytes
///
/// @description
/// `r lifecycle::badge("experimental")`
/// The byte layout carries a magic, a format version and the float width, so
/// reading back a foreign or outdated blob errors rather than decoding into
/// garbage.
///
/// @param ptr External pointer to the Rust graph.
///
/// @returns Raw vector.
///
/// @export
///
/// @keywords internal
#[extendr]
fn rs_graph_to_bytes(ptr: GraphPtr) -> extendr_api::Result<Raw> {
    let bytes = ptr.to_bytes().to_extendr()?;
    Ok(Raw::from_bytes(&bytes))
}

/// Read a Rust neighbour graph back from bytes
///
/// @description
/// `r lifecycle::badge("experimental")`
/// Inverse of `rs_graph_to_bytes()`. The CSR is validated on the way in.
///
/// @param bytes Raw vector written by `rs_graph_to_bytes()`.
///
/// @returns External pointer to the Rust graph.
///
/// @export
///
/// @keywords internal
#[extendr]
fn rs_graph_from_bytes(bytes: Raw) -> extendr_api::Result<GraphPtr> {
    let graph = SparseGraph::<f32>::from_bytes(bytes.as_slice()).to_extendr()?;
    Ok(ExternalPtr::new(graph))
}

/// Flatten a Rust neighbour graph into an edge list
///
/// @description
/// `r lifecycle::badge("experimental")`
/// Each undirected edge comes out once, ready for `igraph::add_edges()`.
///
/// @param ptr External pointer to the Rust graph.
///
/// @returns A list with
/// \itemize{
///   \item edges - Integer vector with the flattened edge pairs (1-indexed!),
///   `c(from_1, to_1, from_2, to_2, ...)`.
///   \item weights - Numeric vector. One weight per pair.
/// }
///
/// @export
///
/// @keywords internal
#[extendr]
fn rs_graph_to_edges(ptr: GraphPtr) -> List {
    let (edges, weights) = ptr.to_edge_pairs();
    let edges: Vec<i32> = edges.iter().map(|&x| x as i32 + 1).collect();
    let weights: Vec<f64> = weights.iter().map(|&w| w as f64).collect();

    list!(edges = edges, weights = weights)
}

/// Build a Rust neighbour graph from an edge list
///
/// @description
/// `r lifecycle::badge("experimental")`
/// Builds an undirected graph. Each edge should be listed once; duplicates are
/// kept as separate entries, as in igraph. Weights are stored as `f32`.
///
/// @param edges Integer vector. Flattened edge pairs (1-indexed!),
/// `c(from_1, to_1, from_2, to_2, ...)`.
/// @param weights Numeric vector. One weight per pair.
/// @param n_nodes Integer. Number of nodes.
///
/// @returns External pointer to the Rust graph.
///
/// @export
///
/// @keywords internal
#[extendr]
fn rs_graph_from_edges(
    edges: &[i32],
    weights: &[f64],
    n_nodes: usize,
) -> extendr_api::Result<GraphPtr> {
    if edges.len() != 2 * weights.len() {
        return Err(Error::Other(format!(
            "edges ({}) must hold two entries per weight ({}).",
            edges.len(),
            weights.len()
        )));
    }

    let mut pairs = Vec::with_capacity(edges.len());
    for &x in edges {
        if x < 1 || x as usize > n_nodes {
            return Err(Error::Other(format!(
                "Node index {} is outside 1..={}.",
                x, n_nodes
            )));
        }
        pairs.push((x - 1) as usize);
    }
    let weights: Vec<f32> = weights.iter().map(|&w| w as f32).collect();

    Ok(ExternalPtr::new(edge_pairs_to_sparse_graph(
        &pairs, &weights, n_nodes,
    )))
}

/// Extract the CSR adjacency of a Rust neighbour graph
///
/// @description
/// `r lifecycle::badge("experimental")`
/// The symmetric adjacency, rows sorted, as ScanPy stores it under
/// `obsp/connectivities`.
///
/// @param ptr External pointer to the Rust graph.
///
/// @returns A list with
/// \itemize{
///   \item indptr - Integer vector. Row pointers (0-indexed).
///   \item indices - Integer vector. Column indices (0-indexed).
///   \item data - Numeric vector. Edge weights.
/// }
///
/// @export
///
/// @keywords internal
#[extendr]
fn rs_graph_to_csr(ptr: GraphPtr) -> List {
    let n = ptr.get_node_number();
    let mut indptr: Vec<i32> = Vec::with_capacity(n + 1);
    let mut indices: Vec<i32> = Vec::new();
    let mut data: Vec<f64> = Vec::new();

    indptr.push(0);
    for i in 0..n {
        let (nbrs, ws) = ptr.get_neighbours_raw(i);
        indices.extend(nbrs.iter().map(|&j| j as i32));
        data.extend(ws.iter().map(|&w| w as f64));
        indptr.push(indices.len() as i32);
    }

    list!(indptr = indptr, indices = indices, data = data)
}

/// Leiden and Louvain clustering on a Rust neighbour graph
///
/// @description
/// `r lifecycle::badge("experimental")`
/// Same engine and parameters as `rs_leiden_clustering()`, but on the graph
/// behind the pointer, so nothing is flattened into R vectors. Louvain is
/// `refine = FALSE`.
///
/// @param ptr External pointer to the Rust graph.
/// @param leiden_params List. See `rs_leiden_clustering()`.
///
/// @returns Integer vector with the 1-based community of each node; `1` is the
/// largest community.
///
/// @export
///
/// @keywords internal
#[extendr]
fn rs_graph_leiden(ptr: GraphPtr, leiden_params: List) -> extendr_api::Result<Vec<i32>> {
    let params = LeidenParams::from_r_list(leiden_params)?;
    let membership = leiden_sparse_graph(&*ptr, Some(params)).to_extendr()?;

    Ok(membership.iter().map(|&x| (x + 1) as i32).collect())
}
