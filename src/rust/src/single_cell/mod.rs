//! Contains the core single cell functionalities that are exposed to R
//! There are dependencies on other parts of the package, specifically
//! around the generation of kNN graphs, sparse matrix methods, etc.

pub mod cell_aggregations;
pub mod dge_pathway_scores;
pub mod fast_ranking;
pub mod methods;
pub mod metrics;
pub mod processing;
pub mod sc_knn_snn;
