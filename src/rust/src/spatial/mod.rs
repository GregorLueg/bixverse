//! Everything and anything related to the Rust <> R interface for spatial
//! transcriptomics.
//!
//! Extends the single cell interface rather than replacing it: the counts sit
//! in the same binary store and the same readers get used, spots are cells
//! with coordinates attached.

pub mod r_sp_analysis;
pub mod r_sp_graph;
pub mod r_sp_image;
pub mod r_sp_processing;
pub mod utils;
