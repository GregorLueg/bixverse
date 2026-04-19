//! Contains the Rust to R interface for bixverse

#![allow(clippy::needless_range_loop)] // I want these loops!

use extendr_api::prelude::*;

pub mod base;
pub mod data;
pub mod enrichment;
pub mod graph;
pub mod meta_cell;
pub mod methods;
pub mod ontology;
pub mod single_cell;

// base
pub use base::r_cors_similarity;
pub use base::r_helpers;
pub use base::r_loess;
pub use base::r_rbf;
pub use base::r_stats;
pub use base::r_svd_pca;

// data
pub use data::r_h5;
pub use data::r_sparse;
pub use data::r_synthetic;

// enrichments
pub use enrichment::r_gsea;
pub use enrichment::r_gsva;
pub use enrichment::r_mitch;
pub use enrichment::r_oea;

// graphs
pub use graph::r_graph_clustering;
pub use graph::r_knn;
pub use graph::r_page_rank;
pub use graph::r_snf;

// methods
pub use methods::r_cistarget;
pub use methods::r_coremo;
pub use methods::r_dgrdl;
pub use methods::r_diffcor;
pub use methods::r_ica;
pub use methods::r_rbh;

// ontology
pub use ontology::r_go_elim;
pub use ontology::r_similiarity;

// single cell
pub use single_cell::r_count_obj;
pub use single_cell::r_sc_analysis;
pub use single_cell::r_sc_batch_corr;
pub use single_cell::r_sc_data;
pub use single_cell::r_sc_metacells;
pub use single_cell::r_sc_plot_extraction;
pub use single_cell::r_sc_processing;
// meta cell
pub use meta_cell::r_mc_processing;

/////////////
// extendR //
/////////////

extendr_module! {
    mod bixverse;
    // base
    use r_cors_similarity;
    use r_helpers;
    use r_rbf;
    use r_stats;
    use r_svd_pca;
    use r_loess;
    // data
    use r_sparse;
    use r_synthetic;
    use r_h5;
    // enrichment
    use r_gsea;
    use r_gsva;
    use r_mitch;
    use r_oea;
    // graphs
    use r_page_rank;
    use r_snf;
    use r_graph_clustering;
    use r_knn;
    // methods
    use r_coremo;
    use r_dgrdl;
    use r_diffcor;
    use r_ica;
    use r_rbh;
    use r_cistarget;
    // ontology
    use r_go_elim;
    use r_similiarity;
    // single cell
    use r_count_obj;
    use r_sc_batch_corr;
    use r_sc_processing;
    use r_sc_analysis;
    use r_sc_metacells;
    use r_sc_plot_extraction;
}
