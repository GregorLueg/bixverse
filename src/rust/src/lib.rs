mod core;
mod r_bindings;
mod utils;

use extendr_api::prelude::*;

// base
pub use r_bindings::r_base::r_cors_similarity;
pub use r_bindings::r_base::r_helpers;
pub use r_bindings::r_base::r_rbf;
pub use r_bindings::r_base::r_stats;
pub use r_bindings::r_base::r_svd_pca;

// data
pub use r_bindings::r_data::r_sparse;
pub use r_bindings::r_data::r_synthetic;

// enrichments
pub use r_bindings::r_enrichment::r_gsea;
pub use r_bindings::r_enrichment::r_gsva;
pub use r_bindings::r_enrichment::r_mitch;
pub use r_bindings::r_enrichment::r_oea;

// graphs
pub use r_bindings::r_graph::r_page_rank;
pub use r_bindings::r_graph::r_snf;

// methods
pub use r_bindings::r_methods::r_coremo;
pub use r_bindings::r_methods::r_dgrdl;
pub use r_bindings::r_methods::r_diffcor;
pub use r_bindings::r_methods::r_ica;
pub use r_bindings::r_methods::r_rbh;

// ontology
pub use r_bindings::r_ontology::r_go_elim;
pub use r_bindings::r_ontology::r_similiarity;

extendr_module! {
    mod bixverse;
    // base
    use r_cors_similarity;
    use r_helpers;
    use r_rbf;
    use r_stats;
    use r_svd_pca;
    // data
    use r_sparse;
    use r_synthetic;
    // enrichment
    use r_gsea;
    use r_gsva;
    use r_mitch;
    use r_oea;
    // graphs
    use r_page_rank;
    use r_snf;
    // methods
    use r_coremo;
    use r_dgrdl;
    use r_diffcor;
    use r_ica;
    use r_rbh;
    // ontology
    use r_go_elim;
    use r_similiarity;
}
