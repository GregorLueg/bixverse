mod helpers;
mod r_function;
mod utils;

use extendr_api::prelude::*;

pub use r_function::r_coremo;
pub use r_function::r_fgsea;
pub use r_function::r_graphs;
pub use r_function::r_gsva;
pub use r_function::r_helpers;
pub use r_function::r_hypergeom;
pub use r_function::r_ica;
pub use r_function::r_linalg;
pub use r_function::r_ontology;
pub use r_function::r_rbh;
pub use r_function::r_single_cell_obj;
pub use r_function::r_sparse_graph_dict;
pub use r_function::r_sparse_struct;
pub use r_function::r_stats;

extendr_module! {
    mod bixverse;
    use r_coremo;
    use r_fgsea;
    use r_graphs;
    use r_gsva;
    use r_helpers;
    use r_hypergeom;
    use r_ica;
    use r_linalg;
    use r_ontology;
    use r_rbh;
    use r_single_cell_obj;
    use r_sparse_graph_dict;
    use r_sparse_struct;
    use r_stats;
}
