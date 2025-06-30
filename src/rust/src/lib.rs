mod helpers;
mod r_function;

mod utils_r_rust;
mod utils_rust;
mod utils_stats;

mod macro_assertions;

use extendr_api::prelude::*;

pub use r_function::r_coremo;
pub use r_function::r_fgsea;
pub use r_function::r_graphs;
pub use r_function::r_helpers;
pub use r_function::r_hypergeom;
pub use r_function::r_ica;
pub use r_function::r_linalg;
pub use r_function::r_ontology;
pub use r_function::r_rbh;
pub use r_function::r_stats;

extendr_module! {
    mod bixverse;
    use r_fgsea;
    use r_hypergeom;
    use r_stats;
    use r_rbh;
    use r_linalg;
    use r_helpers;
    use r_ica;
    use r_ontology;
    use r_coremo;
    use r_graphs;
}
