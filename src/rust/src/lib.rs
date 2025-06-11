mod helpers_fgsea;
mod helpers_geom_elim;
mod helpers_hypergeom;
mod helpers_ica;
mod helpers_linalg;
mod helpers_ontology;
mod helpers_rbh;

mod fun_coremo;
mod fun_fgsea;
mod fun_helpers;
mod fun_hypergeom;
mod fun_ica;
mod fun_linalg;
mod fun_ontology;
mod fun_rbh;
mod fun_stats;

mod utils_r_rust;
mod utils_rust;
mod utils_stats;

mod macro_assertions;

use extendr_api::prelude::*;

extendr_module! {
    mod bixverse;
    use fun_fgsea;
    use fun_hypergeom;
    use fun_stats;
    use fun_rbh;
    use fun_linalg;
    use fun_helpers;
    use fun_ica;
    use fun_ontology;
    use fun_coremo;
}
