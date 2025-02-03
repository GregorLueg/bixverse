mod helpers_hypergeom;
mod helpers_geom_elim;
mod helpers_rbh;

mod fun_hypergeom;
mod fun_stats;
mod fun_rbh;

mod utils_stats;
mod utils_r_rust;
mod utils_rust;

use extendr_api::prelude::*;

extendr_module! {
    mod bixverse;
    use fun_hypergeom;
    use fun_stats;
    use fun_rbh;
}