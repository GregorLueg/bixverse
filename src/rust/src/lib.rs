mod hypergeom_helpers;
mod utils_r_rust;
mod hypergeom_fun;
mod geom_elim_fun;

use extendr_api::prelude::*;

extendr_module! {
    mod BIXverse;
    use hypergeom_fun;
}