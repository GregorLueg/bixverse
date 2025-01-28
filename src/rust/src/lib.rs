mod hypergeom_helpers;
mod geom_elim_helpers;
mod r_rust_utils;
mod hypergeom_fun;
mod stats_fun;
mod stats_utils;

use extendr_api::prelude::*;

extendr_module! {
    mod BIXverse;
    use hypergeom_fun;
    use stats_fun;
}