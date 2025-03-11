# Libraries ----

library(here)
setwd(here())

library(devtools)

# Package-related functions ----

## License ----

use_mit_license()

## Dependencies ----

use_package("data.table")
use_package("checkmate")
use_package("magrittr")
use_package("S7")
use_package("zeallot")
use_package("arrow")
use_package("foreach")
use_package("rextendr")
use_package("igraph")
use_package("rlang")
use_package("dplyr")
use_package("purrr")
use_package("ggplot2")
use_package("R6")

## Functions / R stuff ----

use_r("GSE_hypergeom")
use_r("S7_go_object")

document()

check()


## Rust ----

rextendr::use_extendr()
rextendr::document()
