# Libraries ----

library(here)
setwd(here())

library(devtools)

# Package-related functions ----

## License ----

use_mit_license()

## Dependencies ----

use_package("data.table")
use_package("parallel")
use_package("furrr")
use_package("checkmate")
use_package("magrittr")

## Functions / R stuff ----

use_r("GSE_hypergeom")
