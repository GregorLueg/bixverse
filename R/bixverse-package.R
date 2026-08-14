#' @section Options:
#' \describe{
#'   \item{`bixverse.cache_check`}{Controls how strictly the single cell
#'   classes police stale cached artefacts (PCA, embeddings, kNN, sNN) against
#'   the current cell set. Unset, the default, is two tier: the getters warn,
#'   while [assert_sc_state()] errors inside the functions that hand cached
#'   indices to Rust. `"error"` promotes the getter warnings to errors too,
#'   `"warn"` spells out the default getter behaviour, and `"none"` disables
#'   both tiers. See [get_sc_cache_status()].}
#' }
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom lifecycle deprecated
# Matrix is only ever called as `Matrix::`, which leaves its namespace unloaded
# until the first such call. A `dgRMatrix` deserialised into that session has no
# working S4 dispatch yet, so `nrow()` returns NULL and `[` refuses to subset,
# and the first call after a readRDS/qs2 load fails. This import loads it with
# the package.
#' @importFrom Matrix sparseMatrix
## usethis namespace: end
NULL
