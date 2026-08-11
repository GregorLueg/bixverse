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
## usethis namespace: end
NULL
