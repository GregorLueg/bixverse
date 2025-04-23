# gene set libraries -----------------------------------------------------------

#' Get the Gene Ontology data human
#'
#' @returns data.table with the data ready for [bixverse::gene_ontology_data()].
#'
#' @export
get_go_human_data <- function() {
  path <- system.file("extdata", "go_data_hs.parquet", package = "bixverse")
  if (path != "") {
    go_data_dt <- arrow::read_parquet(path) %>%
      data.table::setDT()

    return(go_data_dt)
  } else {
    error("The expected .parquet file was not found.")
  }
}
