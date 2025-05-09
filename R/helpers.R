# user options -----------------------------------------------------------------

#' Helper function for user option selection
#'
#' @param options Character vector. The options the user has to choose from.
#'
#' @returns The chosen option
select_user_option <- function(options) {
  cat("Please select an option:\n")
  for (i in seq_along(options)) {
    cat(i, ": ", options[i], "\n")
  }

  while (TRUE) {
    cat("\nEnter the number of your choice: ")
    user_choice <- as.integer(readLines(n = 1))

    if (
      !is.na(user_choice) && user_choice > 0 && user_choice <= length(options)
    ) {
      return(options[user_choice])
    } else {
      cat("Invalid selection. Please try again.\n")
    }
  }
}

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
