# classes ----

## bixverse_generic_class ----

bixverse_generic_class <- S7::new_class(
  # Name
  name = "bixverse_generic_class",
  properties = list(params = S7::class_list),

  #' BIXverse generic class
  #'
  #' @description
  #' Generic class that is used for inheritance in certain common methods
  #' across classes.
  #'
  #' @return Returns the S7 object for further operations.
  #'
  #' @export
  constructor = function() {
    S7::new_object(
      S7::S7_object(),
      params = list()
    )
  }
)

# methods ----

## print ----

print <- S7::new_generic("print", "x", function(x, ...) {
  S7::S7_dispatch()
})

## show ----


## getters ----

#' Get the parameters that were used.
#'
#' @description
#' This method accesses the params slot and can return R lists or JSON strings.
#'
#' @export
get_params <- S7::new_generic("get_params", "bixverse_generic_class")

#' @name get_params
#'
#' @description Extracts params from the `bixverse_generic_class` (or child)
#' class and has options to return (pretty) JSONs
#'
#' @usage get_params(
#'  bixverse_generic_class,
#'  to_json = FALSE,
#'  pretty_json = FALSE
#' )
#'
#' @param bixverse_generic_class The underlying `bixverse_generic_class` class.
#' @param to_json Shall the params be returned as a JSON string.
#' @param pretty_json Shall the params be returned as a pretty JSON string.
#'
#' @return Depending on parameters either the R list or a (pretty) JSON string.
#'
#' @method get_params bixverse_generic_class
S7::method(get_params, bixverse_generic_class) <-
  function(bixverse_generic_class,
           to_json = FALSE,
           pretty_json = FALSE) {
    # Checks
    checkmate::assertClass(
      bixverse_generic_class, "BIXverse::bixverse_generic_class"
    )
    checkmate::qassert(to_json, "B1")
    checkmate::qassert(pretty_json, "B1")

    # Body
    to_ret <- S7::prop(bixverse_generic_class, "params")
    if (to_json) {
      to_ret <- jsonlite::toJSON(to_ret)
    }
    if (to_json && pretty_json) {
      to_ret <- jsonlite::prettify(to_ret)
    }

    return(to_ret)
  }
