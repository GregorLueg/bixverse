# classes ----

## bixverse_base_class ----

#' bixverse base class
#'
#' @description
#' Generic base class that is used for inheritance in certain common methods
#' across classes.
#'
#' @section Properties:
#' \describe{
#'   \item{params}{A (nested) list that will store all the parameters of the
#'   applied function.}
#'   \item{final_results}{A data.table that will contain the final results.}
#' }
#'
#' @return Returns the S7 object for further operations.
bixverse_base_class <- S7::new_class(
  # Name
  name = "bixverse_base_class",
  properties = list(params = S7::class_list),

  constructor = function() {
    S7::new_object(
      S7::S7_object(),
      params = list(),
      final_results = data.table()
    )
  }
)

# methods ----

## getters ----

?S7::new_generic

#' Get the parameters that were used.
#'
#' @description Extracts params from the `bixverse_base_class` (or child)
#' class and has options to return (pretty) JSONs
#'
#' @param bixverse_base_class The underlying `bixverse_base_class` class.
#' @param to_json Shall the params be returned as a JSON string.
#' @param pretty_json Shall the params be returned as a pretty JSON string.
#'
#' @return Depending on parameters either the R list or a (pretty) JSON string.
#'
#' @export
get_params <- S7::new_generic(
  name = "get_params",
  dispatch_args = "bixverse_base_class",
  fun = function(bixverse_base_class,
                 to_json = FALSE,
                 pretty_json = FALSE) {
    S7::S7_dispatch()
  })

#' @method get_params bixverse_base_class
#' @export
S7::method(get_params, bixverse_base_class) <-
  function(bixverse_base_class,
           to_json = FALSE,
           pretty_json = FALSE) {
    # Checks
    checkmate::assertClass(
      bixverse_base_class, "bixverse::bixverse_base_class"
    )
    checkmate::qassert(to_json, "B1")
    checkmate::qassert(pretty_json, "B1")

    # Body
    to_ret <- S7::prop(bixverse_base_class, "params")
    if (to_json) {
      to_ret <- jsonlite::toJSON(to_ret)
    }
    if (to_json && pretty_json) {
      to_ret <- jsonlite::prettify(to_ret)
    }

    return(to_ret)
  }


#' Get the final results from the class
#'
#' @description Get the final results from the a
#' [bixverse::bixverse_base_class()] class.
#'
#' @param bixverse_base_class The underlying [bixverse::bixverse_base_class()]
#' class. This class is usually inherited by other classes.
#'
#' @return Returns the final results if any have been stored in the class.
#'
#' @export
get_results <- S7::new_generic(
  name = "get_results",
  dispatch_args = "bixverse_base_class",
  fun = function(bixverse_base_class) {
    S7::S7_dispatch()
  })

#' @export
#' @method get_results bixverse_generic_class
S7::method(get_results, bixverse_base_class) <-
  function(bixverse_base_class) {
    # Checks
    checkmate::assertClass(
      bixverse_base_class, "bixverse::bixverse_base_class"
    )

    # Return
    return(S7::prop(bixverse_base_class, "final_results"))
  }
