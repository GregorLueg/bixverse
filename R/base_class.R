# contains a base class definition and generics for common getters/setters that
# are applied on this base class

# classes ----------------------------------------------------------------------

## BixverseBaseClass -----------------------------------------------------------

#' BixverseBaseClass
#'
#' @description
#' Generic base class that is used for inheritance in certain common methods
#' across other classes.
#'
#' @section Properties:
#' \describe{
#'   \item{params}{A (nested) list that will store all the parameters of the
#'   applied function.}
#'   \item{final_results}{A data.table that will contain the final results.}
#' }
#'
#' @return Returns the S7 object for further operations.
#'
#' @export
BixverseBaseClass <- S7::new_class(
  # Name
  name = "BixverseBaseClass",
  properties = list(params = S7::class_list),
  constructor = function() {
    S7::new_object(
      S7::S7_object(),
      params = list(),
      final_results = data.table()
    )
  }
)

# methods ----------------------------------------------------------------------

## getters ---------------------------------------------------------------------

#' Get the parameters that were used.
#'
#' @description Extracts parameters from the `BixverseBaseClass` class (or
#' child classes) and has options to return (pretty) JSONs. This generic
#' also gets inherited by other classes and can be used to extract parameters.
#' Also, can dispatch to specific methods for certain S3 classes.
#'
#' @param object A class within bixverse that inherits from
#' [bixverse::BixverseBaseClass()] or defined S3 classes.
#' @param to_json Shall the params be returned as a JSON string.
#' @param pretty_json Shall the params be returned as a pretty JSON string.
#'
#' @return Depending on parameters either the R list or a (pretty) JSON string.
#'
#' @export
get_params <- S7::new_generic(
  name = "get_params",
  dispatch_args = "object",
  fun = function(object, to_json = FALSE, pretty_json = FALSE) {
    S7::S7_dispatch()
  }
)

#' @export
#'
#' @method get_params BixverseBaseClass
S7::method(get_params, BixverseBaseClass) <-
  function(object, to_json = FALSE, pretty_json = FALSE) {
    # Checks
    checkmate::assertClass(
      object,
      "bixverse::BixverseBaseClass"
    )
    checkmate::qassert(to_json, "B1")
    checkmate::qassert(pretty_json, "B1")

    # Body
    to_ret <- S7::prop(object, "params")
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
#' @description Get the final results from `BixverseBaseClass` class (or child
#' classes).
#'
#' @param object The underlying [bixverse::BixverseBaseClass()]
#' class. The class functionality is usually inherited by other S7 classes in
#' `bixverse`.
#'
#' @return Returns the final results if any have been stored in the class.
#'
#' @export
get_results <- S7::new_generic(
  name = "get_results",
  dispatch_args = "object",
  fun = function(object) {
    S7::S7_dispatch()
  }
)

#' @export
#'
#' @method get_results BixverseBaseClass
S7::method(get_results, BixverseBaseClass) <-
  function(object) {
    # Checks
    checkmate::assertClass(
      object,
      "bixverse::BixverseBaseClass"
    )

    # Return
    return(S7::prop(object, "final_results"))
  }
