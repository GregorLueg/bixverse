# SpatialSpotDuckDB R6 ---------------------------------------------------------

# Subclass of SingleCellDuckDB that adds spatial-specific helpers:
# - exp_id discovery (the obs column already exists from the multi-file
#   ingest path).
# - spot index lookup per experiment, returning 0-based indices in the
#   ORIGINAL obs space so the Rust streaming reader can use them
#   directly via its `original_index`-keyed lookup map.
# - sidecar table for metaspot assignments (used when
#   `generate_superspot_*` returns a new on-disk SpatialSpot).

#' @title Class for storing SpatialSpot experimental data in DuckDB
#'
#' @description
#' R6 subclass of [bixverse::SingleCellDuckDB()] that adds spatial-
#' specific helpers. The obs table is expected to carry an `exp_id`
#' column identifying which slide/section each spot belongs to.
#'
#' @export
#'
#' @import data.table
#'
#' @keywords internal
SpatialSpotDuckDB <- R6::R6Class(
  classname = "SpatialSpotDuckDB",
  inherit = SingleCellDuckDB,
  public = list(
    #' @description
    #' Returns the distinct experiment identifiers from the obs table.
    #'
    #' @return Character vector of `exp_id` values.
    get_exp_ids = function() {
      private$check_obs_exists()

      con <- private$connect_db()
      on.exit(
        {
          if (exists("con") && !is.null(con)) {
            tryCatch(DBI::dbDisconnect(con), error = function(e) invisible())
          }
        }
      )

      private$check_exp_id_column(con)

      data <- DBI::dbGetQuery(
        conn = con,
        statement = "SELECT DISTINCT exp_id FROM obs ORDER BY exp_id"
      )

      return(as.character(data$exp_id))
    },

    #' @description
    #' Returns the spot indices for a given experiment, expressed in
    #' the original (0-based) obs space - i.e. positions in the
    #' original obs table NOT renumbered after filtering. This matches
    #' the convention used by the Rust streaming reader.
    #'
    #' @param exp_id String. The experiment identifier.
    #' @param filtered Boolean. Honour the `to_keep` filter. Defaults
    #' to `TRUE`.
    #'
    #' @return Integer vector of 0-based original obs indices.
    get_spot_indices_for_exp = function(exp_id, filtered = TRUE) {
      checkmate::qassert(exp_id, "S1")
      checkmate::qassert(filtered, "B1")
      private$check_obs_exists()

      con <- private$connect_db()
      on.exit(
        {
          if (exists("con") && !is.null(con)) {
            tryCatch(DBI::dbDisconnect(con), error = function(e) invisible())
          }
        }
      )

      private$check_exp_id_column(con)

      sql_query <- if (filtered) {
        "SELECT cell_idx FROM obs WHERE exp_id = ? AND to_keep ORDER BY cell_idx"
      } else {
        "SELECT cell_idx FROM obs WHERE exp_id = ? ORDER BY cell_idx"
      }

      data <- DBI::dbGetQuery(
        conn = con,
        statement = sql_query,
        params = list(exp_id)
      )

      # cell_idx is 1-based in the DuckDB; subtract 1 for Rust indexing.
      return(as.integer(data$cell_idx) - 1L)
    },

    #' @description
    #' Create the `metaspot_assignments` sidecar table that links
    #' metaspots in this aggregated SpatialSpot back to the original
    #' spots in the parent. Composite primary key on
    #' `(metaspot_idx, original_spot_idx)`. Used by
    #' `generate_superspot_*` when emitting a new on-disk SpatialSpot.
    #'
    #' @param assignments A `data.frame`/`data.table` with columns
    #' `metaspot_idx` (integer), `original_spot_idx` (integer) and
    #' `exp_id` (character).
    #'
    #' @return Invisibly self.
    create_metaspot_assignments_table = function(assignments) {
      checkmate::assertDataFrame(assignments)
      checkmate::assertNames(
        names(assignments),
        must.include = c("metaspot_idx", "original_spot_idx", "exp_id")
      )
      checkmate::qassert(assignments$metaspot_idx, "I+")
      checkmate::qassert(assignments$original_spot_idx, "I+")
      checkmate::qassert(assignments$exp_id, "S+")

      con <- private$connect_db()
      on.exit(
        {
          if (exists("con") && !is.null(con)) {
            tryCatch(DBI::dbDisconnect(con), error = function(e) invisible())
          }
        }
      )

      DBI::dbExecute(con, "DROP TABLE IF EXISTS metaspot_assignments")
      DBI::dbExecute(
        con,
        paste(
          "CREATE TABLE metaspot_assignments (",
          "  metaspot_idx INTEGER NOT NULL,",
          "  original_spot_idx INTEGER NOT NULL,",
          "  exp_id VARCHAR NOT NULL,",
          "  PRIMARY KEY (metaspot_idx, original_spot_idx)",
          ")"
        )
      )

      DBI::dbWriteTable(
        conn = con,
        name = "metaspot_assignments",
        value = as.data.frame(assignments),
        append = TRUE
      )

      DBI::dbExecute(
        con,
        paste(
          "CREATE INDEX idx_msa_metaspot",
          "ON metaspot_assignments(metaspot_idx)"
        )
      )
      DBI::dbExecute(
        con,
        paste(
          "CREATE INDEX idx_msa_orig_spot",
          "ON metaspot_assignments(original_spot_idx)"
        )
      )

      invisible(self)
    },

    #' @description
    #' Retrieve metaspot assignment rows. Optional filters on
    #' `metaspot_idx`, `original_spot_idx` or `exp_id`.
    #'
    #' @param metaspot_idx Optional integer vector to filter on.
    #' @param original_spot_idx Optional integer vector to filter on.
    #' @param exp_id Optional string to filter on.
    #'
    #' @return A `data.table` with `metaspot_idx`, `original_spot_idx`,
    #' `exp_id` columns.
    get_metaspot_assignments = function(
      metaspot_idx = NULL,
      original_spot_idx = NULL,
      exp_id = NULL
    ) {
      checkmate::qassert(metaspot_idx, c("0", "I+"))
      checkmate::qassert(original_spot_idx, c("0", "I+"))
      checkmate::qassert(exp_id, c("0", "S1"))

      con <- private$connect_db()
      on.exit(
        {
          if (exists("con") && !is.null(con)) {
            tryCatch(DBI::dbDisconnect(con), error = function(e) invisible())
          }
        }
      )

      tables <- DBI::dbGetQuery(con, "SHOW TABLES")[, "name"]
      if (!"metaspot_assignments" %in% tables) {
        stop("metaspot_assignments table not found in the DB.")
      }

      where_clauses <- character()
      params <- list()

      if (!is.null(metaspot_idx)) {
        placeholders <- paste(rep("?", length(metaspot_idx)), collapse = ", ")
        where_clauses <- c(
          where_clauses,
          sprintf("metaspot_idx IN (%s)", placeholders)
        )
        params <- c(params, as.list(metaspot_idx))
      }
      if (!is.null(original_spot_idx)) {
        placeholders <- paste(
          rep("?", length(original_spot_idx)),
          collapse = ", "
        )
        where_clauses <- c(
          where_clauses,
          sprintf("original_spot_idx IN (%s)", placeholders)
        )
        params <- c(params, as.list(original_spot_idx))
      }
      if (!is.null(exp_id)) {
        where_clauses <- c(where_clauses, "exp_id = ?")
        params <- c(params, list(exp_id))
      }

      where_part <- if (length(where_clauses) > 0) {
        paste("WHERE", paste(where_clauses, collapse = " AND "))
      } else {
        ""
      }

      sql_query <- sprintf(
        "SELECT metaspot_idx, original_spot_idx, exp_id FROM metaspot_assignments %s",
        where_part
      )

      res <- data.table::setDT(DBI::dbGetQuery(
        conn = con,
        statement = sql_query,
        params = params
      ))

      return(res)
    }
  ),
  private = list(
    # Helper - confirms the obs table actually carries an exp_id column.
    # Fails fast with a clearer message than a raw SQL error if a user
    # opens a non-spatial directory through this subclass.
    check_exp_id_column = function(con) {
      cols <- DBI::dbGetQuery(
        con,
        "SELECT column_name FROM information_schema.columns WHERE table_name = 'obs'"
      )[["column_name"]]
      if (!"exp_id" %in% cols) {
        stop(
          "obs table is missing the 'exp_id' column required by SpatialSpot."
        )
      }
    }
  )
)
