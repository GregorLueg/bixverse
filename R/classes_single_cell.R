# s7 ---------------------------------------------------------------------------

## single cell class -----------------------------------------------------------

#' @title bixverse single cell class (nightly!)
#'
#' @description
#' This is the `bixverse`-based single cell class. Under the hood it uses a
#' DuckDB for obs and vars storing, and a Rust-based binarised file format to
#' store the raw and normalised counts. In both cases, the idea is not to hold
#' any data that is not needed at a given point of time in memory, but leverage
#' speedy on-disk computations and streaming engines powered by Rust and DuckDB
#' to run the analysis.
#'
#' @param dir_data String. This is the directory in which the experimental files
#' will be
#'
#' @section Properties:
#' \describe{
#'   \item{db_connection}{This contains an R6 class with DuckDB pointers and
#'   wrappers to interact with the table-like data for this experiment.}
#'   \item{count_connection}{This contains an R6-like environment that points
#'   to Rust functions that can work on the counts more specifically.}
#'   \item{dir_data}{Path to the directory in which the data will be saved on
#'   disk.}
#' }
#'
#' @return Returns the `single_cell_exp` class for further operations.
#'
#' @export
single_cell_exp <- S7::new_class(
  name = "single_cell_exp",
  properties = list(
    db_connection = S7::class_any,
    count_connection = S7::class_any,
    dir_data = S7::class_character,
    cache = S7::class_list
  ),
  constructor = function(dir_data) {
    nightly_feature()
    # checks
    checkmate::assertPathForOutput(dir_data)

    # generate the Rust pointer
    count_connection <- SingeCellCountData$new(
      f_path_cells = file.path(dir_data, "counts_cells.bin"),
      f_path_genes = file.path(dir_data, "counts_genes.bin")
    )

    db_connection <- single_cell_duckdb_con$new(
      db_dir = dir_data
    )

    S7::new_object(
      S7::S7_object(),
      db_connection = db_connection,
      count_connection = count_connection,
      dir_data = dir_data,
      cache = list()
    )
  }
)

# r6 ---------------------------------------------------------------------------

## duckdb connector ------------------------------------------------------------

### base class -----------------------------------------------------------------

#' @title Base class for the single cell DuckDB connection
#'
#' @description
#' This is the base class for the single cell experiment DuckDB connection,
#' containing standard functions such as retrievel of tables, connection checks,
#' etc.
#'
#' @export
#'
#' @import data.table
single_cell_duckdb_base <- R6::R6Class(
  # class name
  classname = "single_cell_duckdb_base",
  # public functions, slots
  public = list(
    #' @description
    #' Initialises the Singe Cell DuckDB connection
    #'
    #' @param db_dir String. Path to where to store the db.
    #' @param db_name String. The name of the DB. Defaults to `"sc_duckdb.db"`.
    #'
    #' @return Returns the initialised class
    initialize = function(db_dir, db_name = "sc_duckdb.db") {
      # checks
      checkmate::assertDirectoryExists(db_dir)
      checkmate::qassert(db_name, "S1")

      # define the path
      private$db_path = file.path(db_dir, db_name)
    },

    ###########
    # Getters #
    ###########

    #' @description
    #' Returns the full observation table from the DuckDB
    #'
    #' @param indices Optional cell/obs positions. Not implemented yet.
    #'
    #' @return The observation table (if found) as a data.table
    get_obs_table = function(indices = NULL) {
      # checks
      checkmate::qassert(indices, c("0", "I+"))
      private$check_obs_exists()

      con <- private$connect_db()
      on.exit(
        {
          if (exists("con") && !is.null(con)) {
            tryCatch(DBI::dbDisconnect(con), error = function(e) invisible())
          }
        }
      )

      obs_dt <- data.table::setDT(DBI::dbGetQuery(con, 'SELECT * FROM obs'))

      return(obs_dt)
    },

    #' @description
    #' Returns the full var table from the DuckDB.
    #'
    #' @param indices Optional gene/feature positions. Not implemented yet.
    #'
    #' @return The observation table (if found) as a data.table
    get_vars_table = function(indices = NULL) {
      # checks
      checkmate::qassert(indices, c("0", "I+"))
      private$check_var_exists()

      con <- private$connect_db()
      on.exit(
        {
          if (exists("con") && !is.null(con)) {
            tryCatch(DBI::dbDisconnect(con), error = function(e) invisible())
          }
        }
      )

      var_dt <- data.table::setDT(DBI::dbGetQuery(con, 'SELECT * FROM var'))

      return(var_dt)
    }
  ),
  private = list(
    # Path to the DB
    db_path = NULL,
    # Helper to connect to the DB
    connect_db = function() {
      con <- dbConnect(duckdb::duckdb(), dbdir = private$db_path)
      return(con)
    },
    # Helper to check that obs table exists
    check_obs_exists = function() {
      con <- private$connect_db()
      res <- "obs" %in% DBI::dbGetQuery(con, "SHOW TABLES")[, "name"]
      on.exit(
        {
          if (exists("con") && !is.null(con)) {
            tryCatch(DBI::dbDisconnect(con), error = function(e) invisible())
          }
        }
      )

      if (!res) {
        error("Obs table was not found in the DB.")
      }
    },
    # Helper to check that var table exists
    check_var_exists = function() {
      con <- private$connect_db()
      res <- "var" %in% DBI::dbGetQuery(con, "SHOW TABLES")[, "name"]
      on.exit(
        {
          if (exists("con") && !is.null(con)) {
            tryCatch(DBI::dbDisconnect(con), error = function(e) invisible())
          }
        }
      )

      if (!res) {
        error("Obs table was not found in the DB.")
      }
    }
  )
)

### reader helper --------------------------------------------------------------

#' @title Class for storing single cell experimental data in DuckDB (nightly!)
#'
#' @description
#' This class wraps up the DB connection and methods to interact with the
#' observation/cell metadata and the feature/gene metadata.
#'
#' @export
#'
#' @import data.table
single_cell_duckdb_con <- R6::R6Class(
  # class name
  classname = "single_cell_duckdb_con",
  # public functions, slots
  inherit = single_cell_duckdb_base,
  public = list(
    ##############
    # Readers h5 #
    ##############

    #' @description
    #' This function populates the obs table from an h5 file (if found).
    #'
    #' @param h5_path String. Path to the h5 file from which to load in the
    #' observation table.
    #' @param .verbose Boolean. Controls verbosity of the function.
    #'
    #' @return Returns invisible self. As a side effect, it will load in the
    #' obs data from the h5ad file into the DuckDB.
    populate_obs_from_h5 = function(h5_path, .verbose = TRUE) {
      # checks
      checkmate::assertFileExists(h5_path)
      checkmate::qassert(.verbose, "B1")

      h5_content <- rhdf5::h5ls(
        h5_path
      ) %>%
        data.table::setDT()

      obs <- h5_content[
        group == "/obs" & otype == "H5I_DATASET",
        setNames((paste(group, name, sep = "/")), name)
      ]

      if (length(obs) == 0) {
        warning(
          "No obs data could be found in the h5 file. Nothing was loaded."
        )
        return(invisible(self))
      }

      con <- private$connect_db()

      # close everything no matter what happens
      on.exit(
        {
          if (exists("con") && !is.null(con)) {
            tryCatch(DBI::dbDisconnect(con), error = function(e) invisible())
          }
          tryCatch(rhdf5::h5closeAll(), error = function(e) invisible())
        },
        add = TRUE
      )

      names(obs) <- to_snake_case(obs)

      for (i in seq_along(obs)) {
        col_name <- names(obs)[i]
        col_path <- obs[[i]]
        col_data <- data.table(x = rhdf5::h5read(h5_path, col_path)) %>%
          `names<-`(col_name) %>%
          .[, cell_idx := .I] %>%
          .[, c("cell_idx", col_name), with = FALSE]

        if (i == 1) {
          if (.verbose) {
            message(
              sprintf(
                "Identified %i cells/observations in the h5 file",
                nrow(col_data)
              )
            )
          }
          DBI::dbWriteTable(
            con,
            "obs",
            col_data,
            overwrite = TRUE
          )
        } else {
          DBI::dbWriteTable(
            con,
            "temp_col",
            col_data,
            overwrite = TRUE
          )

          DBI::dbExecute(
            con,
            sprintf(
              'CREATE TABLE obs_new AS
              SELECT obs.*, temp_col.%s
              FROM obs
              JOIN temp_col ON obs.cell_idx = temp_col.cell_idx;
              DROP table obs;
              ALTER TABLE obs_new RENAME TO obs;
              DROP TABLE temp_col',
              col_name
            )
          )
        }
      }

      invisible(self)
    },

    #' @description
    #' This populates the vars table from an h5 file (if found.)
    #'
    #' @param h5_path String. Path to the h5 file from which to load in the
    #' observation table.
    #' @param .verbose Boolean. Controls verbosity of the function.
    #'
    #' @return Returns invisible self. As a side effect, it will load in the
    #' obs data from the h5ad file into the DuckDB.
    populate_vars_from_h5 = function(h5_path, .verbose = TRUE) {
      # checks
      checkmate::assertFileExists(h5_path)
      checkmate::qassert(.verbose, "B1")

      h5_content <- rhdf5::h5ls(
        h5_path
      ) %>%
        data.table::setDT()

      var <- h5_content[
        group == "/var" & otype == "H5I_DATASET",
        setNames((paste(group, name, sep = "/")), name)
      ]

      if (length(var) == 0) {
        warning(
          "No obs data could be found in the h5 file. Nothing was loaded."
        )
        return(invisible(self))
      }

      con <- private$connect_db()

      # close everything no matter what happens
      on.exit(
        {
          if (exists("con") && !is.null(con)) {
            tryCatch(DBI::dbDisconnect(con), error = function(e) invisible())
          }
          tryCatch(rhdf5::h5closeAll(), error = function(e) invisible())
        },
        add = TRUE
      )

      names(var) <- to_snake_case(var)

      for (i in seq_along(var)) {
        col_name <- names(var)[i]
        col_path <- var[[i]]
        col_data <- data.table(x = rhdf5::h5read(h5_path, col_path)) %>%
          `names<-`(col_name) %>%
          .[, gene_idx := .I] %>%
          .[, c("gene_idx", col_name), with = FALSE]

        if (i == 1) {
          if (.verbose) {
            message(
              sprintf(
                "Identified %i genes/features in the h5 file",
                nrow(col_data)
              )
            )
          }
          DBI::dbWriteTable(
            con,
            "var",
            col_data,
            overwrite = TRUE
          )
        } else {
          DBI::dbWriteTable(
            con,
            "temp_col",
            col_data,
            overwrite = TRUE
          )

          DBI::dbExecute(
            con,
            sprintf(
              'CREATE TABLE var_new AS
              SELECT var.*, temp_col.%s
              FROM var
              JOIN temp_col ON obs.cell_idx = temp_col.cell_idx;
              DROP table var;
              ALTER TABLE var_new RENAME TO var;
              DROP TABLE temp_col',
              col_name
            )
          )
        }
      }

      invisible(self)
    }
  )

  ############################
  # Private fields/functions #
  ############################
)
