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
#'   \item{cache}{List with cached data. Future feature, nothing implemented
#'   yet.}
#'   \item{dims}{Dimensions of the original data.}
#'   \item{index_maps}{A list of two named numerics that contains cell id to
#'   cell idx and gene id to gene idx info.}
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
    cache = S7::class_list,
    dims = S7::class_integer,
    index_maps = S7::class_list
  ),
  constructor = function(dir_data) {
    nightly_feature()
    # checks
    checkmate::assertDirectoryExists(dir_data)

    # generate the Rust pointer
    count_connection <- SingeCellCountData$new(
      f_path_cells = file.path(dir_data, "counts_cells.bin"),
      f_path_genes = file.path(dir_data, "counts_genes.bin")
    )

    # generate the DuckDB connector
    db_connection <- single_cell_duckdb_con$new(
      db_dir = dir_data
    )

    S7::new_object(
      S7::S7_object(),
      db_connection = db_connection,
      count_connection = count_connection,
      dir_data = dir_data,
      cache = list(),
      dims = c(0L, 0L),
      index_maps = list()
    )
  }
)

### getters --------------------------------------------------------------------

#### env getters ---------------------------------------------------------------

#' Getter for the single cell DuckDB connection
#'
#' @param object `single_cell_exp` class.
#'
#' @return The DuckDB connector
#'
#' @export
get_sc_duckdb <- S7::new_generic(
  name = "get_sc_duckdb",
  dispatch_args = "object",
  fun = function(
    object
  ) {
    S7::S7_dispatch()
  }
)

#' @method get_sc_duckdb single_cell_exp
#'
#' @export
S7::method(get_sc_duckdb, single_cell_exp) <- function(object) {
  # checks
  checkmate::assertClass(object, "bixverse::single_cell_exp")

  return(S7::prop(object, "db_connection"))
}

#' Getter for the single cell Rust pointer
#'
#' @param object `single_cell_exp` class.
#'
#' @return The Rust structure
#'
#' @export
get_sc_rust_ptr <- S7::new_generic(
  name = "get_sc_rust_ptr",
  dispatch_args = "object",
  fun = function(
    object
  ) {
    S7::S7_dispatch()
  }
)

#' @method get_sc_rust_ptr single_cell_exp
#'
#' @export
S7::method(get_sc_rust_ptr, single_cell_exp) <- function(object) {
  # checks
  checkmate::assertClass(object, "bixverse::single_cell_exp")

  return(S7::prop(object, "count_connection"))
}

#### duckdb getters ------------------------------------------------------------

#' Getter the obs table
#'
#' @param object `single_cell_exp` class.
#' @param indices Optional integer vector. The integer positions of the cells
#' to return.
#' @param cols Optional string vector. The columns from the obs table to return.
#'
#' @return The obs table
#'
#' @export
get_sc_obs <- S7::new_generic(
  name = "get_sc_obs",
  dispatch_args = "object",
  fun = function(
    object,
    indices = NULL,
    cols = NULL
  ) {
    S7::S7_dispatch()
  }
)

#' @method get_sc_obs single_cell_exp
#'
#' @export
S7::method(get_sc_obs, single_cell_exp) <- function(
  object,
  indices = NULL,
  cols = NULL
) {
  # checks
  checkmate::assertClass(object, "bixverse::single_cell_exp")
  checkmate::qassert(indices, c("0", "I+"))

  duckdb_con <- get_sc_duckdb(object)

  obs_table <- duckdb_con$get_obs_table(indices = indices, cols = cols)

  return(obs_table)
}

#' Getter the var table
#'
#' @param object `single_cell_exp` class.
#' @param indices Optional integer vector. The integer positions of the genes
#' to return.
#' @param cols Optional string vector. The columns from the var table to return.
#'
#' @return The vars table
#'
#' @export
get_sc_var <- S7::new_generic(
  name = "get_sc_var",
  dispatch_args = "object",
  fun = function(
    object,
    indices = NULL,
    cols = NULL
  ) {
    S7::S7_dispatch()
  }
)

#' @method get_sc_var single_cell_exp
#'
#' @export
S7::method(get_sc_var, single_cell_exp) <- function(
  object,
  indices = NULL,
  cols = NULL
) {
  # checks
  checkmate::assertClass(object, "bixverse::single_cell_exp")
  checkmate::qassert(indices, c("0", "I+"))

  duckdb_con <- get_sc_duckdb(object)

  var_table <- duckdb_con$get_vars_table(indices = indices, cols = cols)

  return(var_table)
}

#' @method `[[` single_cell_exp
#'
#' @export
S7::method(`[[`, single_cell_exp) <- function(x, i, ...) {
  if (checkmate::qtest(i, "S+")) {
    get_sc_obs(x, cols = i)
  } else if (checkmate::qtest(i, "I+")) {
    get_sc_obs(x, indices = i)
  } else {
    stop("Invalid type")
  }
}

#### count getters -------------------------------------------------------------

#' Getter the counts
#'
#' @param object `single_cell_exp` class.
#' @param assay String. Which slot to return. One of `c("raw", "norm")`.
#' Defaults to `"raw"`.
#' @param return_format String. One of `c("cell", "gene")`. Return data in
#' cell-centric compressed format (CSR) or gene-centric compressed format (CSC).
#' Defaults to `"cell"`.
#' @param cell_indices Optional cell indices.
#' @param gene_indices Optional gene indices.
#'
#' @return The obs table
#'
#' @export
get_sc_counts <- S7::new_generic(
  name = "get_sc_obs",
  dispatch_args = "object",
  fun = function(
    object,
    assay = c("raw", "norm"),
    return_format = c("cell", "gene"),
    cell_indices = NULL,
    gene_indices = NULL
  ) {
    S7::S7_dispatch()
  }
)

#' @method get_sc_counts single_cell_exp
#'
#' @export
S7::method(get_sc_counts, single_cell_exp) <- function(
  object,
  assay = c("raw", "norm"),
  return_format = c("cell", "gene"),
  cell_indices = NULL,
  gene_indices = NULL
) {
  assay <- match.arg(assay)
  return_format <- match.arg(return_format)

  # checks
  checkmate::assertClass(object, "bixverse::single_cell_exp")
  checkmate::assertChoice(assay, c("raw", "norm"))
  checkmate::assertChoice(return_format, c("raw", "norm"))
  checkmate::qassert(cell_indices, c("0", "I+"))
  checkmate::qassert(gene_indices, c("0", "I+"))

  rust_con <- get_sc_rust_ptr(object)

  return(obs_table)
}

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
    #' @param indices Optional cell/obs indices.
    #' @param cols Optional column names to return.
    #'
    #' @return The observation table (if found) as a data.table with optionally
    #' selected indices and/or columns.
    get_obs_table = function(indices = NULL, cols = NULL) {
      # checks
      checkmate::qassert(indices, c("0", "I+"))
      checkmate::qassert(cols, c("0", "S+"))
      private$check_obs_exists()

      con <- private$connect_db()
      on.exit(
        {
          if (exists("con") && !is.null(con)) {
            tryCatch(DBI::dbDisconnect(con), error = function(e) invisible())
          }
        }
      )

      col_part <- if (is.null(cols)) {
        "*"
      } else {
        paste(cols, collapse = ", ")
      }

      sql_query <- if (is.null(indices)) {
        sprintf("SELECT %s FROM obs", col_part)
      } else {
        placeholders <- paste(rep("?", length(indices)), collapse = ", ")
        sprintf(
          "SELECT %s FROM obs WHERE cell_idx IN (%s)",
          col_part,
          placeholders
        )
      }

      obs_dt <- data.table::setDT(DBI::dbGetQuery(
        conn = con,
        statement = sql_query,
        params = as.list(indices)
      ))

      return(obs_dt)
    },

    #' @description
    #' Returns the full var table from the DuckDB.
    #'
    #' @param indices Optional gene/var indices.
    #' @param cols Optional column names to return.
    #'
    #' @return The var table (if found) as a data.table with optionally
    #' selected indices and/or columns.
    get_vars_table = function(indices = NULL, cols = NULL) {
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

      col_part <- if (is.null(cols)) {
        "*"
      } else {
        paste(cols, collapse = ", ")
      }

      sql_query <- if (is.null(indices)) {
        sprintf("SELECT %s FROM var", col_part)
      } else {
        placeholders <- paste(rep("?", length(indices)), collapse = ", ")
        sprintf(
          "SELECT %s FROM var WHERE gene_idx IN (%s)",
          col_part,
          placeholders
        )
      }

      var_dt <- data.table::setDT(DBI::dbGetQuery(
        conn = con,
        statement = sql_query,
        params = as.list(indices)
      ))

      return(var_dt)
    },

    ###########
    # Setters #
    ###########

    #' @description
    #' Add new data to the obs table in the DuckDB
    #'
    #' @param new_data A data.table with additional new columns. The order needs
    #' to be the same as the original in the obs table.
    #'
    #' @return Invisible self while adding the new columns to the obs table
    #' in the DuckDB.
    add_data_obs = function(new_data) {
      # checks
      checkmate::assertDataTable(new_data)
      private$check_obs_exists()

      # add the data
      con <- private$connect_db()
      on.exit(
        {
          if (exists("con") && !is.null(con)) {
            tryCatch(DBI::dbDisconnect(con), error = function(e) invisible())
          }
        }
      )

      new_data[, cell_idx := .I]

      DBI::dbWriteTable(
        con,
        "new_data",
        new_data,
        overwrite = TRUE
      )

      DBI::dbExecute(
        con,
        "
        CREATE OR REPLACE TABLE obs AS
        SELECT
          original.*,
          new_data.* EXCLUDE(cell_idx)
        FROM obs AS original
        JOIN new_data ON original.cell_idx = new_data.cell_idx;
        DROP TABLE new_data
        "
      )

      invisible(self)
    },

    #' @description
    #' Add the information which genes pass threshold to the DuckDB.
    #'
    #' @param new_data A data.table with additional new columns. The order needs
    #' to be the same as the original in the var table.
    #'
    #' @return Invisible self while adding the new columns to the var table
    #' in the DuckDB.
    add_data_var = function(new_data) {
      # checks
      checkmate::assertDataTable(new_data)
      private$check_obs_exists()

      # add the data
      con <- private$connect_db()
      on.exit(
        {
          if (exists("con") && !is.null(con)) {
            tryCatch(DBI::dbDisconnect(con), error = function(e) invisible())
          }
        }
      )

      new_data[, gene_idx := .I]

      DBI::dbWriteTable(
        con,
        "new_data",
        new_data,
        overwrite = TRUE
      )

      DBI::dbExecute(
        con,
        "
        CREATE OR REPLACE TABLE var AS
        SELECT
          original.*,
          new_data.* EXCLUDE(gene_idx)
        FROM var AS original
        JOIN new_data ON original.gene_idx = new_data.gene_idx;
        DROP TABLE new_data
        "
      )

      invisible(self)
    }
  ),
  private = list(
    # Path to the DB
    db_path = NULL,
    # Helper to connect to the DB
    connect_db = function() {
      con <- DBI::dbConnect(duckdb::duckdb(), dbdir = private$db_path)
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
