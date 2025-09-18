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

    #' @description
    #' Returns a mapping between cell index and cell names/barcodes.
    #'
    #' @return A named numeric containing the cell index mapping.
    get_obs_index_map = function() {
      # checks
      private$check_obs_exists()

      con <- private$connect_db()
      on.exit(
        {
          if (exists("con") && !is.null(con)) {
            tryCatch(DBI::dbDisconnect(con), error = function(e) invisible())
          }
        }
      )

      data <- data.table::setDT(DBI::dbGetQuery(
        conn = con,
        statement = 'SELECT cell_idx, cell_id FROM obs'
      )) %>%
        `colnames<-`(c("index", "id"))

      return(setNames(data$index, data$id))
    },

    #' @description
    #' Returns a mapping between variable index and variable names.
    #'
    #' @return A named numeric containing the gene index mapping.
    get_var_index_map = function() {
      # checks
      private$check_var_exists()

      con <- private$connect_db()
      on.exit(
        {
          if (exists("con") && !is.null(con)) {
            tryCatch(DBI::dbDisconnect(con), error = function(e) invisible())
          }
        }
      )

      data <- data.table::setDT(DBI::dbGetQuery(
        conn = con,
        statement = 'SELECT gene_idx, gene_id FROM var'
      )) %>%
        `colnames<-`(c("index", "id"))

      return(setNames(data$index, data$id))
    },

    #' Filter the obs table and reset the cell idx
    #'
    #' @param filter_vec Boolean vector that will be used to filter the obs
    #' table.
    filter_obs_table = function(filter_vec) {
      # checks
      checkmate::qassert(filter_vec, "B1")
      private$check_obs_exists()

      filter_dt <- data.frame(cell_idx = which(filter_vec), keep = TRUE)

      # get the connection
      con <- private$connect_db()
      on.exit(
        {
          if (exists("con") && !is.null(con)) {
            tryCatch(DBI::dbDisconnect(con), error = function(e) invisible())
          }
        }
      )

      DBI::dbWriteTable(
        con,
        "filter_temp",
        filter_dt,
        overwrite = TRUE
      )

      DBI::dbExecute(
        con,
        "
        CREATE OR REPLACE TABLE obs AS
        SELECT 
          ROW_NUMBER() OVER() as cell_idx,
          * EXCLUDE (cell_idx)  -- exclude old cell_idx column
        FROM obs 
        WHERE ROWID IN (SELECT cell_idx FROM temp_filter);
        DROP TABLE filter_temp
        "
      )
    },

    #' Filter the var table and reset the gene idx
    #'
    #' @param filter_vec Boolean vector that will be used to filter the var
    #' table.
    filter_var_table = function(filter_vec) {
      # checks
      checkmate::qassert(filter_vec, "B1")
      private$check_var_exists()

      filter_dt <- data.frame(gene_idx = which(filter_vec), keep = TRUE)

      # get the connection
      con <- private$connect_db()
      on.exit(
        {
          if (exists("con") && !is.null(con)) {
            tryCatch(DBI::dbDisconnect(con), error = function(e) invisible())
          }
        }
      )

      DBI::dbWriteTable(
        con,
        "filter_temp",
        filter_dt,
        overwrite = TRUE
      )

      DBI::dbExecute(
        con,
        "
        CREATE OR REPLACE TABLE var AS
        SELECT 
          ROW_NUMBER() OVER() as gene_idx,
          * EXCLUDE (gene_idx)  -- exclude old gene_idx column
        FROM var 
        WHERE ROWID IN (SELECT gene_idx FROM temp_filter);
        DROP TABLE filter_temp
        "
      )
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
      checkmate::assertTRUE(nrow(new_data) == private$check_obs_row())

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
      checkmate::assertTRUE(nrow(new_data) == private$check_var_row())

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
      on.exit(
        {
          if (exists("con") && !is.null(con)) {
            tryCatch(DBI::dbDisconnect(con), error = function(e) invisible())
          }
        }
      )
      res <- "obs" %in% DBI::dbGetQuery(con, "SHOW TABLES")[, "name"]

      if (!res) {
        error("Obs table was not found in the DB.")
      }
    },
    # Helper to check that var table exists
    check_var_exists = function() {
      con <- private$connect_db()
      on.exit(
        {
          if (exists("con") && !is.null(con)) {
            tryCatch(DBI::dbDisconnect(con), error = function(e) invisible())
          }
        }
      )

      res <- "var" %in% DBI::dbGetQuery(con, "SHOW TABLES")[, "name"]

      if (!res) {
        error("Obs table was not found in the DB.")
      }
    },
    # Check number of rows in the obs table
    check_obs_row = function() {
      private$check_obs_exists()
      con <- private$connect_db()
      on.exit(
        {
          if (exists("con") && !is.null(con)) {
            tryCatch(DBI::dbDisconnect(con), error = function(e) invisible())
          }
        }
      )

      DBI::dbGetQuery(con, "SELECT COUNT(*) FROM obs")
    },
    # Check number of rows in the var table
    check_var_row = function() {
      private$check_obs_exists()
      con <- private$connect_db()
      on.exit(
        {
          if (exists("con") && !is.null(con)) {
            tryCatch(DBI::dbDisconnect(con), error = function(e) invisible())
          }
        }
      )

      DBI::dbGetQuery(con, "SELECT COUNT(*) FROM var")
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
    #' @param filter Optional boolean. If provided, only rows with `TRUE` will
    #' be read in. The length needs to be same as nrow of obs in the h5 object.
    #'
    #' @return Returns invisible self. As a side effect, it will load in the
    #' obs data from the h5ad file into the DuckDB.
    populate_obs_from_h5 = function(h5_path, filter = NULL) {
      # checks
      checkmate::assertFileExists(h5_path)
      checkmate::qassert(filter, c("B+", "0"))

      h5_content <- rhdf5::h5ls(
        h5_path
      ) %>%
        data.table::setDT()

      obs <- h5_content[
        group == "/obs" & otype == "H5I_DATASET",
        setNames((paste(group, name, sep = "/")), name)
      ]

      if (length(obs) == 0) {
        error(
          "No obs data could be found in the h5 file. Nothing was loaded."
        )
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
          `names<-`(col_name)
        if (!is.null(filter)) {
          col_data <- col_data[filter]
        }
        col_data[, cell_idx := .I]
        setcolorder(col_data, c("cell_idx", col_name))

        if (i == 1) {
          colnames(col_data) <- c("cell_idx", "cell_id")
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
    #' @param filter Optional boolean. If provided, only rows with `TRUE` will
    #' be read in. The length needs to be same as nrow of var in the h5 object.
    #'
    #' @return Returns invisible self. As a side effect, it will load in the
    #' obs data from the h5ad file into the DuckDB.
    populate_vars_from_h5 = function(h5_path, filter = NULL) {
      # checks
      checkmate::assertFileExists(h5_path)
      checkmate::qassert(filter, c("B+", "0"))

      h5_content <- rhdf5::h5ls(
        h5_path
      ) %>%
        data.table::setDT()

      var <- h5_content[
        group == "/var" & otype == "H5I_DATASET",
        setNames((paste(group, name, sep = "/")), name)
      ]

      if (length(var) == 0) {
        error(
          "No var data could be found in the h5 file. Nothing was loaded."
        )
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
          `names<-`(col_name)
        if (!is.null(filter)) {
          col_data <- col_data[filter]
        }

        col_data[, gene_idx := .I]
        setcolorder(col_data, c("gene_idx", col_name))

        if (i == 1) {
          colnames(col_data) <- c("gene_idx", "gene_id")
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
    },

    ###############
    # Readers mtx #
    ###############

    #' Function to populate the obs table from plain text files
    #'
    #' @param f_path File path to the plain text file.
    #' @param filter Optional boolean. If provided, only rows with `TRUE` will
    #' be read in. The length needs to be same as nrow.
    #'
    #' @returns Invisible self and populates the internal obs table.
    populate_obs_from_plain_text = function(f_path, filter = NULL) {
      # checks
      checkmate::assertFileExists(f_path)
      checkmate::qassert(filter, c("B+", "0"))

      con <- private$connect_db()
      on.exit(
        {
          if (exists("con") && !is.null(con)) {
            tryCatch(DBI::dbDisconnect(con), error = function(e) invisible())
          }
        }
      )

      header_query <- "SELECT * FROM read_csv(?) LIMIT 0"

      header_df <- DBI::dbGetQuery(
        conn = con,
        statement = header_query,
        params = list(f_path)
      )

      original_col_names <- colnames(header_df)
      sanitised_cols_names <- to_snake_case(original_col_names)
      sanitised_cols_names[1] <- "cell_id"

      col_mapping <- sprintf(
        '"%s" AS %s',
        original_col_names,
        sanitised_cols_names
      ) %>%
        paste(collapse = ",\n    ")

      query <- if (!is.null(filter)) {
        rows_to_keep <- which(filter)
        row_placeholders <- paste(
          rep("?", length(rows_to_keep)),
          collapse = ","
        )
        sprintf(
          "
          CREATE OR REPLACE TABLE obs AS
          SELECT
            ROW_NUMBER() OVER() as cell_idx,
            %s
          FROM (
            SELECT *, ROW_NUMBER() OVER() as original_row_num
            FROM read_csv(?)
          ) t
          WHERE t.original_row_num IN (%s)",
          col_mapping,
          row_placeholders
        )
      } else {
        rows_to_keep <- NULL
        sprintf(
          "
          CREATE OR REPLACE TABLE obs AS
          SELECT 
            ROW_NUMBER() OVER() as cell_idx,
            %s
          FROM (
            SELECT * FROM read_csv(?)
          )",
          col_mapping
        )
      }

      DBI::dbExecute(
        conn = con,
        statement = query,
        params = c(list(f_path), as.list(rows_to_keep))
      )

      invisible(self)
    },

    #' Function to populate the var table from plain text files
    #'
    #' @param f_path File path to the plain text file.
    #' @param filter Optional boolean. If provided, only rows with `TRUE` will
    #' be read in. The length needs to be same as nrow.
    #'
    #' @returns Invisible self and populates the internal var table.
    populate_var_from_plain_text = function(f_path, filter = NULL) {
      # checks
      checkmate::assertFileExists(f_path)
      checkmate::qassert(filter, c("B+", "0"))

      con <- private$connect_db()
      on.exit(
        {
          if (exists("con") && !is.null(con)) {
            tryCatch(DBI::dbDisconnect(con), error = function(e) invisible())
          }
        }
      )

      header_query <- "SELECT * FROM read_csv(?) LIMIT 0"

      header_df <- DBI::dbGetQuery(
        conn = con,
        statement = header_query,
        params = list(f_path)
      )

      original_col_names <- colnames(header_df)
      sanitised_cols_names <- to_snake_case(original_col_names)
      sanitised_cols_names[1] <- "gene_id"

      col_mapping <- sprintf(
        '"%s" AS %s',
        original_col_names,
        sanitised_cols_names
      ) %>%
        paste(collapse = ",\n    ")

      query <- if (!is.null(filter)) {
        rows_to_keep <- which(filter)
        row_placeholders <- paste(
          rep("?", length(rows_to_keep)),
          collapse = ","
        )
        sprintf(
          "
          CREATE OR REPLACE TABLE var AS
          SELECT
            ROW_NUMBER() OVER() AS gene_idx,
            %s
          FROM (
            SELECT *, ROW_NUMBER() OVER() as original_row_num
            FROM read_csv(?)
          ) t
          WHERE t.original_row_num IN (%s)",
          col_mapping,
          row_placeholders
        )
      } else {
        rows_to_keep <- NULL
        sprintf(
          "
          CREATE OR REPLACE TABLE var AS
          SELECT 
            ROW_NUMBER() OVER() AS gene_idx,
            %s
          FROM (
            SELECT * FROM read_csv(?)
          )",
          col_mapping
        )
      }

      DBI::dbExecute(
        conn = con,
        statement = query,
        params = c(list(f_path), as.list(rows_to_keep))
      )

      invisible(self)
    }
  )

  ############################
  # Private fields/functions #
  ############################
)
