# r6 ---------------------------------------------------------------------------

## duckdb connector ------------------------------------------------------------

### base class -----------------------------------------------------------------

#' Base class for the single cell DuckDB connection
#'
#' @description
#' This is the base class for the single cell experiment DuckDB connection,
#' containing standard functions such as retrievel of tables, connection checks,
#' etc.
#'
#' @export
#'
#' @import data.table
#'
#' @keywords internal
SingleCellDuckDBBase <- R6::R6Class(
  # class name
  classname = "SingleCellDuckDBBase",
  # public functions, slots
  public = list(
    #' @description
    #' Initialises the Singe Cell DuckDB connection
    #'
    #' @param db_dir String. Path to where to store the db.
    #' @param db_name String. The name of the DB. Defaults to `"sc_duckdb.db"`.
    #'
    #' @returns Returns the initialised class
    initialize = function(db_dir, db_name = "sc_duckdb.db") {
      # checks
      checkmate::assertDirectoryExists(db_dir)
      checkmate::qassert(db_name, "S1")

      # define the path
      private$db_path <- file.path(db_dir, db_name)
    },

    ###########
    # Getters #
    ###########

    #' @description
    #' Returns the full observation table from the DuckDB
    #'
    #' @param indices Optional cell/obs indices.
    #' @param cols Optional column names to return.
    #' @param filtered Boolean. Whether to return all cells or filtered to to_keep cells.
    #' @param where_equal Optional named list. Only rows where each named column
    #' equals its value are returned; the value is cast to the column's type.
    #' Pushed into the query, so a few rows of a large table do not have to come
    #' through R.
    #'
    #' @returns The observation table (if found) as a data.table with optionally
    #' selected indices and/or columns.
    get_obs_table = function(
      indices = NULL,
      cols = NULL,
      filtered = FALSE,
      where_equal = NULL
    ) {
      # checks
      checkmate::qassert(indices, c("0", "I+"))
      checkmate::qassert(cols, c("0", "S+"))
      checkmate::qassert(filtered, "B1")
      checkmate::assertList(
        where_equal,
        types = "atomic",
        names = "unique",
        null.ok = TRUE
      )
      private$check_obs_exists()

      con <- private$connect_db()
      on.exit(
        {
          if (exists("con") && !is.null(con)) {
            tryCatch(DBI::dbDisconnect(con), error = function(e) invisible())
          }
        }
      )

      if (!is.null(where_equal)) {
        obs_types <- DBI::dbGetQuery(con, "DESCRIBE obs")
        checkmate::assertNames(
          names(where_equal),
          subset.of = obs_types$column_name
        )
      }

      col_part <- if (is.null(cols)) {
        "*"
      } else {
        paste(sprintf('"%s"', cols), collapse = ", ")
      }

      where_clauses <- character()
      params <- list()

      if (filtered) {
        where_clauses <- c(where_clauses, "to_keep = TRUE")
      }

      if (!is.null(indices)) {
        placeholders <- paste(rep("?", length(indices)), collapse = ", ")
        where_clauses <- c(
          where_clauses,
          sprintf("cell_idx IN (%s)", placeholders)
        )
        params <- as.list(indices)
      }

      # casting the value to the column type, not the column to a string, keeps
      # the R semantics of `col == "1"` for a numeric column, whose string form
      # would be "1.0"
      for (col in names(where_equal)) {
        col_type <- obs_types$column_type[obs_types$column_name == col]
        where_clauses <- c(
          where_clauses,
          sprintf('"%s" = TRY_CAST(? AS %s)', col, col_type)
        )
        params <- c(params, list(as.character(where_equal[[col]])))
      }

      where_part <- if (length(where_clauses) > 0) {
        paste("WHERE", paste(where_clauses, collapse = " AND "))
      } else {
        ""
      }

      sql_query <- sprintf("SELECT %s FROM obs %s", col_part, where_part)

      obs_dt <- data.table::setDT(DBI::dbGetQuery(
        conn = con,
        statement = sql_query,
        params = params
      ))

      return(obs_dt)
    },

    #' @description
    #' Returns the var table from the DuckDB.
    #'
    #' @param indices Optional gene/var indices.
    #' @param cols Optional column names to return.
    #'
    #' @returns The var table (if found) as a data.table with optionally
    #' selected indices and/or columns.
    get_vars_table = function(indices = NULL, cols = NULL) {
      # checks
      checkmate::qassert(indices, c("0", "I+"))
      checkmate::qassert(cols, c("S+", "0"))
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
        paste(sprintf('"%s"', cols), collapse = ", ")
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
    #' Returns the var table specific for the ADTs from the DuckDB.
    #'
    #' @param indices Optional ADT/var indices.
    #' @param cols Optional column names to return.
    #'
    #' @returns The var adt table (if found) as a data.table with optionally
    #' selected indices and/or columns.
    get_vars_adt_table = function(indices = NULL, cols = NULL) {
      # checks
      checkmate::qassert(indices, c("0", "I+"))
      checkmate::qassert(cols, c("S+", "0"))
      private$check_var_adt_exists()

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
        paste(sprintf('"%s"', cols), collapse = ", ")
      }

      sql_query <- if (is.null(indices)) {
        sprintf("SELECT %s FROM var_adt", col_part)
      } else {
        placeholders <- paste(rep("?", length(indices)), collapse = ", ")
        sprintf(
          "SELECT %s FROM var_adt WHERE gene_idx IN (%s)",
          col_part,
          placeholders
        )
      }

      var_adt_dt <- data.table::setDT(DBI::dbGetQuery(
        conn = con,
        statement = sql_query,
        params = as.list(indices)
      ))

      return(var_adt_dt)
    },

    #' @description
    #' Returns a mapping between cell index and cell names/barcodes.
    #'
    #' @returns A named numeric containing the cell index mapping.
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
    #' @returns A named numeric containing the gene index mapping.
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

    #' @description
    #' Returns the indices of the cells that have to_keep = TRUE in the DB.
    #'
    #' @returns The index positions (1-index) of the cells to keep.
    get_cells_to_keep = function() {
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
        statement = 'SELECT cell_idx FROM obs WHERE to_keep'
      ))

      return(as.integer(data$cell_idx))
    },

    #' @description
    #' Returns the available column names in the obs table.
    #'
    #' @returns A character vector of column names.
    get_obs_cols = function() {
      private$check_obs_exists()

      con <- private$connect_db()
      on.exit({
        if (exists("con") && !is.null(con)) {
          tryCatch(DBI::dbDisconnect(con), error = function(e) invisible())
        }
      })

      DBI::dbGetQuery(
        con,
        "SELECT column_name FROM information_schema.columns WHERE table_name = 'obs'"
      )$column_name
    },

    #####################
    # Setters / filters #
    #####################

    #' Filter the obs table and reset the cell idx
    #'
    #' @param cell_idx_to_keep Integer vector with the cell indices to keep.
    #' Needs to be 1-indexed!
    #'
    #' @returns Invisible self after updating the to_keep column in the DuckDB.
    set_cells_to_keep = function(cell_idx_to_keep) {
      # checks
      checkmate::qassert(cell_idx_to_keep, "I+")
      private$check_obs_exists()

      filter_dt <- data.frame(cell_idx = cell_idx_to_keep)

      # get the connection
      con <- private$connect_db()
      on.exit(
        {
          if (exists("con") && !is.null(con)) {
            tryCatch(DBI::dbDisconnect(con), error = function(e) invisible())
          }
        }
      )

      DBI::dbWriteTable(con, "temp_filter", filter_dt, overwrite = TRUE)

      DBI::dbExecute(
        con,
        "
        UPDATE obs
        SET to_keep = cell_idx IN (SELECT cell_idx FROM temp_filter);
        DROP TABLE temp_filter
        "
      )

      invisible(self)
    },

    #' Restore every cell in the obs table
    #'
    #' @description
    #' Flips `to_keep` back on for every row. Filtering never deletes rows, it
    #' only flips this flag, so nothing was lost and nothing needs rebuilding.
    #' A single statement rather than a round trip through
    #' `set_cells_to_keep()`, which would write a temp table the size of the
    #' whole obs table just to switch every flag back on.
    #'
    #' @returns Invisible self after resetting the to_keep column.
    reset_cells_to_keep = function() {
      # checks
      private$check_obs_exists()

      # get the connection
      con <- private$connect_db()
      on.exit(
        {
          if (exists("con") && !is.null(con)) {
            tryCatch(DBI::dbDisconnect(con), error = function(e) invisible())
          }
        }
      )

      DBI::dbExecute(con, "UPDATE obs SET to_keep = TRUE")

      invisible(self)
    },

    #' Filter the var table and reset the gene idx
    #'
    #' @param filter_vec Boolean vector that will be used to filter the var
    #' table.
    filter_var_table = function(filter_vec) {
      # TODO needs also updating...

      # checks
      checkmate::qassert(filter_vec, "B+")
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
        "temp_filter",
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
        WHERE gene_idx IN (SELECT gene_idx FROM temp_filter);
        DROP TABLE temp_filter
        "
      )
    },

    #' @description
    #' Add new data to the obs table in the DuckDB
    #'
    #' @param new_data A data.table with additional new columns. The order needs
    #' to be the same as the original in the obs table.
    #'
    #' @returns Invisible self while adding the new columns to the obs table
    #' in the DuckDB.
    add_data_obs = function(new_data) {
      checkmate::assertDataTable(new_data)
      private$check_obs_exists()
      checkmate::assertTRUE(nrow(new_data) == private$check_obs_row())

      con <- private$connect_db()
      on.exit({
        if (exists("con") && !is.null(con)) {
          tryCatch(DBI::dbDisconnect(con), error = function(e) invisible())
        }
      })

      new_data[, cell_idx := .I]
      DBI::dbWriteTable(con, "new_data", new_data, overwrite = TRUE)

      private$write_columns(con, table = "obs", key = "cell_idx")

      invisible(self)
    },

    #' @description
    #' Left join new data to the obs table in the DuckDB by cell_idx
    #'
    #' @param new_data A data.table with a cell_idx column to join on.
    #'
    #' @returns Invisible self while left joining the new data to the obs table
    #' in the DuckDB.
    join_data_obs = function(new_data) {
      checkmate::assertDataTable(new_data)
      private$check_obs_exists()
      checkmate::assertTRUE("cell_idx" %in% names(new_data))

      con <- private$connect_db()
      on.exit({
        if (exists("con") && !is.null(con)) {
          tryCatch(DBI::dbDisconnect(con), error = function(e) invisible())
        }
      })

      DBI::dbWriteTable(con, "new_data", new_data, overwrite = TRUE)

      private$write_columns(con, table = "obs", key = "cell_idx")

      invisible(self)
    },

    #' @description
    #' Independent of the loader, set the to_keep column to `TRUE` initially
    set_to_keep_column = function() {
      private$check_obs_exists()
      con <- private$connect_db()
      on.exit({
        if (exists("con") && !is.null(con)) {
          tryCatch(DBI::dbDisconnect(con), error = function(e) invisible())
        }
      })

      # Ensure to_keep exists and is TRUE by default
      DBI::dbExecute(
        conn = con,
        statement = "
        ALTER TABLE obs ADD COLUMN IF NOT EXISTS to_keep BOOLEAN;
        UPDATE obs SET to_keep = TRUE WHERE to_keep IS NULL;
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
    #' @returns Invisible self while adding the new columns to the var table
    #' in the DuckDB.
    add_data_var = function(new_data) {
      # checks
      checkmate::assertDataTable(new_data)
      private$check_var_exists()
      checkmate::assertTRUE(nrow(new_data) == private$check_var_row())

      con <- private$connect_db()
      on.exit({
        if (exists("con") && !is.null(con)) {
          tryCatch(DBI::dbDisconnect(con), error = function(e) invisible())
        }
      })

      new_data[, gene_idx := .I]
      DBI::dbWriteTable(con, "new_data", new_data, overwrite = TRUE)

      private$write_columns(con, table = "var", key = "gene_idx")

      invisible(self)
    },

    ##########
    # Rename #
    ##########

    #' @description
    #' Rename a column in the obs or var table.
    #'
    #' @param table String. Either "obs" or "var".
    #' @param old String. The current column name.
    #' @param new String. The desired new column name.
    #'
    #' @returns Invisible self.
    rename_column = function(table = c("obs", "var"), old, new) {
      table <- match.arg(table)
      checkmate::qassert(old, "S1")
      checkmate::qassert(new, "S1")

      if (table == "obs") {
        private$check_obs_exists()
      } else {
        private$check_var_exists()
      }

      con <- private$connect_db()
      on.exit({
        if (exists("con") && !is.null(con)) {
          tryCatch(DBI::dbDisconnect(con), error = function(e) invisible())
        }
      })

      existing <- DBI::dbGetQuery(
        con,
        sprintf(
          "SELECT column_name FROM information_schema.columns WHERE table_name = '%s'",
          table
        )
      )$column_name

      if (!old %in% existing) {
        warning(sprintf(
          "Column '%s' not found in table '%s'. No renaming done.",
          old,
          table
        ))
        return(invisible(self))
      }

      DBI::dbExecute(
        con,
        sprintf('ALTER TABLE %s RENAME COLUMN "%s" TO "%s"', table, old, new)
      )

      invisible(self)
    },

    ########
    # Drop #
    ########

    #' @description
    #' Drop columns from the obs or var table.
    #'
    #' Refuses to drop the protected identifier and bookkeeping columns
    #' (`cell_idx`, `cell_id`, `to_keep` for obs; `gene_idx`, `gene_id` for
    #' var). Columns that do not exist trigger a warning and are skipped.
    #'
    #' @param table String. Either "obs" or "var".
    #' @param cols Character vector. The column names to drop.
    #'
    #' @returns Invisible self.
    drop_columns = function(table = c("obs", "var"), cols) {
      table <- match.arg(table)
      checkmate::qassert(cols, "S+")

      if (table == "obs") {
        private$check_obs_exists()
        protected <- c("cell_idx", "cell_id", "to_keep")
      } else {
        private$check_var_exists()
        protected <- c("gene_idx", "gene_id")
      }

      con <- private$connect_db()
      on.exit({
        if (exists("con") && !is.null(con)) {
          tryCatch(DBI::dbDisconnect(con), error = function(e) invisible())
        }
      })

      existing <- DBI::dbGetQuery(
        con,
        sprintf(
          "SELECT column_name FROM information_schema.columns WHERE table_name = '%s'",
          table
        )
      )$column_name

      to_drop_protected <- intersect(cols, protected)
      if (length(to_drop_protected) > 0L) {
        warning(sprintf(
          "Refusing to drop protected column(s) in '%s': %s",
          table,
          paste(to_drop_protected, collapse = ", ")
        ))
      }

      to_drop_missing <- setdiff(cols, c(existing, protected))
      if (length(to_drop_missing) > 0L) {
        warning(sprintf(
          "Column(s) not found in '%s' (skipping): %s",
          table,
          paste(to_drop_missing, collapse = ", ")
        ))
      }

      to_drop <- setdiff(intersect(cols, existing), protected)
      if (length(to_drop) == 0L) {
        return(invisible(self))
      }

      for (col in to_drop) {
        DBI::dbExecute(
          con,
          sprintf('ALTER TABLE %s DROP COLUMN "%s"', table, col)
        )
      }

      invisible(self)
    }
  ),
  private = list(
    # Path to the DB
    db_path = NULL,
    # Whether the obs, var and var_adt tables have been seen to exist
    obs_exists = FALSE,
    var_exists = FALSE,
    var_adt_exists = FALSE,
    # Write the columns of the `new_data` table into `table`, matched on `key`.
    #
    # Adds each column in place and fills it with one UPDATE, rather than
    # rebuilding the whole table around a join. Columns that already exist are
    # replaced, rows without a match in `new_data` are NULL, and `new_data` is
    # dropped afterwards.
    write_columns = function(con, table, key) {
      types <- DBI::dbGetQuery(con, "DESCRIBE new_data")
      types <- types[types$column_name != key, ]
      quoted <- sprintf('"%s"', types$column_name)

      for (i in seq_len(nrow(types))) {
        DBI::dbExecute(
          con,
          sprintf("ALTER TABLE %s DROP COLUMN IF EXISTS %s", table, quoted[i])
        )
        DBI::dbExecute(
          con,
          sprintf(
            "ALTER TABLE %s ADD COLUMN %s %s",
            table,
            quoted[i],
            types$column_type[i]
          )
        )
      }

      if (nrow(types) > 0L) {
        DBI::dbExecute(
          con,
          sprintf(
            "UPDATE %s SET %s FROM new_data WHERE %s.%s = new_data.%s",
            table,
            paste(
              sprintf("%s = new_data.%s", quoted, quoted),
              collapse = ", "
            ),
            table,
            key,
            key
          )
        )
      }

      DBI::dbExecute(con, "DROP TABLE new_data")
    },
    # Helper to connect to the DB
    connect_db = function() {
      # `shared_home` is not about the database file, it is where duckdb keeps
      # downloaded extensions and secrets. We use neither, so there is nothing
      # to gain from writing into the user's ~/.duckdb. Passing it either way
      # also marks the choice as made, which is what stops duckdb announcing
      # its storage location on every single connection.
      con <- DBI::dbConnect(
        duckdb::duckdb(shared_home = FALSE),
        dbdir = private$db_path
      )
      return(con)
    },
    # Helper to check that obs table exists
    check_obs_exists = function() {
      # tables are replaced but never dropped, so once seen this holds and
      # the check does not have to open a connection on every call
      if (isTRUE(private$obs_exists)) {
        return(invisible(TRUE))
      }
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
        stop("Obs table was not found in the DB.")
      }
      private$obs_exists <- TRUE
    },
    # Helper to check that var table exists
    check_var_exists = function() {
      # tables are replaced but never dropped, so once seen this holds and
      # the check does not have to open a connection on every call
      if (isTRUE(private$var_exists)) {
        return(invisible(TRUE))
      }
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
        stop("Var table was not found in the DB.")
      }
      private$var_exists <- TRUE
    },
    # Helper to check that var table exists
    check_var_adt_exists = function() {
      # tables are replaced but never dropped, so once seen this holds and
      # the check does not have to open a connection on every call
      if (isTRUE(private$var_adt_exists)) {
        return(invisible(TRUE))
      }
      con <- private$connect_db()
      on.exit(
        {
          if (exists("con") && !is.null(con)) {
            tryCatch(DBI::dbDisconnect(con), error = function(e) invisible())
          }
        }
      )

      res <- "var_adt" %in% DBI::dbGetQuery(con, "SHOW TABLES")[, "name"]

      if (!res) {
        stop("Var table for ADT was not found in the DB.")
      }
      private$var_adt_exists <- TRUE
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

#' Class for storing single cell experimental data in DuckDB (nightly!)
#'
#' @description
#' This class wraps up the DB connection and methods to interact with the
#' observation/cell metadata and the feature/gene metadata.
#'
#' @export
#'
#' @import data.table
#'
#' @keywords internal
SingleCellDuckDB <- R6::R6Class(
  # class name
  classname = "SingleCellDuckDB",
  # public functions, slots
  inherit = SingleCellDuckDBBase,
  public = list(
    ################
    # Readers h5ad #
    ################

    #' @description
    #' This function populates the obs table from an h5 file (if found).
    #'
    #' @param h5_path String. Path to the h5 file from which to load in the
    #' observation table.
    #' @param filter Optional integer. 1-indexed positions of the obs rows to
    #' keep, in output order.
    #' @param cell_id_col Optional string. If you can see that the first column
    #' in the h5ad obs table is NOT the cell identifier, you can provide
    #' here the correct column name to use as cell_name.
    #'
    #' @returns Returns invisible self. As a side effect, it will load in the
    #' obs data from the h5ad file into the DuckDB.
    populate_obs_from_h5ad = function(
      h5_path,
      filter = NULL,
      cell_id_col = NULL
    ) {
      checkmate::assertFileExists(h5_path)
      checkmate::qassert(filter, c("I+", "0"))
      checkmate::qassert(cell_id_col, c("S1", "0"))

      parquet_path <- tempfile(fileext = ".parquet")
      on.exit(unlink(parquet_path), add = TRUE)

      info <- rs_h5ad_frame_to_parquet(
        h5_path = h5_path,
        frame = "obs",
        keep = if (is.null(filter)) NULL else filter - 1L,
        id_col = cell_id_col,
        out_path = parquet_path
      )
      .warn_skipped_columns(info, "/obs")

      con <- private$connect_db()
      on.exit(
        {
          if (exists("con") && !is.null(con)) {
            tryCatch(DBI::dbDisconnect(con), error = function(e) invisible())
          }
        },
        add = TRUE
      )

      .parquet_to_duckdb(con, "obs", parquet_path, info)

      invisible(self)
    },

    #' @description
    #' This populates the vars table from an h5 file (if found.)
    #'
    #' @param h5_path String. Path to the h5 file from which to load in the
    #' observation table.
    #' @param filter Optional integer. 1-indexed positions of the var rows to
    #' keep, in output order.
    #'
    #' @returns Returns invisible self. As a side effect, it will load in the
    #' var data from the h5ad file into the DuckDB.
    populate_vars_from_h5ad = function(h5_path, filter = NULL) {
      checkmate::assertFileExists(h5_path)
      checkmate::qassert(filter, c("I+", "0"))

      parquet_path <- tempfile(fileext = ".parquet")
      on.exit(unlink(parquet_path), add = TRUE)

      info <- rs_h5ad_frame_to_parquet(
        h5_path = h5_path,
        frame = "var",
        keep = if (is.null(filter)) NULL else filter - 1L,
        id_col = NULL,
        out_path = parquet_path
      )
      .warn_skipped_columns(info, "/var")

      con <- private$connect_db()
      on.exit(
        {
          if (exists("con") && !is.null(con)) {
            tryCatch(DBI::dbDisconnect(con), error = function(e) invisible())
          }
        },
        add = TRUE
      )

      .parquet_to_duckdb(con, "var", parquet_path, info)

      invisible(self)
    },

    #' @description
    #' Populate the obs table from multiple h5ad files.
    #'
    #' @param per_file_info List of lists, each containing: h5_path, exp_id,
    #'   cell_filter (1-indexed integer vector of cells to keep).
    #' @param cell_id_col Optional string. Column name for cell identifiers.
    #'
    #' @returns Invisible self. Populates the obs table in DuckDB.
    populate_obs_from_multi_h5ad = function(per_file_info, cell_id_col = NULL) {
      checkmate::assertList(per_file_info, min.len = 2L)
      checkmate::qassert(cell_id_col, c("S1", "0"))

      parquet_paths <- replicate(
        length(per_file_info),
        tempfile(fileext = ".parquet")
      )
      on.exit(unlink(parquet_paths), add = TRUE)

      parts <- purrr::map2(per_file_info, parquet_paths, \(fi, path) {
        checkmate::assertFileExists(fi$h5_path)
        checkmate::qassert(fi$exp_id, "S1")
        checkmate::qassert(fi$cell_filter, "I+")

        info <- rs_h5ad_frame_to_parquet(
          h5_path = fi$h5_path,
          frame = "obs",
          keep = fi$cell_filter - 1L,
          id_col = cell_id_col,
          out_path = path
        )
        .warn_skipped_columns(info, "/obs")
        list(path = path, info = info, exp_id = fi$exp_id)
      })

      con <- private$connect_db()
      on.exit(
        {
          if (exists("con") && !is.null(con)) {
            tryCatch(DBI::dbDisconnect(con), error = function(e) invisible())
          }
        },
        add = TRUE
      )

      # columns absent from a given file are NULL-filled rather than dropped
      .parquet_parts_to_duckdb(
        con,
        table = "obs",
        parts = parts,
        idx_col = "cell_idx",
        id_col = "cell_id"
      )

      invisible(self)
    },

    #' @description
    #' Populate a minimal var table for multi-file ingestion: gene_idx and
    #' gene_id only. Per-file var annotations are intentionally not merged;
    #' downstream annotation is the caller's responsibility (e.g. via an
    #' external reference and a future enrich_var()).
    #'
    #' @param final_gene_names Character vector. Gene ids in the final order
    #' (matches the Rust binary gene axis).
    #'
    #' @returns Invisible self.
    populate_var_minimal = function(final_gene_names) {
      checkmate::assertCharacter(final_gene_names, min.len = 1L)

      var_dt <- data.table::data.table(
        gene_idx = seq_along(final_gene_names),
        gene_id = final_gene_names
      )

      con <- private$connect_db()
      on.exit(
        {
          if (exists("con") && !is.null(con)) {
            tryCatch(DBI::dbDisconnect(con), error = function(e) invisible())
          }
        },
        add = TRUE
      )

      DBI::dbWriteTable(con, "var", var_dt, overwrite = TRUE)

      invisible(self)
    },

    ###############
    # Readers mtx #
    ###############

    #' Function to populate the obs table from plain text files
    #'
    #' @param f_path File path to the plain text file.
    #' @param has_hdr Boolean. Does the flat file have a header.
    #' @param filter Optional integer. 1-indexed rows of the file to keep.
    #' `cell_idx` follows the file order of the kept rows.
    #'
    #' @returns Invisible self and populates the internal obs table.
    populate_obs_from_plain_text = function(
      f_path,
      has_hdr,
      filter = NULL
    ) {
      # checks
      checkmate::assertFileExists(f_path)
      checkmate::qassert(has_hdr, "B1")
      checkmate::qassert(filter, c("I+", "0"))

      con <- private$connect_db()
      on.exit(
        {
          if (exists("con") && !is.null(con)) {
            tryCatch(DBI::dbDisconnect(con), error = function(e) invisible())
          }
        },
        add = TRUE
      )

      .plain_text_to_duckdb(
        con,
        f_path = f_path,
        has_hdr = has_hdr,
        filter = filter,
        table = "obs",
        idx_col = "cell_idx",
        id_col = "cell_id"
      )

      invisible(self)
    },

    #' Function to populate the var table from plain text files
    #'
    #' @param f_path File path to the plain text file.
    #' @param has_hdr Boolean. Does the flat file have a header.
    #' @param filter Optional integer. 1-indexed rows of the file to keep.
    #' `gene_idx` follows the file order of the kept rows.
    #'
    #' @returns Invisible self and populates the internal var table.
    populate_var_from_plain_text = function(
      f_path,
      has_hdr,
      filter = NULL
    ) {
      # checks
      checkmate::assertFileExists(f_path)
      checkmate::qassert(has_hdr, "B1")
      checkmate::qassert(filter, c("I+", "0"))

      con <- private$connect_db()
      on.exit(
        {
          if (exists("con") && !is.null(con)) {
            tryCatch(DBI::dbDisconnect(con), error = function(e) invisible())
          }
        },
        add = TRUE
      )

      .plain_text_to_duckdb(
        con,
        f_path = f_path,
        has_hdr = has_hdr,
        filter = filter,
        table = "var",
        idx_col = "gene_idx",
        id_col = "gene_id"
      )

      invisible(self)
    },

    #########################
    # Readers R data.tables #
    #########################

    #' Function to populate the obs table from R
    #'
    #' @param obs_dt data.table
    #' @param filter Optional integer. Row indices to keep
    #'
    #' @returns Invisible self and populates the internal obs table.
    populate_obs_from_data.table = function(obs_dt, filter = NULL) {
      # checks
      checkmate::assertDataTable(obs_dt)
      checkmate::qassert(filter, c("I+", "0"))

      obs_dt <- data.table::copy(obs_dt)

      con <- private$connect_db()
      on.exit(
        {
          if (exists("con") && !is.null(con)) {
            tryCatch(DBI::dbDisconnect(con), error = function(e) invisible())
          }
        }
      )
      # seurat has in times cell_id as a column, causing problems

      # deal with column names
      colnames(obs_dt) <- to_snake_case(colnames(obs_dt))
      colnames(obs_dt)[1] <- "cell_id"
      colnames(obs_dt) <- make.unique(colnames(obs_dt))

      if (!is.null(filter)) {
        obs_dt <- obs_dt[filter]
      }
      obs_dt[, cell_idx := .I]
      setcolorder(
        obs_dt,
        c(
          "cell_idx",
          "cell_id",
          setdiff(names(obs_dt), c("cell_idx", "cell_id"))
        )
      )

      DBI::dbWriteTable(
        con,
        "obs",
        obs_dt,
        overwrite = TRUE
      )

      invisible(self)
    },

    #' Function to populate the var table from R
    #'
    #' @param var_dt data.table
    #' @param filter Optional integer. Row indices to keep
    #'
    #' @returns Invisible self and populates the internal obs table.
    populate_var_from_data.table = function(var_dt, filter = NULL) {
      # checks
      checkmate::assertDataTable(var_dt)
      checkmate::qassert(filter, c("I+", "0"))

      var_dt <- data.table::copy(var_dt)

      con <- private$connect_db()
      on.exit(
        {
          if (exists("con") && !is.null(con)) {
            tryCatch(DBI::dbDisconnect(con), error = function(e) invisible())
          }
        }
      )

      # deal with column names
      colnames(var_dt) <- to_snake_case(colnames(var_dt))
      colnames(var_dt)[1] <- "gene_id"

      if (!is.null(filter)) {
        var_dt <- var_dt[filter]
      }

      var_dt[, gene_idx := .I]
      setcolorder(
        var_dt,
        c(
          "gene_idx",
          "gene_id",
          setdiff(names(var_dt), c("gene_idx", "gene_id"))
        )
      )

      DBI::dbWriteTable(
        con,
        "var",
        var_dt,
        overwrite = TRUE
      )

      invisible(self)
    },

    #' Function to populate the var_adt table from R
    #'
    #' @param var_dt data.table with `feature_idx` and `feature_id` columns.
    #' @param filter Optional integer. Row indices to keep.
    #'
    #' @returns Invisible self and populates the internal var_adt table.
    populate_var_adt_from_data.table = function(var_dt, filter = NULL) {
      # checks
      checkmate::assertDataTable(var_dt)
      checkmate::assertSubset(c("feature_idx", "feature_id"), names(var_dt))
      checkmate::qassert(filter, c("I+", "0"))

      var_dt <- data.table::copy(var_dt)

      con <- private$connect_db()
      on.exit(
        {
          if (exists("con") && !is.null(con)) {
            tryCatch(DBI::dbDisconnect(con), error = function(e) invisible())
          }
        }
      )

      if (!is.null(filter)) {
        var_dt <- var_dt[filter]
      }

      data.table::setcolorder(
        var_dt,
        c(
          "feature_idx",
          "feature_id",
          setdiff(names(var_dt), c("feature_idx", "feature_id"))
        )
      )

      DBI::dbWriteTable(
        con,
        "var_adt",
        var_dt,
        overwrite = TRUE
      )

      invisible(self)
    },

    ##################
    # Readers 10x h5 #
    ##################

    #' @description
    #' Populate the obs table from the barcodes in a 10x CellRanger h5 file.
    #'
    #' @param h5_path String. Path to the 10x h5 file.
    #' @param version String. One of `"v2"` or `"v3"`.
    #' @param filter Optional integer. 1-indexed positions of cells to keep.
    #'
    #' @returns Invisible self. Populates the obs table in DuckDB.
    populate_obs_from_tenx_h5 = function(h5_path, version, filter = NULL) {
      checkmate::assertFileExists(h5_path)
      checkmate::assertChoice(version, c("v2", "v3"))
      checkmate::qassert(filter, c("I+", "0"))

      parquet_path <- tempfile(fileext = ".parquet")
      on.exit(unlink(parquet_path), add = TRUE)

      info <- rs_tenx_frame_to_parquet(
        h5_path = h5_path,
        version = version,
        frame = "obs",
        keep = if (is.null(filter)) NULL else filter - 1L,
        out_path = parquet_path
      )

      con <- private$connect_db()
      on.exit(
        {
          if (exists("con") && !is.null(con)) {
            tryCatch(DBI::dbDisconnect(con), error = function(e) invisible())
          }
        },
        add = TRUE
      )

      .parquet_to_duckdb(con, "obs", parquet_path, info)

      invisible(self)
    },

    #' @description
    #' Populate the var table from the features in a 10x CellRanger h5 file.
    #'
    #' Features are read in full file order (all modalities); `filter` selects
    #' the gene-expression features that survived QC.
    #'
    #' @param h5_path String. Path to the 10x h5 file.
    #' @param version String. One of `"v2"` or `"v3"`.
    #' @param filter Optional integer. 1-indexed positions of features to keep.
    #'
    #' @returns Invisible self. Populates the var table in DuckDB.
    populate_vars_from_tenx_h5 = function(h5_path, version, filter = NULL) {
      checkmate::assertFileExists(h5_path)
      checkmate::assertChoice(version, c("v2", "v3"))
      checkmate::qassert(filter, c("I+", "0"))

      parquet_path <- tempfile(fileext = ".parquet")
      on.exit(unlink(parquet_path), add = TRUE)

      info <- rs_tenx_frame_to_parquet(
        h5_path = h5_path,
        version = version,
        frame = "var",
        keep = if (is.null(filter)) NULL else filter - 1L,
        out_path = parquet_path
      )

      con <- private$connect_db()
      on.exit(
        {
          if (exists("con") && !is.null(con)) {
            tryCatch(DBI::dbDisconnect(con), error = function(e) invisible())
          }
        },
        add = TRUE
      )

      .parquet_to_duckdb(con, "var", parquet_path, info)

      invisible(self)
    },

    #' @description
    #' Populate the obs table from multiple 10x CellRanger h5 files.
    #'
    #' Reads each input's barcodes, applies the cell filter, prefixes
    #' `cell_id` with `exp_id`, and rbindlists.
    #'
    #' @param per_file_info List of lists; each must contain `h5_path`
    #' (string), `version` (`"v2"` or `"v3"`), `exp_id` (string) and
    #' `cell_filter` (1-indexed integer vector).
    #'
    #' @returns Invisible self.
    populate_obs_from_multi_tenx_h5 = function(per_file_info) {
      checkmate::assertList(per_file_info, min.len = 2L)

      parquet_paths <- replicate(
        length(per_file_info),
        tempfile(fileext = ".parquet")
      )
      on.exit(unlink(parquet_paths), add = TRUE)

      parts <- purrr::map2(per_file_info, parquet_paths, \(fi, path) {
        checkmate::assertFileExists(fi$h5_path)
        checkmate::assertChoice(fi$version, c("v2", "v3"))
        checkmate::qassert(fi$exp_id, "S1")
        checkmate::qassert(fi$cell_filter, "I+")

        info <- rs_tenx_frame_to_parquet(
          h5_path = fi$h5_path,
          version = fi$version,
          frame = "obs",
          keep = fi$cell_filter - 1L,
          out_path = path
        )
        list(path = path, info = info, exp_id = fi$exp_id)
      })

      con <- private$connect_db()
      on.exit(
        {
          if (exists("con") && !is.null(con)) {
            tryCatch(DBI::dbDisconnect(con), error = function(e) invisible())
          }
        },
        add = TRUE
      )

      .parquet_parts_to_duckdb(
        con,
        table = "obs",
        parts = parts,
        idx_col = "cell_idx",
        id_col = "cell_id"
      )

      invisible(self)
    },

    #########################
    # From multiple DuckDBs #
    #########################

    #' @description
    #' Populate the obs table by merging obs from multiple source DuckDBs.
    #'
    #' Reads obs rows where `to_keep = TRUE` from each source, ordered by
    #' `cell_idx` ascending (matching the order used in the Rust bin merge).
    #' Prefixes `cell_id` with `exp_id`, intersects column names across all
    #' inputs, and rbindlists.
    #'
    #' @param per_file_info List of lists; each must contain `db_path`
    #'   (string) and `exp_id` (string).
    #'
    #' @returns Invisible self.
    populate_obs_from_multi_duckdb = function(per_file_info) {
      checkmate::assertList(per_file_info, min.len = 2L)

      obs_parts <- vector("list", length(per_file_info))

      for (i in seq_along(per_file_info)) {
        fi <- per_file_info[[i]]
        checkmate::assertFileExists(fi$db_path)
        checkmate::qassert(fi$exp_id, "S1")

        src_con <- DBI::dbConnect(
          duckdb::duckdb(shared_home = FALSE),
          dbdir = fi$db_path,
          read_only = TRUE
        )

        obs_dt <- tryCatch(
          data.table::setDT(DBI::dbGetQuery(
            src_con,
            "SELECT * FROM obs WHERE to_keep = TRUE ORDER BY cell_idx"
          )),
          finally = DBI::dbDisconnect(src_con)
        )

        if ("exp_id" %in% names(obs_dt)) {
          stop(sprintf(
            paste(
              "Input '%s' already has an exp_id column in obs.",
              "Merging already-merged objects is not supported."
            ),
            fi$exp_id
          ))
        }

        drop_cols <- intersect(c("cell_idx", "to_keep"), names(obs_dt))
        if (length(drop_cols) > 0L) {
          obs_dt[, (drop_cols) := NULL]
        }

        obs_dt[, cell_id := paste(fi$exp_id, cell_id, sep = "_")]
        obs_dt[, exp_id := fi$exp_id]

        obs_parts[[i]] <- obs_dt
      }

      # keep the union of all obs columns; columns absent from a given file are
      # NA-filled rather than silently dropped (cell_id / exp_id always present)
      obs_combined <- data.table::rbindlist(
        obs_parts,
        use.names = TRUE,
        fill = TRUE
      )
      obs_combined[, cell_idx := .I]
      data.table::setcolorder(
        obs_combined,
        c(
          "cell_idx",
          "cell_id",
          "exp_id",
          setdiff(names(obs_combined), c("cell_idx", "cell_id", "exp_id"))
        )
      )

      con <- private$connect_db()
      on.exit(
        {
          if (exists("con") && !is.null(con)) {
            tryCatch(DBI::dbDisconnect(con), error = function(e) invisible())
          }
        },
        add = TRUE
      )

      DBI::dbWriteTable(con, "obs", obs_combined, overwrite = TRUE)

      invisible(self)
    },

    #' @description
    #' Populate the obs table from a single source DuckDB, restricted to a
    #' subset of cells and reordered to match.
    #'
    #' Needed because `add_data_obs()` requires an obs table to already exist,
    #' so an object written from scratch out of a subset of a parent (as
    #' `cellsweep_sc()` does) has nothing to attach columns to.
    #'
    #' @param source_db_path String. Path to the source DuckDB.
    #' @param cell_idx_to_keep Integer vector. **1-based** source `cell_idx`
    #'   values, in the order the new object stores them.
    #'
    #' @returns Invisible self.
    populate_obs_from_duckdb_subset = function(
      source_db_path,
      cell_idx_to_keep
    ) {
      checkmate::assertFileExists(source_db_path)
      checkmate::qassert(cell_idx_to_keep, "I+[1,)")

      src_con <- DBI::dbConnect(
        duckdb::duckdb(shared_home = FALSE),
        dbdir = source_db_path,
        read_only = TRUE
      )
      obs_dt <- tryCatch(
        # ORDER BY, not insertion order. The positional subset below is only
        # correct if the rows come back in cell_idx order, and DuckDB only
        # guarantees that while its preserve_insertion_order setting is on.
        data.table::setDT(DBI::dbGetQuery(
          src_con,
          "SELECT * FROM obs ORDER BY cell_idx"
        )),
        finally = DBI::dbDisconnect(src_con)
      )

      if (max(cell_idx_to_keep) > nrow(obs_dt)) {
        stop(sprintf(
          "cell_idx_to_keep goes up to %d but the source obs has %d rows.",
          max(cell_idx_to_keep),
          nrow(obs_dt)
        ))
      }

      # Subset by position rather than by joining on cell_idx: the caller's
      # order is the order the counts were written in, and a join would resort.
      obs_subset <- obs_dt[cell_idx_to_keep]

      drop_cols <- intersect(c("cell_idx", "to_keep"), names(obs_subset))
      if (length(drop_cols) > 0L) {
        obs_subset[, (drop_cols) := NULL]
      }
      obs_subset[, cell_idx := .I]
      data.table::setcolorder(
        obs_subset,
        c("cell_idx", setdiff(names(obs_subset), "cell_idx"))
      )

      con <- private$connect_db()
      on.exit(
        {
          if (exists("con") && !is.null(con)) {
            tryCatch(DBI::dbDisconnect(con), error = function(e) invisible())
          }
        },
        add = TRUE
      )

      DBI::dbWriteTable(con, "obs", obs_subset, overwrite = TRUE)

      invisible(self)
    },

    #' @description
    #' Populate the var table from a source DuckDB, filtered and reordered
    #' to match a target gene set.
    #'
    #' @param source_db_path String. Path to the source DuckDB.
    #' @param final_gene_names Character vector. Gene names in the desired
    #'   final order.
    #'
    #' @returns Invisible self.
    populate_vars_from_duckdb_reordered = function(
      source_db_path,
      final_gene_names
    ) {
      checkmate::assertFileExists(source_db_path)
      checkmate::assertCharacter(final_gene_names, min.len = 1L)

      src_con <- DBI::dbConnect(
        duckdb::duckdb(shared_home = FALSE),
        dbdir = source_db_path,
        read_only = TRUE
      )

      var_dt <- tryCatch(
        data.table::setDT(DBI::dbGetQuery(src_con, "SELECT * FROM var")),
        finally = DBI::dbDisconnect(src_con)
      )

      var_dt <- var_dt[match(final_gene_names, gene_id)]
      # gene_universe = "union" can include genes absent from the var reference
      # source; pin gene_id to the canonical universe id so the var <-> counts
      # mapping is never NA-keyed (counts for those genes are ingested correctly)
      n_missing_meta <- sum(is.na(var_dt$gene_id))
      if (n_missing_meta > 0L) {
        warning(sprintf(
          paste(
            "%d gene(s) absent from the var reference source; metadata limited",
            "to gene_id (gene_universe = 'union')."
          ),
          n_missing_meta
        ))
      }
      var_dt[, gene_id := final_gene_names]

      if ("gene_idx" %in% names(var_dt)) {
        var_dt[, gene_idx := NULL]
      }
      var_dt[, gene_idx := .I]
      data.table::setcolorder(
        var_dt,
        c(
          "gene_idx",
          "gene_id",
          setdiff(names(var_dt), c("gene_idx", "gene_id"))
        )
      )

      con <- private$connect_db()
      on.exit(
        {
          if (exists("con") && !is.null(con)) {
            tryCatch(DBI::dbDisconnect(con), error = function(e) invisible())
          }
        },
        add = TRUE
      )

      DBI::dbWriteTable(con, "var", var_dt, overwrite = TRUE)

      invisible(self)
    },

    ##################################
    # From multiple plain text files #
    ##################################

    #' @description
    #' Populate obs from multiple plain-text barcode files.
    #'
    #' Reads each input's barcodes file, filters to the cells that passed QC,
    #' prefixes `cell_id` with `exp_id`, intersects column names across inputs,
    #' and rbindlists.
    #'
    #' @param per_file_info List of lists; each must contain `f_path` (string),
    #' `exp_id` (string), `has_hdr` (boolean), `cell_filter` (1-indexed integer
    #' vector).
    #'
    #' @returns Invisible self.
    populate_obs_from_multi_plain_text = function(per_file_info) {
      checkmate::assertList(per_file_info, min.len = 2L)

      con <- private$connect_db()
      on.exit(
        {
          if (exists("con") && !is.null(con)) {
            tryCatch(DBI::dbDisconnect(con), error = function(e) invisible())
          }
        },
        add = TRUE
      )

      obs_parts <- vector("list", length(per_file_info))

      for (i in seq_along(per_file_info)) {
        fi <- per_file_info[[i]]
        checkmate::assertFileExists(fi$f_path)
        checkmate::qassert(fi$exp_id, "S1")
        checkmate::qassert(fi$has_hdr, "B1")
        checkmate::qassert(fi$cell_filter, "I+")

        delim <- if (grepl("\\.tsv(\\.gz)?$", fi$f_path, ignore.case = TRUE)) {
          "\t"
        } else {
          ","
        }

        dt <- data.table::fread(
          file = fi$f_path,
          sep = delim,
          header = fi$has_hdr
        )

        if (fi$has_hdr) {
          data.table::setnames(dt, to_snake_case(names(dt)))
          data.table::setnames(dt, names(dt)[1L], "cell_id")
        } else {
          data.table::setnames(dt, names(dt)[1L], "cell_id")
        }

        if ("exp_id" %in% names(dt)) {
          stop(sprintf(
            "Input '%s' barcodes file already has an exp_id column.",
            fi$exp_id
          ))
        }

        dt <- dt[fi$cell_filter]
        dt[, cell_id := paste(fi$exp_id, cell_id, sep = "_")]
        dt[, exp_id := fi$exp_id]

        obs_parts[[i]] <- dt
      }

      # keep the union of all obs columns; columns absent from a given file are
      # NA-filled rather than silently dropped (cell_id / exp_id always present)
      combined <- data.table::rbindlist(
        obs_parts,
        use.names = TRUE,
        fill = TRUE
      )
      combined[, cell_idx := .I]
      data.table::setcolorder(
        combined,
        c(
          "cell_idx",
          "cell_id",
          "exp_id",
          setdiff(names(combined), c("cell_idx", "cell_id", "exp_id"))
        )
      )

      DBI::dbWriteTable(con, "obs", combined, overwrite = TRUE)

      invisible(self)
    }
  )
)

# plain text ingestion ---------------------------------------------------------

#' Load a plain-text metadata file into DuckDB
#'
#' @description
#' Reads a CSV / TSV (delimiter from the extension) into `table`, renaming the
#' first column to `id_col` and snake_casing the others when there is a
#' header, and numbering the kept rows in `idx_col`.
#'
#' Row identity has to follow the file, and neither a window function over
#' `read_csv()` nor an `IN` list with one bound parameter per kept row is
#' safe for that at scale. The file is therefore first materialised into a
#' temporary table with insertion order preserved, whose `rowid` is then the
#' 0-based file row. The filter goes in as a registered data.frame and a semi
#' join, and the index is assigned ordered by `rowid`.
#'
#' @param con DuckDB connection.
#' @param f_path String. The file.
#' @param has_hdr Boolean. Whether the file has a header.
#' @param filter Optional integer. 1-indexed rows to keep. `NULL` keeps all.
#' @param table String. Target table, e.g. `"obs"`.
#' @param idx_col String. Index column, e.g. `"cell_idx"`.
#' @param id_col String. Identifier column, e.g. `"cell_id"`.
#'
#' @returns Invisible `NULL`. Writes the table as a side effect.
#'
#' @keywords internal
.plain_text_to_duckdb <- function(
  con,
  f_path,
  has_hdr,
  filter,
  table,
  idx_col,
  id_col
) {
  # checks
  checkmate::assertFileExists(f_path)
  checkmate::qassert(has_hdr, "B1")
  checkmate::qassert(filter, c("I+", "0"))
  checkmate::qassert(table, "S1")
  checkmate::qassert(idx_col, "S1")
  checkmate::qassert(id_col, "S1")

  delimiter <- if (grepl("\\.tsv(\\.gz)?$", f_path, ignore.case = TRUE))
    "\t" else ","
  raw_tbl <- "bixverse_plain_text_raw"
  keep_tbl <- "bixverse_plain_text_keep"

  DBI::dbExecute(con, "SET preserve_insertion_order = true")
  DBI::dbExecute(
    con,
    sprintf(
      "CREATE OR REPLACE TEMP TABLE %s AS SELECT * FROM read_csv(?, delim = %s, header = %s)",
      raw_tbl,
      DBI::dbQuoteString(con, delimiter),
      has_hdr
    ),
    params = list(f_path)
  )
  on.exit(
    DBI::dbExecute(con, sprintf("DROP TABLE IF EXISTS %s", raw_tbl)),
    add = TRUE
  )

  original <- DBI::dbListFields(con, raw_tbl)
  sanitised <- if (has_hdr) to_snake_case(original) else original
  sanitised[1L] <- id_col
  col_mapping <- paste(
    sprintf("t.%s AS %s", DBI::dbQuoteIdentifier(con, original), sanitised),
    collapse = ", "
  )

  where <- ""
  if (!is.null(filter)) {
    duckdb::duckdb_register(con, keep_tbl, data.frame(row_num = filter - 1L))
    on.exit(duckdb::duckdb_unregister(con, keep_tbl), add = TRUE)
    where <- sprintf(
      "WHERE t.rowid IN (SELECT row_num FROM %s)",
      keep_tbl
    )
  }

  DBI::dbExecute(
    con,
    sprintf(
      paste(
        "CREATE OR REPLACE TABLE %s AS",
        "SELECT ROW_NUMBER() OVER (ORDER BY t.rowid) AS %s, %s",
        "FROM %s t %s ORDER BY t.rowid"
      ),
      table,
      idx_col,
      col_mapping,
      raw_tbl,
      where
    )
  )

  invisible(NULL)
}

# parquet ingestion ------------------------------------------------------------

#' Warn about the metadata columns the Parquet writer skipped
#'
#' @param info List. Output of [rs_h5ad_frame_to_parquet()].
#' @param group String. The h5ad group, e.g. `"/obs"`.
#'
#' @returns Invisible `NULL`. Warns once per skipped column.
#'
#' @keywords internal
.warn_skipped_columns <- function(info, group) {
  # checks
  checkmate::assertList(info)
  checkmate::qassert(group, "S1")

  for (col in info$skipped) {
    warning(sprintf("Skipping %s/%s: unsupported h5ad encoding.", group, col))
  }

  invisible(NULL)
}

#' Largest categorical cast to an inline ENUM literal
#'
#' @description
#' Up to this many levels, a categorical is cast with `ENUM('a', 'b', ...)`
#' written into the SQL, which is what `DBI::dbWriteTable()` produces for a
#' factor. Beyond it, parsing the literal dominates the load (a per-cell
#' categorical can carry millions of levels), so the levels go in through a
#' registered data.frame and a named type instead. Not tuned.
#'
#' @keywords internal
SC_ENUM_LITERAL_MAX <- 10000L

#' Drop a table and the named ENUM types its loads created
#'
#' @description
#' The named types of [.parquet_select_col()] are called
#' `<table>_<column>_enum`. Only the names built from the table's current
#' columns are dropped, so a table whose name prefixes another (`var`,
#' `var_adt`) never touches the other's types. A type cannot be dropped while
#' a table uses it, hence the table goes first.
#'
#' @param con DuckDB connection.
#' @param table String. The table, e.g. `"obs"`.
#'
#' @returns Invisible `NULL`.
#'
#' @keywords internal
.drop_table_and_enums <- function(con, table) {
  # checks
  checkmate::qassert(table, "S1")

  cols <- DBI::dbGetQuery(
    con,
    paste(
      "SELECT column_name FROM duckdb_columns()",
      "WHERE database_name = current_database() AND table_name = ?"
    ),
    params = list(table)
  )$column_name

  DBI::dbExecute(
    con,
    sprintf("DROP TABLE IF EXISTS %s", DBI::dbQuoteIdentifier(con, table))
  )
  for (col in cols) {
    DBI::dbExecute(
      con,
      sprintf(
        "DROP TYPE IF EXISTS %s",
        DBI::dbQuoteIdentifier(con, sprintf("%s_%s_enum", table, col))
      )
    )
  }

  invisible(NULL)
}

#' SQL select expression for one Parquet column
#'
#' @description
#' A categorical with at least one level is cast to an ENUM with its levels in
#' order, which is what `DBI::dbWriteTable()` makes of a factor. `ENUM()` with
#' no levels is invalid, so an empty categorical stays VARCHAR. Up to
#' `enum_literal_max` levels the ENUM is written inline; above it a named type
#' `<table>_<column>_enum` is created from a registered data.frame, which is
#' far cheaper for very many levels. Either way the column comes back into R
#' as a factor with the same levels.
#'
#' @param con DuckDB connection.
#' @param table String. The table being built, used to name the type.
#' @param name String. Column name.
#' @param levels Character vector or `NULL`. The ENUM levels.
#' @param enum_literal_max Integer. Most levels written inline.
#'
#' @returns String. The select expression. Creates the named type as a side
#' effect where one is needed.
#'
#' @keywords internal
.parquet_select_col <- function(
  con,
  table,
  name,
  levels,
  enum_literal_max = SC_ENUM_LITERAL_MAX
) {
  # checks
  checkmate::qassert(table, "S1")
  checkmate::qassert(name, "S1")
  checkmate::qassert(levels, c("0", "S*"))
  checkmate::qassert(enum_literal_max, "X1[0,)")

  col <- DBI::dbQuoteIdentifier(con, name)
  if (length(levels) == 0L) {
    return(as.character(col))
  }

  if (length(levels) <= enum_literal_max) {
    return(sprintf(
      "CAST(%s AS ENUM(%s)) AS %s",
      col,
      paste(DBI::dbQuoteString(con, levels), collapse = ", "),
      col
    ))
  }

  type_name <- DBI::dbQuoteIdentifier(con, sprintf("%s_%s_enum", table, name))
  view <- "bixverse_enum_levels_tmp"
  duckdb::duckdb_register(
    con,
    view,
    data.frame(pos = seq_along(levels), lvl = levels)
  )
  on.exit(duckdb::duckdb_unregister(con, view), add = TRUE)
  DBI::dbExecute(con, sprintf("DROP TYPE IF EXISTS %s", type_name))
  DBI::dbExecute(
    con,
    sprintf(
      "CREATE TYPE %s AS ENUM (SELECT lvl FROM %s ORDER BY pos)",
      type_name,
      view
    )
  )

  sprintf("CAST(%s AS %s) AS %s", col, type_name, col)
}

#' Load a Parquet metadata table into DuckDB
#'
#' @description
#' Runs `CREATE OR REPLACE TABLE ... AS SELECT ... FROM read_parquet(...)`,
#' casting the categorical columns to ENUMs.
#'
#' @param con DuckDB connection.
#' @param table String. Target table, e.g. `"obs"`.
#' @param parquet_path String. The Parquet file.
#' @param info List. Output of [rs_h5ad_frame_to_parquet()] or
#' [rs_tenx_frame_to_parquet()] for that file.
#' @param enum_literal_max Integer. Most levels cast with an inline ENUM, see
#' [SC_ENUM_LITERAL_MAX].
#'
#' @returns Invisible `NULL`. Writes the table as a side effect.
#'
#' @keywords internal
.parquet_to_duckdb <- function(
  con,
  table,
  parquet_path,
  info,
  enum_literal_max = SC_ENUM_LITERAL_MAX
) {
  # checks
  checkmate::qassert(table, "S1")
  checkmate::assertFileExists(parquet_path)
  checkmate::assertList(info)
  checkmate::assertNames(
    names(info),
    must.include = c("names", "types", "categories")
  )
  checkmate::qassert(enum_literal_max, "X1[0,)")

  .drop_table_and_enums(con, table)

  cols <- purrr::map_chr(seq_along(info$names), \(i) {
    levels <- if (info$types[[i]] == "categorical") info$categories[[i]]
    .parquet_select_col(con, table, info$names[[i]], levels, enum_literal_max)
  })

  DBI::dbExecute(
    con,
    sprintf(
      "CREATE OR REPLACE TABLE %s AS SELECT %s FROM read_parquet(%s)",
      DBI::dbQuoteIdentifier(con, table),
      paste(cols, collapse = ", "),
      DBI::dbQuoteString(con, parquet_path)
    )
  )

  invisible(NULL)
}

#' ENUM levels of a column across several Parquet parts
#'
#' @description
#' Reproduces what `data.table::rbindlist()` makes of the column when at least
#' one part holds it as a factor: per part in order, a categorical adds its
#' levels in level order and a string column adds its values in row order,
#' each only when not yet seen.
#'
#' @param con DuckDB connection.
#' @param name String. Column name.
#' @param parts List of parts, see [.parquet_parts_to_duckdb()].
#' @param idx_col String. The per-part row index column.
#'
#' @returns Character vector of levels, or `NULL` if the column is not
#' categorical in any part, or is held as something other than a string in
#' some part.
#'
#' @keywords internal
.parquet_parts_levels <- function(con, name, parts, idx_col) {
  # checks
  checkmate::qassert(name, "S1")
  checkmate::assertList(parts)
  checkmate::qassert(idx_col, "S1")

  types <- purrr::map_chr(parts, \(p) {
    i <- match(name, p$info$names)
    if (is.na(i)) NA_character_ else p$info$types[[i]]
  })
  present <- types[!is.na(types)]
  if (
    !any(present == "categorical") ||
      !all(present %in% c("categorical", "string"))
  ) {
    return(NULL)
  }

  levels <- character()
  for (k in seq_along(parts)) {
    p <- parts[[k]]
    new <- switch(
      types[[k]],
      categorical = p$info$categories[[match(name, p$info$names)]],
      string = DBI::dbGetQuery(
        con,
        sprintf(
          paste(
            "SELECT %1$s AS v FROM read_parquet(%2$s) WHERE %1$s IS NOT NULL",
            "GROUP BY %1$s ORDER BY min(%3$s)"
          ),
          DBI::dbQuoteIdentifier(con, name),
          DBI::dbQuoteString(con, p$path),
          DBI::dbQuoteIdentifier(con, idx_col)
        )
      )$v,
      NULL
    )
    levels <- c(levels, setdiff(new, levels))
  }

  levels
}

#' Load several Parquet metadata tables into one DuckDB table
#'
#' @description
#' The multi-file counterpart of [.parquet_to_duckdb()]. Stacks the parts in
#' order with `UNION ALL BY NAME` (columns missing from a part become `NULL`),
#' adds `exp_id`, prefixes the identifier with it and renumbers the index.
#' Column order is the index, the identifier, `exp_id`, then the remaining
#' columns in order of first appearance, as `data.table::rbindlist(fill =
#' TRUE)` produced.
#'
#' @param con DuckDB connection.
#' @param table String. Target table, e.g. `"obs"`.
#' @param parts List. One element per file, each a list with `path` (the
#' Parquet file), `info` (its writer output) and `exp_id`.
#' @param idx_col String. The index column, e.g. `"cell_idx"`.
#' @param id_col String. The identifier column, e.g. `"cell_id"`.
#' @param enum_literal_max Integer. Most levels cast with an inline ENUM, see
#' [SC_ENUM_LITERAL_MAX].
#'
#' @returns Invisible `NULL`. Writes the table as a side effect.
#'
#' @keywords internal
.parquet_parts_to_duckdb <- function(
  con,
  table,
  parts,
  idx_col,
  id_col,
  enum_literal_max = SC_ENUM_LITERAL_MAX
) {
  # checks
  checkmate::qassert(table, "S1")
  checkmate::assertList(parts, min.len = 1L)
  checkmate::qassert(idx_col, "S1")
  checkmate::qassert(id_col, "S1")
  checkmate::qassert(enum_literal_max, "X1[0,)")

  .drop_table_and_enums(con, table)

  q_id <- function(x) as.character(DBI::dbQuoteIdentifier(con, x))
  q_str <- function(x) as.character(DBI::dbQuoteString(con, x))

  other_cols <- setdiff(
    unique(unlist(purrr::map(parts, \(p) p$info$names))),
    c(idx_col, id_col, "exp_id")
  )

  stacked <- paste(
    purrr::imap_chr(parts, \(p, k) {
      sprintf(
        "SELECT *, %i AS _part, %s AS exp_id FROM read_parquet(%s)",
        k,
        q_str(p$exp_id),
        q_str(p$path)
      )
    }),
    collapse = " UNION ALL BY NAME "
  )

  cols <- c(
    sprintf(
      "CAST(row_number() OVER (ORDER BY _part, %s) AS INTEGER) AS %s",
      q_id(idx_col),
      q_id(idx_col)
    ),
    sprintf(
      "exp_id || '_' || coalesce(%s, 'NA') AS %s",
      q_id(id_col),
      q_id(id_col)
    ),
    "exp_id",
    purrr::map_chr(other_cols, \(col) {
      .parquet_select_col(
        con,
        table,
        col,
        .parquet_parts_levels(con, col, parts, idx_col),
        enum_literal_max
      )
    })
  )

  DBI::dbExecute(
    con,
    sprintf(
      "CREATE OR REPLACE TABLE %s AS SELECT %s FROM (%s) ORDER BY _part, %s",
      q_id(table),
      paste(cols, collapse = ", "),
      stacked,
      q_id(idx_col)
    )
  )

  invisible(NULL)
}
