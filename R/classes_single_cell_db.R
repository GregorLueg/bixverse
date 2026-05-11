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
    #' @param filtered Boolean. Whether to return all cells or filtered to to_keep cells.
    #'
    #' @return The observation table (if found) as a data.table with optionally
    #' selected indices and/or columns.
    get_obs_table = function(indices = NULL, cols = NULL, filtered = FALSE) {
      # checks
      checkmate::qassert(indices, c("0", "I+"))
      checkmate::qassert(cols, c("0", "S+"))
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

      col_part <- if (is.null(cols)) {
        "*"
      } else {
        paste(cols, collapse = ", ")
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

    #' @description
    #' Returns the indices of the cells that have to_keep = TRUE in the DB.
    #'
    #' @return The index positions (1-index) of the cells to keep.
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
    #' @return A character vector of column names.
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
    #' @return Invisible self after updating the to_keep column in the DuckDB.
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
    #' @return Invisible self while adding the new columns to the obs table
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

      existing_cols <- DBI::dbGetQuery(
        con,
        "SELECT column_name FROM information_schema.columns WHERE table_name = 'obs'"
      )$column_name
      new_cols <- setdiff(names(new_data), "cell_idx")
      cols_to_exclude <- intersect(existing_cols, new_cols)

      exclude_clause <- if (length(cols_to_exclude) > 0) {
        paste0(" EXCLUDE(", paste(cols_to_exclude, collapse = ", "), ")")
      } else {
        ""
      }

      DBI::dbExecute(
        con,
        sprintf(
          "
          CREATE OR REPLACE TABLE obs AS
          SELECT
            original.*%s,
            new_data.* EXCLUDE(cell_idx)
          FROM obs AS original
          JOIN new_data ON original.cell_idx = new_data.cell_idx
          ORDER BY original.cell_idx;
          DROP TABLE new_data
          ",
          exclude_clause
        )
      )

      invisible(self)
    },

    #' @description
    #' Left join new data to the obs table in the DuckDB by cell_idx
    #'
    #' @param new_data A data.table with a cell_idx column to join on.
    #'
    #' @return Invisible self while left joining the new data to the obs table
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

      existing_cols <- DBI::dbGetQuery(
        con,
        "SELECT column_name FROM information_schema.columns WHERE table_name = 'obs'"
      )$column_name
      new_cols <- setdiff(names(new_data), "cell_idx")
      cols_to_exclude <- intersect(existing_cols, new_cols)

      exclude_clause <- if (length(cols_to_exclude) > 0) {
        paste0(" EXCLUDE(", paste(cols_to_exclude, collapse = ", "), ")")
      } else {
        ""
      }

      DBI::dbExecute(
        con,
        sprintf(
          "
          CREATE OR REPLACE TABLE obs AS
          SELECT
            original.*%s,
            new_data.* EXCLUDE(cell_idx)
          FROM obs AS original
          LEFT JOIN new_data ON original.cell_idx = new_data.cell_idx
          ORDER BY original.cell_idx;
          DROP TABLE new_data
          ",
          exclude_clause
        )
      )

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
    #' @return Invisible self while adding the new columns to the var table
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

      existing_cols <- DBI::dbGetQuery(
        con,
        "SELECT column_name FROM information_schema.columns WHERE table_name = 'var'"
      )$column_name
      new_cols <- setdiff(names(new_data), "gene_idx")
      cols_to_exclude <- intersect(existing_cols, new_cols)

      exclude_clause <- if (length(cols_to_exclude) > 0) {
        paste0(" EXCLUDE(", paste(cols_to_exclude, collapse = ", "), ")")
      } else {
        ""
      }

      DBI::dbExecute(
        con,
        sprintf(
          "
          CREATE OR REPLACE TABLE var AS
          SELECT
            original.*%s,
            new_data.* EXCLUDE(gene_idx)
          FROM var AS original
          JOIN new_data ON original.gene_idx = new_data.gene_idx
          ORDER BY original.gene_idx;
          DROP TABLE new_data
          ",
          exclude_clause
        )
      )

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
    #' @return Invisible self.
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
    #' @return Invisible self.
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
        stop("Obs table was not found in the DB.")
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
        stop("Obs table was not found in the DB.")
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
    #' @param filter Optional integer. Positions of obs to read in from file.
    #' @param cell_id_col Optional string. If you can see that the first column
    #' in the h5ad obs table is NOT the cell identifier, you can provide
    #' here the correct column name to use as cell_name.
    #'
    #' @return Returns invisible self. As a side effect, it will load in the
    #' obs data from the h5ad file into the DuckDB.
    populate_obs_from_h5ad = function(
      h5_path,
      filter = NULL,
      cell_id_col = NULL
    ) {
      checkmate::assertFileExists(h5_path)
      checkmate::qassert(filter, c("I+", "0"))
      checkmate::qassert(cell_id_col, c("S1", "0"))

      h5_content <- rhdf5::h5ls(h5_path) |> data.table::setDT()

      on.exit(
        tryCatch(rhdf5::h5closeAll(), error = function(e) invisible()),
        add = TRUE
      )

      new_cat_names <- h5_content[
        group == "/obs/__categories" & otype == "H5I_DATASET",
        name
      ]
      has_new_cats <- length(new_cat_names) > 0L

      obs_groups <- h5_content[
        group == "/obs" & otype == "H5I_GROUP" & name != "__categories",
        name
      ]

      idx_info <- .resolve_h5_index(h5_path, "/obs", h5_content)

      skip <- c("_index", idx_info$idx_col)
      direct <- h5_content[
        group == "/obs" & otype == "H5I_DATASET" & !name %in% skip,
        name
      ]

      if (
        length(direct) == 0L &&
          length(obs_groups) == 0L &&
          !has_new_cats &&
          is.null(idx_info$idx)
      ) {
        stop("No obs data could be found in the h5 file. Nothing was loaded.")
      }

      cols <- list()

      for (g in obs_groups) {
        sub_path <- paste0("/obs/", g)
        sub_entries <- h5_content[
          group == sub_path & otype == "H5I_DATASET",
          name
        ]
        if (all(c("categories", "codes") %in% sub_entries)) {
          categories <- rhdf5::h5read(h5_path, paste0(sub_path, "/categories"))
          codes <- rhdf5::h5read(h5_path, paste0(sub_path, "/codes"))
          codes[codes < 0L] <- NA_integer_
          cols[[g]] <- factor(categories[codes + 1L], levels = categories)
        }
      }

      for (d in direct) {
        raw <- as.vector(rhdf5::h5read(h5_path, paste0("/obs/", d)))
        if (has_new_cats && d %in% new_cat_names) {
          categories <- as.vector(
            rhdf5::h5read(h5_path, paste0("/obs/__categories/", d))
          )
          raw[raw < 0L] <- NA_integer_
          cols[[d]] <- factor(categories[raw + 1L], levels = categories)
        } else {
          cols[[d]] <- raw
        }
      }

      obs_dt <- data.table::as.data.table(cols)

      if (!is.null(idx_info$idx)) {
        obs_dt[, cell_id := idx_info$idx]
        data.table::setcolorder(
          obs_dt,
          c("cell_id", setdiff(names(obs_dt), "cell_id"))
        )
      } else if (!is.null(cell_id_col)) {
        data.table::setnames(obs_dt, cell_id_col, "cell_id")
      } else {
        data.table::setnames(obs_dt, names(obs_dt)[1L], "cell_id")
      }

      other_cols <- setdiff(names(obs_dt), "cell_id")
      if (length(other_cols) > 0L) {
        data.table::setnames(obs_dt, other_cols, to_snake_case(other_cols))
      }

      if (!is.null(filter)) {
        obs_dt <- obs_dt[filter]
      }

      obs_dt[, cell_idx := .I]
      data.table::setcolorder(
        obs_dt,
        c(
          "cell_idx",
          "cell_id",
          setdiff(names(obs_dt), c("cell_idx", "cell_id"))
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

      DBI::dbWriteTable(con, "obs", obs_dt, overwrite = TRUE)

      invisible(self)
    },

    #' @description
    #' This populates the vars table from an h5 file (if found.)
    #'
    #' @param h5_path String. Path to the h5 file from which to load in the
    #' observation table.
    #' @param filter Optional integer. Positions of obs to read in from file.
    #'
    #' @return Returns invisible self. As a side effect, it will load in the
    #' obs data from the h5ad file into the DuckDB.
    populate_vars_from_h5ad = function(h5_path, filter = NULL) {
      checkmate::assertFileExists(h5_path)
      checkmate::qassert(filter, c("I+", "0"))

      h5_content <- rhdf5::h5ls(h5_path) |> data.table::setDT()

      on.exit(
        tryCatch(rhdf5::h5closeAll(), error = function(e) invisible()),
        add = TRUE
      )

      new_cat_names <- h5_content[
        group == "/var/__categories" & otype == "H5I_DATASET",
        name
      ]
      has_new_cats <- length(new_cat_names) > 0L

      var_groups <- h5_content[
        group == "/var" & otype == "H5I_GROUP" & name != "__categories",
        name
      ]

      idx_info <- .resolve_h5_index(h5_path, "/var", h5_content)

      skip <- c("_index", idx_info$idx_col)
      direct <- h5_content[
        group == "/var" & otype == "H5I_DATASET" & !name %in% skip,
        name
      ]

      if (
        length(direct) == 0L &&
          length(var_groups) == 0L &&
          !has_new_cats &&
          is.null(idx_info$idx)
      ) {
        stop("No var data could be found in the h5 file. Nothing was loaded.")
      }

      cols <- list()

      for (g in var_groups) {
        sub_path <- paste0("/var/", g)
        sub_entries <- h5_content[
          group == sub_path & otype == "H5I_DATASET",
          name
        ]
        if (all(c("categories", "codes") %in% sub_entries)) {
          categories <- rhdf5::h5read(h5_path, paste0(sub_path, "/categories"))
          codes <- rhdf5::h5read(h5_path, paste0(sub_path, "/codes"))
          codes[codes < 0L] <- NA_integer_
          cols[[g]] <- factor(categories[codes + 1L], levels = categories)
        }
      }

      for (d in direct) {
        raw <- as.vector(rhdf5::h5read(h5_path, paste0("/var/", d)))
        if (has_new_cats && d %in% new_cat_names) {
          categories <- as.vector(
            rhdf5::h5read(h5_path, paste0("/var/__categories/", d))
          )
          raw[raw < 0L] <- NA_integer_
          cols[[d]] <- factor(categories[raw + 1L], levels = categories)
        } else {
          cols[[d]] <- raw
        }
      }

      var_dt <- data.table::as.data.table(cols)

      if (!is.null(idx_info$idx)) {
        var_dt[, gene_id := idx_info$idx]
        data.table::setcolorder(
          var_dt,
          c("gene_id", setdiff(names(var_dt), "gene_id"))
        )
      } else {
        data.table::setnames(var_dt, names(var_dt)[1L], "gene_id")
      }

      other_cols <- setdiff(names(var_dt), "gene_id")
      if (length(other_cols) > 0L) {
        data.table::setnames(var_dt, other_cols, to_snake_case(other_cols))
      }

      if (!is.null(filter)) {
        var_dt <- var_dt[filter]
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

    #' @description
    #' Populate the obs table from multiple h5ad files.
    #'
    #' @param per_file_info List of lists, each containing: h5_path, exp_id,
    #'   cell_filter (1-indexed integer vector of cells to keep).
    #' @param cell_id_col Optional string. Column name for cell identifiers.
    #'
    #' @return Invisible self. Populates the obs table in DuckDB.
    populate_obs_from_multi_h5ad = function(per_file_info, cell_id_col = NULL) {
      checkmate::assertList(per_file_info, min.len = 2L)
      checkmate::qassert(cell_id_col, c("S1", "0"))

      obs_parts <- vector("list", length(per_file_info))

      for (i in seq_along(per_file_info)) {
        fi <- per_file_info[[i]]

        h5_content <- rhdf5::h5ls(fi$h5_path) |> data.table::setDT()

        new_cat_names <- h5_content[
          group == "/obs/__categories" & otype == "H5I_DATASET",
          name
        ]
        has_new_cats <- length(new_cat_names) > 0L

        obs_groups <- h5_content[
          group == "/obs" & otype == "H5I_GROUP" & name != "__categories",
          name
        ]

        idx_info <- .resolve_h5_index(fi$h5_path, "/obs", h5_content)

        skip <- c("_index", idx_info$idx_col)
        direct <- h5_content[
          group == "/obs" & otype == "H5I_DATASET" & !name %in% skip,
          name
        ]

        if (
          length(direct) == 0L &&
            length(obs_groups) == 0L &&
            !has_new_cats &&
            is.null(idx_info$idx)
        ) {
          stop(sprintf("No obs data found in %s", fi$h5_path))
        }

        cols <- list()

        for (g in obs_groups) {
          sub_path <- paste0("/obs/", g)
          sub_entries <- h5_content[
            group == sub_path & otype == "H5I_DATASET",
            name
          ]
          if (all(c("categories", "codes") %in% sub_entries)) {
            categories <- rhdf5::h5read(
              fi$h5_path,
              paste0(sub_path, "/categories")
            )
            codes <- rhdf5::h5read(fi$h5_path, paste0(sub_path, "/codes"))
            codes[codes < 0L] <- NA_integer_
            cols[[g]] <- factor(categories[codes + 1L], levels = categories)
          }
        }

        for (d in direct) {
          raw <- as.vector(rhdf5::h5read(fi$h5_path, paste0("/obs/", d)))
          if (has_new_cats && d %in% new_cat_names) {
            categories <- as.vector(
              rhdf5::h5read(fi$h5_path, paste0("/obs/__categories/", d))
            )
            raw[raw < 0L] <- NA_integer_
            cols[[d]] <- factor(categories[raw + 1L], levels = categories)
          } else {
            cols[[d]] <- raw
          }
        }

        rhdf5::h5closeAll()

        obs_dt <- data.table::as.data.table(cols)

        if (!is.null(idx_info$idx)) {
          obs_dt[, cell_id := idx_info$idx]
          data.table::setcolorder(
            obs_dt,
            c("cell_id", setdiff(names(obs_dt), "cell_id"))
          )
        } else if (!is.null(cell_id_col)) {
          data.table::setnames(obs_dt, cell_id_col, "cell_id")
        } else {
          data.table::setnames(obs_dt, names(obs_dt)[1L], "cell_id")
        }

        other_cols <- setdiff(names(obs_dt), "cell_id")
        if (length(other_cols) > 0L) {
          data.table::setnames(obs_dt, other_cols, to_snake_case(other_cols))
        }

        obs_dt <- obs_dt[fi$cell_filter]
        obs_dt[, exp_id := fi$exp_id]
        obs_dt[, cell_id := paste(exp_id, cell_id, sep = "_")]

        obs_parts[[i]] <- obs_dt
      }

      shared_cols <- Reduce(intersect, lapply(obs_parts, names))
      obs_parts <- lapply(obs_parts, function(dt) {
        dt[, .SD, .SDcols = shared_cols]
      })

      obs_combined <- data.table::rbindlist(obs_parts, use.names = TRUE)
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
      on.exit({
        if (exists("con") && !is.null(con)) {
          tryCatch(DBI::dbDisconnect(con), error = function(e) invisible())
        }
      })

      DBI::dbWriteTable(con, "obs", obs_combined, overwrite = TRUE)

      invisible(self)
    },

    #' @description
    #' Populate the var table from an h5ad file, filtered and reordered to
    #' match a target gene set.
    #'
    #' @param h5_path String. Path to the reference h5ad file.
    #' @param final_gene_names Character vector. Gene names in the desired
    #'   final order.
    #'
    #' @return Invisible self. Populates the var table in DuckDB.
    populate_vars_from_h5ad_reordered = function(h5_path, final_gene_names) {
      checkmate::assertFileExists(h5_path)
      checkmate::assertCharacter(final_gene_names, min.len = 1L)

      self$populate_vars_from_h5ad(h5_path = h5_path, filter = NULL)

      var_dt <- self$get_vars_table()
      var_dt <- var_dt[match(final_gene_names, gene_id)]
      var_dt[, gene_idx := .I]

      con <- private$connect_db()
      on.exit({
        if (exists("con") && !is.null(con)) {
          tryCatch(DBI::dbDisconnect(con), error = function(e) invisible())
        }
      })

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
    #' @param filter Optional integer. Positions of obs to read in from file.
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
        }
      )

      # detect delimiter from file extension
      delimiter <- if (grepl("\\.tsv$", f_path, ignore.case = TRUE)) {
        "'\t'"
      } else {
        "','"
      }

      header_query <- sprintf(
        "SELECT * FROM read_csv(?, delim=%s, header=%s) LIMIT 0",
        delimiter,
        has_hdr
      )

      header_df <- DBI::dbGetQuery(
        conn = con,
        statement = header_query,
        params = list(f_path)
      )

      original_col_names <- colnames(header_df)

      if (has_hdr) {
        sanitised_cols_names <- to_snake_case(original_col_names)
        sanitised_cols_names[1] <- "cell_id"
      } else {
        sanitised_cols_names <- original_col_names
        sanitised_cols_names[1] <- "cell_id"
      }

      col_mapping <- sprintf(
        '"%s" AS %s',
        original_col_names,
        sanitised_cols_names
      ) %>%
        paste(collapse = ",\n    ")

      query <- if (!is.null(filter)) {
        row_placeholders <- paste(
          rep("?", length(filter)),
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
            FROM read_csv(?, delim=%s, header=%s)
          ) t
          WHERE t.original_row_num IN (%s)",
          col_mapping,
          delimiter,
          has_hdr,
          row_placeholders
        )
      } else {
        sprintf(
          "
          CREATE OR REPLACE TABLE obs AS
          SELECT 
            ROW_NUMBER() OVER() as cell_idx,
            %s
          FROM (
            SELECT * FROM read_csv(?, delim=%s, header=%s)
          )",
          col_mapping,
          delimiter,
          has_hdr
        )
      }

      DBI::dbExecute(
        conn = con,
        statement = query,
        params = c(list(f_path), as.list(filter))
      )

      invisible(self)
    },

    #' Function to populate the var table from plain text files
    #'
    #' @param f_path File path to the plain text file.
    #' @param has_hdr Boolean. Does the flat file have a header.
    #' @param filter Optional integer. Positions of obs to read in from file.
    #'
    #' @returns Invisible self and populates the internal var table.
    populate_var_from_plain_text = function(f_path, has_hdr, filter = NULL) {
      checkmate::assertFileExists(f_path)
      checkmate::qassert(has_hdr, "B1")
      checkmate::qassert(filter, c("I+", "0"))

      con <- private$connect_db()
      on.exit(
        {
          if (exists("con") && !is.null(con)) {
            tryCatch(DBI::dbDisconnect(con), error = function(e) invisible())
          }
        }
      )

      # Detect delimiter from file extension
      delimiter <- if (grepl("\\.tsv$", f_path, ignore.case = TRUE)) {
        "'\t'"
      } else {
        "','"
      }

      header_query <- sprintf(
        "SELECT * FROM read_csv(?, delim=%s, header=%s) LIMIT 0",
        delimiter,
        has_hdr
      )

      header_df <- DBI::dbGetQuery(
        conn = con,
        statement = header_query,
        params = list(f_path)
      )

      original_col_names <- colnames(header_df)

      if (has_hdr) {
        sanitised_cols_names <- to_snake_case(original_col_names)
        sanitised_cols_names[1] <- "gene_id"
      } else {
        sanitised_cols_names <- original_col_names
        sanitised_cols_names[1] <- "gene_id"
      }

      col_mapping <- sprintf(
        '"%s" AS %s',
        original_col_names,
        sanitised_cols_names
      ) %>%
        paste(collapse = ",\n    ")

      query <- if (!is.null(filter)) {
        row_placeholders <- paste(
          rep("?", length(filter)),
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
            FROM read_csv(?, delim=%s, header=%s)
          ) t
          WHERE t.original_row_num IN (%s)",
          col_mapping,
          delimiter,
          has_hdr,
          row_placeholders
        )
      } else {
        sprintf(
          "
          CREATE OR REPLACE TABLE var AS
          SELECT 
            ROW_NUMBER() OVER() AS gene_idx,
            %s
          FROM (
            SELECT * FROM read_csv(?, delim=%s, header=%s)
          )",
          col_mapping,
          delimiter,
          has_hdr
        )
      }

      DBI::dbExecute(
        conn = con,
        statement = query,
        params = c(list(f_path), as.list(filter))
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

      # regenerate feature_idx from row order for contiguous 1-indexed positions
      var_dt[, feature_idx := .I]
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
    #' @return Invisible self.
    populate_obs_from_multi_duckdb = function(per_file_info) {
      checkmate::assertList(per_file_info, min.len = 2L)

      obs_parts <- vector("list", length(per_file_info))

      for (i in seq_along(per_file_info)) {
        fi <- per_file_info[[i]]
        checkmate::assertFileExists(fi$db_path)
        checkmate::qassert(fi$exp_id, "S1")

        src_con <- DBI::dbConnect(
          duckdb::duckdb(),
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

      shared_cols <- Reduce(intersect, lapply(obs_parts, names))
      obs_parts <- lapply(obs_parts, function(dt) {
        dt[, .SD, .SDcols = shared_cols]
      })

      obs_combined <- data.table::rbindlist(obs_parts, use.names = TRUE)
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
    #' Populate the var table from a source DuckDB, filtered and reordered
    #' to match a target gene set.
    #'
    #' @param source_db_path String. Path to the source DuckDB.
    #' @param final_gene_names Character vector. Gene names in the desired
    #'   final order.
    #'
    #' @return Invisible self.
    populate_vars_from_duckdb_reordered = function(
      source_db_path,
      final_gene_names
    ) {
      checkmate::assertFileExists(source_db_path)
      checkmate::assertCharacter(final_gene_names, min.len = 1L)

      src_con <- DBI::dbConnect(
        duckdb::duckdb(),
        dbdir = source_db_path,
        read_only = TRUE
      )

      var_dt <- tryCatch(
        data.table::setDT(DBI::dbGetQuery(src_con, "SELECT * FROM var")),
        finally = DBI::dbDisconnect(src_con)
      )

      var_dt <- var_dt[match(final_gene_names, gene_id)]

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
    #' @return Invisible self.
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

      shared_cols <- Reduce(intersect, lapply(obs_parts, names))
      obs_parts <- lapply(obs_parts, function(d) {
        d[, .SD, .SDcols = shared_cols]
      })

      combined <- data.table::rbindlist(obs_parts, use.names = TRUE)
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
    },

    #' @description
    #' Populate the var table from a single plain-text features file,
    #' filtered and reordered to match a target gene set.
    #'
    #' @param f_path String. Path to one input's features file.
    #' @param has_hdr Boolean.
    #' @param final_gene_names Character vector. Gene IDs in the desired final
    #' order.
    #'
    #' @return Invisible self.
    populate_vars_from_plain_text_reordered = function(
      f_path,
      has_hdr,
      final_gene_names
    ) {
      checkmate::assertFileExists(f_path)
      checkmate::qassert(has_hdr, "B1")
      checkmate::assertCharacter(final_gene_names, min.len = 1L)

      delim <- if (grepl("\\.tsv(\\.gz)?$", f_path, ignore.case = TRUE)) {
        "\t"
      } else {
        ","
      }

      dt <- data.table::fread(
        file = f_path,
        sep = delim,
        header = has_hdr
      )

      if (has_hdr) {
        data.table::setnames(dt, to_snake_case(names(dt)))
        data.table::setnames(dt, names(dt)[1L], "gene_id")
      } else {
        data.table::setnames(dt, names(dt)[1L], "gene_id")
      }

      dt <- dt[match(final_gene_names, gene_id)]

      if ("gene_idx" %in% names(dt)) {
        dt[, gene_idx := NULL]
      }
      dt[, gene_idx := .I]
      data.table::setcolorder(
        dt,
        c("gene_idx", "gene_id", setdiff(names(dt), c("gene_idx", "gene_id")))
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

      DBI::dbWriteTable(con, "var", dt, overwrite = TRUE)

      invisible(self)
    }
  )
)
