# obs / var ingestion via Parquet ----------------------------------------------

# The obs and var DuckDB tables are written from Rust through a Parquet file.
# They must match what the previous R path (rhdf5 + data.table +
# dbWriteTable) produced: column names and order, DuckDB types (ENUM levels
# included) and values. The old path is rebuilt below as the reference.

source("helper_sc.R", local = TRUE)

library(magrittr)

test_temp_dir <- sc_test_dir("io_obs_var")

## fixtures --------------------------------------------------------------------

fixture <- sc_test_fixture()

f_syn <- file.path(test_temp_dir, "syn.h5ad")
write_h5ad_sc(
  f_path = f_syn,
  counts = fixture$counts,
  obs = fixture$obs,
  fixture$var,
  .verbose = FALSE
)

f_v2 <- file.path(test_temp_dir, "tenx_v2.h5")
write_tenx_h5_sc(
  f_path = f_v2,
  counts = fixture$counts,
  barcodes = fixture$obs$cell_id,
  features = data.table::data.table(
    id = fixture$var$gene_id,
    name = fixture$var$ensembl_id
  ),
  version = "v2"
)

f_v3 <- file.path(test_temp_dir, "tenx_v3.h5")
write_tenx_h5_sc(
  f_path = f_v3,
  counts = fixture$counts,
  barcodes = fixture$obs$cell_id,
  features = data.table::data.table(
    id = fixture$var$gene_id,
    name = fixture$var$ensembl_id,
    feature_type = "Gene Expression"
  ),
  version = "v3"
)

# anndata 0.13: nullable strings, booleans and integers, enum booleans
f_null <- "synthetic_data/anndata_0_13_nullable.h5ad"
# anndata 0.13: categoricals with an unused level and a missing value, int64,
# float with NaN, mixed case / dotted column names
f_rich <- "synthetic_data/anndata_rich.h5ad"

filt <- function(n) as.integer(sort(c(seq(1L, n, by = 3L), n)))

## reference: the previous R path ----------------------------------------------

old_h5ad_frame <- function(h5_path, group, filter, id_col, idx_name, id_name) {
  h5_content <- data.table::setDT(rhdf5::h5ls(h5_path))
  on.exit(rhdf5::h5closeAll())
  frame <- bixverse:::.read_h5ad_frame(h5_path, group, h5_content)
  dt <- data.table::as.data.table(frame$cols)
  if (!is.null(frame$idx)) {
    dt[, (id_name) := frame$idx]
  } else if (!is.null(id_col)) {
    data.table::setnames(dt, id_col, id_name)
  } else {
    data.table::setnames(dt, names(dt)[1L], id_name)
  }
  others <- setdiff(names(dt), id_name)
  if (length(others) > 0L) {
    data.table::setnames(dt, others, bixverse:::to_snake_case(others))
  }
  if (!is.null(filter)) {
    dt <- dt[filter]
  }
  dt[, (idx_name) := .I]
  data.table::setcolorder(dt, c(idx_name, id_name))
  dt
}

old_tenx_frame <- function(h5_path, version, frame, filter) {
  on.exit(rhdf5::h5closeAll())
  rd <- function(p) as.character(rhdf5::h5read(h5_path, p))
  dt <- if (frame == "obs") {
    data.table::data.table(
      cell_id = rd(if (version == "v3") "matrix/barcodes" else "barcodes")
    )
  } else if (version == "v3") {
    data.table::data.table(
      gene_id = rd("matrix/features/id"),
      gene_name = rd("matrix/features/name"),
      feature_type = rd("matrix/features/feature_type")
    )
  } else {
    data.table::data.table(gene_id = rd("genes"), gene_name = rd("gene_names"))
  }
  if (!is.null(filter)) {
    dt <- dt[filter]
  }
  idx_name <- if (frame == "obs") "cell_idx" else "gene_idx"
  dt[, (idx_name) := .I]
  data.table::setcolorder(dt, idx_name)
  dt
}

old_multi <- function(parts) {
  combined <- data.table::rbindlist(
    purrr::map(parts, \(p) {
      dt <- p$dt[, !"cell_idx"]
      dt[, exp_id := p$exp_id]
      dt[, cell_id := paste(exp_id, cell_id, sep = "_")]
      dt
    }),
    use.names = TRUE,
    fill = TRUE
  )
  combined[, cell_idx := .I]
  data.table::setcolorder(combined, c("cell_idx", "cell_id", "exp_id"))
  combined
}

## comparison helpers ----------------------------------------------------------

# a fresh DuckDB per case; returns list(types, data) of `table`
read_back <- function(table, populate = NULL, dt = NULL) {
  dir <- tempfile("obs_var_db", tmpdir = test_temp_dir)
  dir.create(dir)
  db <- bixverse:::SingleCellDuckDB$new(db_dir = dir)
  if (!is.null(populate)) {
    populate(db)
  }
  con <- DBI::dbConnect(duckdb::duckdb(), file.path(dir, "sc_duckdb.db"))
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE))
  if (!is.null(dt)) {
    DBI::dbWriteTable(con, table, dt, overwrite = TRUE)
  }
  list(
    types = DBI::dbGetQuery(
      con,
      sprintf("SELECT column_name, column_type FROM (DESCRIBE %s)", table)
    ),
    data = DBI::dbGetQuery(con, sprintf("SELECT * FROM %s", table))
  )
}

expect_same_table <- function(new, old, label) {
  expect_equal(
    current = new$types,
    target = old$types,
    info = sprintf("%s - column names and DuckDB types", label)
  )
  expect_equal(
    current = new$data,
    target = old$data,
    info = sprintf("%s - values", label)
  )
}

## h5ad ------------------------------------------------------------------------

h5ad_cases <- list(
  list(f = f_syn, frame = "obs", n = 1000L, filter = TRUE),
  list(f = f_syn, frame = "obs", n = 1000L, filter = FALSE),
  list(f = f_syn, frame = "var", n = 100L, filter = TRUE),
  list(f = f_null, frame = "obs", n = 40L, filter = TRUE),
  list(f = f_null, frame = "var", n = 5L, filter = FALSE),
  list(f = f_rich, frame = "obs", n = 50L, filter = TRUE),
  list(f = f_rich, frame = "var", n = 8L, filter = TRUE)
)

for (case in h5ad_cases) {
  filter <- if (case$filter) filt(case$n) else NULL
  label <- sprintf(
    "h5ad %s %s%s",
    basename(case$f),
    case$frame,
    if (case$filter) " (filtered)" else ""
  )
  table <- case$frame

  new <- read_back(table, populate = \(db) {
    if (case$frame == "obs") {
      db$populate_obs_from_h5ad(case$f, filter = filter)
    } else {
      db$populate_vars_from_h5ad(case$f, filter = filter)
    }
  })
  old <- read_back(
    table,
    dt = old_h5ad_frame(
      case$f,
      group = paste0("/", case$frame),
      filter = filter,
      id_col = NULL,
      idx_name = if (case$frame == "obs") "cell_idx" else "gene_idx",
      id_name = if (case$frame == "obs") "cell_id" else "gene_id"
    )
  )
  expect_same_table(new, old, label)
}

expect_true(
  current = any(grepl(
    "^ENUM",
    read_back("obs", populate = \(db) {
      db$populate_obs_from_h5ad(f_rich)
    })$types$column_type
  )),
  info = "h5ad - categoricals land as ENUM"
)

### named ENUM types ----------------------------------------------------------

# categoricals above the literal threshold get a named type; force that path
# with a threshold of one level and check it against the inline literal
named_dir <- tempfile("obs_var_named", tmpdir = test_temp_dir)
dir.create(named_dir)
pq <- file.path(named_dir, "obs.parquet")
info <- rs_h5ad_frame_to_parquet(
  h5_path = f_rich,
  frame = "obs",
  keep = filt(50L) - 1L,
  id_col = NULL,
  out_path = pq
)
con <- DBI::dbConnect(duckdb::duckdb(), file.path(named_dir, "named.db"))

bixverse:::.parquet_to_duckdb(con, "obs_literal", pq, info)
# twice, so the stale table and type are dropped on a re-populate
for (rep in 1:2) {
  bixverse:::.parquet_to_duckdb(con, "obs", pq, info, enum_literal_max = 1L)
}

named_types <- DBI::dbGetQuery(
  con,
  "SELECT column_name, column_type FROM (DESCRIBE obs)"
)
literal <- DBI::dbGetQuery(con, "SELECT * FROM obs_literal")
named <- DBI::dbGetQuery(con, "SELECT * FROM obs")

expect_equal(
  current = DBI::dbGetQuery(
    con,
    paste(
      "SELECT type_name FROM duckdb_types()",
      "WHERE logical_type = 'ENUM' AND NOT internal ORDER BY type_name"
    )
  )$type_name,
  target = c("obs_cell_type_enum", "obs_donor_enum"),
  info = "named ENUM - one named type per categorical, none left over"
)

expect_equal(
  current = named_types,
  target = DBI::dbGetQuery(
    con,
    "SELECT column_name, column_type FROM (DESCRIBE obs_literal)"
  ),
  info = "named ENUM - DESCRIBE reports the same type as the inline literal"
)

expect_identical(
  current = levels(named$cell_type),
  target = c("t cell", "b cell", "unused"),
  info = "named ENUM - level order kept, unused level included"
)

expect_identical(
  current = named,
  target = literal,
  info = "named ENUM - same factors and values as the inline literal"
)

expect_equal(
  current = DBI::dbGetQuery(
    con,
    paste(
      "SELECT count(*) AS n FROM obs",
      "WHERE cell_type = TRY_CAST('b cell' AS obs_cell_type_enum)"
    )
  )$n,
  target = sum(named$cell_type == "b cell", na.rm = TRUE),
  info = "named ENUM - the DESCRIBE type works in a TRY_CAST filter"
)

DBI::dbDisconnect(con, shutdown = TRUE)

## 10x -------------------------------------------------------------------------

for (version in c("v2", "v3")) {
  f <- if (version == "v2") f_v2 else f_v3
  for (frame in c("obs", "var")) {
    n <- if (frame == "obs") 1000L else 100L
    new <- read_back(frame, populate = \(db) {
      if (frame == "obs") {
        db$populate_obs_from_tenx_h5(f, version, filter = filt(n))
      } else {
        db$populate_vars_from_tenx_h5(f, version, filter = filt(n))
      }
    })
    old <- read_back(
      frame,
      dt = old_tenx_frame(f, version, frame, filter = filt(n))
    )
    expect_same_table(new, old, sprintf("10x %s %s", version, frame))
  }
}

## multi file ------------------------------------------------------------------

multi_h5ad_info <- list(
  list(h5_path = f_null, exp_id = "a", cell_filter = filt(40L)),
  list(h5_path = f_rich, exp_id = "b", cell_filter = filt(50L)),
  list(h5_path = f_null, exp_id = "c", cell_filter = seq_len(40L))
)

new <- read_back("obs", populate = \(db) {
  db$populate_obs_from_multi_h5ad(multi_h5ad_info)
})
old <- read_back(
  "obs",
  dt = old_multi(purrr::map(multi_h5ad_info, \(fi) {
    list(
      exp_id = fi$exp_id,
      dt = old_h5ad_frame(
        fi$h5_path,
        "/obs",
        fi$cell_filter,
        NULL,
        "cell_idx",
        "cell_id"
      )
    )
  }))
)
expect_same_table(new, old, "multi h5ad obs")

multi_tenx_info <- list(
  list(h5_path = f_v2, version = "v2", exp_id = "a", cell_filter = filt(1000L)),
  list(h5_path = f_v3, version = "v3", exp_id = "b", cell_filter = 1:1000)
)

new <- read_back("obs", populate = \(db) {
  db$populate_obs_from_multi_tenx_h5(multi_tenx_info)
})
old <- read_back(
  "obs",
  dt = old_multi(purrr::map(multi_tenx_info, \(fi) {
    list(
      exp_id = fi$exp_id,
      dt = old_tenx_frame(fi$h5_path, fi$version, "obs", fi$cell_filter)
    )
  }))
)
expect_same_table(new, old, "multi 10x obs")

## plain text ------------------------------------------------------------------

# the kept rows must come out in file order with a dense index, whatever order
# a parallel CSV scan produces; the large file is there to get a parallel scan.
# The filter is a set of rows here (a repeated row is kept once), so it is
# deduplicated before building the expectation.
pt_dir <- sc_test_dir("plain_text", test_temp_dir)
write_cellranger_output(
  f_path = pt_dir,
  counts = fixture$counts,
  obs = fixture$obs,
  var = fixture$var,
  rows = "cells",
  format_type = "csv",
  .verbose = FALSE
)

set.seed(1L)
n_big <- 200000L
big_dt <- data.table::data.table(
  `Cell Barcode` = sprintf("cell_%07i", seq_len(n_big)),
  nCount = sample.int(5000L, n_big, replace = TRUE),
  score = round(stats::runif(n_big), 6)
)
f_big <- file.path(pt_dir, "big.csv")
data.table::fwrite(big_dt, f_big)
big_filter <- sort(sample.int(n_big, 180000L))

pt_cases <- list(
  list(
    label = "plain text obs (filtered)",
    f = file.path(pt_dir, "barcodes.csv"),
    frame = "obs",
    filter = unique(filt(1000L))
  ),
  list(
    label = "plain text var (filtered)",
    f = file.path(pt_dir, "features.csv"),
    frame = "var",
    filter = unique(filt(100L))
  ),
  list(
    label = "plain text obs (unfiltered)",
    f = file.path(pt_dir, "barcodes.csv"),
    frame = "obs",
    filter = NULL
  ),
  list(
    label = "plain text large obs",
    f = f_big,
    frame = "obs",
    filter = big_filter
  )
)

for (case in pt_cases) {
  new <- read_back(case$frame, populate = \(db) {
    if (case$frame == "obs") {
      db$populate_obs_from_plain_text(case$f, TRUE, filter = case$filter)
    } else {
      db$populate_var_from_plain_text(case$f, TRUE, filter = case$filter)
    }
  })
  idx_col <- if (case$frame == "obs") "cell_idx" else "gene_idx"
  id_col <- if (case$frame == "obs") "cell_id" else "gene_id"

  expected <- data.table::fread(case$f)
  data.table::setnames(expected, bixverse:::to_snake_case(names(expected)))
  data.table::setnames(expected, 1L, id_col)
  if (!is.null(case$filter)) {
    expected <- expected[case$filter]
  }

  # rows as stored, no ORDER BY: the table itself must be in file order
  data <- new$data
  expect_equal(
    current = data[[idx_col]],
    target = seq_len(nrow(expected)),
    info = sprintf("%s - dense index in stored order", case$label)
  )
  expect_equal(
    current = data[[id_col]],
    target = expected[[id_col]],
    info = sprintf("%s - kept rows in file order", case$label)
  )
  expect_equal(
    current = names(data),
    target = c(idx_col, names(expected)),
    info = sprintf("%s - column names", case$label)
  )
}

## cleanup ---------------------------------------------------------------------

sc_test_cleanup(test_temp_dir)
