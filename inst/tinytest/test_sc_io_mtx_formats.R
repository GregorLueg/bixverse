# mtx formats: gzip, real and pattern fields -----------------------------------

source("helper_sc.R", local = TRUE)

test_temp_dir <- sc_test_dir("io_mtx_formats")

## fixtures --------------------------------------------------------------------

fixture <- sc_test_fixture()
sc_qc_param <- sc_test_qc_params(fixture)

dir_plain <- sc_test_dir("plain", test_temp_dir)
write_cellranger_output(
  f_path = dir_plain,
  counts = fixture$counts,
  obs = fixture$obs,
  var = fixture$var,
  rows = "cells",
  format_type = "csv",
  .verbose = FALSE
)
mtx_lines <- readLines(file.path(dir_plain, "matrix.mtx"))

# write `lines` to `path`, gzipped when the path ends in .gz
write_lines <- function(lines, path) {
  con <- if (grepl("\\.gz$", path)) gzfile(path, "w") else file(path, "w")
  on.exit(close(con))
  writeLines(lines, con)
}

# same trio in a new directory with the given matrix lines
mtx_dir <- function(name, lines, mtx_name = "matrix.mtx", gz_obs = FALSE) {
  d <- sc_test_dir(name, test_temp_dir)
  write_lines(lines, file.path(d, mtx_name))
  obs_name <- if (gz_obs) "barcodes.csv.gz" else "barcodes.csv"
  write_lines(
    readLines(file.path(dir_plain, "barcodes.csv")),
    file.path(d, obs_name)
  )
  file.copy(file.path(dir_plain, "features.csv"), d)
  params_sc_mtx_io(
    path_mtx = file.path(d, mtx_name),
    path_obs = file.path(d, obs_name),
    path_var = file.path(d, "features.csv"),
    cells_as_rows = TRUE,
    has_hdr = TRUE
  )
}

load_dir <- function(name, io_params, qc = sc_qc_param) {
  dir_data <- sc_test_dir(paste0("sc_", name), test_temp_dir)
  load_mtx(
    object = SingleCells(dir_data = dir_data),
    sc_mtx_io_param = io_params,
    sc_qc_param = qc,
    .verbose = FALSE
  )
}

# the data lines are everything after the banner, comments and shape line
is_data <- !grepl("^%", mtx_lines)
is_data[which(is_data)[1L]] <- FALSE
with_field <- function(lines, field) {
  sub("coordinate integer", paste("coordinate", field), lines, fixed = TRUE)
}

sc_plain <- load_dir(
  "plain",
  params_sc_mtx_io(
    path_mtx = file.path(dir_plain, "matrix.mtx"),
    path_obs = file.path(dir_plain, "barcodes.csv"),
    path_var = file.path(dir_plain, "features.csv"),
    cells_as_rows = TRUE,
    has_hdr = TRUE
  )
)
counts_plain <- sc_plain[,, return_format = "cell"]

## gzip ------------------------------------------------------------------------

sc_gz <- load_dir(
  "gz",
  mtx_dir("gz", mtx_lines, mtx_name = "matrix.mtx.gz", gz_obs = TRUE)
)

expect_equal(
  current = sc_gz[[]],
  target = sc_plain[[]],
  info = "mtx.gz - same obs table as the plain file"
)

expect_equal(
  current = get_sc_var(sc_gz),
  target = get_sc_var(sc_plain),
  info = "mtx.gz - same var table as the plain file"
)

expect_equal(
  current = sc_gz[,, return_format = "cell"],
  target = counts_plain,
  info = "mtx.gz - same counts as the plain file"
)

## real field ------------------------------------------------------------------

real_lines <- with_field(mtx_lines, "real")
parts <- strsplit(real_lines[is_data], " ", fixed = TRUE)
real_lines[is_data] <- vapply(
  parts,
  \(p) sprintf("%s %s %.1f", p[1L], p[2L], as.numeric(p[3L])),
  character(1L)
)

sc_real <- load_dir("real", mtx_dir("real", real_lines))

expect_equal(
  current = sc_real[,, return_format = "cell"],
  target = counts_plain,
  info = "real mtx - integral values load as counts"
)

bad_lines <- real_lines
first_data <- which(is_data)[1L]
p <- strsplit(bad_lines[first_data], " ", fixed = TRUE)[[1L]]
bad_lines[first_data] <- sprintf("%s %s 1.5", p[1L], p[2L])

expect_error(
  current = load_dir("real_bad", mtx_dir("real_bad", bad_lines)),
  pattern = "not a non-negative integral count",
  info = "real mtx - a fractional value is refused"
)

## pattern field ---------------------------------------------------------------

pattern_lines <- with_field(mtx_lines, "pattern")
pattern_lines[is_data] <- sub(" [^ ]+$", "", pattern_lines[is_data])

sc_pattern <- load_dir(
  "pattern",
  mtx_dir("pattern", pattern_lines),
  qc = params_sc_min_quality(
    min_unique_genes = 1L,
    min_lib_size = 1L,
    min_cells = 1L
  )
)
counts_pattern <- sc_pattern[,, return_format = "cell"]
kept <- match(sc_pattern[[]]$cell_id, fixture$obs$cell_id)

expect_true(
  current = length(counts_pattern@x) > 0L && all(counts_pattern@x == 1),
  info = "pattern mtx - every stored count is 1"
)

expect_equivalent(
  current = Matrix::rowSums(counts_pattern),
  target = Matrix::rowSums(fixture$counts[kept, ] != 0),
  info = "pattern mtx - one count per stored entry"
)

## deprecated argument ---------------------------------------------------------

expect_warning(
  current = load_mtx(
    object = SingleCells(
      dir_data = sc_test_dir("sc_deprecated", test_temp_dir)
    ),
    sc_mtx_io_param = mtx_dir("deprecated", mtx_lines),
    sc_qc_param = sc_qc_param,
    mtx_streaming = FALSE,
    .verbose = FALSE
  ),
  pattern = "mtx_streaming",
  info = "load_mtx - mtx_streaming is deprecated"
)

## cleanup ---------------------------------------------------------------------

sc_test_cleanup(test_temp_dir)
