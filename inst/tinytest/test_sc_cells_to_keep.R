# cells to keep: unique, ascending, DuckDB authoritative ----------------------

sc <- demo_single_cells(prepped = FALSE)
dir_data <- S7::prop(sc, "dir_data")

## canonicalisation ------------------------------------------------------------

unsorted <- c(50L, 3L, 17L, 3L, 200L, 1L)
sc <- suppressWarnings(set_cells_to_keep(sc, unsorted))

expect_equal(
  current = get_cells_to_keep(sc),
  target = sort(unique(unsorted)) - 1L,
  info = "set_cells_to_keep - stores unique ascending 0-based indices"
)

expect_equal(
  current = get_sc_duckdb(sc)$get_cells_to_keep(),
  target = sort(unique(unsorted)),
  info = "set_cells_to_keep - DuckDB agrees with the ScMap"
)

cell_names <- get_cell_names(sc)
sc_by_name <- suppressWarnings(
  set_cells_to_keep(sc, cell_names[unsorted])
)

expect_equal(
  current = get_cells_to_keep(sc_by_name),
  target = get_cells_to_keep(sc),
  info = "set_cells_to_keep - names and indices give the same filter"
)

expect_error(
  current = set_cells_to_keep(sc, c(cell_names[1:3], "not_a_cell")),
  pattern = "not in the object",
  info = "set_cells_to_keep - unknown cell names error"
)

## obs columns follow the kept cells -------------------------------------------

new_vals <- sprintf("v_%i", get_cells_to_keep(sc) + 1L)
sc <- set_sc_new_obs_col(sc, col_name = "probe", new_data = new_vals)
obs <- get_sc_obs(sc, filtered = TRUE)

expect_equal(
  current = obs$probe,
  target = sprintf("v_%i", obs$cell_idx),
  info = "set_sc_new_obs_col - values land on the right cells"
)

## extractors join obs by cell id ----------------------------------------------

gene <- get_gene_names(sc)[1]
expr_dt <- extract_gene_expression(sc, features = gene, obs_cols = "probe")

expect_equal(
  current = expr_dt$probe,
  target = obs$probe[match(expr_dt$cell_id, obs$cell_id)],
  info = "extract_gene_expression - obs joined by cell id"
)

## load_existing takes the filter from DuckDB ----------------------------------

save_sc_exp_to_disk(sc, type = "rds")

# another front end changes the filter in the DB after the memory file was
# written
get_sc_duckdb(sc)$set_cells_to_keep(cell_idx_to_keep = c(2L, 4L, 6L))

expect_warning(
  current = sc_reloaded <- load_existing(
    SingleCells(dir_data = dir_data),
    .verbose = FALSE
  ),
  pattern = "differs from the one in the memory file",
  info = "load_existing - warns when DuckDB and memory file disagree"
)

expect_equal(
  current = get_cells_to_keep(sc_reloaded),
  target = c(1L, 3L, 5L),
  info = "load_existing - DuckDB filter wins"
)

get_sc_duckdb(sc)$set_cells_to_keep(
  cell_idx_to_keep = get_cells_to_keep(sc) + 1L
)

expect_silent(
  current = sc_quiet <- load_existing(
    SingleCells(dir_data = dir_data),
    .verbose = FALSE
  ),
  info = "load_existing - no warning when both agree"
)

expect_equal(
  current = get_cells_to_keep(sc_quiet),
  target = get_cells_to_keep(sc),
  info = "load_existing - agreeing filter round trips"
)

unlink(dir_data, recursive = TRUE, force = TRUE)
