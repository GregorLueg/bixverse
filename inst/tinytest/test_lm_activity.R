# ulm / mlm activity -----------------------------------------------------------

## data ------------------------------------------------------------------------

# generated with decoupler 2.3.0 (dc.mt.ulm / dc.mt.mlm, tmin = 5): a
# genes x samples matrix and a net with mixed signed and continuous weights,
# targets absent from the data and one source that drops below tmin
fixture <- qs2::qs_read("./test_data/lm_activity_decoupler.qs")

mat <- fixture$mat
net <- fixture$net

## parity with decoupler -------------------------------------------------------

ulm_res <- calc_lm_activity(
  mat = mat,
  net = net,
  lm_params = params_lm_activity(p_values = TRUE)
)

mlm_res <- calc_lm_activity(
  mat = mat,
  net = net,
  lm_params = params_lm_activity(method = "mlm", p_values = TRUE)
)

expect_equal(
  current = colnames(ulm_res$es),
  target = colnames(fixture$ulm_es),
  info = "ulm: same sources survive the pruning as in decoupler"
)

expect_equal(
  current = ulm_res$es,
  target = fixture$ulm_es,
  tolerance = 1e-10,
  info = "ulm: t-values match decoupler"
)

# decoupler adjusts the ULM p-values per sample across sources
expect_equal(
  current = ulm_res$fdr,
  target = fixture$ulm_padj,
  tolerance = 1e-6,
  info = "ulm: BH-adjusted p-values match decoupler"
)

expect_equal(
  current = mlm_res$es,
  target = fixture$mlm_es,
  tolerance = 1e-10,
  info = "mlm: t-values match decoupler"
)

# decoupler leaves the MLM p-values unadjusted
expect_equal(
  current = mlm_res$pvals,
  target = fixture$mlm_pv,
  tolerance = 1e-8,
  info = "mlm: p-values match decoupler"
)

## output shapes ---------------------------------------------------------------

coef_res <- calc_lm_activity(
  mat = mat,
  net = net,
  lm_params = params_lm_activity(stat = "coef")
)

expect_true(
  current = checkmate::testMatrix(
    coef_res$es,
    nrows = ncol(mat),
    ncols = ncol(fixture$ulm_es)
  ),
  info = "coef: samples x sources matrix"
)

expect_null(
  current = coef_res$pvals,
  info = "no p-values unless requested"
)

expect_true(
  current = all(sign(coef_res$es) == sign(ulm_res$es)),
  info = "coef and t share their sign"
)

vec_res <- calc_lm_activity(
  mat = stats::setNames(mat[, 1], rownames(mat)),
  net = net,
  lm_params = params_lm_activity(p_values = TRUE)
)

expect_true(
  current = checkmate::testDataTable(
    vec_res,
    nrows = ncol(fixture$ulm_es)
  ),
  info = "vector input returns a data.table"
)

expect_equal(
  current = vec_res$es,
  target = unname(ulm_res$es[1, ]),
  info = "vector input equals the first column of the matrix"
)

## input handling --------------------------------------------------------------

mat_na <- mat
mat_na[1, 2] <- NA

expect_warning(
  current = calc_lm_activity(mat = mat_na, net = net),
  pattern = "non-finite",
  info = "rows with NA are dropped with a warning"
)

net_dup <- rbind(net, net[1, ])

expect_warning(
  current = calc_lm_activity(mat = mat, net = net_dup),
  pattern = "duplicated",
  info = "duplicated source-target pairs give a warning"
)

expect_equal(
  current = suppressWarnings(calc_lm_activity(mat = mat, net = net_dup))$es,
  target = ulm_res$es,
  info = "deduplicated net gives the same scores"
)

expect_warning(
  current = no_res <- calc_lm_activity(
    mat = mat,
    net = net,
    lm_params = params_lm_activity(tmin = 1000L)
  ),
  pattern = "tmin",
  info = "no surviving source warns"
)

expect_null(current = no_res, info = "no surviving source returns NULL")

net_unweighted <- net[, c("source", "target")]
unweighted_res <- calc_lm_activity(mat = mat, net = net_unweighted)

expect_equal(
  current = dim(unweighted_res$es),
  target = dim(ulm_res$es),
  info = "absent weight column is filled with 1"
)

# a source duplicated under another name is perfectly collinear
net_collinear <- rbind(
  net,
  transform(net[net$source == "tf_1", ], source = "tf_1_copy")
)

expect_error(
  current = calc_lm_activity(
    mat = mat,
    net = net_collinear,
    lm_params = params_lm_activity(method = "mlm")
  ),
  pattern = "rank deficient",
  info = "mlm errors on collinear sources"
)

expect_error(
  current = params_lm_activity(method = "wmean"),
  pattern = "should be one of",
  info = "unknown method is refused"
)
