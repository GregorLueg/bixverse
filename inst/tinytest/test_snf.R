# snf test ---------------------------------------------------------------------

## synthetic data --------------------------------------------------------------

set.seed(123)
n_samples <- 12
n_features <- 8

# Dataset 1: Continuous only
# Create some structure with correlated features
continuous_data <- matrix(
  rnorm(n_samples * n_features),
  nrow = n_samples,
  ncol = n_features
)

# Add some correlation structure
continuous_data[, 2] <- continuous_data[, 1] + rnorm(n_samples, sd = 0.3)
continuous_data[, 4] <- continuous_data[, 3] + rnorm(n_samples, sd = 0.3)
continuous_data[, 6] <- continuous_data[, 5] + rnorm(n_samples, sd = 0.3)

# Scale to different ranges for more realistic data
continuous_data[, 1:2] <- continuous_data[, 1:2] * 10 + 50 # Age-like
continuous_data[, 3:4] <- continuous_data[, 3:4] * 5 + 20 # BMI-like
continuous_data[, 5:6] <- continuous_data[, 5:6] * 0.5 + 5 # Lab value-like
continuous_data[, 7:8] <- abs(continuous_data[, 7:8]) * 100 # Count-like

rownames(continuous_data) <- sprintf("sample_%02i", 1:12)
colnames(continuous_data) <- sprintf("cont_feat_%i", 1:8)

continuous_df <- as.data.table(continuous_data, keep.rownames = "sample_id")

# Dataset 2: Mixed categorical and continuous
mixed_data <- continuous_df[, 1:4] # Keep first 4 as continuous

# Add categorical features
mixed_data$sex <- factor(sample(c("M", "F"), n_samples, replace = TRUE))
mixed_data$smoking <- factor(sample(
  c("never", "former", "current"),
  n_samples,
  replace = TRUE
))
mixed_data$stage <- factor(sample(
  c("I", "II", "III", "IV"),
  n_samples,
  replace = TRUE
))
mixed_data$status <- factor(sample(
  c("healthy", "diseased"),
  n_samples,
  replace = TRUE
))

# Dataset 3: Categorical only
categorical_df <- data.table(
  sample_ids = mixed_data$sample_id,
  sex = factor(sample(c("M", "F"), n_samples, replace = TRUE)),
  ethnicity = factor(sample(c("A", "B", "C"), n_samples, replace = TRUE)),
  smoking = factor(sample(
    c("never", "former", "current"),
    n_samples,
    replace = TRUE
  )),
  education = factor(sample(
    c("high_school", "bachelors", "masters", "phd"),
    n_samples,
    replace = TRUE
  )),
  stage = factor(sample(c("I", "II", "III", "IV"), n_samples, replace = TRUE)),
  mutation = factor(sample(c("WT", "MUT"), n_samples, replace = TRUE)),
  status = factor(sample(c("healthy", "diseased"), n_samples, replace = TRUE)),
  response = factor(sample(
    c("responder", "non_responder"),
    n_samples,
    replace = TRUE
  ))
)

# tests ------------------------------------------------------------------------

## class generation ------------------------------------------------------------

snf_obj <- snf()

# getter behaving
expect_true(
  current = class(get_snf_params(snf_obj)) == "list",
  info = paste("snf param getter working")
)

expect_warning(
  current = get_snf_adjcacency_mat(snf_obj, name = "x"),
  info = paste("warning if no adjacency matrix can be found")
)

expect_warning(
  current = run_snf(snf_obj),
  info = paste("warning if no adjacency matrices can be found")
)


### generation with already a matrix -------------------------------------------

# should error
expect_error(
  current = snf(data = continuous_data),
  info = paste("error without modality name")
)

snf_obj <- snf(
  data = continuous_data,
  data_name = "continous",
  snf_params = params_snf(k = 3L)
)

retrieved_mat_continuous <- get_snf_adjcacency_mat(snf_obj, name = "continous")

expect_true(
  current = checkmate::testMatrix(
    retrieved_mat_continuous,
    mode = "numeric",
    row.names = "named",
    col.names = "named",
    nrows = n_samples,
    ncols = n_samples
  ),
  info = paste("continuous adjacency matrix has expected type, names and dim")
)

snf_obj <- snf(
  data = mixed_data,
  data_name = "mixed",
  snf_params = params_snf(k = 3L)
)

retrieved_mat_mixed <- get_snf_adjcacency_mat(snf_obj, name = "mixed")

expect_true(
  current = checkmate::testMatrix(
    retrieved_mat_mixed,
    mode = "numeric",
    row.names = "named",
    col.names = "named",
    nrows = n_samples,
    ncols = n_samples
  ),
  info = paste("mixed adjacency matrix has expected type, names and dim")
)

snf_obj <- snf(
  data = categorical_df,
  data_name = "categorical",
  snf_params = params_snf(k = 3L)
)

retrieved_mat_cat <- get_snf_adjcacency_mat(snf_obj, name = "categorical")

expect_true(
  current = checkmate::testMatrix(
    retrieved_mat_cat,
    mode = "numeric",
    row.names = "named",
    col.names = "named",
    nrows = n_samples,
    ncols = n_samples
  ),
  info = paste("categorical adjacency matrix has expected type, names and dim")
)

## methods ---------------------------------------------------------------------

### addition of data -----------------------------------------------------------

snf_obj <- snf(
  data = continuous_data,
  data_name = "continous",
  snf_params = params_snf(k = 3L)
)

expect_error(
  current = add_snf_data_modality(
    object = snf_obj,
    data = mixed_data[-1, ],
    data_name = "mixed"
  ),
  info = paste("error with too few samples")
)


snf_obj <- add_snf_data_modality(
  object = snf_obj,
  data = mixed_data,
  data_name = "mixed"
)

expect_true(
  current = checkmate::testMatrix(
    get_snf_adjcacency_mat(snf_obj, name = "mixed"),
    mode = "numeric",
    row.names = "named",
    col.names = "named",
    nrows = n_samples,
    ncols = n_samples
  ),
  info = paste(
    "s7 method - mixed adjacency matrix has expected type, names and dim"
  )
)

snf_obj <- add_snf_data_modality(
  object = snf_obj,
  data = categorical_df,
  data_name = "categorical"
)

expect_true(
  current = checkmate::testMatrix(
    get_snf_adjcacency_mat(snf_obj, name = "categorical"),
    mode = "numeric",
    row.names = "named",
    col.names = "named",
    nrows = n_samples,
    ncols = n_samples
  ),
  info = paste(
    "s7 method - categorical adjacency matrix has expected type, names and dim"
  )
)

### snf algorithm --------------------------------------------------------------

# are all of the adjacency matrices in the object ... ?
expect_true(
  current = checkmate::testList(
    S7::prop(snf_obj, "adj_matrices"),
    types = "matrix",
    names = "named",
    len = 3L
  ),
  info = "addition of all data modalities works as expected"
)

snf_obj <- run_snf(snf_obj)

matrix_1 <- get_snf_final_mat(snf_obj)

expect_true(
  current = checkmate::testMatrix(
    matrix_1,
    mode = "numeric",
    row.names = "named",
    col.names = "named",
    nrows = n_samples,
    ncols = n_samples
  ),
  info = paste(
    "SNF matrix generation works"
  )
)

# version where I only keep the mixed and continuous one
snf_obj <- run_snf(snf_obj, to_include = c("continous", "mixed"))

matrix_2 <- get_snf_final_mat(snf_obj)

# they should not be the same
expect_true(
  current = any(matrix_1 != matrix_2),
  info = paste(
    "subselection of adjacency matrices works"
  )
)
