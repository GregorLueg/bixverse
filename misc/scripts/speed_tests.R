set.seed(123)

nrows <- 100
ncols <- 100

mat <- matrix(data = rnorm(nrows * ncols), nrow = nrows, ncol = ncols)
rownames(mat) <- sprintf("sample_%i", 1:nrows)
colnames(mat) <- sprintf("feature_%i", 1:ncols)


rust_result <- rs_mutual_info(mat, NULL)

# ensure that the same discretisation is used
infotheo_res <- infotheo::mutinformation(infotheo::discretize(
  random_data,
  disc = "equalwidth",
  nbins = sqrt(nrow(random_data))
))
