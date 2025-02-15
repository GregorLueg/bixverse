# Check fastICA implementation

?create_synthetic_signal_matrix()

rextendr::document()

?rs_contrastive_pca

synthetic_data = create_synthetic_signal_matrix()

plot_synthetic_GEX_HT(synthetic_data)

data = synthetic_data$mat

dim(data)

set.seed(10101)
random_ord = sample(seq(1, 1000), 1000)

data = data[random_ord, ]

data_whitened <- rs_whiten_matrix(t(data))

dim(data_whitened)

n_ica = 200

set.seed(10101)
w.init <- matrix(rnorm(n_ica ^ 2), n_ica, n_ica)

dim(data_whitened)

data_whitened_red = data_whitened[1:n_ica, ]

dim(data_whitened_red)

dim(w.init)

?rs_fast_ica

a = rs_fast_ica(
  data_whitened_red,
  w.init,
  maxit = 200L,
  alpha = 1,
  tol = 0.0001,
  ica_type = "exp",
  verbose = FALSE
)

S <- matrix(runif(15000), 5000, 3)
A <- matrix(c(1, 1, 1, -1, 3, 2, 0.5, -1, 2), 3, 3, byrow = TRUE)
X <- S %*% A

X_white = rs_whiten_matrix(X)

n_ica = dim(X_white)[1]

set.seed(10101)
w.init <- matrix(rnorm(n_ica ^ 2), n_ica, n_ica)

a = rs_fast_ica(
  X_white,
  w.init,
  maxit = 200L,
  alpha = 1,
  tol = 0.0001,
  ica_type = "logcosh",
  verbose = FALSE
)
