# Check fastICA implementation

library(zeallot)

devtools::document()
rextendr::document()

devtools::check()

syn_data = synthetic_signal_matrix()



meta_data = data.table::data.table(
  sample_id = names(syn_data$group),
  case_control = 'case',
  grp = syn_data$group
)


ica_test = bulk_coexp(X, meta_data)

ica_test = ica_processing(ica_test)

c(X_norm, K) %<-% rs_prepare_whitening(X)

dim(X_norm)

dim(K)

ica_res <- fast_ica_rust(
  X_norm,
  K,
  n_icas = 10L,
  ica_fun = "logcosh",
  seed = 246L,
  .verbose = FALSE
)

# Implement a Rust version that does the approach from the MSTD paper

rextendr::document()

syn_data = synthetic_signal_matrix()

X = t(syn_data$mat) # Artificial count matrix; transposed for rows = samples, cols = features

c(X_norm, K) %<-% rs_prepare_whitening(X)

test <- rs_ica_iters(
  x_whiten = X_norm,
  k = K,
  no_comp = 2L,
  no_iters = 100L,
  maxit = 200L,
  alpha = 1,
  tol = 1e-04,
  ica_type = "logcosh",
  random_seed = 123L,
  verbose = FALSE
)

pearson_r <- abs(rs_cor(test$s_combined, spearman = FALSE))

pearson_r[1:5, 1:5]

dim(test$s_combined)

test$converged

length(test$converged)

n_icas = 2L

ica_fun = "logcosh"
seed = NULL
maxit = 200L
alpha = 1
tol =  1e-04
.verbose = FALSE

K_red <- matrix(K[1:n_icas, ], n_icas, dim(X_norm)[1])

all(matrix(K[1:n_icas, ], n_icas, dim(X_norm)[1]) == K[1:n_icas, ])

dim(K[1:n_icas, ])

dim(K_red)

dim(K_red)

seed = 245

dim(X_norm)

dim(X_norm)

X1 <- K_red %*% X_norm

dim(X1)

set.seed(seed)
w_init <- matrix(rnorm(n_icas * n_icas), nrow = n_icas, ncol = n_icas)

c(a, converged) %<-% rs_fast_ica(
  X1,
  w_init,
  maxit = maxit,
  alpha = alpha,
  tol = tol,
  ica_type = ica_fun,
  verbose = FALSE
)

dim(a)

w <- a %*% K_red
S <- w %*% X_norm

to_solve = w %*% t(w)

solved <- solve(w %*% t(w))

A <- t(w) %*% solved

cbind(A, A1)

 A1 <- A

dim(A)

dim(A)
