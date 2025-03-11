# Check fastICA implementation

library(zeallot)

devtools::document()
rextendr::document()


syn_data = synthetic_signal_matrix()

X = t(syn_data$mat) # Artificial count matrix; transposed for rows = samples, cols = features

meta_data = data.table::data.table(
  sample_id = names(syn_data$group),
  case_control = 'case',
  grp = syn_data$group
)


ica_test = bulk_coexp(X, meta_data)

ica_test = ica_processing(ica_test)

c(X_norm, K) %<-% rs_prepare_whitening(X)

dim(X_norm)

ica_res <- fast_ica_rust(
  X_norm,
  K,
  n_icas = 90L,
  ica_fun = "logcosh",
  seed = 246L,
  .verbose = FALSE
)

results <- ica_res$A

ica_res$A[1:10, ]

results[1:10, ]


ica_res

dim(X_norm)

nrow(X_norm)

n_ica = 15

K_2 <- matrix(K[1:n_ica, ], n_ica, dim(X_norm)[1])

X1 <- K_2 %*% X_norm

set.seed(10101)
w.init <- matrix(rnorm(n_ica ^ 2), n_ica, n_ica)

a = rs_fast_ica(
  X1,
  w.init,
  maxit = 200L,
  alpha = 1,
  tol = 1e-04,
  ica_type = "logcosh",
  verbose = TRUE
)

w <- a$mixing %*% K_2
S <- w %*% X_norm
A <- t(w) %*% solve(w %*% t(w))

A

dim(w)

dim(S)

dim(A)


t(A)

res$A[1:5, 1:5]

# [,1]       [,2]       [,3]
# [1,]  0.7677161 -0.3683827 -0.5243150
# [2,] -0.1558387 -0.9010015  0.4048587
# [3,] -0.6215515 -0.2291080 -0.7491216
# Iteration 1 tol = 0.2508784
# [,1]       [,2]        [,3]
# [1,]  0.6673540 -0.7220799 -0.18231637
# [2,] -0.7372638 -0.6751538 -0.02468729
# [3,] -0.1052654  0.1508904 -0.98292995
# Iteration 2 tol = 0.332646
# [,1]        [,2]       [,3]
# [1,]  0.98664947 -0.07880524 -0.1425221
# [2,] -0.04167781 -0.96816819  0.2468062
# [3,] -0.15743498 -0.23757122 -0.9585271
# Iteration 3 tol = 0.04147293
