# Check fastICA implementation

library(zeallot)

?create_synthetic_signal_matrix()

?fastICA::fastICA()

rextendr::document()

features = 10000
signals = 3

S <- matrix(runif(3 * features), features, signals)
A <- matrix(c(1, 1, 1, -1, 3, 2, 0.5, -1, 2), 3, 3, byrow = TRUE)
X <- S %*% A

tictoc::tic()
res = fastICA::fastICA(X, n.comp = 3, method = "R")
tictoc::toc()

tictoc::tic()
c(X_white, K) %<-% rs_whiten_matrix(X)

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

w <- a$mixing %*% K
S <- w %*% X_white
A <- t(w) %*% solve(w %*% t(w))
tictoc::toc()

w

S

A
