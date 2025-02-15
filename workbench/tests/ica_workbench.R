# Check fastICA implementation

?fastICA::fastICA()

set.seed(10101)
S <- matrix(runif(10000), 5000, 2)
A <- matrix(c(1, 1, -1, 3), 2, 2, byrow = TRUE)
X <- S %*% A

a <- fastICA::fastICA(
  X,
  2,
  alg.typ = "parallel",
  fun = "logcosh",
  alpha = 1,
  method = "C",
  row.norm = FALSE,
  maxit = 200,
  tol = 0.0001,
  verbose = TRUE
)

rextendr::document()

X <- S %*% A

X <- S %*% A
X1_rs <- rs_whiten_matrix(X)

n.comp = 2
maxit = 200
alpha = 1
tol = 0.0001
verbose = TRUE

Diag <- function(d) {
  if (length(d) > 1L) {
    diag(d)
  } else {
    as.matrix(d)
  }
}

set.seed(10101)
w.init <- matrix(rnorm(n.comp ^ 2), n.comp, n.comp)

X = X1_rs
p <- ncol(X)
W <- w.init
sW <- La.svd(W)
W <- sW$u %*% Diag(1/sW$d) %*% t(sW$u) %*% W
lim <- rep(1000, maxit)
it <- 1
while (lim[it] > tol && it < maxit) {
  wx <- W %*% X
  gwx <- tanh(alpha * wx)
  v1 <- gwx %*% t(X) / p
  g.wx <- alpha * (1 - (gwx)^2)
  v2 <- Diag(apply(g.wx, 1, FUN = mean)) %*% W
  W1 <- v1 - v2
  sW1 <- La.svd(W1)
  W1 <- sW1$u %*% Diag(1 / sW1$d) %*% t(sW1$u) %*% W1
  lim[it + 1] <- max(abs(diag(W1 %*% t(W)) - 1))
  W <- W1
  if (verbose)
    message("Iteration ", it, " tol = ", format(lim[it + 1]))
  it <- it + 1
}


rextendr::document()

X <- S %*% A
X1_rs <- rs_whiten_matrix(X)
n.comp = 2

set.seed(10101)
w.init <- matrix(rnorm(n.comp ^ 2), n.comp, n.comp)


a = rs_fast_ica(
  X1_rs,
  w.init,
  maxit = 200L,
  alpha = 1,
  tol = 0.0001,
  verbose = FALSE
)


X = X1_rs
p <- ncol(X)
W <- w.init
sW <- La.svd(W)
W <- sW$u %*% Diag(1/sW$d) %*% t(sW$u) %*% W
lim <- rep(1000, maxit)
wx <- W %*% X
gwx <- tanh(alpha * wx)
v1 <- gwx %*% t(X) / p
v1
g.wx <- alpha * (1 - (gwx)^2)
v2 <- Diag(apply(g.wx, 1, FUN = mean)) %*% W
v2
W1 <- v1 - v2
sW1 <- La.svd(W1)
W1 <- sW1$u %*% Diag(1 / sW1$d) %*% t(sW1$u) %*% W1

W1

abs(diag(W1 %*% t(W)) - 1)

W1
