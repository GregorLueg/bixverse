# Check fastICA implementation

?fastICA::fastICA()





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

S <- matrix(runif(10000), 5000, 2)
A <- matrix(c(1, 1, -1, 3), 2, 2, byrow = TRUE)
X <- S %*% A

X <- S %*% A
rs_whiten_matrix(X)


n.comp = 2

X <- S %*% A

n <- nrow(X)
p <- ncol(X)

X <- scale(X, scale = FALSE)

X <- t(X)

V <- X %*% t(X)/n

V

?La.svd

1 / sqrt(0.12760457422372956)

s <- La.svd(V)

D <- diag(c(1/sqrt(s$d)))

D

K <- D %*% t(s$u)

dim(K)

K <- matrix(K[1:n.comp, ], n.comp, p)
X1 <- K %*% X

X1[1:2, 1:5]

dim(X)
dim(X1)
