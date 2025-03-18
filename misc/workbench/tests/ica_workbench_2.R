# Random SVD

A = X <- ica_test@processed_data$processed_data

rank = 50

m <- nrow(A)
n <- ncol(A)

P <- matrix(rnorm(n * rank), nrow = n, ncol = rank)

dim(P)

Z <- A %*% P

dim(Z)

QR <- qr(Z)
Q <- qr.Q(QR)

Q

Y <- t(Q) %*% A

dim(Y)

dim(A)

svd_res <- La.svd(Y)

U_tilde <- svd_res$u
Sigma <- svd_res$d
Vt <- svd_res$v

U <- Q %*% U_tilde

dim(U)
dim(Vt)

randomized_svd <- function(A, rank) {
  # Get matrix dimensions
  m <- nrow(A)
  n <- ncol(A)

  # Generate a random Gaussian matrix
  P <- matrix(rnorm(n * rank), nrow = n, ncol = rank)

  # Form the sample matrix Z, which is m x k
  Z <- A %*% P

  print("Matrix Z:")
  print(dim(Z))

  # Orthonormalize Z using QR decomposition
  QR <- qr(Z)
  Q <- qr.Q(QR)


  print("Matrix Q:")
  print(dim(Q))

  # Obtain the low-rank approximation of A
  Y <- t(Q) %*% A

  print("Matrix Y:")
  print(dim(Y))

  # Perform SVD on the low-rank approximation
  SVD <- La.svd(Y)

  U_tilde <- SVD$u

  print("Matrix U_tilde:")
  print(dim(U_tilde))

  Sigma <- SVD$d


  Vt <- SVD$v

  print("Matrix Vt:")
  print(dim(Vt))
  # Obtain the final singular vectors
  U <- Q %*% U_tilde

  # Return the components
  return(list(U = U, Sigma = Sigma, Vt = t(Vt)))
}

A <- ica_test@processed_data$processed_data[, 1:2500]

# n <- nrow(A)
n <- ncol(A)

centered <- scale(A)

centered <- t(A)

v = centered %*% t(centered) / n

rextendr::document()

tictoc::tic()
result_random <- randomized_svd(v, rank = 10)
tictoc::toc()

dim(result_random$U)

length(results_rs$s)

tictoc::tic()
results_rs <- rs_random_svd(v, rank = 10, seed = 123L, oversampling = 10L, n_power_iter = 2L)
tictoc::toc()


results_rs$v[1:5, 1:5]

result_random$Vt[1:5, 1:5]

results_rs$u[1:5, 1:5]

result_random$U[1:5, 1:5]

dim(results_rs$u)

dim(result_random$U)

dim(result_random$Vt)

plot(sqrt(results_rs$s), sqrt(result_random$Sigma))

plot(sqrt(results_rs$s), sqrt(result_random$Sigma))

x <- FALSE
switch(as.integer(x) + 1, "No", "Yes")

tictoc::tic()
result_real <- La.svd(v)
tictoc::toc()


D_real <- diag(c(1/sqrt(result_real$d)))
K_real <- D_real %*% t(result_real$u)

K_real[1:5, 1:5]

K[1, ]

D_random <- diag(c(1/sqrt(result_random$Sigma)))
K_random <- D_random %*% t(result_random$U)

K_random[1:5, 1:5]

D_rust <- diag(c(1/sqrt(results_rs$s)))
K_rust <- D_rust %*% t(results_rs$u)

dim(K_rust)

K_rust[1:2, 1:5]

dim(K_rust)

plot(K_rust[1, ], K_random[1, ])

plot(sqrt(results_rs$s), sqrt(result_real$d[1:10]))

x_whiten <- K_rust[1:10, ] %*% v

x_whiten_original <- K_real[1:10, ] %*% v


cor(t(x_whiten), t(x_whiten_original))

x_whiten[1:2, 1:5]

x_whiten_original[1:2, 1:5]

plot(x_whiten[3, ], x_whiten_original[3, ])

cor(x_whiten[3, ], x_whiten_original[3, ])

dim(K_real)

dim(K_random)

result_real$u

plot(sqrt(result_random$Sigma), sqrt(result_real$d[1:n]))


a <- "A"

params = list(b = "B", c = "C", d = "D")

with(params, paste0(a, b, c, d))

