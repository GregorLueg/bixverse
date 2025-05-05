# tests ica --------------------------------------------------------------------

S <- cbind(sin((1:1000) / 20), rep((((1:200) - 100) / 100), 5))
A <- matrix(c(0.291, 0.6557, -0.5439, 0.5572), 2, 2)
X <- S %*% A

## ica -------------------------------------------------------------------------

expected_k <- matrix(c(-0.8243566, -2.2245938, -1.7226241, 1.0645727), 2, 2)

# prepare the whitening
c(X_norm, K) %<-% rs_prepare_whitening(X, TRUE, 123L, NULL, NULL, NULL)

expect_equal(
  current = K,
  target = expected_k,
  tolerance = 10e-6,
  info = paste(
    "ICA pre-whitening matrix k"
  )
)

### logcosh implementation -----------------------------------------------------

# conscious complete randomisation
random_seed <- as.integer(sample(1:1000, 1))

# logcosh implementation
rs_ica_res_logcosh <- fast_ica_rust(
  X_norm,
  K,
  n_icas = 2L,
  ica_fun = "logcosh",
  seed = random_seed
)

# pending on seed the sign and the position of the extracted signal
# changes. However, in one case the correlation should be close to abs(1),
# indicating that ICA correctly identified the two sources of the signal

logcosh_correlation_signal_11 <- abs(cor(S[, 1], rs_ica_res_logcosh$S[1, ]))
logcosh_correlation_signal_22 <- abs(cor(S[, 2], rs_ica_res_logcosh$S[2, ]))
logcosh_correlation_signal_12 <- abs(cor(S[, 1], rs_ica_res_logcosh$S[2, ]))
logcosh_correlation_signal_21 <- abs(cor(S[, 2], rs_ica_res_logcosh$S[1, ]))

logcosh_correctly_reconstructed <- purrr::map_lgl(
  c(
    logcosh_correlation_signal_11,
    logcosh_correlation_signal_22,
    logcosh_correlation_signal_12,
    logcosh_correlation_signal_21
  ),
  ~ {
    .x > 0.99
  }
)

expect_true(
  sum(logcosh_correctly_reconstructed) == 2,
  info = paste(
    "ICA signal identification with logcosh"
  )
)

#### exp -----------------------------------------------------------------------

# exp implementation
rs_ica_res_exp <- fast_ica_rust(
  X_norm,
  K,
  n_icas = 2L,
  ica_fun = "exp",
  seed = random_seed
)

exp_correlation_signal_11 <- abs(cor(S[, 1], rs_ica_res_exp$S[1, ]))
exp_correlation_signal_22 <- abs(cor(S[, 2], rs_ica_res_exp$S[2, ]))
exp_correlation_signal_12 <- abs(cor(S[, 1], rs_ica_res_exp$S[2, ]))
exp_correlation_signal_21 <- abs(cor(S[, 2], rs_ica_res_exp$S[1, ]))

exp_correctly_reconstructed <- purrr::map_lgl(
  c(
    logcosh_correlation_signal_11,
    logcosh_correlation_signal_22,
    logcosh_correlation_signal_12,
    logcosh_correlation_signal_21
  ),
  ~ {
    .x > 0.99
  }
)

expect_true(
  sum(exp_correctly_reconstructed) == 2,
  info = paste(
    "ICA signal identification with exp"
  )
)
