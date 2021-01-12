library(lme4)
# Settings are (n, T, sigma_u1, sigma_u2, sigma_e, gamma),
# where gamma is scale factor for lambda entries sigma_u1 and sigma_u2
do_one_sim <- function(seed, settings){
  set.seed(seed)
  # Settings
  n_i <- settings[1]
  n_t <- settings[2]
  sigma_e <- settings[3] 
  sigma_u1 <-settings[4] * settings[6]
  sigma_u2 <- settings[5] * settings[6]
 
  # Generate data
  X <- matrix(runif(n_i * n_t, -1, 2) , nrow = n_i, ncol = n_t)
  U1 <- rnorm(n_i, sd = sigma_u1)
  U2 <- rnorm(n_i, sd = sigma_u2)
  E <- matrix(rnorm(n_i * n_t, sd = sigma_e), nrow = n_i, ncol = n_t)
  Y <- X + 1 + E
  for(ii in 1:n_i){
    Y[ii, ] <- Y[ii, ] + U1[ii] + U2[ii] * X[ii, ]
  }
  D <- data.frame(c(t(Y)), c(t(X)))
  D$time <- rep(1:n_t, n_i)
  D$obs <- as.factor(rep(1:n_i, each = n_t))
  names(D)[1:2] <- c("y", "x")
  
  # Fit model
  fit <- lme4::lmer(y ~ x + (1|obs) + (0 + x|obs), data = D, REML = FALSE)
  VC <- as.data.frame(lme4::VarCorr(fit))
  
  # Fit null model and get lrt stat
  Z_tall <- as.matrix(lme4::getME(fit, "Z"))
  X_tall <- lme4::getME(fit, "X")
  Lam_null <- diag(rep(c(sigma_u1, sigma_u2), each = n_i), n_i * 2)
  Sigma_null <- diag(sigma_e^2, n_i * n_t) + Z_tall %*% Lam_null^2 %*% t(Z_tall)
  beta_null <- c(qr.solve(crossprod(X_tall, qr.solve(Sigma_null, X_tall)),
                          crossprod(X_tall, qr.solve(Sigma_null, D$y))))
  ll_null <- lmmstest::log_lik(y = D$y,
                               X = X_tall,
                               Z = Z_tall,
                               Beta = beta_null,
                               sigma = sigma_e,
                               lambda = c(sigma_u1, sigma_u2),
                               lam_idx = rep(1:2, each = n_i),
                               diffs = 0)[[1]]
  ll_alt <- lmmstest::log_lik(y = D$y,
                               X = X_tall,
                               Z = Z_tall,
                               Beta = unname(fixef(fit)),
                               sigma = VC[3, 5],
                               lambda = VC[1:2, 5],
                               lam_idx = rep(1:2, each = n_i),
                               diffs = 0)[[1]]
  lrt_stat <- 2 * (ll_alt - ll_null)
  lrt_p_val <- pchisq(lrt_stat, df = 3, lower.tail = FALSE)
  
  # Wald test
  finf <- lmmstest:::fish_inf(y = D$y,
                      X = X_tall,
                      Z = Z_tall,
                      Beta = beta_null,
                      sigma = sigma_e,
                      lambda = c(sigma_u1, sigma_u2),
                      lam_idx = rep(1:2, each = n_i))
  e <- VC[c(3, 1, 2), 5] - c(sigma_e, sigma_u1, sigma_u2)
  wald_stat <- c(crossprod(e, finf[3:5, 3:5] %*% e))
  wald_p_val <- pchisq(wald_stat, df = 3, lower.tail = FALSE)
  
  # Score test for random effects
  c(lmmstest::score_test(
    y = D$y,
    X = X_tall,
    Z = Z_tall,
    Beta = unname(fixef(fit)),
    sigma = sigma_e,
    lambda = c(sigma_u1, sigma_u2),
    lam_idx = rep(1:2, each = n_i),
    test_idx = 2 + 1:3 # All variance parameters
  ),
  "lrt_chi_sq" = lrt_stat, "lrt_p_val" =lrt_p_val,
  "wald_chi_sq" = wald_stat, "wald_p_val" = wald_p_val)
}

# Do simulation
library(doParallel)
library(doRNG)
cl <- makeCluster(8)
registerDoParallel(cl)
n_i <- c(20, 40, 80)
gamma <- c(0, 0.01, 0.05, seq(0.1, 0.5, length.out = 5))
n_sims <- 1e4
results <- list()
idx <- 1
for(ii in 1:length(n_i)){
  for(jj in 1:length(gamma)){
    # Settings are (n, T, sigma_e, sigma_u1, sigma_u2, gamma)
    settings <- c(n_i[ii], 10, 1, 1, 1, gamma[jj])
    res_mat <- foreach(kk = 1:n_sims, .combine = rbind,
                       .errorhandling = "remove",
                       .packages = c("lmmstest", "lme4")) %dorng%{
                         do_one_sim((idx - 1) * n_sims + kk, settings)
                       }
    results[[idx]] <- res_mat
    idx <- idx + 1
  }
}
stopCluster(cl)


names(results) <- paste(rep(paste("n", n_i, sep = "_"), length(gamma)),
                        paste("gamma", gamma, sep = "_"), sep = "_")
saveRDS(results, "~/GitHub/int-est/sims/sims_big.Rds")
