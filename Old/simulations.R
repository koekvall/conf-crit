library(lme4)
# Settings are (n, T, lambda_1, lambda_2, sigma, gamma),
# where gamma is scale factor for lambda
do_one_sim <- function(seed, settings){
  
  set.seed(seed)
  
  # Unpack settings
  n_i <- settings[1]
  n_t <- settings[2]
  sig <- settings[3] 
  lam_1 <-settings[4] * settings[6]
  lam_2 <- settings[5] * settings[6]
  do_pow <- settings[7]
 
  # Generate data
  X <- matrix(runif(n_i * n_t, -1, 2) , nrow = n_i, ncol = n_t)
  U1 <- rnorm(n_i, sd = lam_1)
  U2 <- rnorm(n_i, sd = lam_2)
  E <- matrix(rnorm(n_i * n_t, sd = sig), nrow = n_i, ncol = n_t)
  Y <- X + 1 + E
  for(ii in 1:n_i){
    Y[ii, ] <- Y[ii, ] + U1[ii] + U2[ii] * X[ii, ]
  }
  D <- data.frame(c(t(Y)), c(t(X)))
  D$time <- rep(1:n_t, n_i)
  D$unit <- as.factor(rep(1:n_i, each = n_t))
  names(D)[1:2] <- c("y", "x")
  
  # Fit model with sigma estimated
  fit <- lme4::lmer(y ~ x + (1|unit) + (0 + x|unit), data = D, REML = FALSE)
  VC <- as.data.frame(lme4::VarCorr(fit))
  
  Z_tall <- as.matrix(lme4::getME(fit, "Z"))
  X_tall <- lme4::getME(fit, "X")
  
  # Re-fit model with sigma fixed at truth
  obj <- function(theta){
    -lmmstest::log_lik(y = D$y,
                           X = X_tall,
                           Z = Z_tall,
                           Beta = theta[1:2],
                           sigma = sig,
                           lambda = theta[3:4],
                           lam_idx = rep(1:2, each = n_i),
                           diffs = 0)[[1]]
  }
  obj_grad <- function(theta){
    -lmmstest:::score(y = D$y,
                      X = X_tall,
                      Z = Z_tall,
                      Beta = theta[1:2],
                      sigma = sig,
                      lambda = theta[3:4],
                      lam_idx = rep(1:2, each = n_i))[c(1:2, 4:5)]
  }
  
  theta_start <- c(fixef(fit), VC[1:2, 5])
  opt <- optim(theta_start, obj, obj_grad, method = "L-BFGS-B", lower = c(-Inf, -Inf, 0, 0), upper = c(Inf, Inf, Inf, Inf),
               control = list(factr = 1e3))
  
  if(do_pow){ # power against (lam_1, lam_2) = (1e-12, 1e-12)
    lam_null <- c(1e-12, 1e-12)
  } else{
    lam_null <- c(lam_1, lam_2)
  }
  L_null <- diag(rep(lam_null, each = n_i), n_i * 2)
  Sigma_null <- diag(sig^2, n_i * n_t) + Z_tall %*%
    L_null^2 %*% t(Z_tall)

  beta_null <- c(qr.solve(crossprod(X_tall, qr.solve(Sigma_null, X_tall)),
                          crossprod(X_tall, qr.solve(Sigma_null, D$y))))
  
  ll_null <- lmmstest::log_lik(y = D$y,
                               X = X_tall,
                               Z = Z_tall,
                               Beta = beta_null,
                               sigma = sig,
                               lambda = lam_null,
                               lam_idx = rep(1:2, each = n_i),
                               diffs = 0)[[1]]
  ll_alt <- lmmstest::log_lik(y = D$y,
                               X = X_tall,
                               Z = Z_tall,
                               Beta = opt$par[1:2],
                               sigma = sig,
                               lambda = opt$par[3:4],
                               lam_idx = rep(1:2, each = n_i),
                               diffs = 0)[[1]]
  lrt_stat <- 2 * (ll_alt - ll_null)
  lrt_p_val <- pchisq(lrt_stat, df = 2, lower.tail = FALSE)
  
  # Wald test
  finf <- lmmstest:::fish_inf(y = D$y,
                      X = X_tall,
                      Z = Z_tall,
                      Beta = opt$par[1:2],
                      sigma = sig,
                      lambda = opt$par[3:4],
                      lam_idx = rep(1:2, each = n_i))
  e <- opt$par[3:4] - lam_null
  wald_stat <- c(crossprod(e, finf[4:5, 4:5] %*% e))
  wald_p_val <- pchisq(wald_stat, df = 2, lower.tail = FALSE)
  
  # Score test for random effects
  c(lmmstest::score_test(
    y = D$y,
    X = X_tall,
    Z = Z_tall,
    Beta = beta_null,
    sigma = sig,
    lambda = lam_null,
    lam_idx = rep(1:2, each = n_i),
    test_idx = 4:5, # lambda parameters
    fix_idx = 3 # Sigma is known
  ),
  "lrt_chi_sq" = lrt_stat, "lrt_p_val" =lrt_p_val,
  "wald_chi_sq" = wald_stat, "wald_p_val" = wald_p_val)
}

# Do simulation
library(doParallel)
library(doRNG)
cl <- makeCluster(8)
registerDoParallel(cl)
n_i <- c(20, 80)
gamma <- c(1e-12, 0.01, 0.05, seq(0.1, 0.5, length.out = 5))
do_pow <- FALSE
n_sims <- 1e4
results <- list()
idx <- 1
for(ii in 1:length(n_i)){
  for(jj in 1:length(gamma)){
    # Settings are (n, T, sigma,lambda_1, lambda_2, gamma, do_pow)
    settings <- c(n_i[ii], 10, 1, 1, 1, gamma[jj], do_pow)
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
if(do_pow){
  saveRDS(results, "~/GitHub/conf-crit/sims/sims_power_new.Rds")
} else{
  saveRDS(results, "~/GitHub/conf-crit/sims/sims_cover_new.Rds")
}

