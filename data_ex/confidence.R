library(foreign)
library(lme4)
library(lmmstest)
library(merDeriv)
library(doParallel)
library(doRNG)
library(tictoc)


fev1_dat <- read.dta("http://www.hsph.harvard.edu/fitzmaur/ala2e/fev1.dta")

# Fit model
fit <- lmer(exp(logfev1) ~ age + ht + baseage + baseht + (1|id) +
                    (0 + age|id), data = fev1_dat, REML = FALSE)

# Get confidence regions ------------------------------------------------------

X <- getME(fit, "X")
Z <- as.matrix(getME(fit, "Z"))
y <- getME(fit, "y")
Beta_hat <- unname(fixef(fit))
VC <- as.data.frame(VarCorr(fit))
sig_hat <- VC[3, 5]
lam_hat <- VC[1:2, 5]
theta_hat <- c(Beta_hat, sig_hat, lam_hat)


# lambda_1
get_score_stat_1 <- function(lam_null){
  lmmstest::score_test(y = y,
                       X = X,
                       Z = Z,
                       Beta = Beta_hat,
                       sigma = sig_hat,
                       lambda = c(lam_null, lam_hat[2]),
                       lam_idx = rep(1:2, each = 300),
                       test_idx = 7,
                       efficient = TRUE)["chi_sq"]
}
lam_vec_1 <- seq(0, 0.08, length.out = 50)

cl <- makeCluster(8)
registerDoParallel(cl)
tic()
test_vec_1 <- foreach(kk = 1:length(lam_vec_1), .combine = c,
        .errorhandling = "remove",
        .packages = c("lmmstest")) %dorng%{
        get_score_stat_1(lam_vec_1[kk])
}
toc()
# lambda_2
get_score_stat_2 <- function(lam_null){
  lmmstest::score_test(y = y,
                       X = X,
                       Z = Z,
                       Beta = Beta_hat,
                       sigma = sig_hat,
                       lambda = c(lam_hat[1], lam_null),
                       lam_idx = rep(1:2, each = 300),
                       test_idx = 8, efficient = TRUE)["chi_sq"]
}

lam_vec_2 <- seq(0.018, 0.023, length.out = 50)
tic()
test_vec_2 <- foreach(kk = 1:length(lam_vec_1), .combine = c,
                      .errorhandling = "remove",
                      .packages = c("lmmstest")) %dorng%{
                        get_score_stat_2(lam_vec_2[kk])
                      }
toc()
# Joint lambda
get_score_stat_3 <- function(lam_null){
  lmmstest::score_test(y = y,
                       X = X,
                       Z = Z,
                       Beta = Beta_hat,
                       sigma = sig_hat,
                       lambda = lam_null,
                       lam_idx = rep(1:2, each = 300),
                       test_idx = 7:8)["chi_sq"]
}

test_mat <- matrix(NA, nrow = length(lam_vec_1), ncol = length(lam_vec_2))
tic()
for(ii in 1:length(lam_vec_1)){
  cat("ii = ", ii, "\n")
  out <- foreach(jj = 1:length(lam_vec_2), .combine = c) %dopar% {
    res <- try(get_score_stat_3(c(lam_vec_1[ii], lam_vec_2[jj])), TRUE)
    if(class(res) == "try-error"){
      res <- NA
    }
    res
  }
  test_mat[ii, ] <- out
}
toc()

# sigma
get_score_stat_4 <- function(sig_null){
  lmmstest::score_test(y = y,
                       X = X,
                       Z = Z,
                       Beta = Beta_hat,
                       sigma = sig_null,
                       lambda = lam_hat,
                       lam_idx = rep(1:2, each = 300),
                       test_idx = 6)["chi_sq"]
}

sig_vec <- seq(0.15, 0.17, length.out = 50)
tic()
test_vec_3 <- foreach(kk = 1:length(lam_vec_1), .combine = c,
                      .errorhandling = "remove",
                      .packages = c("lmmstest")) %dorng%{
                        get_score_stat_4(sig_vec[kk])
                      }
toc()
stopCluster(cl)

data_ex_conf <- list(lam_vec_1, lam_vec_2, sig_vec, test_vec_1, test_vec_2,
                     test_mat, test_vec_3, fit)
saveRDS(data_ex_conf, "~/GitHub/conf-crit/data_ex/data_ex_R1.Rds")
