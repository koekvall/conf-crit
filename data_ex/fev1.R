library(foreign)
library(lme4)
library(lmmstest)
library(var)

fev1_dat <- read.dta("http://www.hsph.harvard.edu/fitzmaur/ala2e/fev1.dta")

# Fit model
fit <- lmer(exp(logfev1) ~ age + ht + baseage + baseht + (1|id) +
                    (0 + age|id), data = fev1_dat, REML = FALSE)

# Fit model with sigma "known" (not used)
X <- getME(fit, "X")
Z <- as.matrix(getME(fit, "Z"))
y <- getME(fit, "y")
Beta_hat <- unname(fixef(fit))
VC <- as.data.frame(VarCorr(fit))
lam_hat <- VC[1:2, 5]
# obj <- function(theta){
#   -lmmstest::log_lik(y = y,
#                      X = X,
#                      Z = Z,
#                      Beta = theta[1:5],
#                      sigma = 0.1,
#                      lambda = theta[6:7],
#                      lam_idx = rep(1:2, each = 300),
#                      diffs = 0)[[1]]
# }
# obj_grad <- function(theta){
#   -lmmstest:::score(y = y,
#                      X = X,
#                      Z = Z,
#                      Beta = theta[1:5],
#                      sigma = 0.1,
#                      lambda = theta[6:7],
#                      lam_idx = rep(1:2, each = 300))[c(1:5, 7:8)]
# }
#
# opt <- optim(par = c(Beta_hat, lam_hat),
#              fn = obj,
#              gr = obj_grad,
#              method = "L-BFGS-B",
#              lower = c(rep(-Inf, 5), 0, 0),
#              upper = rep(Inf, 7))

## Gives this
Beta_hat_fix <- c(-1.46427992,  0.07675971,  2.86959668, -0.02103172, -0.96130005,  0.00000000,  0.02237456)


# CI LRT
confint(fit)

# CI Wald
# Use midpoint CI estimate to avoid singularity
inv_finf <- solve(lmmstest:::fish_inf(y = y, X = X, Z = Z, Beta = Beta_hat,
                            sigma = sig_hat, lambda = c(0.02408163, lam_hat[2]),
                            lam_idx = rep(1:2, each = 300))[6:8, 6:8]) # By block-diagonality
c(sig_hat, lam_hat) + 1.96 * sqrt(diag(inv_finf))
c(sig_hat, lam_hat) - 1.96 * sqrt(diag(inv_finf))

# Proposed CI
ex_dat <- readRDS("~/GitHub/conf-crit/data_ex/data_ex_R1.Rds")

# Confidence interval for random intercept standard deviation
lam_val <- ex_dat[[1]]
test_val <- ex_dat[[4]] - qchisq(0.95, 1)
test_val
# Take mean where sign switches
test_val[30:31]
c(0, mean(lam_val[30:31]))

# Confidence interval for random slope standard deviation
lam_val <- ex_dat[[2]]
test_val <- ex_dat[[5]] - qchisq(0.95, 1)
test_val
# Take mean where sign switches
test_val[c(8:9, 39:40)]
c(mean(lam_val[8:9]), mean(lam_val[39:40]))

# Confidence interval for error term standard deviation
sig_val <- ex_dat[[3]]
test_val <- ex_dat[[7]] - qchisq(0.95, 1)
test_val
# Take mean where sign switches
test_val[c(5:6, 30:31)]
c(mean(sig_val[5:6]), mean(sig_val[30:31]))

# Testing
fit_A <- lmer(exp(logfev1) ~ age + ht + baseage + baseht +
                (0 + age|id), data = fev1_dat, REML = FALSE)
fit_A2 <- lmer(exp(logfev1) ~ age + ht + baseage + baseht + (1|id),
               data = fev1_dat, REML = FALSE)
varCompTest(fit, fit_A)
varCompTest(fit, fit_A2)
