out_dir <- "/Users/karekv/Dropbox/Apps/Overleaf/conf_crit/" # end in "/"
PDF <- TRUE


### Introduction ###
get_statistics <- function(y, lam)
{
  n <- length(y)
  ll_fun <- function(y, lam){
    sum(- 0.5 * log(1 + lam^2) - 0.5 * y^2 / (1 + lam^2))
  }
  ll <- ll_fun(y, lam)
  s <- sum(-lam / (1 + lam^2) + y^2 * lam / (1 + lam^2)^2)
  h <- sum(y^2 * (1 - 3 * lam^2) / (1 + lam^2)^3 - (1 - lam^2) / (1 + lam^2)^2)
  finf <- n * 2 * lam^2 / (1 + lam^2)^2
  our_stat <- sum(-1 + y^2 / (1 + lam^2))^2 / (2 * n)

  score_obs_stat <- s^2 / (-h)

  mle <- sqrt(max(c(mean(y^2) - 1, 0)))
  lrt_stat <- 2 * (ll_fun(y, mle) - ll)
  wald_stat <- (mle - lam)^2 * finf
  wald_obs_stat <-  (mle - lam)^2 * (-h)
  return(c(ll = ll, s = s, h = h, finf = finf, our_stat = our_stat,
              score_obs_stat = score_obs_stat, lrt_stat = lrt_stat,
              wald_stat = wald_stat, wald_obs_stat = wald_obs_stat))
}

set.seed(3)
y1 <- rnorm(100, sd = 1.1)
y2 <- rnorm(100, sd = 1.1)
lam_vals <- seq(-0.1, 1, length.out = 1000)
dat1 <- t(sapply(lam_vals, function(x)get_statistics(y1, x)))
dat2 <-  t(sapply(lam_vals, function(x)get_statistics(y2, x)))


# Create figure
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442",
"#0072B2", "#D55E00", "#CC79A7")

if(PDF) pdf(paste(out_dir, "fig_loglik.pdf", sep = ""), width = 6, height = 4)
par(cex.axis = 1.1, cex.lab = 1.1, cex = 1.1)
par(mgp=c(2.2,0.45,0), tcl=-0.4, mar=c(3.3,3.6,1.5,1.1))

# Likelihood
plot(lam_vals, dat1[, "ll"], type = "l", lwd = 2,
     ylab = "log-likelihood", xlab = expression(lambda) ,
     ylim = c(min(dat2[, "ll"]), max(dat1[, "ll"])),
     col = cbbPalette[1])
lines(lam_vals, dat2[, "ll"], type = "l", lwd = 2, lty = 1,
      col = cbbPalette[1])

if(PDF) dev.off()

### Cover figures ####
sim_data <- readRDS("~/GitHub/conf-crit/sims/full_simsR1_1e4.Rds")
n_sims <- nrow(sim_data[[1]])
# Verify this using names(sim_data)
gamma <- c(1e-12, 0.01, 0.05, seq(0.1, 0.5, length.out = 5))
n_vec <- c(80, 20)

# Cover curves
cover <- matrix(0, nrow = length(sim_data), ncol = 8)
for(ii in 1:nrow(cover)){
  cover[ii, 1:3] <- 1 - colMeans(sim_data[[ii]][, c("our_p_val_cov", "lrt_p_val_cov",
                                                    "wald_p_val_cov")] < 0.05)
  cover[ii, 4:6] <- apply(sim_data[[ii]][, c("our_p_val_cov", "lrt_p_val_cov",
                                             "wald_p_val_cov")] < 0.05,
                          2, sd) / sqrt(n_sims)
}
cover[, 7] <- rep(n_vec, each = length(gamma))
cover[, 8] <- rep(gamma, length(n_vec))
colnames(cover) <- c("our", "lrt", "wald", "se_our", "se_lrt", "se_wald",
                     "n", "gamma")

if(PDF) pdf(paste(out_dir, "fig_cover_R1.pdf", sep = ""), width = 10, height = 4)
par(cex.axis = 1.1, cex.lab = 1.1, cex = 1.1)
par(mgp=c(2.2,0.45,0), tcl=-0.4, mar=c(3.3,3.6,1.5,1.1))
par(mfrow = c(1, 2))

for(jj in c(20, 80)){
  plot_dat <- cover[cover[, "n"] == jj, ]
  plot(x = gamma,
       y = plot_dat[, "our"],
       type = "l",
       ylim = c(0.87, 1),
       lwd = 2,
       ylab = "estimated cover probability",
       xlab = expression(lambda),
       main = paste0("N = ", jj),
       col = cbbPalette[1])
  arrows(gamma,
         plot_dat[, "our"],
         gamma,
         plot_dat[, "our"] + 2 * plot_dat[, "se_our"],
         length = 0.05,
         angle = 90,
         col = cbbPalette[1])
  arrows(gamma,
         plot_dat[, "our"],
         gamma,
         plot_dat[, "our"] - 2 * plot_dat[, "se_our"],
         length = 0.05,
         angle = 90,
         col = cbbPalette[1])
  abline(h = 0.95, lwd = 1)
  lines(x = gamma,
        y = plot_dat[, "lrt"],
        lty = 2,
        lwd = 2,
        col = cbbPalette[2])
  arrows(gamma,
         plot_dat[, "lrt"],
         gamma,
         plot_dat[, "lrt"] + 2 * plot_dat[, "se_lrt"],
         length = 0.05,
         angle = 90,
         col = cbbPalette[2])
  arrows(gamma,
         plot_dat[, "lrt"],
         gamma,
         plot_dat[, "lrt"] - 2 * plot_dat[, "se_lrt"],
         length = 0.05,
         angle = 90,
         col = cbbPalette[2])
  lines(x = gamma,
        y = plot_dat[, "wald"],
        lty = 3,
        lwd = 2,
        col = cbbPalette[3])
  arrows(gamma,
         plot_dat[, "wald"],
         gamma,
         plot_dat[, "wald"] + 2 * plot_dat[, "se_wald"],
         length = 0.05,
         angle = 90,
         col = cbbPalette[3])
  arrows(gamma,
         plot_dat[, "wald"],
         gamma,
         plot_dat[, "wald"] - 2 * plot_dat[, "se_wald"],
         length = 0.05,
         angle = 90,
         col = cbbPalette[3])
}
if(PDF) dev.off()
# Distribution
pp <- ppoints(n_sims)
if(PDF) pdf(paste(out_dir, "fig_qq_R1.pdf", sep = ""), width = 10, height = 4)
par(cex.axis = 1.1, cex.lab = 1.1, cex = 1.1)
par(mgp=c(2.2,0.45,0), tcl=-0.4, mar=c(3.3,3.6,1.5,1.1))
par(mfrow = c(1, 2))

# Small but non-zero scale parameters
plot_dat <- sim_data[[2]]
plot(x = qchisq(pp, df = 2),
     y = quantile(plot_dat[, "our_stat_cov"], pp),
     xlab = "theoretical quantiles",
     ylab = "sample quantiles",
     main = bquote(~ lambda == 0.01),
     col = cbbPalette[1])
points(x = qchisq(pp, df = 3),
       y = quantile(plot_dat[, "lrt_stat_cov"], pp),
       pch = 2,
       col = cbbPalette[2])
points(x = qchisq(pp, df = 3),
       y = quantile(plot_dat[, "wald_stat_cov"], pp),
       pch = 3,
       col = cbbPalette[3])
abline(a = 0, b = 1, col = "red")

# Large scale parameters
plot_dat <- sim_data[[8]]
plot(x = qchisq(pp, df = 2),
     y = quantile(plot_dat[, "our_stat_cov"], pp),
     xlab = "theoretical quantiles",
     ylab = "sample quantiles",
     main = bquote(~ lambda == 0.5),
     col = cbbPalette[1])
points(x = qchisq(pp, df = 3),
       y = quantile(plot_dat[, "lrt_stat_cov"], pp),
       pch = 2,
       col = cbbPalette[2])
points(x = qchisq(pp, df = 3),
       y = quantile(plot_dat[, "wald_stat_cov"], pp),
       pch = 3,
       col = cbbPalette[3])
abline(a = 0, b = 1, col = "red")
if(PDF) dev.off()

### Power figures ####
n_sims <- nrow(sim_data[[1]])
# Verify this using names(sim_data)
gamma <- c(1e-12, 0.01, 0.05, seq(0.1, 0.5, length.out = 5))
n_vec <- c(80, 20)

# Power curves
power <- matrix(0, nrow = length(sim_data), ncol = 8)
for(ii in 1:nrow(cover)){
  power[ii, 1:3] <- colMeans(sim_data[[ii]][, c("our_p_val_pow", "lrt_p_val_pow",
                                                    "wald_p_val_pow")] < 0.05)
  power[ii, 4:6] <- apply(sim_data[[ii]][, c("our_p_val_pow", "lrt_p_val_pow",
                                             "wald_p_val_pow")] < 0.05,
                          2, sd) / sqrt(n_sims)
}
power[, 7] <- rep(n_vec, each = length(gamma))
power[, 8] <- rep(gamma, length(n_vec))
colnames(power) <- c("our", "lrt", "wald", "se_our", "se_lrt", "se_wald",
                     "n", "gamma")

if(PDF) pdf(paste(out_dir, "fig_power_R1.pdf", sep = ""), width = 10, height = 4)
par(cex.axis = 1.1, cex.lab = 1.1, cex = 1.1)
par(mgp=c(2.2,0.45,0), tcl=-0.4, mar=c(3.3,3.6,1.5,1.1))
par(mfrow = c(1, 2))

for(jj in c(20, 80)){
  plot_dat <- power[power[, "n"] == jj, ]
  plot(x = gamma,
       y = plot_dat[, "our"],
       type = "l",
       ylim = c(0, 1),
       lwd = 2,
       ylab = "estimated rejection probability",
       xlab = expression(lambda),
       main = paste0("N = ", jj),
       col = cbbPalette[1])
  arrows(gamma,
         plot_dat[, "our"],
         gamma,
         plot_dat[, "our"] + 2 * plot_dat[, "se_our"],
         length = 0.05,
         angle = 90,
         col = cbbPalette[1])
  arrows(gamma,
         plot_dat[, "our"],
         gamma,
         plot_dat[, "our"] - 2 * plot_dat[, "se_our"],
         length = 0.05,
         angle = 90,
         col = cbbPalette[1])
  abline(h = 0.05, lwd = 1)
  lines(x = gamma,
        y = plot_dat[, "lrt"],
        lty = 2,
        lwd = 2,
        col = cbbPalette[2])
  arrows(gamma,
         plot_dat[, "lrt"],
         gamma,
         plot_dat[, "lrt"] + 2 * plot_dat[, "se_lrt"],
         length = 0.05,
         angle = 90,
         col = cbbPalette[2])
  arrows(gamma,
         plot_dat[, "lrt"],
         gamma,
         plot_dat[, "lrt"] - 2 * plot_dat[, "se_lrt"],
         length = 0.05,
         angle = 90,
         col = cbbPalette[2])
  lines(x = gamma,
        y = plot_dat[, "wald"],
        lty = 3,
        lwd = 2,
        col = cbbPalette[3])
  arrows(gamma,
         plot_dat[, "wald"],
         gamma,
         plot_dat[, "wald"] + 2 * plot_dat[, "se_wald"],
         length = 0.05,
         angle = 90,
         col = cbbPalette[3])
  arrows(gamma,
         plot_dat[, "wald"],
         gamma,
         plot_dat[, "wald"] - 2 * plot_dat[, "se_wald"],
         length = 0.05,
         angle = 90,
         col = cbbPalette[3])
}
if(PDF) dev.off()
# Data example ----------------------------------------------------------------
ex_dat <- readRDS("~/GitHub/conf-crit/data_ex/data_ex_R1.Rds")
if(PDF) pdf(paste(out_dir, "fig_data_ex_R1.pdf", sep = ""), width = 9, height = 3)
par(cex.axis = 1.6, cex.lab = 1.6, cex = 1.6)
par(mgp=c(3,1,0), tcl=-0.4, mar=c(4,5,1.5,1.1))
par(mfrow = c(1, 3))

# Confidence interval for first component
plot(x = ex_dat[[1]], y = ex_dat[[4]], type = "l", xlab = bquote(~ lambda[1]),
     ylab = "test statistic", ylim = c(0, 8))
abline(h = qchisq(0.95, 1), lty = 2)

plot(x = ex_dat[[2]], y = ex_dat[[5]], type = "l", xlab = bquote(~ lambda[2]),
     ylab = "test statistic", ylim = c(0, 8))
abline(h = qchisq(0.95, 1), lty = 2)

contour(x = ex_dat[[1]],
        y = ex_dat[[2]],
        z = ex_dat[[6]],
        nlevels = 2,
        levels = c(0, round(qchisq(c(0.8, 0.9, 0.95, 0.99), df = 2), digits = 2)),
        labels = c(0, 0.8, 0.9, 0.95, 0.99),
        labcex = 1,
        xlab = bquote(~ lambda[1]),
        ylab = bquote(~ lambda[2]))
#VC <- as.data.frame(lme4::VarCorr(ex_dat[[6]]))
abline(v = VC[1, 5], lty = 2)
abline(h = VC[2, 5], lty = 2)
dev.off()

# Supplementary Material -------------------------------------------------------
### Compare vover figures with psi known ####
sim_data <- readRDS("~/GitHub/conf-crit/sims/full_simsR1_1e4.Rds")
n_sims <- nrow(sim_data[[1]])

# Verify this using names(sim_data)
gamma <- c(1e-12, 0.01, 0.05, seq(0.1, 0.5, length.out = 5))
n_vec <- c(80, 20)

# Cover curves
cover <- matrix(0, nrow = length(sim_data), ncol = 6)
for(ii in 1:nrow(cover)){
  cover[ii, 1:2] <- 1 - colMeans(sim_data[[ii]][, c("our_p_val_cov",
                                                    "our_p_val_true")] < 0.05)
  cover[ii, 3:4] <- apply(sim_data[[ii]][, c("our_p_val_cov",
                                             "our_p_val_true")] < 0.05,
                          2, sd) / sqrt(n_sims)
}
cover[, 5] <- rep(n_vec, each = length(gamma))
cover[, 6] <- rep(gamma, length(n_vec))
colnames(cover) <- c("our", "our_true", "se_our", "se_our_true",
                     "n", "gamma")

if(PDF) pdf(paste(out_dir, "fig_cover_true.pdf", sep = ""), width = 10, height = 4)
par(cex.axis = 1.1, cex.lab = 1.1, cex = 1.1)
par(mgp=c(2.2,0.45,0), tcl=-0.4, mar=c(3.3,3.6,1.5,1.1))
par(mfrow = c(1, 2))

for(jj in c(20, 80)){
  plot_dat <- cover[cover[, "n"] == jj, ]
  plot(x = gamma,
       y = plot_dat[, "our"],
       type = "l",
       ylim = c(0.87, 1),
       lwd = 2,
       ylab = "estimated cover probability",
       xlab = expression(lambda),
       main = paste0("N = ", jj),
       col = cbbPalette[1])
  arrows(gamma,
         plot_dat[, "our"],
         gamma,
         plot_dat[, "our"] + 2 * plot_dat[, "se_our"],
         length = 0.05,
         angle = 90,
         col = cbbPalette[1])
  arrows(gamma,
         plot_dat[, "our"],
         gamma,
         plot_dat[, "our"] - 2 * plot_dat[, "se_our"],
         length = 0.05,
         angle = 90,
         col = cbbPalette[1])
  abline(h = 0.95, lwd = 1)
  lines(x = gamma,
        y = plot_dat[, "our_true"],
        lty = 2,
        lwd = 2,
        col = cbbPalette[2])
  arrows(gamma,
         plot_dat[, "our_true"],
         gamma,
         plot_dat[, "our_true"] + 2 * plot_dat[, "se_our_true"],
         length = 0.05,
         angle = 90,
         col = cbbPalette[2])
  arrows(gamma,
         plot_dat[, "our_true"],
         gamma,
         plot_dat[, "our_true"] - 2 * plot_dat[, "se_our_true"],
         length = 0.05,
         angle = 90,
         col = cbbPalette[2])
}
if(PDF) dev.off()
