###############################################################################
#  Stochastic Lotka–Volterra SDE fitting + diagnostics for lynx–hare
#  – Reads Lynx–Hare CSV (user must download from repo)
#  – Stage 1: drift estimation via linear regression
#  – Stage 2: diffusion estimation via residual variance
#  – Outputs: time‐series & simulation plots, fit results CSV
#
#  Requirements: readr   dplyr   ggplot2   tibble
#  Author:       Arjun Nair        Date: 2024-05-27
###############################################################################

# 1) Load / Install Required Packages ─────────────────────────────────────────
needed <- c("readr", "dplyr", "ggplot2", "tibble")
for (pkg in needed) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}
library(readr)
library(dplyr)
library(ggplot2)
library(tibble)

# 2) Read & preprocess data ───────────────────────────────────────────────────
#    Make sure “Leigh1968_harelynx.csv” is downloaded into your working directory
data <- read_csv("Leigh1968_harelynx.csv")

# 2.1 Filter non-positive counts then log-transform
data2 <- data %>%
  filter(hare > 0, lynx > 0) %>%
  mutate(
    log_hare = log(hare),
    log_lynx = log(lynx)
  )

# 2.2 Compute one-step log‐increments and rescale
data3_log2 <- data2 %>%
  mutate(
    log_hare_next = lead(log_hare),
    log_lynx_next = lead(log_lynx)
  ) %>%
  filter(!is.na(log_hare_next)) %>%
  mutate(
    delta_log_hare = log_hare_next - log_hare,
    delta_log_lynx  = log_lynx_next - log_lynx,
    hare_s = hare / 1e4,
    lynx_s = lynx / 1e4
  )

# 3) Drift estimation (Stage 1) ────────────────────────────────────────────────
fit_log_hare2 <- lm(delta_log_hare ~ 1 + lynx_s, data = data3_log2)
a_hat <-  coef(fit_log_hare2)["(Intercept)"]
b_hat <- -coef(fit_log_hare2)["lynx_s"]

fit_log_lynx2 <- lm(delta_log_lynx ~ 1 + hare_s, data = data3_log2)
c_hat <- -coef(fit_log_lynx2)["(Intercept)"]
d_hat <-  coef(fit_log_lynx2)["hare_s"]

# 4) Diffusion estimation (Stage 2) ───────────────────────────────────────────
resid_hare <- data3_log2$delta_log_hare -
  (coef(fit_log_hare2)["(Intercept)"] +
     coef(fit_log_hare2)["lynx_s"] * data3_log2$lynx_s)
resid_lynx <- data3_log2$delta_log_lynx -
  (coef(fit_log_lynx2)["(Intercept)"] +
     coef(fit_log_lynx2)["hare_s"] * data3_log2$hare_s)

sigma1_hat <- sqrt(mean(resid_hare^2))
sigma2_hat <- sqrt(mean(resid_lynx^2))

# Report drift & diffusion estimates
cat("Drift estimates:\n")
cat(" a =", round(a_hat,4), " b =", signif(b_hat,3), "\n")
cat(" c =", round(c_hat,4), " d =", signif(d_hat,3), "\n")
cat("Noise estimates:\n")
cat(" sigma1 =", round(sigma1_hat,4), " sigma2 =", round(sigma2_hat,4), "\n\n")

# 5) Simulate sLV trajectories & diagnostics ─────────────────────────────────
set.seed(123)
a <- a_hat; b <- b_hat; c <- c_hat; d <- d_hat
s1 <- sigma1_hat; s2 <- sigma2_hat

n    <- nrow(data3_log2)
M    <- 20
simH <- matrix(NA, n+1, M)
simL <- matrix(NA, n+1, M)
simH[1,] <- data2$log_hare[1]
simL[1,] <- data2$log_lynx[1]

for (j in 1:M) {
  for (k in 1:n) {
    Hk <- simH[k,j]; Lk <- simL[k,j]
    dH <- (a - b*Lk) + s1 * rnorm(1)
    dL <- (-c + d*Hk) + s2 * rnorm(1)
    simH[k+1,j] <- Hk + dH
    simL[k+1,j] <- Lk + dL
  }
}

sim_df <- tibble(
  year         = rep(data2$year[1:(n+1)], times = M),
  sim_log_hare = as.vector(simH),
  sim_log_lynx = as.vector(simL),
  traj         = rep(1:M, each = n+1)
)

# Simulation plot
p_sim <- ggplot() +
  geom_line(data = sim_df,
            aes(year, sim_log_hare, group = traj),
            colour = "forestgreen", alpha = .3) +
  geom_line(data = data2,
            aes(year, log_hare),
            colour = "darkgreen", size = 1) +
  geom_line(data = sim_df,
            aes(year, sim_log_lynx, group = traj),
            colour = "firebrick", alpha = .3) +
  geom_line(data = data2,
            aes(year, log_lynx),
            colour = "darkred", size = 1) +
  scale_x_continuous(limits = c(min(data2$year), max(data2$year)),
                     breaks = seq(min(data2$year), max(data2$year), by = 10),
                     expand = c(0,0)) +
  labs(x = "Year", y = "Log Count",
       title = "sLV Simulations (light) vs. Data (dark)") +
  theme_minimal()

print(p_sim)

# 6) Export plots & results table ────────────────────────────────────────────
# Explicit time‐series plot
p_ts <- ggplot(data2, aes(year)) +
  geom_line(aes(y = log_hare), colour = "darkgreen") +
  geom_line(aes(y = log_lynx), colour = "darkred") +
  labs(x = "Year", y = "Log Count", title = "Lynx–Hare Time Series") +
  theme_minimal()

ggsave("lynx_hare_timeseries.png", p_ts,  width = 6, height = 4, dpi = 300)
ggsave("sLV_simulations.png",       p_sim, width = 6, height = 4, dpi = 300)

results <- tibble(
  parameter = c("a","b","c","d","sigma1","sigma2",
                "rmse_hare","rmse_lynx",
                "obs_cycle_mean","obs_cycle_sd",
                "sim_cycle_mean","sim_cycle_sd"),
  estimate  = c(a_hat, b_hat, c_hat, d_hat, sigma1_hat, sigma2_hat,
                sqrt(mean((sim_df %>%
                             group_by(year) %>%
                             summarize(mean_log_hare = mean(sim_log_hare)) %>%
                             left_join(data2 %>% select(year, log_hare), by="year") %>%
                             pull(mean_log_hare) - data2$log_hare)^2)),
                sqrt(mean((sim_df %>%
                             group_by(year) %>%
                             summarize(mean_log_lynx = mean(sim_log_lynx)) %>%
                             left_join(data2 %>% select(year, log_lynx), by="year") %>%
                             pull(mean_log_lynx) - data2$log_lynx)^2)),
                mean(diff(find_peaks(data2$log_hare, data2$year))),
                sd(diff(find_peaks(data2$log_hare, data2$year))),
                mean(unlist(lapply(split(sim_df, sim_df$traj),
                                   function(df) diff(find_peaks(df$sim_log_hare, df$year))))),
                sd(unlist(lapply(split(sim_df, sim_df$traj),
                                 function(df) diff(find_peaks(df$sim_log_hare, df$year))))))
)
write_csv(results, "sLV_fit_results.csv")

message("Plots and results saved to working directory.")
