# ─────────────────────────────────────────────────────────────────────────────
# 1) Load / Install Required Packages
# ─────────────────────────────────────────────────────────────────────────────
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

# ─────────────────────────────────────────────────────────────────────────────
# 2) Read and preprocess data
# ─────────────────────────────────────────────────────────────────────────────
# Ensure the file Leigh1968_harelynx.csv is downloaded from the repo into the working directory
data <- read_csv("Leigh1968_harelynx.csv")

# 2.1 Filter out non-positive counts then Log-transform
data2 <- data %>%
  filter(hare > 0, lynx > 0) %>%
  mutate(
    log_hare = log(hare),
    log_lynx = log(lynx)
  )

# 2.2 Compute one-step log-increments and rescale
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

# ─────────────────────────────────────────────────────────────────────────────
# 3) Drift estimation (Stage 1)
# ─────────────────────────────────────────────────────────────────────────────
fit_log_hare2 <- lm(delta_log_hare ~ 1 + lynx_s, data = data3_log2)
a_hat <-  coef(fit_log_hare2)["(Intercept)"]
b_hat <- -coef(fit_log_hare2)["lynx_s"]

fit_log_lynx2 <- lm(delta_log_lynx ~ 1 + hare_s, data = data3_log2)
c_hat <- -coef(fit_log_lynx2)["(Intercept)"]
d_hat <-  coef(fit_log_lynx2)["hare_s"]

# ─────────────────────────────────────────────────────────────────────────────
# 4) Diffusion estimation (Stage 2)
# ─────────────────────────────────────────────────────────────────────────────
resid_hare <- data3_log2$delta_log_hare -
  ( coef(fit_log_hare2)["(Intercept)"]
    + coef(fit_log_hare2)["lynx_s"] * data3_log2$lynx_s )
resid_lynx <- data3_log2$delta_log_lynx -
  ( coef(fit_log_lynx2)["(Intercept)"]
    + coef(fit_log_lynx2)["hare_s"] * data3_log2$hare_s )

sigma1_hat <- sqrt(mean(resid_hare^2))
sigma2_hat <- sqrt(mean(resid_lynx^2))

# Report drift & diffusion estimates
cat("Drift estimates:\n")
cat(" a =", round(a_hat,4),   " b =", signif(b_hat,3), "\n")
cat(" c =", round(c_hat,4),   " d =", signif(d_hat,3), "\n")
cat("Noise estimates:\n")
cat(" sigma1 =", round(sigma1_hat,4), " sigma2 =", round(sigma2_hat,4), "\n\n")

# ─────────────────────────────────────────────────────────────────────────────
# 5) Simulate sLV trajectories and diagnostics
# ─────────────────────────────────────────────────────────────────────────────
set.seed(123)
a <- a_hat; b <- b_hat; c <- c_hat; d <- d_hat
s1 <- sigma1_hat; s2 <- sigma2_hat

T   <- nrow(data3_log2)
M   <- 20
sim_logH <- matrix(NA, T+1, M)
sim_logL <- matrix(NA, T+1, M)
sim_logH[1,] <- data2$log_hare[1]
sim_logL[1,] <- data2$log_lynx[1]

for (j in 1:M) {
  for (k in 1:T) {
    Hk <- sim_logH[k,j]; Lk <- sim_logL[k,j]
    dH <- (a - b*Lk)*1 + s1*rnorm(1)
    dL <- (-c + d*Hk)*1 + s2*rnorm(1)
    sim_logH[k+1,j] <- Hk + dH
    sim_logL[k+1,j] <- Lk + dL
  }
}

sim_df <- tibble(
  year         = rep(data2$year[1:(T+1)], times = M),
  sim_log_hare = as.vector(sim_logH),
  sim_log_lynx = as.vector(sim_logL),
  traj         = rep(1:M, each = T+1)
)

# Simulation plot
p_sim <- ggplot() +
  geom_line(data = sim_df,
            aes(year, sim_log_hare, group = traj),
            color = "forestgreen", alpha = .3) +
  geom_line(data = data2,
            aes(year, log_hare),
            color = "darkgreen", size = 1) +
  geom_line(data = sim_df,
            aes(year, sim_log_lynx, group = traj),
            color = "firebrick", alpha = .3) +
  geom_line(data = data2,
            aes(year, log_lynx),
            color = "darkred", size = 1) +
  scale_x_continuous(limits = c(min(data2$year), max(data2$year)),
                     breaks = seq(min(data2$year), max(data2$year), by = 10),
                     expand = c(0,0)) +
  labs(x="Year", y="Log Count",
       title="sLV Simulations (light) vs. Data (dark)") +
  theme_minimal()

print(p_sim)

# Diagnostics: RMSE
sim_mean <- sim_df %>%
  group_by(year) %>%
  summarize(mean_log_hare = mean(sim_log_hare),
            mean_log_lynx = mean(sim_log_lynx)) %>%
  left_join(data2 %>% select(year, log_hare, log_lynx), by="year")
rmse_hare <- sqrt(mean((sim_mean$mean_log_hare - sim_mean$log_hare)^2))
rmse_lynx <- sqrt(mean((sim_mean$mean_log_lynx - sim_mean$log_lynx)^2))
cat("RMSE: hare =", round(rmse_hare,3),
    " lynx =", round(rmse_lynx,3), "\n\n")

# Helper for cycle lengths
find_peaks <- function(x,y) {
  idx <- which(x[-c(1,length(x))] > x[-c(length(x)-1,length(x))] &
                 x[-c(1,length(x))] > x[-c(1,2)]) + 1
  y[idx]
}
obs_peaks   <- find_peaks(data2$log_hare, data2$year)
obs_periods <- diff(obs_peaks)
sim_periods <- unlist(lapply(split(sim_df, sim_df$traj),
                             function(df) diff(find_peaks(df$sim_log_hare, df$year))))
cat("Observed cycle mean =", round(mean(obs_periods),1),
    "sd =", round(sd(obs_periods),1), "\n")
cat("Simulated cycle mean =", round(mean(sim_periods),1),
    "sd =", round(sd(sim_periods),1), "\n\n")

# ─────────────────────────────────────────────────────────────────────────────
# 6) Export plots and results table
# ─────────────────────────────────────────────────────────────────────────────

p_ts <- ggplot(data2, aes(year)) +
  geom_line(aes(y = log_hare), color = "darkgreen") +
  geom_line(aes(y = log_lynx), color = "darkred") +
  labs(x = "Year", y = "Log Count", title = "Lynx–Hare Time Series") +
  theme_minimal()

# Save plots
ggsave("lynx_hare_timeseries.png", p_ts,  width = 6, height = 4, dpi = 300)
ggsave("sLV_simulations.png",    p_sim, width = 6, height = 4, dpi = 300)

# Assemble results
results <- tibble(
  parameter       = c("a","b","c","d","sigma1","sigma2",
                      "rmse_hare","rmse_lynx",
                      "obs_cycle_mean","obs_cycle_sd",
                      "sim_cycle_mean","sim_cycle_sd"),
  estimate        = c(a_hat, b_hat, c_hat, d_hat, sigma1_hat, sigma2_hat,
                      rmse_hare, rmse_lynx,
                      mean(obs_periods), sd(obs_periods),
                      mean(sim_periods), sd(sim_periods))
)
write_csv(results, "sLV_fit_results.csv")

message("Plots and results saved to working directory.")
