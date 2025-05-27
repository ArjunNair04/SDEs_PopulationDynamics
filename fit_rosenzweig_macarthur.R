###############################################################################
#  Stochastic Rosenzweig–MacArthur SDE fitting + diagnostics
#  – Reads yearly lynx–hare counts
#  – Stage 1: OLS drift (fixed K,h,α from literature)
#  – Stage 2: variance‐ratio diffusion
#  – EM‐simulation overlay & diagnostics
#
#  Requirements: readr   dplyr   ggplot2   tibble
#  Author:       Arjun Nair         Date: 2024-05-27
###############################################################################

# 1) Load / Install Required Packages ─────────────────────────────────────────
pkgs <- c("readr", "dplyr", "ggplot2", "tibble")
new  <- pkgs[!pkgs %in% installed.packages()[,"Package"]]
if (length(new)) install.packages(new, quietly = TRUE)
invisible(lapply(pkgs, library, character.only = TRUE))

# 2) Data import & preprocessing ──────────────────────────────────────────────
#    Make sure “Leigh1968_harelynx.csv” is in your working directory
data <- read_csv("Leigh1968_harelynx.csv", show_col_types = FALSE)
data2 <- data %>%
  filter(hare > 0, lynx > 0) %>%             # avoid log(0)
  mutate(
    log_hare = log(hare),
    log_lynx = log(lynx),
    log_hare_next = lead(log_hare),
    log_lynx_next = lead(log_lynx)
  ) %>%
  filter(!is.na(log_hare_next)) %>%
  mutate(
    dlog_hare = log_hare_next - log_hare,
    dlog_lynx = log_lynx_next - log_lynx
  )
T_steps <- nrow(data2)
cat("N obs =", T_steps,
    "| Δlog_hare in", round(range(data2$dlog_hare), 2),
    "| Δlog_lynx in", round(range(data2$dlog_lynx), 2), "\n\n")

# 3) Fixed Rosenzweig–MacArthur constants ─────────────────────────────────────
scale_fac <- 1e4           # 10^4 pelts
K_fix     <- 6             # 6 × 10^4 pelts
h_fix     <- 0.0022        # yr/prey
alpha_fix <- 0.02          # dimensionless

# 4) Stage 1 – OLS drift estimates ────────────────────────────────────────────
data4 <- data2 %>%
  mutate(
    H      = hare / scale_fac,
    P      = lynx / scale_fac,
    logis_H= 1 - H / K_fix,
    fr_H   = H / (1 + h_fix * H)
  )

fit_H <- lm(dlog_hare ~ 1 + logis_H + I(fr_H * P), data = data4)
r_hat     <- coef(fit_H)["logis_H"]
alpha_hat <- coef(fit_H)["I(fr_H * P)"]

fit_P <- lm(dlog_lynx ~ 1 + I(fr_H * H), data = data4)
beta_alpha_hat <- coef(fit_P)["I(fr_H * H)"]
m_hat          <- -coef(fit_P)["(Intercept)"]
beta_hat       <- beta_alpha_hat / alpha_hat

# 5) Stage 2 – variance‐ratio diffusion ───────────────────────────────────────
sigma_H_hat <- sqrt(mean(resid(fit_H)^2))
sigma_P_hat <- sqrt(mean(resid(fit_P)^2))

# 6) Parameter table ─────────────────────────────────────────────────────────
rma_params <- tibble(
  Parameter = c("r", "K (fixed)", "α",   "β",   "h (fixed)", "m",
                "σ_H", "σ_P"),
  Estimate  = c(r_hat, K_fix,    alpha_hat, beta_hat,
                h_fix, m_hat,    sigma_H_hat, sigma_P_hat),
  Units     = c("per yr", "×10⁴ pelts", "dimless", "dimless",
                "yr/prey", "per yr", "per √yr", "per √yr")
)
print(rma_params, n = Inf)

# 7) Simulation overlay ───────────────────────────────────────────────────────
set.seed(2025)
M <- 20
logH_mat <- matrix(NA_real_, T_steps + 1, M)
logP_mat <- matrix(NA_real_, T_steps + 1, M)
logH_mat[1, ] <- data2$log_hare[1]
logP_mat[1, ] <- data2$log_lynx[1]

for (j in seq_len(M)) {
  for (k in seq_len(T_steps)) {
    Hk <- exp(logH_mat[k, j]) / scale_fac
    Pk <- exp(logP_mat[k, j]) / scale_fac
    
    mu_H <- r_hat * (1 - Hk / K_fix) -
      alpha_hat * Hk / (1 + h_fix * Hk) * Pk
    mu_P <- beta_hat * alpha_hat * Hk / (1 + h_fix * Hk) - m_hat
    
    logH_mat[k + 1, j] <- logH_mat[k, j] +
      mu_H + sigma_H_hat * rnorm(1)
    logP_mat[k + 1, j] <- logP_mat[k, j] +
      mu_P + sigma_P_hat * rnorm(1)
  }
}

sim_df <- tibble(
  year = rep(data$year[1:(T_steps + 1)], times = M),
  logH = as.vector(logH_mat),
  logP = as.vector(logP_mat),
  path = rep(seq_len(M), each = T_steps + 1)
)

# assign the ggplot to a variable so ggsave() has an explicit object
p_rma <- ggplot() +
  geom_line(aes(year, logH, group = path), data = sim_df,
            colour = "forestgreen", alpha = 0.25) +
  geom_line(aes(year, logP, group = path), data = sim_df,
            colour = "firebrick", alpha = 0.25) +
  geom_line(aes(year, log_hare), data = data2,
            colour = "darkgreen", size = 1.1) +
  geom_line(aes(year, log_lynx), data = data2,
            colour = "darkred",   size = 1.1) +
  labs(title = "Stochastic Rosenzweig–MacArthur: simulations (pale) vs data",
       x = "Year", y = "Log(count)") +
  scale_x_continuous(breaks = seq(min(data2$year), max(data2$year), by = 10),
                     expand = c(0, 0)) +
  theme_minimal()

print(p_rma)

# 8) Diagnostics & save outputs ───────────────────────────────────────────────
sim_mean <- sim_df %>%
  group_by(year) %>%
  summarize(mean_logH = mean(logH),
            mean_logP = mean(logP),
            .groups = "drop") %>%
  left_join(data2 %>% select(year, log_hare, log_lynx), by = "year")

rmse_H <- sqrt(mean((sim_mean$mean_logH - sim_mean$log_hare)^2))
rmse_P <- sqrt(mean((sim_mean$mean_logP - sim_mean$log_lynx)^2))

cat("RMSE of mean simulated vs observed:\n")
cat("  Hare RMSE =", round(rmse_H,4), "log-units\n")
cat("  Lynx RMSE =", round(rmse_P,4), "log-units\n\n")

find_peaks <- function(x, yrs) {
  idx <- which(
    x[-c(1,length(x))] > x[-c(length(x)-1,length(x))] &
      x[-c(1,length(x))] > x[-c(1,2)]
  ) + 1
  yrs[idx]
}

obs_peaks   <- find_peaks(data2$log_hare, data2$year)
obs_periods <- diff(obs_peaks)
all_sim_periods <- unlist(
  lapply(split(sim_df, sim_df$path), function(df)
    diff(find_peaks(df$logH, df$year)))
)

cat("Observed hare cycle lengths (years):\n")
print(obs_periods)
cat("  Mean =", round(mean(obs_periods),1),
    "SD =", round(sd(obs_periods),1), "\n\n")
cat("Simulated hare cycle lengths (all paths):\n")
cat("  Mean =", round(mean(all_sim_periods),1),
    "SD =", round(sd(all_sim_periods),1),
    "(n =", length(all_sim_periods), ")\n\n")

write_csv(rma_params,      "rma_fit_results.csv")
write_csv(sim_mean,        "rma_sim_mean_vs_data.csv")
write_csv(
  tibble(obs_periods = obs_periods,
         sim_periods_all = all_sim_periods[1:length(obs_periods)]),
  "rma_cycle_lengths.csv"
)

# fatal flaw fixed: ggsave now references an explicit plot object
ggsave("rma_simulations.png", p_rma, width = 8, height = 4, dpi = 300)
