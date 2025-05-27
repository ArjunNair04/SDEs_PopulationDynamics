###############################################################################
#  Stochastic–Logistic SDE fitting + diagnostics for three media
#  – Reads one worksheet per medium (user-specified below)
#  – Fits r & K (via NLS) and σ (method of moments)
#  – Saves drift, QQ and EM-simulations overlays
#
#  Requirements: readxl   data.table   ggplot2
#  Author:       Arjun Nair         Date: 2024-05-25
#
#  * Ensure the files “KHK growth curves_*.xlsx” are downloaded from the
#    repository into your working directory before running this script.
###############################################################################

# 0) Prepare output directory ────────────────────────────────────────────────
if (!dir.exists("plots")) dir.create("plots")

# 1) Load / Install Required Packages ─────────────────────────────────────────
pkgs <- c("readxl", "data.table", "ggplot2")
new  <- pkgs[!pkgs %in% installed.packages()[,"Package"]]
if (length(new)) install.packages(new, quiet = TRUE)
lapply(pkgs, library, character.only = TRUE)

# 2) User choices ─────────────────────────────────────────────────────────────
TARGETS <- data.table(
  medium = c("LB", "MAA", "M63"),
  file   = c("KHK growth curves_LB.xlsx",
             "KHK growth curves_MAA.xlsx",
             "KHK growth curves_M63.xlsx"),
  sheet  = c(32, 32, 32)
)

# 3) Helpers ─────────────────────────────────────────────────────────────────
## 3.1 Read & average logOD
read_mean_logOD <- function(xlsx, sht) {
  dt <- as.data.table(read_excel(xlsx, sheet = sht))
  setnames(dt, 1, "Time")
  dt <- melt(dt, id.vars = "Time", variable.name = "Well", value.name = "OD")[OD > 0]
  dt[, .(logOD = mean(log(OD))), by = Time]
}

## 3.2 Fit r, K, σ
fit_logistic <- function(mt) {
  df    <- as.data.frame(mt)
  df$x  <- exp(df$logOD)
  df$dt <- c(NA, diff(df$Time))
  df$der<- c(NA, diff(df$x)) / df$dt
  
  fit <- nls(
    der ~ r * x * (1 - x / K),
    data    = df[df$dt > 0, ],
    start   = list(r = 1, K = max(df$x)),
    control = nls.control(maxiter = 100, warnOnly = TRUE)
  )
  pars <- coef(fit)
  
  drift <- pars["r"] * df$x * (1 - df$x / pars["K"])
  res   <- diff(df$x) - drift[-length(drift)] * df$dt[-1]
  sigma <- sqrt(sum(res^2) / sum(df$x[-1]^2 * df$dt[-1]))
  
  list(r = pars["r"], K = pars["K"], sigma = sigma, df = df)
}

## 3.3 EM sampler
em_paths <- function(t, x0, r, K, sig, M = 20) {
  dt <- diff(t); n <- length(t)
  X  <- matrix(NA, n, M); X[1, ] <- x0
  for (k in 1:(n - 1)) {
    dW <- rnorm(M, 0, sqrt(dt[k]))
    X[k+1, ] <- X[k, ] + r * X[k, ] * (1 - X[k, ] / K) * dt[k] +
      sig * X[k, ] * dW
  }
  X
}

## 3.4 Plotting
plot_drift <- function(df, r, K, tag) {
  drift <- r * df$x * (1 - df$x / K)
  p <- ggplot(data.frame(
    x   = df$x[-nrow(df)],
    emp = diff(df$x)      / df$dt[-1],
    fit = drift[-length(drift)]
  ), aes(x)) +
    geom_point(aes(y = emp), colour = "grey40", size = 1.6) +
    geom_line(aes(y = fit), colour = "steelblue", size = 1) +
    labs(title = sprintf("Drift–Fit (%s)", tag),
         x = "Abundance x", y = expression(dX/dt)) +
    theme_minimal(base_size = 13)
  print(p); p
}

plot_qq <- function(df, r, K, tag) {
  drift <- r * df$x * (1 - df$x / K)
  std   <- (diff(df$x) - drift[-length(drift)] * df$dt[-1]) /
    (df$x[-1] * sqrt(df$dt[-1]))
  p <- ggplot(data.frame(std = std), aes(sample = std)) +
    stat_qq(size = 1.6, colour = "grey30") +
    stat_qq_line(colour = "steelblue") +
    labs(title = sprintf("QQ-Plot (%s)", tag)) +
    theme_minimal(base_size = 13)
  print(p); p
}

plot_em <- function(df, r, K, sig, tag) {
  paths <- em_paths(df$Time, df$x[1], r, K, sig, M = 20)
  sim   <- data.frame(
    Time = rep(df$Time, ncol(paths)),
    X    = c(paths),
    Run  = factor(rep(1:ncol(paths), each = nrow(paths)))
  )
  p <- ggplot() +
    geom_line(data = sim, aes(Time, X, group = Run),
              colour = "grey80", alpha = 0.6) +
    geom_line(data = df,  aes(Time, x), colour = "steelblue", size = 1.2) +
    labs(title = sprintf("EM paths (%s)", tag),
         x = "Time (h)", y = "Abundance") +
    theme_minimal(base_size = 13)
  print(p); p
}

# 4) Main loop ────────────────────────────────────────────────────────────────
for (i in seq_len(nrow(TARGETS))) {
  med  <- TARGETS[i, medium]
  file <- TARGETS[i, file]
  sht  <- TARGETS[i, sheet]
  tag  <- sprintf("%s, strain %s", med, sht)
  
  cat("\n---", tag, "---\n")
  mean_tab <- read_mean_logOD(file, sht)
  fit      <- fit_logistic(mean_tab)
  cat("  r =", round(fit$r,4),
      "  K =", round(fit$K,4),
      "  σ =", round(fit$sigma,5), "\n")
  
  # 4.1 Drift diagnostic
  pd <- plot_drift(fit$df, fit$r, fit$K, tag)
  ggsave(sprintf("plots/drift_%s.png", med), pd,
         width = 5, height = 4, dpi = 300)
  
  # 4.2 QQ diagnostic
  pq <- plot_qq(fit$df, fit$r, fit$K, tag)
  ggsave(sprintf("plots/qq_%s.png", med), pq,
         width = 5, height = 4, dpi = 300)
  
  # 4.3 EM sample paths
  pe <- plot_em(fit$df, fit$r, fit$K, fit$sigma, tag)
  ggsave(sprintf("plots/em_%s.png", med), pe,
         width = 6, height = 4, dpi = 300)
}

###############################################################################
#  End of script
###############################################################################
