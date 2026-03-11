rm(list = ls())
graphics.off()

# File paths

gee_file <- "YOUR PATH HERE"
q_file   <- "YOUR PATH HERE"

# Plot colors

COL_OBS <- "black"
COL_RAW <- "red"
COL_CAL <- "blue"

# Read source tables
met  <- read.csv(gee_file, stringsAsFactors = FALSE)
qobs <- read.csv(q_file, stringsAsFactors = FALSE)

# build standardized columns from the raw file headers
met$Date  <- as.Date(met$date)
qobs$Date <- as.Date(qobs$time)
qobs$Q    <- as.numeric(qobs$value)

# keep only needed columns
keep_met <- c("Date", "tmean_c", "prcp_mm", "srad_wm2", "swe_mm", "swe_b1_mm", "swe_b2_mm", "swe_b3_mm")
met  <- met[, keep_met]
qobs <- qobs[, c("Date", "Q")]

# drop bad dates if any
met  <- met[!is.na(met$Date), ]
qobs <- qobs[!is.na(qobs$Date), ]

# merge and sort
q <- merge(met, qobs, by = "Date", all = FALSE)
q <- q[order(q$Date), ]

cat("Met rows:", nrow(met), "\n")
cat("Q rows:", nrow(qobs), "\n")
cat("Merged rows:", nrow(q), "\n")
cat("Met dates not in Q:", sum(!(met$Date %in% qobs$Date)), "\n")
cat("Q dates not in met:", sum(!(qobs$Date %in% met$Date)), "\n")

# Coerce numeric columns

numcols <- c("Q", "tmean_c", "prcp_mm", "srad_wm2", "swe_mm", "swe_b1_mm", "swe_b2_mm", "swe_b3_mm")
q[numcols] <- lapply(q[numcols], as.numeric)

# Basic date fields

q$YEAR  <- as.integer(format(q$Date, "%Y"))
q$Month <- as.integer(format(q$Date, "%m"))
q$Day   <- as.integer(format(q$Date, "%d"))
q$DOY   <- as.integer(format(q$Date, "%j"))

# Water year

q$WY <- ifelse(q$Month >= 10, q$YEAR + 1, q$YEAR)
q$WY <- as.integer(q$WY)

# Water-year day

wy_start <- as.Date(paste0(q$WY - 1, "-10-01"))
q$WYDOY <- as.integer(q$Date - wy_start) + 1
q$WYDOY <- ifelse(q$WYDOY >= 1 & q$WYDOY <= 366, q$WYDOY, NA_integer_)

# Seasonal sine term

q$SIN_DOY <- sin(2 * pi * q$DOY / 365.25)

# Positive degree days from mean temperature

q$pdd <- pmax(q$tmean_c, 0)

# Rolling sums

roll_sum <- function(x, k) {
  out <- rep(NA_real_, length(x))
  for (i in seq_along(x)) {
    i0 <- max(1, i - k + 1)
    out[i] <- sum(x[i0:i], na.rm = TRUE)
  }
  out
}

q$pdd_7  <- roll_sum(q$pdd, 7)
q$pdd_30 <- roll_sum(q$pdd, 30)

# Rain / snow split using mean temperature

q$rain_mm <- ifelse(is.finite(q$tmean_c) & q$tmean_c > 0, q$prcp_mm, 0)
q$snow_mm <- ifelse(is.finite(q$tmean_c) & q$tmean_c <= 0, q$prcp_mm, 0)

# Antecedent precipitation index

compute_api <- function(prcp, decay = 0.9) {
  api <- rep(NA_real_, length(prcp))
  api[1] <- ifelse(is.finite(prcp[1]), prcp[1], 0)
  for (i in 2:length(prcp)) {
    p <- ifelse(is.finite(prcp[i]), prcp[i], 0)
    api[i] <- decay * api[i - 1] + p
  }
  api
}

q$api <- compute_api(q$prcp_mm, decay = 0.9)

# Rain-on-snow flag

q$ROSflag <- ifelse(
  is.finite(q$rain_mm) & is.finite(q$swe_mm) &
    q$rain_mm > 0 & q$swe_mm > 1,
  1, 0
)

# 7-day dSWE by band
# Negative values indicate snow loss

lag_k <- function(x, k) {
  c(rep(NA_real_, k), x[1:(length(x) - k)])
}

q$dswe17 <- q$swe_b1_mm - lag_k(q$swe_b1_mm, 7)
q$dswe27 <- q$swe_b2_mm - lag_k(q$swe_b2_mm, 7)
q$dswe37 <- q$swe_b3_mm - lag_k(q$swe_b3_mm, 7)

# Clean copies for modeling

q$dswe17_clean <- q$dswe17
q$dswe27_clean <- q$dswe27
q$dswe37_clean <- q$dswe37

month_num <- as.integer(format(q$Date, "%m"))
winter_period <- !(month_num %in% 3:9)

precount <- sum(!complete.cases(q[, c("dswe17_clean", "dswe27_clean", "dswe37_clean")]))

for (v in c("dswe17_clean", "dswe27_clean", "dswe37_clean")) {
  winter_vals <- q[[v]][winter_period]
  mu  <- mean(winter_vals, na.rm = TRUE)
  sdv <- sd(winter_vals, na.rm = TRUE)
  z <- (q[[v]] - mu) / sdv
  bad <- winter_period & abs(z) > 3
  q[[v]][bad] <- NA
}

postcount <- sum(!complete.cases(q[, c("dswe17_clean", "dswe27_clean", "dswe37_clean")]))

cat("\ndSWE cleaning summary\n")
cat("New rows with any NA introduced by filter:", postcount - precount, "\n")
cat("Added NAs in dswe17_clean:", sum(is.na(q$dswe17_clean)) - sum(is.na(q$dswe17)), "\n")
cat("Added NAs in dswe27_clean:", sum(is.na(q$dswe27_clean)) - sum(is.na(q$dswe27)), "\n")
cat("Added NAs in dswe37_clean:", sum(is.na(q$dswe37_clean)) - sum(is.na(q$dswe37)), "\n")

cat("\nRows in merged dataset:", nrow(q), "\n")
cat("Water years present:\n")
print(sort(unique(q$WY)))

# Metrics

NSE <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  1 - sum((x - y)^2) / sum((x - mean(x))^2)
}

RMSE <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  sqrt(mean((x - y)^2))
}

PBIAS <- function(obs, pred) {
  ok <- is.finite(obs) & is.finite(pred)
  obs <- obs[ok]
  pred <- pred[ok]
  100 * sum(pred - obs) / sum(obs)
}

KGE <- function(obs, pred) {
  ok <- is.finite(obs) & is.finite(pred)
  obs <- obs[ok]
  pred <- pred[ok]
  
  r <- cor(obs, pred)
  alpha <- sd(pred) / sd(obs)
  beta  <- mean(pred) / mean(obs)
  
  1 - sqrt((r - 1)^2 + (alpha - 1)^2 + (beta - 1)^2)
}

# Recompute WYDOY for subsets

add_WYDOY <- function(df) {
  wy_start <- as.Date(paste0(df$WY - 1, "-10-01"))
  df$WYDOY <- as.integer(df$Date - wy_start) + 1
  df$WYDOY <- ifelse(df$WYDOY >= 1 & df$WYDOY <= 366, df$WYDOY, NA_integer_)
  df
}

# Training-only daily Q climatology

compute_doy_mean_q <- function(df_model, df_all) {
  df_model <- add_WYDOY(df_model)
  df_all   <- add_WYDOY(df_all)
  
  clim <- aggregate(Q ~ WYDOY, data = df_model, FUN = mean, na.rm = TRUE)
  names(clim)[2] <- "DOY_mean_Q"
  
  if ("DOY_mean_Q" %in% names(df_all)) df_all$DOY_mean_Q <- NULL
  
  df_all <- merge(df_all, clim, by = "WYDOY", all.x = TRUE, sort = FALSE)
  
  fallback <- mean(df_model$Q, na.rm = TRUE)
  df_all$DOY_mean_Q[is.na(df_all$DOY_mean_Q)] <- fallback
  
  df_all
}

# Holdout runner

run_holdout <- function(test_wy, q_all) {
  
  cat("\n-------------------------------------\n")
  cat("Running holdout for WY", test_wy, "\n")
  cat("-------------------------------------\n")
  
  q_test  <- subset(q_all, WY == test_wy)
  q_train <- subset(q_all, WY != test_wy)
  
  cat("Training rows:", nrow(q_train), "\n")
  cat("Testing rows :", nrow(q_test), "\n")
  
  q_all2 <- compute_doy_mean_q(df_model = q_train, df_all = q_all)
  
  q_test  <- subset(q_all2, WY == test_wy)
  q_train <- subset(q_all2, WY != test_wy)
  
  q_train$logq <- ifelse(is.finite(q_train$Q) & q_train$Q > 0, log(q_train$Q), NA_real_)
  q_test$logq  <- ifelse(is.finite(q_test$Q)  & q_test$Q  > 0, log(q_test$Q),  NA_real_)
  
  e <- lm(
    logq ~ prcp_mm + pdd + api + rain_mm + ROSflag +
      SIN_DOY + srad_wm2 +
      dswe17_clean + dswe27_clean + dswe37_clean +
      DOY_mean_Q +
      pdd:dswe17_clean + pdd:dswe27_clean + pdd:dswe37_clean,
    data = q_train
  )
  
  pred_log_train <- predict(e, newdata = q_train)
  pred_log_test  <- predict(e, newdata = q_test)
  
  pred_Q_train <- exp(pred_log_train)
  pred_Q_test  <- exp(pred_log_test)
  
  ok_train <- is.finite(q_train$Q) & is.finite(pred_Q_train)
  k <- sum(q_train$Q[ok_train]) / sum(pred_Q_train[ok_train])
  
  pred_Q_test_cal   <- k * pred_Q_test
  pred_log_test_cal <- log(pred_Q_test_cal)
  
  keep <- !is.na(q_test$Date) &
    is.finite(q_test$logq) &
    is.finite(q_test$Q) &
    is.finite(pred_log_test) &
    is.finite(pred_log_test_cal)
  
  d <- q_test[keep, ]
  d <- d[order(d$Date), ]
  
  d$pred_logq_raw <- pred_log_test[keep]
  d$pred_Q_raw    <- pred_Q_test[keep]
  d$pred_logq_cal <- pred_log_test_cal[keep]
  d$pred_Q_cal    <- pred_Q_test_cal[keep]
  
  metrics <- data.frame(
    WY = test_wy,
    n_days = nrow(d),
    NSE_log_raw    = NSE(d$logq, d$pred_logq_raw),
    NSE_log_cal    = NSE(d$logq, d$pred_logq_cal),
    RMSE_log_raw   = RMSE(d$logq, d$pred_logq_raw),
    RMSE_log_cal   = RMSE(d$logq, d$pred_logq_cal),
    KGE_Q_raw      = KGE(d$Q, d$pred_Q_raw),
    KGE_Q_cal      = KGE(d$Q, d$pred_Q_cal),
    PBIAS_Q_raw    = PBIAS(d$Q, d$pred_Q_raw),
    PBIAS_Q_cal    = PBIAS(d$Q, d$pred_Q_cal),
    ratio_raw      = sum(d$pred_Q_raw) / sum(d$Q),
    ratio_cal      = sum(d$pred_Q_cal) / sum(d$Q),
    k_volume_scale = k
  )
  
  print(metrics)
  
  list(
    test_wy = test_wy,
    model = e,
    k = k,
    d = d,
    metrics = metrics
  )
}

# Plotting

plot_holdout <- function(res) {
  d  <- res$d
  wy <- res$test_wy
  k  <- res$k
  
  nse_log_cal <- NSE(d$logq, d$pred_logq_cal)
  pbias_cal   <- PBIAS(d$Q, d$pred_Q_cal)
  kge_cal     <- KGE(d$Q, d$pred_Q_cal)
  
  yl <- range(c(d$logq, d$pred_logq_raw, d$pred_logq_cal), na.rm = TRUE)
  
  plot(
    d$Date, d$logq, type = "l", ylim = yl,
    col = COL_OBS,
    xlab = "Date", ylab = "ln(Q)",
    main = paste0("WY", wy, " Hydrograph: Observed vs Predicted (log space)")
  )
  lines(d$Date, d$pred_logq_raw, lty = 2, col = COL_RAW)
  lines(d$Date, d$pred_logq_cal, lty = 3, col = COL_CAL)
  
  legend(
    "topright",
    legend = c("Observed", "Predicted (raw)", "Predicted (cal)"),
    lty = c(1, 2, 3),
    col = c(COL_OBS, COL_RAW, COL_CAL),
    bty = "n"
  )
  
  mtext(
    paste0(
      "Log NSE cal=", round(nse_log_cal, 3),
      " | KGE cal=", round(kge_cal, 3),
      " | PBIAS cal=", round(pbias_cal, 1),
      "% | k=", signif(k, 4)
    ),
    side = 3, line = 0.2
  )
  
  plot(
    d$logq, d$pred_logq_cal,
    xlab = "Observed log(Q)",
    ylab = "Predicted log(Q) (cal)",
    main = paste0("WY", wy, " Observed vs Predicted (log space) - calibrated")
  )
  abline(0, 1, lty = 2)
  
  mtext(
    paste0(
      "Log NSE cal=", round(nse_log_cal, 3),
      " | KGE cal=", round(kge_cal, 3),
      " | PBIAS cal=", round(pbias_cal, 1),
      "% | k=", signif(k, 4)
    ),
    side = 3, line = 0.2
  )
  
  plot(
    d$Date, cumsum(d$Q), type = "l",
    col = COL_OBS,
    xlab = "Date", ylab = "Cumulative Q",
    main = paste0("WY", wy, " Cumulative Volume Check")
  )
  lines(d$Date, cumsum(d$pred_Q_raw), lty = 2, col = COL_RAW)
  lines(d$Date, cumsum(d$pred_Q_cal), lty = 3, col = COL_CAL)
  
  legend(
    "topleft",
    legend = c("Observed", "Predicted (raw)", "Predicted (cal)"),
    lty = c(1, 2, 3),
    col = c(COL_OBS, COL_RAW, COL_CAL),
    bty = "n"
  )
}

# Holdouts

res_2010 <- run_holdout(2010, q)
res_2012 <- run_holdout(2012, q)
res_2023 <- run_holdout(2023, q)
res_2024 <- run_holdout(2024, q)

# Comparison table

comparison <- rbind(
  res_2010$metrics,
  res_2012$metrics,
  res_2023$metrics,
  res_2024$metrics
)

cat("\nHoldout comparison table\n")
print(comparison)

# Plots

plot_holdout(res_2010)
plot_holdout(res_2012)
plot_holdout(res_2023)
plot_holdout(res_2024)

