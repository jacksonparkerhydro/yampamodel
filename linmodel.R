library(dplyr)


dat <- read.csv("C:\\Users\\jacks\\Documents\\YAMPAFLOW\\analysisCSV.csv")

hist(dat$Q)
hist(dat$SWEAvg)
hist(dat$dSWE)
hist(dat$precip)
hist(dat$T)

#need to transform Q variables (many zeroes, right skew), SWE variables not dSWE (same), precip

datlogs <- dat %>%
  mutate(
    Q_log        = log1p(Q), 
    Q_lag1_log = log1p(QMinus1),
    Q_lag7_log = log1p(QMinus7),
    SWE_log      = log1p(SWEAvg),
    SWE_lag1_log = log1p(SWEMinus1),
    SWE_lag7_log = log1p(SWEMinus7),
    dSWE7_asinh    = asinh(dSWE),
    P7_log       = log1p(precip)
  )

#check again

hist(datlogs$Q_log)
hist(datlogs$SWE_log)
hist(datlogs$SWE_lag1_log)
hist(datlogs$SWE_lag7_log)
hist(datlogs$dSWE7_asinh)
hist(datlogs$P7_log)

#first model run with no transform

lmnotrans <- lm(Q ~ SWEAvg + SWEMinus1 + SWEMinus7 + dSWE + QMinus1 + QMinus7 + T + precip, data = dat)
summary(lmnotrans)
plot(lmnotrans)

#second model with transforms

lmtrans <- lm(Q_log ~ SWEAvg + SWEMinus1 + SWEMinus7 + dSWE + Q_lag1_log + Q_lag7_log + T + P7_log, data = datlogs)
summary(lmtrans)
plot(lmtrans)

#third model throwing out streamflow lag

lmpredict <- lm(Q_log ~ SWEAvg + SWEMinus1 + SWEMinus7 + dSWE + T + P7_log, data = datlogs)
summary(lmpredict)
plot(lmpredict)

#final model removing 1 day SWE lag
lmfinal <- lm(Q_log ~ SWEAvg + SWEMinus7 + T + P7_log, data = datlogs)
summary(lmfinal)
plot(lmfinal)

#cut out WY2024 so i can retrain on prior data and validate
class(datlogs$DateSpine)
install.packages("stringr")
library(stringr)
#converting date string to usable water year column
datlogs <- datlogs %>%
  mutate(
    m = as.integer(str_extract(DateSpine, "^[0-9]{1,2}")),
    y = as.integer(str_extract(DateSpine, "[0-9]{4}$")),
    WY = ifelse(m >= 10, y + 1, y)
  ) %>%
  select(-m, -y)

#new dataframes for training/testing model
train <- datlogs %>% filter(WY < 2024)
test  <- datlogs %>% filter(WY == 2024)

#train test model (dropped dSWE due to irrelevance at this point)
lmtrained <- lm(Q_log ~ SWEAvg + SWEMinus7 + P7_log + T, data = train)
summary(lmtrained)
plot(lmtrained)

#test test model!
test$Q_log_pred <- predict(lmtrained, newdata = test)
test$Q_pred <- exp(test$Q_log_pred) - 1
# predict on WY2024
test$Q_log_pred <- as.numeric(predict(lmtrained, newdata = test))
test$Q_pred <- exp(test$Q_log_pred) - 1

# dates + validity mask
x <- as.Date(test$DateSpine, format = "%m/%d/%Y")

ok <- complete.cases(x, test$Q_log, test$Q_log_pred, test$Q, test$Q_pred) &
  test$Q > 0 & test$Q_pred > 0

if (sum(ok) < 2) {
  stop(paste0(
    "No valid points to plot. sum(ok)=", sum(ok),
    "\nCheck DateSpine parsing and missing predictors in WY2024."
  ))
}

# RMSE in log scale (log1p space)
rmse_log <- sqrt(mean((test$Q_log_pred[ok] - test$Q_log[ok])^2))

# model R^2 (from training fit)
model_r2_adj <- summary(lmtrained)$adj.r.squared

# prediction R^2 (out-of-sample, WY2024, log1p)
pred_r2_log <- 1 - sum((test$Q_log_pred[ok] - test$Q_log[ok])^2) /
  sum((test$Q_log[ok] - mean(test$Q_log[ok]))^2)

metrics_txt <- paste0(
  "Adj R² (model, train) = ", sprintf("%.3f", model_r2_adj),
  "   R² (pred, WY2024, log1p) = ", sprintf("%.3f", pred_r2_log),
  "\nRMSE (WY2024, log1p) = ", sprintf("%.3f", rmse_log)
)

ymin <- min(c(test$Q[ok], test$Q_pred[ok]))
ymax <- max(c(test$Q[ok], test$Q_pred[ok]))
ylim <- c(ymin, ymax * 1.3)

yticks <- c(50, 100, 250, 500, 1000, 2500, 5000, 10000, 20000)
yticks <- yticks[yticks >= ylim[1] & yticks <= ylim[2]]

op <- par(no.readonly = TRUE)
par(mar = c(5, 4, 4, 2) + 0.1)

plot(
  x[ok], test$Q[ok],
  type = "l",
  col = "black", lwd = 2,
  log = "y",
  ylim = ylim,
  xlab = "Date",
  ylab = "",
  main = "WY 2024: Observed vs Predicted Discharge on the Yampa at Deerlodge Park",
  yaxt = "n"
)

axis(2, at = yticks, labels = format(yticks, big.mark = ","))
mtext("Discharge (CFS, log scale)", side = 2, line = 3)

lines(x[ok], test$Q_pred[ok], col = "red", lwd = 2)

legend(
  "topright",
  c("Observed", "Predicted"),
  col = c("black", "red"),
  lwd = 2,
  bty = "n"
)
metrics_lines <- c(
  sprintf("Adj. R² (model) = %.3f", model_r2_adj),
  sprintf("Adj. R² (prediction) = %.3f", pred_r2_log),
  sprintf("RMSE = %.3f", rmse_log)
)

legend(
  x = "topleft",
  inset = c(0.01, 0.01),
  legend = metrics_lines,
  bty = "n",
  cex = 0.9
)



