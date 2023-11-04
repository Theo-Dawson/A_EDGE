#Importing packages and data

install.packages("FinTS")
library(tseries)
library(forecast)
library(lubridate)
library(dplyr)
library(lmtest)
library(splines)
library(rugarch)
library(FinTS)
library(tidyverse)

bond <- read.csv("C:/Users/tanne/Downloads/BOND_5yr_daily.csv")

#Cleaning data

##Dummy variables for months
bond$Date <- as.Date(bond$Date, format="%Y-%m-%d")
bond <- bond %>%
  mutate(m01 = as.integer(month(Date) == 1),
         m02 = as.integer(month(Date) == 2),
         m03 = as.integer(month(Date) == 3),
         m04 = as.integer(month(Date) == 4),
         m05 = as.integer(month(Date) == 5),
         m06 = as.integer(month(Date) == 6),
         m07 = as.integer(month(Date) == 7),
         m08 = as.integer(month(Date) == 8),
         m09 = as.integer(month(Date) == 9),
         m10 = as.integer(month(Date) == 10),
         m11 = as.integer(month(Date) == 11),
         m12 = as.integer(month(Date) == 12))

bond$y = bond$Close/100
bond$ly = log(bond$y)
bond$times = 1:nrow(bond)

##Training (80%) and test (20%) sets
cutoff_index = floor(0.95* nrow(bond))
train = bond[1:cutoff_index, ]
test = bond[(cutoff_index + 1):nrow(bond), ]

##Remove unused variables
train = train[, -c(1:7,20)]
test = test[, -c(1:7,20)]

##Output variable, ly is log transformed
y.train = exp(train$ly)
ly.train = train$ly

y.test = exp(test$ly)
ly.test = test$ly

ts.plot(y.train)

##Times variable, just an integer sequence
times = 1:length(bond$times)
times.train = 1:length(ly.train)
times.test = (length(ly.train)+1):(length(ly.train)+length(ly.test))

##Matrix of everything
X.train <- train
X.test <- test

#Initial fitting (cubic splines LR with no autocorrelation structures)

knots <- quantile(times.train, probs = c(0.25, 0.5, 0.75))
reg.model<- lm(ly ~poly(times, degree=3)*(.-times), data = X.train)
#reg.model<- lm(ly ~ poly(times, degree = 3) + ., data = X.train)
plot(times.train, ly.train, type = "l", col = "black", xlab = "Date", 
     ylab = "Observed values", main = "Time Series with Cubic Splines Fit")
lines(times.train, fitted(reg.model), type = "l", col = "red")

#Residual analysis

res = reg.model$residuals

##Time series of residuals
ts.plot(res)

##ACH plot to check if there is autocorrelation
acf(res, lag.max = 100)

##Ljung-Box test to check autocorrelation (low p value means there is autocorrelation)
Box.test(res, type = "Ljung-Box")

##White test to check heterosexuality (low p value means heterosexual)
white_test = bptest(reg.model, ~ fitted.values(reg.model) + I(fitted.values(reg.model)^2))
print(white_test)

##Summary plus normalcy check
checkresiduals(res)




#Fitting ARMA-GARCH

##auto.arima to automatically select optimal ARMA parameters
res = diff(res, lag =21)
auto.arima(res)


##ARMA-GARCH
spec <- ugarchspec(
  variance.model = list(garchOrder = c(5, 5)),
  mean.model = list(armaOrder = c(25,25), include.mean = FALSE),
  distribution.model = "sstd")
garch.model <- ugarchfit(spec = spec, data = res, solver = "hybrid")

#GARCH Residual Analysis
res.garch <- residuals(garch.model, standardize = TRUE)
res.garch = as.vector(coredata(res.garch))
ts.plot(res.garch)
acf(res.garch^2, main = "ACF of Squared GARCH Residuals", lag.max = 100)
pacf(res.garch^2, lag.max = 100)
ArchTest(res.garch^2, lag = 22)

# Perform Ljung-Box test on the squared residuals
Box.test(res.garch^2, lag = 22, type = "Ljung-Box")

checkresiduals(res.garch^2)

#Training set performance
reg.predict = predict(reg.model, newdata = X.train) 
res.predict = fitted(garch.model)
ly.predict <- reg.predict + res.predict
y.predict = exp(ly.predict)

interval = 1
plot(times.train, y.train, type = 'l', col = 'black', xlab = 'Time', ylab = 'y', main = 'Actual vs. Predicted')
lines(times.train, y.predict, type = 'l', col = 'red')

#Test set performance
reg.predict = predict(reg.model, newdata = X.test) ##WE NEED TO CHECK FOR SE SO WE CAN GET CONF INT FROM JUST USING ARIMA
res.predict = ugarchforecast(garch.model, n.ahead = length(times.test))
res.predict <- fitted(res.predict)
undifferenced_fitted_values = c()
for (i in 1:length(res.predict)) {
  undifferenced_fitted_values[i] <- res.predict[i] + res[i + (length(res) - length(res.predict)) - 21]
}
ly.predict = reg.predict + res.predict
forecasted_sigma <- as.numeric(res.predict@forecast$sigmaFor)
upper_bound = exp(ly.predict + forecasted_sigma^0.5)
lower_bound = exp(ly.predict - forecasted_sigma^0.5)
y.predict = coredata(exp(ly.predict))

plot(y.test[1:14]*100, type = 'l', col = 'black', xlab = 'Time', ylab = 'y', main = 'Test Set: Actual vs. Predicted', ylim=c(80, 100))
lines(y.predict[1:14]*100, type = 'l', col = 'red')
lines(upper_bound[1:14]*100, type = 'l', col = 'blue')
lines(lower_bound[1:14]*100, type = 'l', col = 'blue')

reg.predict$se

z = diff(ly.train,1)
ts.plot(z)
acf(z)
pacf(z)

arima.model = auto.arima(ly.train, seasonal = TRUE, D = 12)
res = residuals(arima.model)
ts.plot(res)
spec = ugarchspec(
  variance.model = list(garchOrder = c(2, 2)),
  mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),
  distribution.model = "sstd")
garch.model <- ugarchfit(spec = spec, data = res, solver = "hybrid")
forecast_object <- forecast(arima.model, h = 12)  # 'h' is the horizon of the forecast
arima_forecasts <- as.numeric(forecast_object$mean)

garch_forecast <- ugarchforecast(garch.model, n.ahead = 12)
garch_volatility <- as.numeric((sigma(garch_forecast))^0.5)
upper_bound <- arima_forecasts + garch_volatility  # 95% prediction interval upper bound
lower_bound <- arima_forecasts - garch_volatility
plot(bond$y, type = "l", col = "blue", main = "ARIMA-GARCH Forecasts vs Actual", xlab = "Time", ylab = "Values", xlim = c(0,80), ylim = c(2.6, 4))
lines(67:78, exp(arima_forecasts), col = "green", type = "l")

lines(exp(upper_bound)[1:7], col = "green", type = "l")
lines(exp(lower_bound)[1:7], col = "green", lty = 2)
abline(v = 2, col = "purple", lwd = 2)  # lwd is the line width
