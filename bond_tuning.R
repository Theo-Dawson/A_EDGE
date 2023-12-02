#Importing packages and data

library(tseries)
library(forecast)
library(lubridate)
library(dplyr)
library(lmtest)
library(splines)
library(rugarch)
library(FinTS)
library(tidyverse)
library(beepr)

################  FUNCTIONS  ####################
getData <- function(file_path, test_window, bool_plots) {
    bond <- read.csv(file_path)

    # Dummy vars for months
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

    # Dummy vars for days of month
    bond <- bond %>%
        mutate(DayOfMonth = day(Date))
    for (day in 1:31) {
        col_name <- sprintf("d%02d", day)
        bond <- bond %>%
            mutate(!!col_name := as.integer(DayOfMonth == day))
    }
    bond$DayOfMonth <- NULL

    y <- bond$Close/100
    bond$ly <- log(y)
    bond$times <- 1:nrow(bond)

    # Train on everything except test_window (2 week intervals)
    cutoff_index <- floor(nrow(bond) - test_window)  #floor to ensure int
    train <- bond[1:cutoff_index, ]
    test <- bond[(cutoff_index+1):nrow(bond), ]

    # Ex: nrow(bond) = 1257; test_window = 100
    # Train on 0 -> 1157
    # Test on 1158 -> nrow(bond) [length 100]

    # Ex: nrow(bond) = 1257; test_window = 2
    # Train on 0 -> 1255
    # Test on 1256 -> nrow(bond) [length 2]

    # NOTE: STILL ONLY PREDICT 2 WEEKS (14 DAYS) INTO TESTING WINDOW

    # Remove additional variables
    train <- train[, -c(1:7)]
    test <- test[, -c(1:7)]

    # Output variable, `ly` is log transformed
    y.train <- exp(train$ly)
    ly.train <- train$ly

    y.test <- exp(test$ly)
    ly.test <- test$ly

    if (bool_plots == TRUE) {
        ts.plot(ly.train)
    }

    # Times variable, just an integer sequence
    times <- 1:length(bond$times)
    times.train <- 1:length(ly.train)
    times.test <- (length(ly.train)+1):(length(ly.train) + length(ly.test))

    # Ex: nrow(bond) = 1257; nrow(train) = 1157 [given test_window = 1000]
    # times.train = 1:1157
    # times.test = 1158:1257

    # Convert to Matrix
    X.train <- data.frame(train)
    X.test <- data.frame(test)

    return(list(X.train, y.train, ly.train, times.train,
                X.test, y.test, ly.test, times.test,
                bond, times, test, train, y))
}

fitRegression <- function(degree,
                          X.train, times.train, ly.train,
                          bool_plots) {
    # Fit initial regression model to time-series data

    reg.model <- lm(ly ~ poly(times, degree=degree)*(.-times)+0, data=X.train)

    if (bool_plots == TRUE) {
        plot(times.train, ly.train, type="l", col="black", xlab="Date",
             ylab="Observed Values", main="Time Series with Regression Fit")
        lines(times.train, fitted(reg.model), type="l", col="red")
    }

    return(reg.model)
}

residualAnalysis <- function(reg.model, res, bool_plots, bool_output) {
    # Perform analysis on residuals from `reg.model` performance
    # Do not have to run for time-series analysis [just checking assumptions]

    if (bool_plots == TRUE) {
        ts.plot(res)  #time series of residuals
        acf(res, lag.max=14)  #ACH plot to check presence of autocorrelation
    }

    box_test <- Box.test(res, type="Ljung-Box")  #Ljung-Box test to check autocorrelation
    box_p_value <- box_test$p.value[[1]]
    # Note: low p-value means autocorrelation present

    white_test <- bptest(reg.model, ~ fitted.values(reg.model) + I(fitted.values(reg.model)^2))
    # Check heteroscedasticity (residuals have constant variance)
    # Note: low p-value means heteroscedasticity
    white_p_value <- white_test$p.value[[1]]

    if (bool_output == TRUE) {
        checkresiduals(res)  #summary plus normality check
        print(box_test)
        print(white_test)
    }

    return(list(box_p_value, white_p_value))
}

fitARMAGARCH <- function(res, reg.model,
                         times.train, y.train,
                         X.test, y.test,
                         garch1, garch2,
                         bool_plots, bool_output) {

    auto_fit <- auto.arima(res)  #auto.arima to select optimal params
    auto_order <- arimaorder(auto_fit)  #returns list(P, D, Q) of auto.arima()
    auto_P <- auto_order[[1]]
    auto_D <- auto_order[[2]]  #need more difference, should be 0
    auto_Q <- auto_order[[3]]

    # ARMA with GARCH
    spec <- ugarchspec(
        variance.model = list(garchOrder = c(garch1, garch2)),
        mean.model = list(armaOrder = c(auto_P, auto_Q), include.mean = FALSE),
        distribution.model = "sstd")
    garch.model <- ugarchfit(spec = spec, data = res, solver = "hybrid")

    # GARCH Residual Analysis
    res.garch <- residuals(garch.model, standardize=TRUE)
    res.garch <- as.vector(coredata(res.garch))

    if (bool_plots == TRUE) {
        ts.plot(res.garch^2)
        acf(res.garch^2, main="ACF of Squared GARCH Residuals", lag.max=14)
        pacf(res.garch^2, lag.max=14)
    }


    white_test <- ArchTest(res.garch^2, lag=14)
    white_p_value <- white_test$p.value[[1]]

    box_test <- Box.test(res.garch^2, lag = 14, type = "Ljung-Box")
    box_p_value <- box_test$p.value[[1]]


    if (bool_output == TRUE) {
        print(white_test)
        print(box_test)
    }

    return(list(garch.model,
                auto_P, auto_D, auto_Q,
                box_p_value, white_p_value))
}

getMAE <- function(y.actual, y.predicted, n) {

    mae = 0

    for (i in 1:n) {
        error = abs(y.actual[i] - y.predicted[i])
        mae = mae + ((1/n)*error)
    }
    return(mae)
}

getRMSE <- function(y.actual, y.predicted, n) {

    mse = 0
    rmse = 0

    for (i in 1:n) {
        error = y.actual[i] - y.predicted[i]
        squared_error = error^2
        mse = mse + ((1/n)*squared_error)
    }

    rmse = sqrt(mse)
    return(rmse)
}

getMDA <- function(y.actual, y.predicted, n) {
    total_right = 0

    for (i in 2:n) {
        delta_actual = y.actual[i] - y.actual[i-1]
        delta_predicted = y.predicted[i] - y.predicted[i-1]

        if ((delta_actual < 0) & (delta_predicted < 0)){
            total_right = total_right + 1
        } else if ((delta_actual >= 0) & (delta_predicted >= 0)) {
            total_right = total_right + 1
        } else {
            total_right = total_right
        }
    }

    return(total_right/(n-1))
}

train_performance <- function(reg.model, garch.model,
                              X.train, y.train, times.train,
                              bool_plots) {

    reg.predict <- predict(reg.model, newdata = X.train)
    res.predict <- fitted(garch.model)

    ly.predict <- reg.predict + res.predict
    y.predict <- coredata(exp(ly.predict))

    if (bool_plots == TRUE) {
        plot(times.train, y.train, type = 'l', col = 'black', xlab = 'Time', ylab = 'y', main = 'Training: Actual vs. Predicted')
        lines(times.train, y.predict, type = 'l', col = 'red')
    }

    y.predict <- as.vector(y.predict)
    n = length(y.train)

    #// TODO
    mae <- getMAE(y.train, y.predict, n)
    rmse <- getRMSE(y.train, y.predict, n)
    mda <- getMDA(y.train, y.predict, n)
    # coverage_prob <- getCoverage(y.train, y.predict, n)
    # interval_score <- getInterval(y.train, y.predict, n)

    return(list(y.predict, mae, rmse, mda))
}

test_performance <- function(reg.model, garch.model,
                             X.test, y.test, times.test,
                             test_window,
                             z_score,
                             bool_plots) {

    reg.predict <- predict(reg.model, newdata=X.test)
    garch.predict <- ugarchforecast(garch.model, n.ahead=test_window)
    res.predict <- fitted(garch.predict)
    ly.predict <- reg.predict + res.predict

    forecasted_sigma <- as.numeric(garch.predict@forecast$sigmaFor)
    upper_bound <- exp(ly.predict + z_score*forecasted_sigma)
    lower_bound <- exp(ly.predict - z_score*forecasted_sigma)

    y.predict <- coredata(exp(ly.predict))

    if (bool_plots == TRUE) {
        plot(y.test[1:test_window]*100, type = 'l', col = 'black', xlab = 'Time', ylab = 'y', main = 'Test Set: Actual vs. Predicted')
        lines(y.predict[1:test_window]*100, type = 'l', col = 'red')
        lines(upper_bound[1:14]*100, type = 'l', col = 'blue')
        lines(lower_bound[1:14]*100, type = 'l', col = 'blue')
    }

    y.predict <- as.vector(y.predict)
    n = min(length(y.test), 14)

    #// TODO
    mae <- getMAE(y.test, y.predict, n)
    rmse <- getRMSE(y.test, y.predict, n)
    mda <- getMDA(y.test, y.predict, n)
    # coverage_prob <- getCoverage(y.test, y.predict, n)
    # interval_score <- getInterval(y.test, y.predict, n)

    return(list(y.predict, mae, rmse, mda))
}
#################################################

degree_list <- c()
arima_p_list <- c()
arima_d_list <- c()
arima_q_list <- c()
garch_1_list <- c()
garch_2_list <- c()
z_stat_list <- c()
test_window_list <- c()
is_aic_list <- c()
is_mae_list <- c()
is_rmse_list <- c()
is_mda_list <- c()
is_coverage_list <- c()
is_interval_list <- c()
os_ais_list <- c()
os_mae_list <- c()
os_rmse_list <- c()
os_mda_list <- c()
os_coverage_list <- c()
os_interval_list <- c()


setwd("~/Documents/MIT/15.072_Advanced_Analytics_Edge/Project/A_EDGE")

bond_path <- "Data/BOND_5yr_daily.csv"

## Critical Functions: (must run)
# 1. getData(bond_path, test_window, bool_plots)
# 2. fitRegression(degree, X.train, times.train, ly.train, bool_plots)
# 3. residualAnalysis(reg.model, res, bool_plots, bool_output)
# 4. gitARMAGARCH(res, reg.model, times.train, y.train, X.test, y.test, bool_plots, bool_output)
# 5. train_performance(reg.model, garch.model, X.train, y.train, times.train, bool_plots))
# 6. test_performance(reg.model, garch.model, X.test, y.test, times.test, test_window, z_score, bool_plots)


for (p in c(1, 2, 3, 4)) {
    for (g1 in c(1, 2, 3)) {
        for (g2 in c(1, 2, 3)) {
            for (t_window in seq(100, 14, by=-1)) {
                test_window = t_window
                degree_p = p
                plots = FALSE
                output = FALSE
                garch1 = g1
                garch2 = g2
                z_stat = 1.96

                data_return <- getData(bond_path, test_window, plots)
                # Args: file_path, test_window, bool_plots
                # Return: list(X.train, y.train, ly.train, times.train,
                #              X.test, y.test, ly.test, times.test,
                #              bond, times))

                X.train <- data_return[[1]]
                y.train <- data_return[[2]]
                ly.train <- data_return[[3]]
                times.train <- data_return[[4]]

                X.test <- data_return[[5]]
                y.test <- data_return[[6]]
                ly.test <- data_return[[7]]
                times.test <- data_return[[8]]

                bond <- data_return[[9]]
                times <- data_return[[10]]
                test <- data_return[[11]]
                train <- data_return[[12]]
                y <- data_return[[13]]

                remove(data_return)



                reg.model <- fitRegression(degree_p, X.train, times.train, ly.train, plots)
                # Args: degree, X.train, times.train, ly.train, bool_plots
                # Return: reg.model

                res <- reg.model$residuals



                residual_return <- residualAnalysis(reg.model, res, plots, output)
                # Args: reg.model, res, bool_plots, bool_output
                # Return: list(box_p_value, white_p_value)

                box_p_value_initial <- residual_return[[1]]
                white_p_value_initial <- residual_return[[2]]

                remove(residual_return)



                garch_return <- fitARMAGARCH(res, reg.model,
                                             times.train, y.train,
                                             X.test, y.test,
                                             garch1, garch2,
                                             plots, output)
                # Args: res, reg.model, times.train, y.train, X.test, y.test, garch1, garch2, bool_plots, bool_output
                # Return: list(garch.model, auto_P, auto_D, auto_Q, box_p_value, white_p_value)

                garch.model <- garch_return[[1]]

                arima_P <- garch_return[[2]]
                arima_D <- garch_return[[3]]
                arima_Q <- garch_return[[4]]

                box_p_value_final <- garch_return[[5]]
                white_p_value_final <- garch_return[[6]]

                remove(garch_return)



                train_return <- train_performance(reg.model, garch.model,
                                                  X.train, y.train, times.train,
                                                  plots)

                # Args: reg.model, garch.model, X.train, y.train, times.train, bool_plots
                # Return: y.predict, mae, rmse, mda
                y.train.predict <- train_return[[1]]
                is_mae <- train_return[[2]]
                is_rmse <- train_return[[3]]
                is_mda <- train_return[[4]]

                remove(train_return)


                test_return <- test_performance(reg.model, garch.model,
                                                X.test, y.test, times.test,
                                                test_window,
                                                z_stat,
                                                plots)
                # Agrs: reg.model, garch.model, X.test, y.test, times.test, test_window, z_score, bool_plots
                # Return: y.predict mae, rmse, mda
                y.test.predict <- test_return[[1]]
                os_mae <- test_return[[2]]
                os_rmse <- test_return[[3]]
                os_mda <- test_return[[4]]

                remove(test_return)


                # //TODO
                is_aic = 0
                is_coverage = 0
                is_interval = 0
                os_aic = 0
                os_coverage = 0
                os_interval = 0

                degree_list <- append(degree_list, degree_p)
                arima_p_list <- append(arima_p_list, arima_P)
                arima_d_list <- append(arima_d_list, arima_D)
                arima_q_list <- append(arima_q_list, arima_Q)
                garch_1_list <- append(garch_1_list, garch1)
                garch_2_list <- append(garch_2_list, garch2)
                z_stat_list <- append(z_stat_list, z_stat)
                test_window_list <- append(test_window_list, test_window)
                is_aic_list <- append(is_aic_list, is_aic)
                is_mae_list <- append(is_mae_list, is_mae)
                is_rmse_list <- append(is_rmse_list, is_rmse)
                is_mda_list <- append(is_mda_list, is_mda)
                is_coverage_list <- append(is_coverage_list, is_coverage)
                is_interval_list <- append(is_interval_list, is_interval)
                os_ais_list <- append(os_ais_list, os_aic)
                os_mae_list <- append(os_mae_list, os_mae)
                os_rmse_list <- append(os_rmse_list, os_rmse)
                os_mda_list <- append(os_mda_list, os_mda)
                os_coverage_list <- append(os_coverage_list, os_coverage)
                os_interval_list <- append(os_interval_list, os_interval)
            }
        }
    }
}

# Write to CSV:
# 1. Type of regression (degree_p)
# 2. ARMIA Values (arima_P, arima_D, arima_Q)
# 3. GARCH Values (garch1, garch2)
# 4. Testing Window (test_window)
# 4. In-Sample Evaluation (AIC, MAE, RMSE, MDA, Coverage Prob., Interval Score)
# 5. Out-of-Sample Evaluation (AIC, MAE, RMSE, MDA, Coverage Prob., Interval Score)


params_df = data.frame(`Degree Polynomial` = degree_list,
                       `ARIMA P` = arima_p_list,
                       `ARIMA D` = arima_d_list,
                       `ARIMA Q` = arima_q_list,
                       `GARCH 1` = garch_1_list,
                       `GARCH 2` = garch_2_list,
                       `CI Z Stat` = z_stat_list,
                       `Test Window` = test_window_list,
                       `IS AIC` = is_aic_list,
                       `IS MAE` = is_mae_list,
                       `IS RMSE` = is_rmse_list,
                       `IS MDA` = is_mda_list,
                       `IS Coverage Prob.` = is_coverage_list,
                       `IS Interval Score` = is_interval_list,
                       `OS AIC` = os_ais_list,
                       `OS MAE` = os_mae_list,
                       `OS RMSE` = os_rmse_list,
                       `OS MDA` = os_mda_list,
                       `OS Coverage Prob.` = os_coverage_list,
                       `OS Interval Score` = os_interval_list)
write.csv(params_df, "params_test.csv", row.names=TRUE)

beep()
beep()
beep()
