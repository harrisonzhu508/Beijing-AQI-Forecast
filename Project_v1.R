library(readr)
library(car)
library(tseries)
library(astsa)


############################################################STEP 1
#Read dataframe
Beijing <- suppressMessages(read_csv("~/Desktop/Semestre 2/Time Series/Project/beijing23 - Sheet1.csv"))

View(Beijing)
par(mfrow = c(1, 3),pty = "s")
#time series conversion
Beijing_TS <- ts(Beijing$AQI)
plot(Beijing_TS,main = "Original Dataset: Plot", xlab = "Time", ylab = "AQI")
acf(Beijing_TS, main = "Original Dataset: ACF")
pacf(Beijing_TS, main = "Original Dataset: PACF")

#kpss test
kpss.test(Beijing_TS)

############################################################
############logdiff transformation
log_Beijing_TS <- log(Beijing_TS)
diff_log_beijing <- diff(log_Beijing_TS)

#plots
par(mfrow = c(1, 3),lheight = 1,pty = "s")

#diff
diff_log_beijing <- diff(log_Beijing_TS)
plot(diff_log_beijing,main = "Logdiff Dataset: Plot")

#acf
acf(diff_log_beijing, main = "Logdiff Dataset: ACF")

#pacf
pacf(diff_log_beijing, main = "Logdiff Dataset: PACF")
#kpss
kpss.test(diff_log_beijing)

############################################################
#####Box-Cox Transformations

Beijing_root <- Beijing_TS^(0.25)

#polynomial of order 10 try

par(mfrow = c(1, 3),lheight = 1,pty = "s")


plot(Beijing_root,main = "1/4 Dataset: Plot")

#acf
acf(Beijing_root, main = "1/4 Dataset: ACF")

#pacf
pacf(Beijing_root, main = "1/4 Dataset: PACF")

#kpss
kpss.test(Beijing_root)

############################################################
#######Trying out STL Decomposition

Beijing_TS1 <- ts(Beijing$AQI^(0.25),frequency = 90)
plot(stl(Beijing_TS1,s.window = "per",t.window = 101))


############################################################STEP2

###MA 4
model_MA4 <- sarima(diff_log_beijing, p = 0, d = 0, q = 4) 
fit4 <- model_MA4$fit


###MA 3
model_MA1 <- sarima(diff_log_beijing, p = 0, d = 0, q = 3) 
fit1 <- model_MA1$fit


###MA 2
model_MA2 <- sarima(diff_log_beijing, p = 0, d = 0, q = 2) 
fit2 <- model_MA2$fit


###ARMA 12
model_ARMA12 <- sarima(diff_log_beijing, p = 1, d = 0, q = 2) 
fit_ARMA12 <- model_ARMA12$fit

###ARMA 13
model_ARMA13 <- sarima(diff_log_beijing, p = 1, d = 0, q = 3) 
fit_ARMA13 <- model_ARMA13$fit

###ARMA 21
model_ARMA21 <- sarima(diff_log_beijing, p = 2, d = 0, q = 1) 
fit_ARMA21 <- model_ARMA21$fit

###ARMA 22
model_ARMA22 <- sarima(diff_log_beijing, p = 2, d = 0, q = 2) 
fit_ARMA22 <- model_ARMA22$fit

###ARMA 22
model_ARMA23 <- sarima(diff_log_beijing, p = 2, d = 0, q = 3) 
fit_ARMA23 <- model_ARMA23$fit

###ARMA 31
model_ARMA31 <- sarima(diff_log_beijing, p = 3, d = 0, q = 1) 
fit_ARMA31 <- model_ARMA31$fit

###ARMA 32
model_ARMA32 <- sarima(diff_log_beijing, p = 3, d = 0, q = 2) 
fit_ARMA32 <- model_ARMA32$fit

###ARMA 33
model_ARMA33 <- sarima(diff_log_beijing, p = 3, d = 0, q = 3) 
fit_ARMA33 <- model_ARMA33$fit

#summary
c("04", fit4$aic,fit4$loglik)
c("03", fit1$aic,fit1$loglik)
c("02", fit2$aic,fit2$loglik)
c("11", fit_ARMA11$aic,fit_ARMA11$loglik)
c("12", fit_ARMA12$aic,fit_ARMA12$loglik)
c("13", fit_ARMA13$aic,fit_ARMA13$loglik)
c("21", fit_ARMA21$aic,fit_ARMA21$loglik)
c("22", fit_ARMA22$aic,fit_ARMA22$loglik)
c("23", fit_ARMA23$aic,fit_ARMA23$loglik)
c("31", fit_ARMA31$aic,fit_ARMA31$loglik)
c("32", fit_ARMA32$aic,fit_ARMA32$loglik)
c("33", fit_ARMA33$aic,fit_ARMA33$loglik)

