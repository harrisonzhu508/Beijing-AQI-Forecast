library(readr)
library(car)
library(tseries)
library(astsa)
library(smooth)
library(Mcomp)
library(fGarch)

#square values of residuals
#monte carlo simulate t noise
#

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

#periodogram
cpgram(Beijing_TS, main = "Cumulative periodogram")


#kpss test
kpss.test(Beijing_TS)

############################################################
############diff log transformation
log_Beijing_TS <- log(Beijing_TS)
diff_log_beijing <- diff(log_Beijing_TS)

#plots
par(mfrow = c(1, 3),lheight = 1,pty = "s")

#diff
diff_log_beijing <- diff(log_Beijing_TS)
plot(diff_log_beijing,main = "Logdiff Dataset: Plot")

#acf
acf(diff_log_beijing, main = "Logdiff Dataset: ACF", lag.max = 1000)
points(c(365, 730), c(0.1,0.1), pch = 22, col = 'red')


#pacf
pacf(diff_log_beijing, main = "Logdiff Dataset: PACF", lag.max = 1000)
#kpss
kpss.test(diff_log_beijing)

############################################################
#######diff diff^365 smooth log




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

############################################################kernel estimation

plot(Beijing_TS, type = "l")
lines(ksmooth(time(Beijing_TS), Beijing_TS, "normal", bandwidth=3), lwd=2, col=2)
par(fig = c(.65, 1, .65, 1), new = TRUE) # the insert


smoothed_diff <- Beijing_TS - ksmooth(time(Beijing_TS), Beijing_TS, "normal", bandwidth=3)$y

model_test <- sarima(smoothed_diff, p = 0, d = 0, q = 3 ) 

model_test$fit$aic

############################################################STEP2


#########ARMA
###MA 4
model_MA4 <- sarima(diff_log_beijing, p = 0, d = 0, q = 4) 
fit4 <- model_MA4$fit


###MA 3
model_MA1 <- sarima(diff_log_beijing, p = 0, d = 0, q = 3) 
fit1 <- model_MA1$fit

plot(model_MA1$fit$residuals)

acf(model_MA1$fit$residuals^2, lag.max = 500)

###MA 2
model_MA2 <- sarima(diff_log_beijing, p = 0, d = 0, q = 2, no.constant = TRUE) 
fit2 <- model_MA2$fit


###ARMA 12
model_ARMA12 <- sarima(diff_log_beijing, p = 1, d = 0, q = 2, no.constant = TRUE) 
fit_ARMA12 <- model_ARMA12$fit

###ARMA 13
model_ARMA13 <- sarima(diff_log_beijing, p = 1, d = 0, q = 3, no.constant = TRUE) 
fit_ARMA13 <- model_ARMA13$fit

###ARMA 21
model_ARMA21 <- sarima(diff_log_beijing, p = 2, d = 0, q = 1 ,no.constant = TRUE) 
fit_ARMA21 <- model_ARMA21$fit

###ARMA 22
model_ARMA22 <- sarima(diff_log_beijing, p = 2, d = 0, q = 2, no.constant = TRUE) 
fit_ARMA22 <- model_ARMA22$fit

qqPlot(model_ARMA22$fit$residuals)

###ARMA 22
model_ARMA23 <- sarima(diff_log_beijing, p = 2, d = 0, q = 3, no.constant = TRUE) 
fit_ARMA23 <- model_ARMA23$fit

###ARMA 31
model_ARMA31 <- sarima(diff_log_beijing, p = 3, d = 0, q = 1, no.constant = TRUE) 
fit_ARMA31 <- model_ARMA31$fit



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


##########S#########S#########S#########S#########S#########S#########SLBT



#########SARMA
diff2 <- diff(diff_log_beijing, lag = 365)

par(mfrow = c(1,2))
acf(diff2)
pacf(diff2)

model_sarima <- sarima(diff2, p = 3, d = 0, q = 3, P = 0, D = 0, Q = 0, S = 0)
qqPlot(model_sarima$fit$residuals)

##########S#########S#########S#########S#########S#########S#########SLBT

lbtsarima <- c(); 
for (h in 6:20) lbtsarima[h] <- Box.test(model_sarima$residuals,lag=h,type='Ljung-Box',fitdf=30 - 2 - 2)$p.value

plot(lbtsarima, ylim=c(0,1)); abline(h=0.05,col='blue',lty='dotted')

##########S#########S#########S#########S#########S#########S#########Forecasting

plot(forecast(model_sarima, h = 100, level = .95))

############################################################
##########PLAYGROUND########################################
############################################################


#get ARMA(2,2)
s <- auto.arima(diff_log_beijing , max.q = 10, max.Q = 10, 
                max.d = 1, allowdrift = TRUE, max.D = 1, max.P = 10, 
                max.p = 10, max.order = 25, stepwise = FALSE, ic = "aic")
roots1 <- polyroot(c(1, -coef(s)[grep("^ar",names(s$coef))]))
for_std <- forecast(s)
autoplot(for_std, include = 5)

#
s <- auto.arima(diff_log_beijing , max.q = 10, max.Q = 10, 
                max.d = 1, allowdrift = TRUE, max.D = 1, max.P = 10, 
                max.p = 10, max.order = 25, stepwise = FALSE, ic = "aic")
#GARCH
g <- garchFit(~ arma(2,0) + garch(1, 0), data = diff_log_beijing, cond.dist='norm')
stdg <- g@residuals/g@sigma.t
n <- length(stdg)

par(mfrow = c(1,2))
qqplot(x=qt(p=(1:n)/(n+1),df=10)*sqrt((10-2)/10), y=stdg )
qqPlot(stdg, distribution = "t", df = 10)





