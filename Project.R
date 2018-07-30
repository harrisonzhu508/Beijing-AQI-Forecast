#libraries used
library(astsa)
library(car)
library(fGarch)
library(forecast)
library(ggplot2)
library(KernSmooth)
library(lubridate)
library(readr)
library(rugarch)
library(tseries)
library(TSA)
library(xtable)
library(xts)
################################################################################Data
################################################################################################################################################################
################################################################################################################################################################

Beijing <-read_csv("~/Desktop/Semestre 2/7. Time Series/Project/beijing23 - Sheet1.csv")

########################################Test data for 01/02/2018 to end of March
BEIJING02_03 <- read_csv("~/Desktop/BEIJING02-03.csv")

######################################convert data to time series object
Beijing.ts <- ts(Beijing$AQI)
#####################################transformations
Beijing.log <- log(Beijing.ts)

date <- seq(as.Date("2015/01/01"), by="1 day", length.out=1125)
plot(date, Beijing$AQI, xlab = "Date", ylab = (ylab <- "AQI Index"), main = "Beijing AQI January 2015 to February 2018", type = "l")


#####################################plot the Y_{t} data
par(mfrow = c(1, 1),pty = "s")
plot(Beijing.log,main = "Original Dataset: Plot", xlab = "Time", ylab = "AQI")
acf(Beijing.ts, main = "Original Dataset: ACF")
pacf(Beijing.ts, main = "Original Dataset: PACF")

#kpss test
kpss.test(Beijing.ts)

################################################################################Data
################################################################################################################################################################
################################################################################################################################################################
#functions to estimate the variance
##function to approximate the variance using a window of length 5
##assuming that the variance in the variance is over 365 days
window_func <- function(m, data) {
  Y <- rep(0, 365 / m)
  Z  <- rep(0 , 365)
  y1 <- data[1:365]^2
  y2 <- data[366:730]^2
  y3 <- data[731:1095]^2
  
  for (i in 1:(365 / m)){
    
    pivot <- (i - 1) * m + 1
    Y[i] <- sum(y1[pivot:(pivot + m - 1)]) + sum(y2[pivot:(pivot + m - 1)]) + sum(y3[pivot:(pivot + m - 1)])
    Z[pivot:(pivot + m - 1)] <- rep(Y[i],m)
  }
  
  Z / (3 * m - 1)
  
}

mean(Beijing.log)
#Box-Jung Dr. Thibaud's code modified
box.jung <- function(model, p, q){
  
  e <- model$residuals
  lbt <- c()
  m <- p + q + 1
  for (h in m:50) {
    lbt[h] <- Box.test(e, lag=h, type='Ljung-Box', fitdf=(p + q))$p.value
  }
  plot(m:50, lbt[m:50], ylim=c(0,1), xlab='DF', ylab='p value', main = "Ljung-Box")
  abline(h=0.05,col = 'blue',lty = 'dotted')
  
}

#aic_table from Leo Belzile
aic_table <- function(data, P, Q) {
  table <- matrix(NA, (P + 1), (Q + 1))
  for (p in 0:P) {
    for (q in 0:Q) {
      table[p + 1, q + 1] <- Arima(data, order = c(p, 1, q))$aic
    }
  }
  dimnames(table) <- list(paste("<b> AR", 0:P, "</b>", sep = ""), paste("MA", 
                                                                        0:Q, sep = ""))
  table
}
################################################################################################################################################################
################################################################################################################################################################
#m = 5
st_dev.2 <- sqrt(window_func(5, Beijing.log))
stabseries.log <- Beijing.log / c(rep(st_dev.2, 3), st_dev.2[1:30])

#truncation of series

par(mfrow = c(1, 3))
plot(date[2:366], diff(stabseries.log)[1:365], type = 'l')
plot(date[367:731], diff(stabseries.log)[366:730], type = 'l')
plot(date[732:1096], diff(stabseries.log)[731:1095], type = 'l')

par(mfrow = c(3, 1))
plot(date[1:366], Beijing.log[1:366], type = 'l', ylab = expression(LY[t]), xlab = "Date in 2015")
plot(date[367:731], Beijing.log[367:731], type = 'l', ylab = expression(LY[t]), xlab = "Date in 2016")
plot(date[732:1096], Beijing.log[732:1096], type = 'l', ylab = expression(LY[t]), xlab = "Date in 2017")

#assessment plots
plot(date[-1], diff(stabseries.log), main = expression(paste("Plot of D[log(", Y[t]/s[t]^(1),")]")), ylab = "AQI", xlab = "Date" ,type ="l")
acf(diff(stabseries.log), main = "ACF")
pacf(diff(stabseries.log), main = "PACF")

model.stabseries.log <- Arima(stabseries.log, c(2, 1, 2), method = "ML")

layout(matrix(c(2,3,1,1), 2, 2, byrow = TRUE))
qqPlot(model.stabseries.log$residuals, main = "Q-Q Plot", ylab = "Empirical Quantiles")
box.jung(model.stabseries.log, p = 2, q = 2)
cpgram(model.stabseries.log$residuals, main = "Cumulative Periodogram")

par(mfrow = c(1,2), pty = 's')
acf(model.stabseries.log$residuals, main = "ACF")
pacf(model.stabseries.log$residuals, main = "PACF")

Beijing_aic_table.2 <- aic_table(stabseries.log, 4, 5)
dimnames(Beijing_aic_table.2) <- list(paste0('AR', 0:4), paste0('MA', 0:5))
latex_tab <- xtable::xtable(Beijing_aic_table.2, booktabs = TRUE, caption = 'AIC values for ARMA models for the Lake Huron dataset')
print(latex_tab,booktabs = TRUE, caption.placement = 'top')

#log likelihoods
logliks2 <- matrix(0L, 6, 6)
for (p in 0:5){
  for (q in 0:5){
    logliks2[p+1,q+1] <- Arima(stabseries.log, c(p, 1, q), method = "ML")$loglik
    
  }
  
}

rbind(logliks2[1,] , Beijing_aic_table.2[1,])

xtable::xtable(rbind(logliks2[1,] , Beijing_aic_table.2[1,]), booktabs = TRUE, caption = 'AIC values for ARMA models')


#extract coefficients
model.stabseries.log$coef
model.stabseries.log$sigma2
#ar
polyroot(c(1, -coef(model.stabseries.log)[grep("^ar", names(model.stabseries.log$coef))]))
#ma
polyroot(c(1, coef(model.stabseries.log)[grep("^ma", names(model.stabseries.log$coef))]))

#prediction
h <- 300
forecast.log <- forecast(model.stabseries.log, h = h) 

f1 <- c(sapply(forecast.log$lower[,2], as.numeric))
f2 <- c(sapply(forecast.log$mean, as.numeric))
f3 <- c(sapply(forecast.log$upper[,2], as.numeric))


forecast <- exp(f2 * st_dev.2[30:329])
forecast.upper <- exp(f3 * st_dev.2[30:329])

Prediction <- data.frame(c(c(Beijing.ts[800:1125]), forecast))
Prediction$newcolumn <- c(rep(NA, 326),forecast.upper)
colnames(Prediction) <- c("AQI", "U 95%")

Date <- seq(as.Date("2017/03/11"), by="1 day", length.out= 326 + 300)
AQI <- Prediction$AQI
U <- Prediction$`U 95%`

#forecast plot
ggplot(Prediction, aes(Date,AQI)) +
  ggtitle("Beijing AQI Forecast") + 
  xlab("Date") + 
  ylab("AQI") +  
  geom_path() + 
  geom_ribbon(data = Prediction, aes(x = Date, y = AQI, ymin = 0, ymax = U), 
              color= "blue", alpha = .3) +
  geom_line(data = BEIJING02_03, mapping = aes(x = Date[327:385], y = AOTI), color = "red")

#calculate points lying outside CI
nb_error <- 0
for (i in 1:50){
  
  if (i != 15 ){
    
    if (BEIJING02_03$AOTI[i] > forecast.upper[i]){
      nb_error <- nb_error + 1
      
    } 
    
  }
  
}

for (i in 52:59){
  
  if (BEIJING02_03$AOTI[i] > forecast.upper[i]){
    nb_error <- nb_error + 1
    
  }
  
}

#error percentage
nb_error / 57 * 100
