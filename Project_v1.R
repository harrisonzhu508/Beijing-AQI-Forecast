library(readr)
library(car)
library(tseries)

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
