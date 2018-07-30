# Beijing-AQI-Forecast
Time Series forecasting of AQI values in Beijing

In collaboration with 林欣诺 (Lin, Xinnuo), we were able to give a prediction interval of the air quality index (AQI) in Beijing from March 2018 to August 2018. 

## Methods

We used ARIMA($p,d,q$) models to model the AQI. We also tested GARCH, SARIMA and ARMA-GARCH models but a relatively simple ARIMA($p,d,q$) model was sufficient. For more details, please refer to ``Project Report.pdf``.
