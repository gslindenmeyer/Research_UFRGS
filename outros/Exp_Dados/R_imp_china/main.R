# Requer:

source('functions\\func_sarima.R')
source('functions\\func_boost-ic.R')
source('functions\\func_lasso-cv.R')
source('functions\\func_ridge-cv.R')
source('functions\\func_elast-cv.R')
source("R\\functions.R")
library(forecast)
library(mboost)
library(readr)
library(glmnet)
library(ggplot2)

#### Entrada de dados ####

X = head(read_csv('X.csv', col_names = T, progress = T)[,-c(1,2)],-1)
Y_or = read.csv('Y_or.csv')[,2]
date = as.Date(read.csv('Y_or.csv')[,4])
Y = read.csv('Y_or.csv')[-1,3]

X_lag = add_lags(X,Y='n')
Y_lag = tail(Y,-11)
Y_or_lag = tail(Y_or,-11)
names = colnames(X_lag)
colnames(X_lag) = 1:(ncol(X_lag))

#### Pré-processamento ####

ratio_start_in = .757
ratio_start_lag = 0.75

n_tot_lag <- length(Y_lag)
n_out_lag <- ceiling(n_tot_lag - ratio_start_lag*n_tot_lag)
ind_out_lag <- seq(to = n_tot_lag, by = 1, length = n_out_lag)

n_tot <- length(Y)
n_out <- ceiling(n_tot - ratio_start_in*n_tot)
ind_out <- seq(to = n_tot, by = 1, length = n_out)

lambdas <- 10^seq(2, -3, by = -.1)


# Série temporal China

date[floor((length(Y_or)-1)*ratio_start_in)]

df_y = data.frame(Date = date,Y = Y_or)

ggplot(df_y, aes(Date, Y)) +
  geom_line(size = 0.2, color = 'black') +
  geom_vline(xintercept = date[floor((length(Y_or)-1)*ratio_start_in)],
             linetype="dashed",
             color = "red", size=0.6)  + 
  labs(title="Exportações Brasil para China",
       x ="Data", y = "Log das Exportações (U$)")

##### Forecasting h=1 ####

sarima_1 = sarima_forecast(Y_or, Y, h=1, ratio_start = ratio_start_in)

lasso_1 = lassocv_forecast(Y_or = Y_or_lag, Y = Y_lag,
                          X=X_lag,lambdas = lambdas, h=1, ratio_start = 0.75)

ridge_1 = ridgecv_forecast(Y_or = Y_or_lag, Y = Y_lag,
                          X=X_lag,lambdas = lambdas, h=1, ratio_start = 0.75)

elast_1 = elastcv_forecast(Y_or = Y_or_lag, Y = Y_lag,
                           X=X_lag,lambdas = lambdas, h=1, ratio_start = 0.75,
                           alpha = 0.5)

boost_1 = boostic_forecast(Y_or = Y_or_lag,Y = Y_lag, X_lag, v = v_in,
                           h = 1, ratio_start = 0.75, Mstop = 50)

#### plots and evaluation ####

df_forecast_1 = data.frame(Date = tail(date,74), Y = tail(Y_or,74),
                           sarima = sarima_1$benchmark,
                           lasso = lasso_1$forecast, ridge = ridge_1$forecast,
                           boost = boost_1$forecast, elast = elast_1$forecast)

ggplot(df_forecast_1, aes(Date)) +
  geom_line(aes(y = Y, colour = "Y"), size = 1, color = 'black') + 
  geom_line(aes(y = sarima, colour = "SARIMA")) +
  geom_line(aes(y = lasso, colour = 'lasso')) +
  geom_line(aes(y = ridge, colour = 'ridge')) +
  geom_line(aes(y = boost, colour = 'boost')) +
  geom_line(aes(y = elast, colour = 'elast')) + 
  labs(title="Previsão Exportações Brasil para China com h = 1",
       x ="Data", y = "Log das Exportações (U$)",
       colour = 'Modelo')

sarima_e_1 = evaluation(sarima_1$benchmark, Y_or, ind_out, "sarima h=1")
lasso_e_1 = evaluation(lasso_1$forecast, Y_or_lag, ind_out_lag, "lasso h=1")
ridge_e_1 = evaluation(ridge_1$forecast, Y_or_lag, ind_out_lag, "ridge h=1")
boost_e_1 = evaluation(boost_1$forecast, Y_or_lag, ind_out_lag, "boost h=1")
elast_e_1 = evaluation(elast_1$forecast, Y_or_lag, ind_out_lag, "elast h=1")

##### Forecasting h=2 ####

sarima_2 = sarima_forecast(Y_or, Y, h=2, ratio_start = ratio_start_in)

lasso_2 = lassocv_forecast(Y_or = Y_or_lag, Y = Y_lag,
                           X=X_lag,lambdas = lambdas, h=2, ratio_start = 0.75)

ridge_2 = ridgecv_forecast(Y_or = Y_or_lag, Y = Y_lag,
                           X=X_lag,lambdas = lambdas, h=2, ratio_start = 0.75)

elast_2 = elastcv_forecast(Y_or = Y_or_lag, Y = Y_lag,
                           X=X_lag,lambdas = lambdas, h=2, ratio_start = 0.75,
                           alpha = 0.5)

boost_2 = boostic_forecast(Y_or = Y_or_lag,Y = Y_lag, X_lag, v = v_in,
                           h = 2, ratio_start = 0.75, Mstop = 50)

#### plots and evaluation ####

df_forecast_2 = data.frame(Date = tail(date,74), Y = tail(Y_or,74),
                           sarima = sarima_2$benchmark,
                           lasso = lasso_2$forecast, ridge = ridge_2$forecast,
                           boost = boost_2$forecast, elast = elast_2$forecast)

ggplot(df_forecast_1, aes(Date)) +
  geom_line(aes(y = Y, colour = "Y"), size = 1, color = 'black') + 
  geom_line(aes(y = sarima, colour = "SARIMA")) +
  geom_line(aes(y = lasso, colour = 'lasso')) +
  geom_line(aes(y = ridge, colour = 'ridge')) +
  geom_line(aes(y = boost, colour = 'boost')) +
  geom_line(aes(y = elast, colour = 'elast')) + 
  labs(title="Previsão Exportações Brasil para China com h = 2",
       x ="Data", y = "Log das Exportações (U$)",
       colour = 'Modelo')

sarima_e_2 = evaluation(sarima_2$benchmark, Y_or, ind_out, "sarima h=2")
lasso_e_2 = evaluation(lasso_2$forecast, Y_or_lag, ind_out_lag, "lasso h=2")
ridge_e_2 = evaluation(ridge_2$forecast, Y_or_lag, ind_out_lag, "ridge h=2")
boost_e_2 = evaluation(boost_2$forecast, Y_or_lag, ind_out_lag, "boost h=2")
elast_e_2 = evaluation(elast_2$forecast, Y_or_lag, ind_out_lag, "elast h=2")

#### Forecasting h=3 ####

sarima_3 = sarima_forecast(Y_or, Y, h=3, ratio_start = ratio_start_in)

lasso_3 = lassocv_forecast(Y_or = Y_or_lag, Y = Y_lag,
                           X=X_lag,lambdas = lambdas, h=3, ratio_start = 0.75)

ridge_3 = ridgecv_forecast(Y_or = Y_or_lag, Y = Y_lag,
                           X=X_lag,lambdas = lambdas, h=3, ratio_start = 0.75)

elast_3 = elastcv_forecast(Y_or = Y_or_lag, Y = Y_lag,
                           X=X_lag,lambdas = lambdas, h=3, ratio_start = 0.75,
                           alpha = 0.5)

boost_3 = boostic_forecast(Y_or = Y_or_lag,Y = Y_lag, X_lag, v = v_in,
                           h = 3, ratio_start = 0.75, Mstop = 50)

#### plots and evaluation ####

df_forecast_3 = data.frame(Date = tail(date,74), Y = tail(Y_or,74),
                           sarima = sarima_3$benchmark,
                           lasso = lasso_3$forecast, ridge = ridge_3$forecast,
                           boost = boost_3$forecast, elast = elast_3$forecast)

ggplot(df_forecast_1, aes(Date)) +
  geom_line(aes(y = Y, colour = "Y"), size = 1, color = 'black') + 
  geom_line(aes(y = sarima, colour = "SARIMA")) +
  geom_line(aes(y = lasso, colour = 'lasso')) +
  geom_line(aes(y = ridge, colour = 'ridge')) +
  geom_line(aes(y = boost, colour = 'boost')) +
  geom_line(aes(y = elast, colour = 'elast')) + 
  labs(title="Previsão Exportações Brasil para China com h = 3",
       x ="Data", y = "Log das Exportações (U$)",
       colour = 'Modelo')

sarima_e_3 = evaluation(sarima_3$benchmark, Y_or, ind_out, "sarima h=3")
lasso_e_3 = evaluation(lasso_3$forecast, Y_or_lag, ind_out_lag, "lasso h=3")
ridge_e_3 = evaluation(ridge_3$forecast, Y_or_lag, ind_out_lag, "ridge h=3")
boost_e_3 = evaluation(boost_3$forecast, Y_or_lag, ind_out_lag, "boost h=3")
elast_e_3 = evaluation(elast_3$forecast, Y_or_lag, ind_out_lag, "elast h=3")

