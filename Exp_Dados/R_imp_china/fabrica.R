# Indústria de Forecastings

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

X = read_csv('X.csv', col_names = T, progress = T)[,-c(1,2)]
Y_or = read.csv('Y_or.csv')[,3]
date = as.Date(read.csv('Y_or.csv')[,1])
Y = read.csv('Y_or.csv')[-1,4]

X_lag = add_lags(X,Y)
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

reg_names = c('china_imp', 'china_exp', 'asia_imp', 'asia_exp',
              'euro_imp', 'euro_exp', 'imp', 'exp')

for (i in reg_names) {
  
  #### Entrada de dados ####
  
  Y_or = read.csv(paste('Y_or_',i,'.csv',sep=''))[,3]
  date = as.Date(read.csv(paste('Y_or_',i,'.csv',sep=''))[,1])
  Y = read.csv(paste('Y_or_',i,'.csv',sep=''))[-1,4]
  
  X_lag = add_lags(X,Y)
  Y_lag = tail(Y,-11)
  Y_or_lag = tail(Y_or,-11)
  names = colnames(X_lag)
  colnames(X_lag) = 1:(ncol(X_lag))
  
  ##### Forecasting h=1 ####
  
  sarima_1 = sarima_forecast(Y_or, Y, h=1, ratio_start = ratio_start_in)
  saveRDS(sarima_1, paste('sarima_1_', i, '.RData', sep=''))
  
  lasso_1 = lassocv_forecast(Y_or = Y_or_lag, Y = Y_lag,
                             X=X_lag,lambdas = lambdas, h=1, ratio_start = 0.75)
  saveRDS(lasso_1, paste('lasso_1_', i, '.RData', sep=''))
  
  ridge_1 = ridgecv_forecast(Y_or = Y_or_lag, Y = Y_lag,
                             X=X_lag,lambdas = lambdas, h=1, ratio_start = 0.75)
  saveRDS(ridge_1, paste('ridge_1_', i, '.RData', sep=''))
  
  elast_1 = elastcv_forecast(Y_or = Y_or_lag, Y = Y_lag,
                             X=X_lag,lambdas = lambdas, h=1, ratio_start = 0.75,
                             alpha = 0.5)
  saveRDS(elast_1, paste('elast_1_', i, '.RData', sep=''))
  
  boost_1 = boostic_forecast(Y_or = Y_or_lag,Y = Y_lag, X_lag, v = v_in,
                             h = 1, ratio_start = 0.75, Mstop = 50)
  saveRDS(boost_1, paste('boost_1_', i, '.RData', sep=''))
  
  ##### Forecasting h=2 ####
  
  sarima_2 = sarima_forecast(Y_or, Y, h=2, ratio_start = ratio_start_in)
  saveRDS(sarima_2, paste('sarima_2_', i, '.RData', sep=''))
  
  lasso_2 = lassocv_forecast(Y_or = Y_or_lag, Y = Y_lag,
                             X=X_lag,lambdas = lambdas, h=2, ratio_start = 0.75)
  saveRDS(lasso_2, paste('lasso_2_', i, '.RData', sep=''))
  
  ridge_2 = ridgecv_forecast(Y_or = Y_or_lag, Y = Y_lag,
                             X=X_lag,lambdas = lambdas, h=2, ratio_start = 0.75)
  saveRDS(ridge_2, paste('ridge_2_', i, '.RData', sep=''))
  
  elast_2 = elastcv_forecast(Y_or = Y_or_lag, Y = Y_lag,
                             X=X_lag,lambdas = lambdas, h=2, ratio_start = 0.75,
                             alpha = 0.5)
  saveRDS(elast_2, paste('elast_2_', i, '.RData', sep=''))
  
  boost_2 = boostic_forecast(Y_or = Y_or_lag,Y = Y_lag, X_lag, v = v_in,
                             h = 2, ratio_start = 0.75, Mstop = 50)
  saveRDS(boost_2, paste('boost_2_', i, '.RData', sep=''))
  
  #### Forecasting h=3 ####
  
  sarima_3 = sarima_forecast(Y_or, Y, h=3, ratio_start = ratio_start_in)
  saveRDS(sarima_3, paste('sarima_3_', i, '.RData', sep=''))
  
  lasso_3 = lassocv_forecast(Y_or = Y_or_lag, Y = Y_lag,
                             X=X_lag,lambdas = lambdas, h=3, ratio_start = 0.75)
  saveRDS(lasso_3, paste('lasso_3_', i, '.RData', sep=''))
  
  ridge_3 = ridgecv_forecast(Y_or = Y_or_lag, Y = Y_lag,
                             X=X_lag,lambdas = lambdas, h=3, ratio_start = 0.75)
  saveRDS(ridge_3, paste('ridge_3_', i, '.RData', sep=''))
  
  elast_3 = elastcv_forecast(Y_or = Y_or_lag, Y = Y_lag,
                             X=X_lag,lambdas = lambdas, h=3, ratio_start = 0.75,
                             alpha = 0.5)
  saveRDS(elast_3, paste('elast_3_', i, '.RData', sep=''))
  
  boost_3 = boostic_forecast(Y_or = Y_or_lag,Y = Y_lag, X_lag, v = v_in,
                             h = 3, ratio_start = 0.75, Mstop = 50)
  saveRDS(boost_3, paste('boost_3_', i, '.RData', sep=''))
  
}
