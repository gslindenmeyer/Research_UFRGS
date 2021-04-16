# SARIMA Forecast

sarima_forecast = function(Y_or, Y, h, ratio_start = 0.8, arima = c(1,1,0),
                           seasonal = c(1,1,0)) {
  n_tot <- length(Y)
  n_out <- ceiling(n_tot - ratio_start*n_tot)
  ind_out <- seq(to = n_tot, by = 1, length = n_out)
  Y_arima = c(Y_or[ind_out[1]])
  for(i in 1:n_out){
    ind_in <- seq(from = 1, to = ind_out[i] - h, by = 1)
    bench = arima(exp(Y_or[ 1:(ind_out[i] - h + 1) ]), arima
                  , seasonal = list(order = seasonal, period = 12)
    )
    forecast_bench = forecast(bench, h)
    y_predicted_bench = forecast_bench$mean[h]
    y_predicted_arima = log(y_predicted_bench)
    Y_arima = append(Y_arima, (y_predicted_arima))
    print(i/n_out)
  }
  results <- list(benchmark = Y_arima)
  return(results)
}