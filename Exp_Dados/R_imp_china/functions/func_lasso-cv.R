# Lasso functions

lassocv_forecast = function(Y_or, Y, X, lambdas, h, ratio_start = 0.8) {
  n_tot <- length(Y)
  n_out <- ceiling(n_tot - ratio_start*n_tot)
  ind_out <- seq(to = n_tot, by = 1, length = n_out)
  Y_predicted = c(Y_or[ind_out[1]])
  varimp_df = data.frame(rep(0, (ncol(X)+1)))
  selected_var = c()
  for(i in 1:n_out){
    ind_in <- seq(from = 1, to = ind_out[i] - h, by = 1)
    y_extra = c()
    x_reg <- as.matrix(as.data.frame(X[head(ind_in,-1),])) # x independent t = 1, ..., T.in-h
    x0_reg <- as.matrix(X[tail(ind_in,1),])
    for(j in 1:h) {
      # expanding window
      y_dep <- append(Y[tail(ind_in,-j)], y_extra)
      y_reg <- as.matrix(y_dep)
      # finding m*
      lasso_reg <- cv.glmnet(x_reg, y_reg,
                             alpha = 1, lambda = lambdas,
                             standardize = TRUE, nfolds = 5)
      # Best 
      lambda_best <- lasso_reg$lambda.min
      lasso_model <- glmnet(x_reg, y_reg, alpha = 1, lambda = lambda_best,
                            standardize = TRUE)
      predictions_test <- predict(lasso_model, s = lambda_best,
                                  newx = x0_reg)
      y_predicted = unname(predictions_test)
      cat("Selected lambda is: ", lambda_best, "\n")
      # visualizing selected predictors varimp
      #varimp_df_partial = data.frame(varimp(model_1))
      #sum_reduction = sum(varimp_df_partial[,1])
      #varimp_partial = varimp_df_partial[,1]/sum_reduction
      #varimp_df = cbind(varimp_df, varimp_partial)
      # visualizing selected predictors frequency
      #selected_var = append(selected_var, list(model_1$xselect()))
      # output
      y_extra = append(y_extra, y_predicted)
    }
    Y_predicted = append(Y_predicted, Y_or[(ind_out[1]+i-(h))] + sum(y_extra))
    print(i/n_out)
  }
  results <- list(forecast = Y_predicted
                  #varimp = varimp_df[,-1],
                  #selected = selected_var
  )
  return(results)
}