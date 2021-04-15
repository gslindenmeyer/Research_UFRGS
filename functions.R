linear_model = function(yt, xt, h, p, n_out=200){
  ## xt nao faz nada ainda
  
  df = data.frame(lag_1 = NA, lag_2=NA, lag_3 = NA, lag_4 = NA, lag_5 = NA, lag_6=NA, 
                  lag_7 = NA, lag_8 = NA,lag_9 = NA, lag_10=NA, lag_11 = NA, lag_12 = NA)
  df = df[,1:p]
  for(i in 1:p){ 
    #  cat((i/p)*100, '%\n')
    if(p>1){
      df[i:(length(yt)-1),i] <- head(yt,-i)
    } else {
      df[i:(length(yt)-1)] <- head(yt,-i)
    }
  }
  
  df <- na.exclude(df, )
  rownames(df) <- NULL
  if(p==1){
    df = data.frame(lag_1=df)
  }
  
  ytp <- tail(yt,-p)
  
  n_tot <- length(ytp)
  #n_out <- ceiling(n_tot - ratio*n_tot)
  ind_out <- seq(to = n_tot, by = 1, length = n_out)
  Y_predicted = c()
  data <- cbind(y_reg = ytp, df)
  for(i in 1:n_out){
    ind_in <- seq(from = 1, to = ind_out[i] - h, by = 1)
    if(p>1){
      x_reg <- data[ind_in,colnames(df)] # x independent t = 1, ..., T.in-h
      x0_reg <- data[(ind_out[i]-h+1),colnames(df)]
    } 
    else { 
      x_reg <- data.frame(lag_1=data[ind_in,colnames(df)])
      x0_reg <- data.frame(lag_1=data[(ind_out[i]-h+1),colnames(df)])
      row.names(x0_reg) <- (ind_out[i]-h+1)
    }
    for(j in 1:1) {
      # expanding window
      if(h>1){
        y_dep <- ytp[tail((h-1):(length(ind_in)+h-1),-j)] #ytp[ind_in] #[tail(ind_in,-(h))]  
      }else{
        y_dep <- ytp[ind_in] #ytp[ind_in] #[tail(ind_in,-(h))]
      }
      
      data2 <- as.data.frame(cbind(y_reg = y_dep, x_reg)) 
      #y_reg <- as.matrix(y_dep)
      # finding m*
      model_1 = lm(y_reg ~., data=data2)
      y_predicted = predict(model_1, newdata=x0_reg)
      # visualizing selected predictors varimp
      
      #y_extra = append(y_extra, y_predicted)
    }
    Y_predicted = append(Y_predicted, y_predicted)
    #print(i/n_out)
  }
  mse <- mean((tail(ytp,length(Y_predicted))-Y_predicted)^2)
  varxt <- var(tail(ytp,length(Y_predicted)))
  
  results_list <- list('y_hat' = Y_predicted,
                       'y_filtered' = tail(ytp,length(Y_predicted)),
                       'mse' = mse,
                       'varxt'= varxt,
                       'r2' = 1-mse/varxt)
  return(results_list)
}

