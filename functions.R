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
  aic_predicted = c()
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
      aic_predicted = c(aic_predicted, AIC(model_1))
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
                       'r2' = 1-mse/varxt,
                       'aic' = aic_predicted)
  return(results_list)
}

best_linear_model = function(yt, xt, n_out_out){
  pb <- txtProgressBar(min = 0, max = 12, style = 3)
  
  df_linear_estimate = data.frame(lag_1 = NA, lag_2=NA, lag_3 = NA, lag_4 = NA, lag_5 = NA, lag_6=NA, 
                                  lag_7 = NA, lag_8 = NA,lag_9 = NA, lag_10=NA, lag_11 = NA, lag_12 = NA)
  df_best_linear = data.frame(lag_1 = NA, lag_2=NA, lag_3 = NA, lag_4 = NA, lag_5 = NA, lag_6=NA, 
                              lag_7 = NA, lag_8 = NA,lag_9 = NA, lag_10=NA, lag_11 = NA, lag_12 = NA)
  for(j in 1:12){
    #progress(j, max.value = 12, progress.bar = T)
    for(i in 1:12){
      model = linear_model(yt=yt,h=i, p=j, n_out = n_out_out)
      df_linear_estimate[i,j] <- mean(model$aic)
      df_best_linear[i,j] <- mean(model$r2)}
    setTxtProgressBar(pb, j)
    if (j==12) cat(" Done!\n")
  }
  result = list('p_linear_estimate' = as.integer(which.min(lapply(df_linear_estimate, mean))),
                'p_best_linear' = as.integer(which.max(lapply(df_best_linear, mean))))
  return(result)
}

catch_r2 = function(yt, h_max = 12, p_out=12, n_out_out=200,name = NULL, best=F){
  vetor_h = c()
  if(best){
    best = best_linear_model(yt, n_out_out = n_out_out)
    linear_estimate = catch_r2(yt, h_max = 12, p_out=best$p_linear_estimate, n_out_out=n_out_out,name = paste(name,"linear_estimate", sep = ""), best=F)
    best_linear = catch_r2(yt, h_max = 12, p_out=best$p_best_linear, n_out_out=n_out_out,name = paste(name,"best_linear", sep = ""), best=F)
    a = list('linear_estimate' = linear_estimate, 'best_linear'= best_linear)
  } 
  else{
    pb <- txtProgressBar(min = 0, max = h_max, style = 3)
    
     for(i in 1:h_max){
       model <- linear_model(yt, h=i,p=p_out, n_out = n_out_out)
       x = ifelse(model$r2 > 0, model$r2,0)
       vetor_h = c(vetor_h, x)
       setTxtProgressBar(pb, i)
       if (i==h_max) cat(" Done!\n")
       
     }
     a = list('vetor_h' = vetor_h, 'h'=1:h_max, 'lag' = p_out, 'n_predictions' = n_out_out, 'name' = name)
     
   }
  return(a)}


