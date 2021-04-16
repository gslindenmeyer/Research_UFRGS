require(mboost)
require(forecast)
require(zoo)

p = 12 #lags from 1 to 12
yt = array(dim = 500) ## Criação do vetor 1 linear
t = 500
yt[c(1,2,3)] = rnorm(3, mean = 0, sd = 0.01)
for(i in 4:t){
yt[i] = 0.21*yt[i-1] + 0.35*yt[i-2] + 0.17*yt[i-3] + 0.1*rnorm(1, mean = 0, sd = 0.1)
    
}
plot.ts(yt)

yt2 = array(dim = 500) ## Criação do vetor 2 não linear
t = 500
yt2[c(1)] = rnorm(1, mean = 0, sd = 0.25)
for(i in 2:t){
  yt2[i] = 0.4*(5-yt2[i-1]^2)/(1+yt2[i-1]^2)+rnorm(1, mean = 0, sd = 0.5)

}
plot.ts(yt2)


  
  

df = data.frame(lag_1 = NA, lag_2=NA, lag_3 = NA, lag_4 = NA, lag_5 = NA, lag_6=NA, 
               lag_7 = NA, lag_8 = NA,lag_9 = NA, lag_10=NA, lag_11 = NA, lag_12 = NA)
df = df[,1:p]
df2 = df

for(i in 1:p){ 
#  cat((i/p)*100, '%\n')
  if(p>1){
  df[i:(length(yt)-1),i] <- head(yt,-i)
  } else {
    df[i:(length(yt)-1)] <- head(yt,-i)
  }
}
#
#df = df[,1]
#
df <- na.exclude(df, )
rownames(df) <- NULL
for(i in 1:p){ 
  #  cat((i/p)*100, '%\n')
  if(p>1){
    df2[i:(length(yt2)-1),i] <- head(yt2,-i)
  } else {
    df2[i:(length(yt2)-1)] <- head(yt2,-i)
  }
}
df2 <- na.exclude(df2)
rownames(df2) <- NULL

ytp <- tail(yt,-p)
ytp2 <- tail(yt2,-p)

ratio = 0.603

n_tot <- length(ytp)
n_out <- ceiling(n_tot - ratio*n_tot)
ind_out <- seq(to = n_tot, by = 1, length = n_out)
Y_predicted = c()

data <- cbind(y_reg = ytp, df)

for(i in 1:n_out){
  ind_in <- seq(from = 1, to = ind_out[i] - h, by = 1)
  if(p>1){
  x_reg <- data[ind_in,colnames(df)] # x independent t = 1, ..., T.in-h
  x0_reg <- data[ind_out[i],colnames(df)]
  } 
  else { 
    x_reg <- data.frame(lag_1=df[head(ind_in,-1)]) # x independent t = 1, ..., T.in-h
    x0_reg <- data.frame(lag_1=df[tail(ind_in,1)])  
    }
  for(j in 1:h) {
    # expanding window
    y_dep <- ytp[ind_in] #[tail(ind_in,-(h))]
    data2 <- cbind(y_reg = y_dep, x_reg) 
    #y_reg <- as.matrix(y_dep)
    # finding m*
    model_1 = lm(y_reg ~., data=data2)
    y_predicted = predict(model_1, newdata=x0_reg)
    # visualizing selected predictors varimp
    
    #y_extra = append(y_extra, y_predicted)
  }
  Y_predicted = append(Y_predicted, y_predicted)
  print(i/n_out)
}

plot.ts(tail(ytp,length(Y_predicted)))
lines(Y_predicted, col = 'red')
mse <- mean((tail(ytp,length(Y_predicted))-Y_predicted)^2)
varxt <- var(tail(ytp,length(Y_predicted)))
1-mse/varxt

###
Y_predicted = c()

data <- cbind(y_reg = ytp2, df2)

for(i in 1:n_out){
  ind_in <- seq(from = 1, to = ind_out[i] - h, by = 1)
  if(p>1){
    x_reg <- data[ind_in,colnames(df2)] # x independent t = 1, ..., T.in-h
    x0_reg <- data[ind_out[i],colnames(df2)]
  } 
  else { 
    x_reg <- data.frame(lag_1=df2[head(ind_in,-1)]) # x independent t = 1, ..., T.in-h
    x0_reg <- data.frame(lag_1=df2[tail(ind_in,1)])  
  }
  for(j in 1:h) {
    # expanding window
    y_dep <- ytp2[ind_in] #[tail(ind_in,-(h))]
    data2 <- cbind(y_reg = y_dep, x_reg) 
    #y_reg <- as.matrix(y_dep)
    # finding m*
    model_1 = lm(y_reg ~., data=data2)
    y_predicted = predict(model_1, newdata=x0_reg)
    # visualizing selected predictors varimp
    
    #y_extra = append(y_extra, y_predicted)
  }
  Y_predicted = append(Y_predicted, y_predicted)
  print(i/n_out)
}

plot.ts(tail(ytp2,length(Y_predicted)))
lines(Y_predicted, col = 'red')
mse <- mean((tail(ytp2,length(Y_predicted))-Y_predicted)^2)
varxt <- var(tail(ytp2,length(Y_predicted)))
1-mse/varxt




