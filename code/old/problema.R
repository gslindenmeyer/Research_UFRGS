
### FOR --> Lag function <-- install Hmisc
library(Hmisc)
library(R.matlab)

###############

mse <- function(X, Y) { ## X is the forecasted vector, Y is the realized vector
  seq_an = !is.na(X)
  ind = mean((Y[seq_an] - X[seq_an])^2)
  return(ind)
}

r_squared <- function(X, Y){
  seq_an = !is.na(X)
  mse2 = mse(X,Y)
  var2 = var(Y[seq_an])
  
  ind = 1 - (mse2/var2)
  return(ifelse(ind < 0, 0,ind ))
}


###############
h=1 #### you will fix an H for the predictions. You will receive one R² only per H.
#### You can change this to obtain the others R² OR do a for loop
data_test = readMat('others/FREDTests/FRED_MSW_recursive_stationary_1959_1999_204_1(RPI).mat') ### !!!! you need to fix this path

yt = data_test$true.Ytph[h,]
or_t = data_test$Xt

xt = data_test$bic.pred[h,]



### Initial tests
yt=Lag(yt,h)

yt[1:h] = or_t[1:h]
yt[(h+1):length(yt)] = or_t[1:(length(yt)-h)]*exp(yt[(h+1):length(yt)])
plot.ts(yt)
plot.ts(or_t)
## you can see that they are very identical

#### We want to find the same R² curve from RPI. We propose four methods

#### 1) considering always a true or_t

xt = data_test$bic.pred[h,]
start=481
xt[(start-h):(start-1)] = or_t[(start-h):(start-1)]
xt[(start):length(xt)] = or_t[(start-h):(length(xt)-h)]*exp(xt[(start):length(xt)])

plot.ts(or_t)
lines(xt, col = 'red')

r_squared(xt,or_t) ### R_squared is too good

#### considering always a predicted or_t

xt = data_test$bic.pred[h,]
start=481
xt[(start-h):(start-1)] = or_t[(start-h):(start-1)]
for(i in (start):length(xt)) {
xt[i] = xt[(i-h)]*exp(xt[i])
}
plot.ts(or_t)
lines(xt, col = 'red')

r_squared(xt,or_t) ### R_squared is tooo bad


#### considering rounds
xt = data_test$bic.pred[h,]

start = data_test$test.data.start
end = data_test$test.data.end
rounds=17
for (round in 1:rounds){
  for (i in start[round]:(start[round] + h - 1)) {
  
  xt[i] = or_t[i-h]*exp(xt[i])
  }
  if(!h==12){
    for (i in (start[round] + h):end[round]) {
      xt[i] = xt[i-h]*exp(xt[i])
      
    }
  }
}

plot.ts(or_t)
lines(xt, col = 'red')

1 - mse(xt,or_t)/var(or_t[start[1]:end[17]]) ## r_squared is too good
r_squared(xt,or_t)

#### considering true ytph (no level) and lagging the ytph

tsboost = c()
for(h in 1:12){
xt = data_test$tsboost.pred[h,]
tyt = Lag(data_test$true.Ytph[h,],h)

# plot.ts(tyt)
# lines(xt, col = 'red')

tsboost = c(tsboost, r_squared(xt,tyt))
}

bic = c()
for(h in 1:12){
  xt = data_test$bic.pred[h,]
  tyt = Lag(data_test$true.Ytph[h,],h)
  
  # plot.ts(tyt)
  # lines(xt, col = 'red')
  
  bic = c(bic, r_squared(xt,tyt))
}

cbind(tsboost,bic, tsboost>=bic)

