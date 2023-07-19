
### FOR --> Lag function <-- install Hmisc
library(Hmisc)
library(R.matlab)

###############

mse <- function(X, Y) { ## X is the forecasted vector, Y is the realized vector
  seq_an = !is.na(X)
  ind = mean((Y[seq_an] - X[seq_an])^2)
  return(ind)
}



###############
#### You can change this to obtain the others R² OR do a for loop
data_test = readMat('FREDTests/FRED_MSW_recursive_stationary_1959_1999_216_3(DPCERA3M086SBEA).mat') ### !!!! you need to fix this path
####### a transformação realizada para essa série é a 5: diferença do log
############### ############### ############### ############### ############### 
or_t = log(data_test$Xt) ## como o artigo diz que as séries estão geralmente em log, vamos manter a original desse jeito

######## de acordo com o artigo, Yt+h = xt+h + Yt
## portanto, se o fim do training é t = 480, temos que Y480+1 = x480+1 + Y480
## Y480+2 = x480+2 + Y480
## Y469+12 = x469+12 + Y469
######## Vamos seguir essa lógica agora! 
### MÉTODO 1


start=data_test$test.data.start[1] # primeira predicão é no t = 481 
end= tail(data_test$test.data.end,1)
model = c()
model2 = c()

for(h in 1:12){ 
xt = data_test$bic.pred[h,] ## modelo linear
var = data_test$const.pred[h,] ## modelo constante (que consideramos ser a variancia()
xt2 = data_test$tsboost.pred[h,] ## modelo two stage boosting

xt[start:end] = or_t[(start-h):(end-h)]+(xt[start:end])
xt2[start:end] = or_t[(start-h):(end-h)]+(xt2[start:end])


var[start:end] = or_t[(start-h):(end-h)]+(var[start:end])

# ind =(1-mse(xt,or_t)/mse(var,or_t)) # const pred ################### descomentar esta linha e comentar a debaixo
########### caso queiras trocar para calcular o indicador R² usando o constant prediction em vez da variancia. O mesmo vale para os outros métodos
ind =(1-mse(xt,or_t)/ var(or_t[481:696])) # variance

model = c(model,   ifelse(ind<0,0,ind) )

# ind =(1-mse(xt2,or_t)/mse(var,or_t)) # const pred  ################### descomentar esta linha e comentar a debaixo para trocar
ind =(1-mse(xt2,or_t)/ var(or_t[481:696])) # variance

model2 = c(model2,   ifelse(ind<0,0,ind) )

}

cbind(model,model2, model2>=model)

### MÉTODO 2

#### considering always a predicted or_t apart from t-h until t = 480
model = c()
model2 = c()
for(h in 1:12){ 
  
xt = data_test$bic.pred[h,]
xt2 = data_test$tsboost.pred[h,]

var = data_test$const.pred[h,]

for(i in (start-h):(start-1)) {
  
  xt[i+h] = or_t[(i)] + (xt[i+h])
  xt2[i+h] = or_t[(i)] + (xt2[i+h])
  
  var[i+h] = or_t[(i)] + (var[i+h])
  
}
for(i in (start):(end-h)) {
  
  xt[i+h] = xt[(i)] + (xt[i+h])
  xt2[i+h] = xt2[(i)] + (xt2[i+h])
  
  var[i+h] = var[(i)] + (var[i+h])
  
}

ind =(1-mse(xt,or_t)/mse(var,or_t))
# ind =(1-mse(xt,or_t)/ var(or_t[481:696])) # variance   ################### descomentar esta linha e comentar a de cima para trocar

model = c(model,   ifelse(ind<0,0,ind) )
ind =(1-mse(xt2,or_t)/mse(var,or_t))
# ind =(1-mse(xt2,or_t)/ var(or_t[481:696])) # variance    ################### descomentar esta linha e comentar a de cima para trocar

model2 = c(model2,   ifelse(ind<0,0,ind) )

}
cbind(model,model2, model2>=model)

### MÉTODO 3

#### considering rounds

start = data_test$test.data.start
end = data_test$test.data.end
rounds=18
model = c()
model2 = c()


for(h in 1:12){
  xt = data_test$bic.pred[h,]
  var = data_test$const.pred[h,]
  xt2 = data_test$tsboost.pred[h,]
for (round in 1:rounds){
  for (i in start[round]:(start[round] + h - 1)) {
    
    xt[i] = or_t[i-h]+(xt[i])
    xt2[i] = or_t[i-h]+(xt2[i])
    
    var[i] = or_t[i-h]+(var[i])
    
  }
  if(!h==12){
    for (i in (start[round] + h):end[round]) {
      xt[i] = xt[i-h]+(xt[i])
      xt2[i] = xt[i-h]+(xt2[i])
      
      var[i] = var[i-h]+(var[i])
      
    }
  }
}
  # ind =(1-mse(xt,or_t)/ var(or_t[481:696]))# variance ##     ################### descomentar esta linha e comentar a debaixo  para trocar
  ind =(1-mse(xt,or_t)/ mse(var,or_t)) # const pred
  
  model = c(model,   ifelse(ind<0,0,ind) )
  # ind =(1-mse(xt2,or_t)/var(or_t[481:696])) #variance ##     ################### descomentar esta linha e comentar a debaixo  para trocar
  ind =(1-mse(xt2,or_t)/mse(var,or_t)) # const pred
  
  model2 = c(model2,   ifelse(ind<0,0,ind) )
  
}

cbind(model,model2, model2>=model)

### MÉTODO 4
#### considering true ytph (no level) and lagging the ytph

model = c()
model2 = c()

for(h in 1:12){
  xt = data_test$bic.pred[h,]
  xt2 = data_test$bspline.pred[h,]
  
  var = data_test$const.pred[h,]
  tyt = data_test$true.Ytph[h,]
  # ind =(1-mse(xt,tyt)/var(tyt[492:696]))# variance  ##     ################### descomentar esta linha e comentar a debaixo  para trocar
  ind =1-(mse(xt[481:695],tyt[481:695])/mse(var[481:695],tyt[481:695])) # constant pred
  
  model = c(model,   ifelse(ind<0,0,ind) )
  # ind =(1-mse(xt2,tyt)/ var(tyt[600:696]))# variance  ##     ################### descomentar esta linha e comentar a debaixo  para trocar
  ind =1-(mse(xt2[481:695],tyt[481:695])/mse(var[481:695],tyt[481:695])) #const pred
  
  model2 = c(model2,   ifelse(ind<0,0,ind) )
  
}

cbind(model,model2, model2>model)


rm(list = ls())
require(R.matlab)
dados <- readMat("FREDTests/FRED_MSW_recursive_stationary_1959_1999_216_3(DPCERA3M086SBEA).mat")
attach(dados)

hor <- 1:12
r2 <- rep(NA, 12)
ind <- seq(from = 481, to = 695, by = 1)

for (i in seq_along(hor)) {
  r2[i] <- 1 - mean((ar.pred[i, ind] - true.Ytph[i, ind])^2) /
    mean((true.Ytph[i, ind] - const.pred[i, ind])^2)
}
print(r2)
detach(dados)








x <- function(i){
  if (i < 10) warning("A warning")
  i
}

tt <- try(x(5))
tt <- try(x(5))
