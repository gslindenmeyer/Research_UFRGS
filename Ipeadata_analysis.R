#############################################################
## Here we do the analysis
## gslindenmeyer@gmail.com - Nov 2021
##
## -based on H. Kauppi and T. Virtanen code
#############################################################

##############################
# 0 - Load libraries
##############################
rm(list = ls())
set.seed(123)

require(mboost)
require(R.matlab)
require(forecast)
require(sandwich)
require(gridExtra)
library(ggpubr)
library(afmtools)
library(readxl)
library(dplyr)
library(readr)
# data <- readMat("FRED.mat")
source("helper_functions_analysis.R")

##
library(ggplot2)
library(ggthemes)
library(paletteer)

## Dataset

load("data/dataset.RData")


M <- 300

filenames <- list.files(paste("IPEAtests2", M, "/", sep = ""))
x <- c()

tsboost <- c()
tsboost_aic <- c()

tsboost_noextra <- c()
tsboost_noextra_aic <- c()

bspline_aic <- c()
bspline <- c()

bspline_noextra <- c()
bspline_noextra_aic <- c()

bols <- c()
bols_aic <- c()

tree <- c()
tree2 <- c()

for (name in filenames) {
  data_test <- readMat(paste("IPEAtests2", M, "/", name, sep = ""))
  r_benchmark <- c()
  r_nonlinear <- c()
  tsboost_in <- c()
  tsboost_noextra_in <- c()
  bspline_in <- c()
  bspline_noextra_in <- c()
  bols_in <- c()
  tree_in <- c()

  tsboost_in2 <- c()
  tsboost_noextra_in2 <- c()
  bspline_in2 <- c()
  bspline_noextra_in2 <- c()
  bols_in2 <- c()
  tree_in2 <- c()

  start <- data_test$test.data.start[1]
  end <- tail(data_test$test.data.end, 1)


  for (h in 1:12) {
    xt <- data_test$ar.pred[h, ]
    xt2 <- data_test$tsboost.pred[h, ]
    xt3 <- data_test$tsboost.noextra.pred[h, ]
    xt4 <- data_test$bspline.pred[h, ]
    xt5 <- data_test$bspline.noextra.pred[h, ]
    xt6 <- data_test$bols.pred[h, ]
    xt7 <- data_test$tree.pred[h, ]

    xt2a <- data_test$tsboost.pred.aic[h, ]
    xt3a <- data_test$tsboost.noextra.pred.aic[h, ]
    xt4a <- data_test$bspline.pred.aic[h, ]
    xt5a <- data_test$bspline.noextra.pred.aic[h, ]
    xt6a <- data_test$bols.pred.aic[h, ]
    xt7a <- data_test$tree2.pred[h, ]

    var <- data_test$const.pred[h, ]
    yt <- data_test$true.Ytph[h, ]

    tl <- start:(end - h)


    # ind <- 1 - (mse(xt[tl], yt[tl]) / mse(var[tl], yt[tl])) # constant pred

    # r_benchmark <- c(r_benchmark, ifelse(ind < 0, 0, ind))
    # ind <- 1 - (mse(xt2[tl], yt[tl]) / mse(var[tl], yt[tl])) # const pred

    # r_nonlinear <- c(r_nonlinear, ifelse(ind < 0, 0, ind))
    tsboost_in <- c(tsboost_in, mse(xt[tl], yt[tl]) / mse(xt2[tl], yt[tl]))
    tsboost_noextra_in <- c(tsboost_noextra_in, mse(xt[tl], yt[tl]) / mse(xt3[tl], yt[tl]))
    bspline_in <- c(bspline_in, mse(xt[tl], yt[tl]) / mse(xt4[tl], yt[tl]))
    bspline_noextra_in <- c(bspline_noextra_in, mse(xt[tl], yt[tl]) / mse(xt5[tl], yt[tl]))
    bols_in <- c(bols_in, mse(xt[tl], yt[tl]) / mse(xt6[tl], yt[tl]))
    tree_in <- c(tree_in, mse(xt[tl], yt[tl]) / mse(xt7[tl], yt[tl]))

    tsboost_in2 <- c(tsboost_in, mse(xt[tl], yt[tl]) / mse(xt2a[tl], yt[tl]))
    tsboost_noextra_in2 <- c(tsboost_noextra_in2, mse(xt[tl], yt[tl]) / mse(xt3a[tl], yt[tl]))
    bspline_in2 <- c(bspline_in2, mse(xt[tl], yt[tl]) / mse(xt4a[tl], yt[tl]))
    bspline_noextra_in2 <- c(bspline_noextra_in2, mse(xt[tl], yt[tl]) / mse(xt5a[tl], yt[tl]))
    bols_in2 <- c(bols_in2, mse(xt[tl], yt[tl]) / mse(xt6a[tl], yt[tl]))
    tree_in2 <- c(tree_in2, mse(xt[tl], yt[tl]) / mse(xt7a[tl], yt[tl]))
  }
  tsboost <- c(tsboost, mean(tsboost_in))
  tsboost_noextra <- c(tsboost_noextra, mean(tsboost_noextra_in))
  bspline <- c(bspline, mean(bspline_in))
  bspline_noextra <- c(bspline_noextra, mean(bspline_noextra_in))
  bols <- c(bols, mean(bols_in))
  tree <- c(tree, mean(tree_in))

  tsboost_aic <- c(tsboost_aic, mean(tsboost_in2))
  tsboost_noextra_aic <- c(tsboost_noextra_aic, mean(tsboost_noextra_in2))
  bspline_aic <- c(bspline_aic, mean(bspline_in2))
  bspline_noextra_aic <- c(bspline_noextra_aic, mean(bspline_noextra_in2))
  bols_aic <- c(bols_aic, mean(bols_in2))
  tree2 <- c(tree2, mean(tree_in2))
}
sum(tsboost < 1) / length(tsboost)
sum(tsboost_noextra < 1) / length(tsboost_noextra)
sum(bspline < 1) / length(bspline)
sum(bspline_noextra < 1) / length(bspline_noextra)
sum(bols < 1) / length(bols)
sum(tree < 1) / length(tree)

sum(tsboost_aic < 1) / length(tsboost_aic)
sum(tsboost_noextra_aic < 1) / length(tsboost_noextra_aic)
sum(bspline_aic < 1) / length(bspline_aic)
sum(bspline_noextra_aic < 1) / length(bspline_noextra_aic)
sum(bols_aic < 1) / length(bols_aic)
sum(tree2 < 1) / length(tree2)


## Is aic better ?
sum(tsboost_aic < tsboost) / length(tsboost_aic)
sum(tsboost_noextra_aic < tsboost_noextra) / length(tsboost_noextra_aic)
sum(bspline_aic < bspline) / length(bspline_aic)
sum(bspline_noextra_aic < bspline_noextra) / length(bspline_noextra_aic)
sum(bols_aic < bols) / length(bols_aic)
sum(tree2 < tree) / length(tree2)



df <- data.frame(
  "h" = 1:12, "R2_ar" = r_linear,
  "R2_tsboost" = r_nonlinear
)
plot <- ggplot(data = df, mapping = aes(x = h)) +
  geom_line(aes(y = R2_ar, color = "darkred")) +
  geom_point(aes(y = R2_ar), shape = 21, colour = "black") +
  geom_line(aes(y = R2_tsboost, color = "steelblue")) +
  geom_point(aes(y = R2_tsboost), shape = 21, colour = "black") +
  ylab("R²") +
  # ggtitle(paste("model_", substring(models[[i]]$best_linear$name, 1, 1), sep = "")) +
  ylim(0, 1) + # scale_y_continuous(labels=as.character(seq(0,1, by=0.1)),breaks=seq(0,1, by=0.1)) +
  scale_x_continuous(labels = as.character(1:12), breaks = 1:12) +
  scale_color_discrete(name = "Models", labels = c("AR", "TS boost"))
#+theme(legend.position='bottom')
plot


### Average R² for out of sample forecasts. # 3 different cases, R² min = 0 and min = .1
M <- 300

filenames <- list.files(paste("IPEAtests2", M, "/", sep = ""))

linear <- rep(0, 140)
bols <- rep(0, 140)
bspline <- rep(0, 140)
tsboost <- rep(0, 140)
bspline2 <- rep(0, 140)
tsboost2 <- rep(0, 140)
tree  <- rep(0, 140)
all_cases <- list()

for (i in 1:12) {
  all_cases[[i]] <- data.frame(linear, bols, bspline, tsboost, bspline2, tsboost2, tree)
}

for (h in 1:12) {
  for (i in 1:length(filenames)) {
    data_test <- readMat(paste("IPEAtests2", M, "/", filenames[i], sep = ""))


    start <- data_test$test.data.start[1]
    end <- tail(data_test$test.data.end, 1)


    xt <- data_test$ar.pred[h, ] # linear
    xt2 <- data_test$tsboost.pred[h, ] # tsboost
    xt3 <- data_test$bspline.pred[h, ] # bspline
    xt4 <- data_test$bols.pred[h, ] # bols
    
    xt5 <- data_test$bspline.noextra.pred[h,]
    xt6 <- data_test$tsboost.noextra.pred[h,]
    xt7 <- data_test$tree2.pred[h,]

    var <- data_test$const.pred[h, ] # for var
    yt <- data_test$true.Ytph[h, ] # true value

    tl <- start:(end - h)


    ind <- 1 - mean((xt4[tl] - yt[tl])^2) /
      mean((yt[tl] - var[tl])^2)
    all_cases[[h]][i, "bols"] <- ifelse(ind < 0, 0, ind)

    ind <- 1 - mean((xt[tl] - yt[tl])^2) /
      mean((yt[tl] - var[tl])^2)
    all_cases[[h]][i, "linear"] <- ifelse(ind < 0, 0, ind)

    ind <- 1 - mean((xt3[tl] - yt[tl])^2) /
      mean((yt[tl] - var[tl])^2)
    all_cases[[h]][i, "bspline"] <- ifelse(ind < 0, 0, ind)

    ind <- 1 - mean((xt2[tl] - yt[tl])^2) /
      mean((yt[tl] - var[tl])^2)
    all_cases[[h]][i, "tsboost"] <- ifelse(ind < 0, 0, ind)
    
    ind <- 1 - mean((xt5[tl] - yt[tl])^2) /
      mean((yt[tl] - var[tl])^2)
    all_cases[[h]][i, "bspline2"] <- ifelse(ind < 0, 0, ind)
    
    ind <- 1 - mean((xt6[tl] - yt[tl])^2) /
      mean((yt[tl] - var[tl])^2)
    all_cases[[h]][i, "tsboost2"] <- ifelse(ind < 0, 0, ind)
    
    ind <- 1 - mean((xt7[tl] - yt[tl])^2) /
      mean((yt[tl] - var[tl])^2)
    all_cases[[h]][i, "tree"] <- ifelse(ind < 0, 0, ind)
    
    
  }
}

all <- data.frame(
  "linear" = rep(0, 12),
  "bols" = rep(0, 12),
  "bspline" = rep(0, 12),
  "tsboost" = rep(0, 12),
  "bspline2" = rep(0, 12),
  "tsboost2" = rep(0, 12),
  "tree" = rep(0, 12)
  
)

all1 <- all

spl <- all
spl1 <- all

lin <- all
lin1 <- all



N0 = rep(140,12)
N1 = c()

for (i in 1:12) {
  base = all_cases[[i]]
  
  all[i, ] <- sapply(base, mean)
  all1[i, ] <- sapply(base, mean_01)
  N1 = c(N1, sum(base['tsboost']>.1)) 
  
}


N2 = c()
N3 = c()
N4 = c()
N5 = c()
for (i in 1:12) {
  ## Any of the splines method is best of linear cases
  base = all_cases[[i]]
  cut = base['linear'] < base['tsboost'] | base['linear'] < base['bspline']
  
    spl[i,'linear'] = mean(base[cut,'linear'])
    spl[i,'bols'] = mean(base[cut,'bols'])
    spl[i,'bspline'] = mean(base[cut,'bspline'])
    spl[i,'tsboost'] = mean(base[cut,'tsboost'])
    spl[i,'bspline2'] = mean(base[cut,'bspline2'])
    spl[i,'tsboost2'] = mean(base[cut,'tsboost2'])
    spl[i,'tree'] = mean(base[cut,'tree'])
    
  N2 = c(N2, sum(cut))
  
  cut = (base['linear'] < base['tsboost'] | base['linear'] < base['bspline']) & base['linear'] >=0.1 & base['bols'] >=0.1 & base['bspline'] >=0.1 & base['tsboost'] >=0.1
  
  spl1[i,'linear'] = mean(base[cut,'linear'])
  spl1[i,'bols'] = mean(base[cut,'bols'])
  spl1[i,'bspline'] = mean(base[cut,'bspline'])
  spl1[i,'tsboost'] = mean(base[cut,'tsboost'])
  spl1[i,'bspline2'] = mean(base[cut,'bspline2'])
  spl1[i,'tsboost2'] = mean(base[cut,'tsboost2'])
  spl1[i,'tree'] = mean(base[cut,'tree'])
  
  N3 = c(N3, sum(cut))
  
  ## Linear methods are better 
  
  cut = !(base['linear'] < base['tsboost'] | base['linear'] < base['bspline'])
  
  lin[i,'linear'] = mean(base[cut,'linear'])
  lin[i,'bols'] = mean(base[cut,'bols'])
  lin[i,'bspline'] = mean(base[cut,'bspline'])
  lin[i,'tsboost'] = mean(base[cut,'tsboost'])
  lin[i,'bspline2'] = mean(base[cut,'bspline2'])
  lin[i,'tsboost2'] = mean(base[cut,'tsboost2'])
  lin[i,'tree'] = mean(base[cut,'tree'])
  
  
  N4 = c(N4, sum(cut))
  
  cut = !(base['linear'] < base['tsboost'] | base['linear'] < base['bspline']) & base['linear'] >=0.1 & base['bols'] >=0.1 & base['bspline'] >=0.1 & base['tsboost'] >=0.1
  
  lin1[i,'linear'] = mean(base[cut,'linear'])
  lin1[i,'bols'] = mean(base[cut,'bols'])
  lin1[i,'bspline'] = mean(base[cut,'bspline'])
  lin1[i,'tsboost'] = mean(base[cut,'tsboost'])
  lin1[i,'bspline2'] = mean(base[cut,'bspline2'])
  lin1[i,'tsboost2'] = mean(base[cut,'tsboost2'])
  lin1[i,'tree'] = mean(base[cut,'tree'])
  
  N5 = c(N5, sum(cut))
  
  
}

sapply(lin, mean)
sapply(lin1, mean)

clipr::write_clip(t(all))
clipr::write_clip(t(N0))

clipr::write_clip(t(all1))
clipr::write_clip(t(N1))

clipr::write_clip(t(spl))
clipr::write_clip(t(N2))

clipr::write_clip(t(spl1))
clipr::write_clip(t(N3))

clipr::write_clip(t(lin))
clipr::write_clip(t(N4))

clipr::write_clip(t(lin1))
clipr::write_clip(t(N5))


### rRMSE for out of sample forecasts. # 3 different cases

M <- 300

filenames <- list.files(paste("IPEAtests2", M, "/", sep = ""))

linear <- rep(0, 140)
bols <- rep(0, 140)
bspline <- rep(0, 140)
tsboost <- rep(0, 140)
bspline2 <- rep(0, 140)
tsboost2 <- rep(0, 140)
tree  <- rep(0, 140)
all_cases <- list()

for (i in 1:12) {
  all_cases[[i]] <- data.frame(linear, bols, bspline, tsboost, bspline2, tsboost2, tree)
}

for (h in 1:12) {
  for (i in 1:length(filenames)) {
    data_test <- readMat(paste("IPEAtests2", M, "/", filenames[i], sep = ""))
    
    
    start <- data_test$test.data.start[1]
    end <- tail(data_test$test.data.end, 1)
    
    
    xt <- data_test$ar.pred[h, ] # linear
    xt2 <- data_test$tsboost.pred[h, ] # tsboost
    xt3 <- data_test$bspline.pred[h, ] # bspline
    xt4 <- data_test$bols.pred[h, ] # bols
    
    xt5 <- data_test$bspline.noextra.pred[h,]
    xt6 <- data_test$tsboost.noextra.pred[h,]
    xt7 <- data_test$tree2.pred[h,]
    
    var <- data_test$const.pred[h, ] # for var
    yt <- data_test$true.Ytph[h, ] # true value
    
    tl <- start:(end - h)
    
    
    ind <- sqrt(mean((xt4[tl] - yt[tl])^2) )/
      sqrt(mean((yt[tl] - xt[tl])^2))
    all_cases[[h]][i, "bols"] <- ifelse(ind < 0, 0, ind)
    
    ind <- sqrt(mean((xt[tl] - yt[tl])^2)) /
      sqrt(mean((yt[tl] - xt[tl])^2))
    all_cases[[h]][i, "linear"] <- ifelse(ind < 0, 0, ind)
    
    ind <- sqrt(mean((xt3[tl] - yt[tl])^2) )/
      sqrt(mean((yt[tl] - xt[tl])^2))
    all_cases[[h]][i, "bspline"] <- ifelse(ind < 0, 0, ind)
    
    ind <- sqrt(mean((xt2[tl] - yt[tl])^2)) /
      sqrt( mean((yt[tl] - xt[tl])^2))
    all_cases[[h]][i, "tsboost"] <- ifelse(ind < 0, 0, ind)
    
    ind <- sqrt(mean((xt5[tl] - yt[tl])^2)) /
      sqrt( mean((yt[tl] - xt[tl])^2))
    all_cases[[h]][i, "bspline2"] <- ifelse(ind < 0, 0, ind)
    
    ind <- sqrt(mean((xt6[tl] - yt[tl])^2) )/
      sqrt(mean((yt[tl] - xt[tl])^2))
    
    all_cases[[h]][i, "tsboost2"] <- ifelse(ind < 0, 0, ind)
    
    ind <- sqrt(mean((xt7[tl] - yt[tl])^2) )/
      sqrt(  mean((yt[tl] - xt[tl])^2))
    
    all_cases[[h]][i, "tree"] <- ifelse(ind < 0, 0, ind)
    
    
  }
}

all <- data.frame(
  "linear" = rep(0, 12),
  "bols" = rep(0, 12),
  "bspline" = rep(0, 12),
  "tsboost" = rep(0, 12),
  "bspline2" = rep(0, 12),
  "tsboost2" = rep(0, 12),
  "tree" = rep(0, 12)
  
)

all1 <- all
spl <- all
lin <- all

for (i in 1:12) {
  base = all_cases[[i]]
  
  all[i, ] <- sapply(base, mean)
  
}


N1 = c()
N2 = c()
for (i in 1:12) {
  ## Any of the splines method is best of linear cases
  base = all_cases[[i]]
  cut = base['linear'] > base['tsboost'] | base['linear'] > base['bspline']
  
  spl[i,'linear'] = mean(base[cut,'linear'])
  spl[i,'bols'] = mean(base[cut,'bols'])
  spl[i,'bspline'] = mean(base[cut,'bspline'])
  spl[i,'tsboost'] = mean(base[cut,'tsboost'])
  spl[i,'bspline2'] = mean(base[cut,'bspline2'])
  spl[i,'tsboost2'] = mean(base[cut,'tsboost2'])
  spl[i,'tree'] = mean(base[cut,'tree'])
  
  N1 = c(N1, sum(cut))
  
  ## Linear methods are better 
  
  cut = !(base['linear'] > base['tsboost'] | base['linear'] > base['bspline'])
  
  lin[i,'linear'] = mean(base[cut,'linear'])
  lin[i,'bols'] = mean(base[cut,'bols'])
  lin[i,'bspline'] = mean(base[cut,'bspline'])
  lin[i,'tsboost'] = mean(base[cut,'tsboost'])
  lin[i,'bspline2'] = mean(base[cut,'bspline2'])
  lin[i,'tsboost2'] = mean(base[cut,'tsboost2'])
  lin[i,'tree'] = mean(base[cut,'tree'])
  
  N2 = c(N2, sum(cut))
  
  
  
}

sapply(lin, mean)
sapply(lin1, mean)

clipr::write_clip(t(all))

clipr::write_clip(t(spl))
clipr::write_clip(t(N1))

clipr::write_clip(t(lin))
clipr::write_clip(t(N2))



### DM / GW tests

M <- 300

filenames <- list.files(paste("tests/IPEAtests2", M, "/", sep = ""))

linear <- rep(0, 140)
bols <- rep(0, 140)
bspline <- rep(0, 140)
tsboost <- rep(0, 140)
bspline2 <- rep(0, 140)
tsboost2 <- rep(0, 140)
tree  <- rep(0, 140)
all_05 <- list()

for (i in 1:12) {
  all_05[[i]] <- data.frame(linear, bols, bspline, tsboost, bspline2, tsboost2, tree)
}

all_10 = all_05

for (h in 1:12) {
  for (i in c(1:32, 34:44, 47:65, 67:90, 92:107, 109:115, 117:125, 127:140)) {
    data_test <- readMat(paste("tests/IPEAtests2", M, "/", filenames[i], sep = ""))
    
    
    start <- data_test$test.data.start[1]
    end <- tail(data_test$test.data.end, 1)
    
    yt <- data_test$true.Ytph[h, ] # true value
    
    xt <- data_test$ar.pred[h, ] # linear
    
    xt2 <- data_test$tsboost.pred[h, ] # tsboost
    xt3 <- data_test$bspline.pred[h, ] # bspline
    xt4 <- data_test$bols.pred[h, ] # bols
    
    xt5 <- data_test$bspline.noextra.pred[h,]
    xt6 <- data_test$tsboost.noextra.pred[h,]
    xt7 <- data_test$tree2.pred[h,]
    
    var <- data_test$const.pred[h, ] # for var
    
    tl <- start:(end - h)
    
    
    ind <- gw.test(xt4[tl],xt[tl],yt[tl], tau = h, T=288, alternative = 'less', method = 'HAC')$p.value
    all_10[[h]][i, "bols"] <- ifelse(ind <= 0.1, 1, 0)
    all_05[[h]][i, "bols"] <- ifelse(ind <= 0.05, 1, 0)
    
    ind <- gw.test(xt3[tl],xt[tl],yt[tl], tau = h, T=288, alternative = 'less', method = 'HAC')$p.value
    all_10[[h]][i, "bspline"] <- ifelse(ind <= 0.1, 1, 0)
    all_05[[h]][i, "bspline"] <- ifelse(ind <= 0.05, 1, 0)
    
    ind <- gw.test(xt2[tl],xt[tl],yt[tl], tau = h, T=288, alternative = 'less', method = 'HAC')$p.value
    all_10[[h]][i, "tsboost"] <- ifelse(ind <= 0.1, 1, 0)
    all_05[[h]][i, "tsboost"] <- ifelse(ind <= 0.05, 1, 0)
    
    ind <- gw.test(xt5[tl],xt[tl],yt[tl], tau = h, T=288, alternative = 'less', method = 'HAC')$p.value
    all_10[[h]][i, "bspline2"] <- ifelse(ind <= 0.1, 1, 0)
    all_05[[h]][i, "bspline2"] <- ifelse(ind <= 0.05, 1, 0)
    
    ind <- gw.test(xt6[tl],xt[tl],yt[tl], tau = h, T=288, alternative = 'less', method = 'HAC')$p.value
    all_10[[h]][i, "tsboost2"] <- ifelse(ind <= 0.1, 1, 0)
    all_05[[h]][i, "tsboost2"] <- ifelse(ind <= 0.05, 1, 0)
    
    ind <- gw.test(xt7[tl],xt[tl],yt[tl], tau = h, T=288, alternative = 'greater', method = 'HAC')$p.value
    all_10[[h]][i, "tree"] <- ifelse(ind <= 0.1, 1, 0)
    all_05[[h]][i, "tree"] <- ifelse(ind <= 0.05, 1, 0)
    
    
  }
}

all10 <- data.frame(
  "linear" = rep(0, 12),
  "bols" = rep(0, 12),
  "bspline" = rep(0, 12),
  "tsboost" = rep(0, 12),
  "bspline2" = rep(0, 12),
  "tsboost2" = rep(0, 12),
  "tree" = rep(0, 12)
  
)

all05 <- all10

for (i in 1:12) {
  base = all_10[[i]]
  
  all10[i, ] <- sapply(base, sum)
  
  base = all_05[[i]]
  
  all05[i, ] <- sapply(base, sum)
  
}



clipr::write_clip(t(all10))
clipr::write_clip(t(all05))


### rRMSE as advisor requested

M <- 300

filenames <- list.files(paste("IPEAtests2", M, "/", sep = ""))


linear <- rep(0, 140)
bols <- rep(0, 140)
bspline <- rep(0, 140)
tsboost <- rep(0, 140)
bspline2 <- rep(0, 140)
tsboost2 <- rep(0, 140)
tree  <- rep(0, 140)
all_cases <- list()

for (i in 1:12) {
  all_cases[[i]] <- data.frame(linear, bols, bspline, tsboost, bspline2, tsboost2, tree)
}
for (h in 12) {
  for (i in 1:length(filenames)) {
    data_test <- readMat(paste("IPEAtests2", M, "/", filenames[i], sep = ""))
    
    
    start <- data_test$test.data.start[1]
    end <- tail(data_test$test.data.end, 1)
    
    
    xt <- data_test$ar.pred[h, ] # linear
    xt2 <- data_test$tsboost.pred[h, ] # tsboost
    xt3 <- data_test$bspline.pred[h, ] # bspline
    xt4 <- data_test$bols.pred[h, ] # bols
    
    xt5 <- data_test$bspline.noextra.pred[h,]
    xt6 <- data_test$tsboost.noextra.pred[h,]
    xt7 <- data_test$tree2.pred[h,]
    
    var <- data_test$const.pred[h, ] # for var
    yt <- data_test$true.Ytph[h, ] # true value
    
    tl <- start:(end - h)
    
    
    ind <- sqrt(mean((xt4[tl] - yt[tl])^2) )/
      sqrt(mean((yt[tl] - xt[tl])^2))
    all_cases[[h]][i, "bols"] <- ifelse(ind < 0, 0, ind)
    
    ind <- sqrt(mean((xt[tl] - yt[tl])^2)) /
      sqrt(mean((yt[tl] - xt[tl])^2))
    all_cases[[h]][i, "linear"] <- ifelse(ind < 0, 0, ind)
    
    ind <- sqrt(mean((xt3[tl] - yt[tl])^2) )/
      sqrt(mean((yt[tl] - xt[tl])^2))
    all_cases[[h]][i, "bspline"] <- ifelse(ind < 0, 0, ind)
    
    ind <- sqrt(mean((xt2[tl] - yt[tl])^2)) /
      sqrt( mean((yt[tl] - xt[tl])^2))
    all_cases[[h]][i, "tsboost"] <- ifelse(ind < 0, 0, ind)
    
    ind <- sqrt(mean((xt5[tl] - yt[tl])^2)) /
      sqrt( mean((yt[tl] - xt[tl])^2))
    all_cases[[h]][i, "bspline2"] <- ifelse(ind < 0, 0, ind)
    
    ind <- sqrt(mean((xt6[tl] - yt[tl])^2) )/
      sqrt(mean((yt[tl] - xt[tl])^2))
    
    all_cases[[h]][i, "tsboost2"] <- ifelse(ind < 0, 0, ind)
    
    ind <- sqrt(mean((xt7[tl] - yt[tl])^2) )/
      sqrt(  mean((yt[tl] - xt[tl])^2))
    
    all_cases[[h]][i, "tree"] <- ifelse(ind < 0, 0, ind)
    
  }
}

names <- pull(read_excel("dados.xlsx", sheet = "Planilha1", 
                    col_names = FALSE))
order <-  as.numeric(gsub("([0-9]+).*$", "\\1", (substr(filenames, 40, 42))))

Series <- names[order]

clipr::write_clip(cbind(Series,round(all_cases[[12]],3)))


#### All the plots
create_df = function (data, index, h){
  result = data[[h]][index,]
  
  if(h>1){
    for( i in 2:h){
      result = rbind(result, data[[i]][index,])
    }
  }
  
  h = 1:h
  return(cbind(result, h))  
}

filenames
dataset$metadados$name


plot_data = function (index_out, data, h =12) {
  df = create_df(data = data, index = index_out, h =12)
  col = paletteer_d("RColorBrewer::Paired", 7)
  
  return(plot <- ggplot(data = df, mapping = aes(x = h)) +
    geom_line(aes(y = linear, color = col[1]), ) +
    geom_point(aes(y = linear), shape = 21, colour = "black") +
      
    geom_line(aes(y = bols, color = col[2])) +
    geom_point(aes(y = bols), shape = 21, colour = "black") +
      
    geom_line(aes(y = bspline, color = col[3])) +
    geom_point(aes(y = bspline), shape = 21, colour = "black") +
      
    geom_line(aes(y = bspline2, color = col[4])) +
    geom_point(aes(y = bspline2), shape = 21, colour = "black") +
      
    geom_line(aes(y = tsboost, color = col[5])) +
    geom_point(aes(y = tsboost), shape = 21, colour = "black") +
      
    geom_line(aes(y = tsboost2, color = col[6])) +
    geom_point(aes(y = tsboost2), shape = 21, colour = "black") +
      
    geom_line(aes(y = tree, color = col[8])) +
    geom_point(aes(y = tree), shape = 21, colour = "black") +
      
      
    ylab("R²") +
    ggtitle(Series[index_out]) +
    ylim(0, 1) + # scale_y_continuous(labels=as.character(seq(0,1, by=0.1)),breaks=seq(0,1, by=0.1)) +
    scale_x_continuous(labels = as.character(df$h), breaks = df$h) +
    theme_calc() + 
    scale_colour_manual(values = col,name = "Models", labels = c("Linear", "BOLS",
                                                     'BSpline', 'BSpline*', "TSBoost", 'TSBoost*','Tree')))

}


plot_data = function (index_out, data, h =12) {
  df = create_df(data = data, index = index_out, h =12)
  col = paletteer_d("RColorBrewer::Paired", 3)
  
  return(plot <- ggplot(data = df, mapping = aes(x = h)) +
           geom_line(aes(y = linear, color = col[1]), ) +
           geom_point(aes(y = linear), shape = 21, colour = "black") +
           
           geom_line(aes(y = bols, color = col[2])) +
           geom_point(aes(y = bols), shape = 21, colour = "black") +
           
           
         #  geom_line(aes(y = bspline, color = col[3])) +
        #   geom_point(aes(y = bspline), shape = 21, colour = "black") +
           
        #   geom_line(aes(y = tsboost, color = col[4])) +
        #   geom_point(aes(y = tsboost), shape = 21, colour = "black") +
           
           geom_line(aes(y = tree, color = col[3])) +
           geom_point(aes(y = tree), shape = 21, colour = "black") +
           
           ggtitle(Series[index_out]) +
           ylim(0, 1) + # scale_y_continuous(labels=as.character(seq(0,1, by=0.1)),breaks=seq(0,1, by=0.1)) +
           scale_x_continuous(labels = as.character(df$h), breaks = df$h) +
           theme_excel_new() + xlab("h") + ylab("R²") + theme(axis.title=element_text(size=10)) +
           scale_colour_manual(values = col,name = "Models", labels = c("Linear", "BOLS",
                                                                        #'BSpline', "TSBoost",
                                                                        "Tree")))
  
}
cut = which_splines_superior(all_cases)

myplots <- vector('list', length(cut))
for(i in cut){
  temp_plot <- plot_data(index_out = i, all_cases)
  myplots[[i]] <- temp_plot
  ggsave(temp_plot, path="plots4", filename = paste0('plot_', i, '.tiff'), width = 6, height = 4, units = "in")
}

which_splines_superior <- function(data_out){
  vector = c()
  for(i in 1:140){
    temp <- create_df(data = data_out, i, 12)
    condition = sum( ((temp[,'linear'] < temp[,'tree']) & (temp[,'bols'] < temp[,'tree']))
                     & ((temp[,'linear'] < temp[,'tree']) & (temp[,'bols'] < temp[,'tree'])) )
    if(condition>=9){
      vector = c(vector, i)
    }
  }
  return(vector)
}




cut = which_splines_superior(all_cases)
myplots <- vector('list', length(cut))
k=1
for(i in cut){
  temp_plot <- plot_data(index_out = i, all_cases)
  myplots[[k]] <- temp_plot
  k=k+1
  ggsave(temp_plot, path="plots\\plots2", filename = paste0('plot_', i, '.tiff'), width = 6, height = 4, units = "in")
}
ggarrange(myplots[[3]],myplots[[6]], myplots[[7]],myplots[[8]], 
          ncol = 2, nrow=2, common.legend = TRUE, legend="bottom")

ggarrange(myplots[[1]],myplots[[2]], myplots[[3]],myplots[[4]], myplots[[5]],myplots[[6]], 
          ncol = 2, nrow=3, common.legend = TRUE, legend="bottom")


ggarrange(myplots[[7]],myplots[[8]], myplots[[11]],myplots[[9]], myplots[[13]],myplots[[14]], 
          ncol = 2, nrow=3, common.legend = TRUE, legend="bottom")


## analysis with AIC method

### calculation of all cases for AIC methods
M <- 300

filenames <- list.files(paste("tests/IPEAtests2", M, "/", sep = ""))

linear <- rep(0, 140)
bols <- rep(0, 140)
bspline <- rep(0, 140)
tsboost <- rep(0, 140)
bspline2 <- rep(0, 140)
tsboost2 <- rep(0, 140)
all_cases_aic <- list()

for (i in 1:12) {
  all_cases_aic[[i]] <- data.frame(linear, bols, bspline, tsboost, bspline2, tsboost2)
}

for (h in 1:12) {
  for (i in 1:length(filenames)) {
    data_test <- readMat(paste("tests/IPEAtests2", M, "/", filenames[i], sep = ""))
    
    
    start <- data_test$test.data.start[1]
    end <- tail(data_test$test.data.end, 1)
    
    
    xt <- data_test$ar.pred[h, ] # linear
    xt2 <- data_test$tsboost.pred.aic[h, ] # tsboost
    xt3 <- data_test$bspline.pred.aic[h, ] # bspline
    xt4 <- data_test$bols.pred.aic[h, ] # bols
    
    xt5 <- data_test$bspline.noextra.pred.aic[h,]
    xt6 <- data_test$tsboost.noextra.pred.aic[h,]
    
    var <- data_test$const.pred[h, ] # for var
    yt <- data_test$true.Ytph[h, ] # true value
    
    tl <- start:(end - h)
    
    
    ind <- 1 - mean((xt4[tl] - yt[tl])^2) /
      mean((yt[tl] - var[tl])^2)
    all_cases_aic[[h]][i, "bols"] <- ifelse(ind < 0, 0, ind)
    
    ind <- 1 - mean((xt[tl] - yt[tl])^2) /
      mean((yt[tl] - var[tl])^2)
    all_cases_aic[[h]][i, "linear"] <- ifelse(ind < 0, 0, ind)
    
    ind <- 1 - mean((xt3[tl] - yt[tl])^2) /
      mean((yt[tl] - var[tl])^2)
    all_cases_aic[[h]][i, "bspline"] <- ifelse(ind < 0, 0, ind)
    
    ind <- 1 - mean((xt2[tl] - yt[tl])^2) /
      mean((yt[tl] - var[tl])^2)
    all_cases_aic[[h]][i, "tsboost"] <- ifelse(ind < 0, 0, ind)
    
    ind <- 1 - mean((xt5[tl] - yt[tl])^2) /
      mean((yt[tl] - var[tl])^2)
    all_cases_aic[[h]][i, "bspline2"] <- ifelse(ind < 0, 0, ind)
    
    ind <- 1 - mean((xt6[tl] - yt[tl])^2) /
      mean((yt[tl] - var[tl])^2)
    all_cases_aic[[h]][i, "tsboost2"] <- ifelse(ind < 0, 0, ind)
    
  }
}

all <- data.frame(
  "linear" = rep(0, 12),
  "bols" = rep(0, 12),
  "bspline" = rep(0, 12),
  "tsboost" = rep(0, 12),
  "bspline2" = rep(0, 12),
  "tsboost2" = rep(0, 12)
  
)

all1 <- all

spl <- all
spl1 <- all

lin <- all
lin1 <- all

N0 = rep(140,12)
N1 = c()

for (i in 1:12) {
  base = all_cases_aic[[i]]
  
  all[i, ] <- sapply(base, mean)
  all1[i, ] <- sapply(base, mean_01)
  N1 = c(N1, sum(base['tsboost']>.1)) 
  
}


N2 = c()
N3 = c()
N4 = c()
N5 = c()
for (i in 1:12) {
  ## Any of the splines method is best of linear cases
  base = all_cases_aic[[i]]
  cut = base['linear'] < base['tsboost'] | base['linear'] < base['bspline']
  
  spl[i,'linear'] = mean(base[cut,'linear'])
  spl[i,'bols'] = mean(base[cut,'bols'])
  spl[i,'bspline'] = mean(base[cut,'bspline'])
  spl[i,'tsboost'] = mean(base[cut,'tsboost'])
  spl[i,'bspline2'] = mean(base[cut,'bspline2'])
  spl[i,'tsboost2'] = mean(base[cut,'tsboost2'])
  
  N2 = c(N2, sum(cut))
  
  cut = (base['linear'] < base['tsboost'] | base['linear'] < base['bspline']) & base['linear'] >=0.1 & base['bols'] >=0.1 & base['bspline'] >=0.1 & base['tsboost'] >=0.1
  
  spl1[i,'linear'] = mean(base[cut,'linear'])
  spl1[i,'bols'] = mean(base[cut,'bols'])
  spl1[i,'bspline'] = mean(base[cut,'bspline'])
  spl1[i,'tsboost'] = mean(base[cut,'tsboost'])
  spl1[i,'bspline2'] = mean(base[cut,'bspline2'])
  spl1[i,'tsboost2'] = mean(base[cut,'tsboost2'])
  
  N3 = c(N3, sum(cut))
  
  ## Linear methods are better 
  
  cut = !(base['linear'] < base['tsboost'] | base['linear'] < base['bspline'])
  
  lin[i,'linear'] = mean(base[cut,'linear'])
  lin[i,'bols'] = mean(base[cut,'bols'])
  lin[i,'bspline'] = mean(base[cut,'bspline'])
  lin[i,'tsboost'] = mean(base[cut,'tsboost'])
  lin[i,'bspline2'] = mean(base[cut,'bspline2'])
  lin[i,'tsboost2'] = mean(base[cut,'tsboost2'])
  
  
  N4 = c(N4, sum(cut))
  
  cut = !(base['linear'] < base['tsboost'] | base['linear'] < base['bspline']) & base['linear'] >=0.1 & base['bols'] >=0.1 & base['bspline'] >=0.1 & base['tsboost'] >=0.1
  
  lin1[i,'linear'] = mean(base[cut,'linear'])
  lin1[i,'bols'] = mean(base[cut,'bols'])
  lin1[i,'bspline'] = mean(base[cut,'bspline'])
  lin1[i,'tsboost'] = mean(base[cut,'tsboost'])
  lin1[i,'bspline2'] = mean(base[cut,'bspline2'])
  lin1[i,'tsboost2'] = mean(base[cut,'tsboost2'])
  
  N5 = c(N5, sum(cut))
  
  
}

sapply(lin, mean)
sapply(lin1, mean)

clipr::write_clip(t(all))
clipr::write_clip(t(N0))

clipr::write_clip(t(all1))
clipr::write_clip(t(N1))

clipr::write_clip(t(spl))
clipr::write_clip(t(N2))

clipr::write_clip(t(spl1))
clipr::write_clip(t(N3))

clipr::write_clip(t(lin))
clipr::write_clip(t(N4))

clipr::write_clip(t(lin1))
clipr::write_clip(t(N5))


plot_data = function (index_out, data, h =12) {
  df = create_df(data = data, index = index_out, h =12)
  col = paletteer_d("RColorBrewer::Paired", 4)
  
  return(plot <- ggplot(data = df, mapping = aes(x = h)) +
           geom_line(aes(y = linear, color = col[1]), ) +
           geom_point(aes(y = linear), shape = 21, colour = "black") +
           
           geom_line(aes(y = bols, color = col[2])) +
           geom_point(aes(y = bols), shape = 21, colour = "black") +
           
           
            geom_line(aes(y = bspline, color = col[3])) +
              geom_point(aes(y = bspline), shape = 21, colour = "black") +
           
             geom_line(aes(y = tsboost, color = col[4])) +
              geom_point(aes(y = tsboost), shape = 21, colour = "black") +
           
          # geom_line(aes(y = tree, color = col[3])) +
          # geom_point(aes(y = tree), shape = 21, colour = "black") +
           
           ggtitle(Series[index_out]) +
           ylim(0, 1) + # scale_y_continuous(labels=as.character(seq(0,1, by=0.1)),breaks=seq(0,1, by=0.1)) +
           scale_x_continuous(labels = as.character(df$h), breaks = df$h) +
           theme_excel_new() + xlab("h") + ylab("R²") + theme(axis.title=element_text(size=10)) +
           scale_colour_manual(values = col,name = "Models", labels = c("Linear", "BOLS",
                                                                        'BSpline', "TSBoost")))
                                                                        #"Tree"
  
}

which_splines_superior <- function(data_out){
  vector = c()
  for(i in 1:140){
    temp <- create_df(data = data_out, i, 12)
    condition = sum( ((temp[,'linear'] < temp[,'bspline']) & (temp[,'bols'] < temp[,'bspline']))
                     & ((temp[,'linear'] < temp[,'tsboost']) & (temp[,'bols'] < temp[,'tsboost'])) )
    if(condition>=9){
      vector = c(vector, i)
    }
  }
  return(vector)
}

cut = which_splines_superior(all_cases_aic)
myplots <- vector('list', length(cut))
k=1
for(i in cut){
  temp_plot <- plot_data(index_out = i, all_cases_aic)
  myplots[[k]] <- temp_plot
  k=k+1
  ggsave(temp_plot, path="plots\\plots3", filename = paste0('plot_', i, '.tiff'), width = 6, height = 4, units = "in")
}
