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
# data <- readMat("FRED.mat")
source("helper_functions_analysis.R")

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
x <- c()

linear <- c()
tsboost <- c()
bspline <- c()
bols <- c()

linear1 <- c()
tsboost1 <- c()
bspline1 <- c()
bols1 <- c()
N1 <- c()

linearb <- c()
tsboostb <- c()
bsplineb <- c()
bolsb <- c()
N2 <- c()

linear1b <- c()
tsboost1b <- c()
bspline1b <- c()
bols1b <- c()
N3 <- c()

linearl <- c()
tsboostl <- c()
bsplinel <- c()
bolsl <- c()
N4 <- c()

linear1l <- c()
tsboost1l <- c()
bspline1l <- c()
bols1l <- c()
N5 <- c()

for (h in 1:12) {
  linear_in <- c()
  tsboost_in <- c()
  bspline_in <- c()
  bols_in <- c()

  for (name in filenames) {
    data_test <- readMat(paste("IPEAtests2", M, "/", name, sep = ""))


    start <- data_test$test.data.start[1]
    end <- tail(data_test$test.data.end, 1)


    xt <- data_test$bic.pred[h, ] # linear
    xt2 <- data_test$tsboost.pred[h, ] # tsboost
    xt3 <- data_test$bspline.pred[h, ] # bspline
    xt4 <- data_test$bols.pred[h, ] # bols


    var <- data_test$const.pred[h, ] # for var
    yt <- data_test$true.Ytph[h, ] # true value

    tl <- start:(end - h)


    ind <- 1 - (mse(xt[tl], yt[tl]) / mse(var[tl], yt[tl])) # constant pred
    linear_in <- c(linear_in, ifelse(ind < 0, 0, ind))

    ind <- 1 - (mse(xt2[tl], yt[tl]) / mse(var[tl], yt[tl])) # const pred
    tsboost_in <- c(tsboost_in, ifelse(ind < 0, 0, ind))

    ind <- 1 - (mse(xt3[tl], yt[tl]) / mse(var[tl], yt[tl])) # constant pred
    bspline_in <- c(bspline_in, ifelse(ind < 0, 0, ind))

    ind <- 1 - (mse(xt4[tl], yt[tl]) / mse(var[tl], yt[tl])) # const pred
    bols_in <- c(bols_in, ifelse(ind < 0, 0, ind))
  }

  ### All cases
  linear <- c(linear, mean(linear_in))
  tsboost <- c(tsboost, mean(tsboost_in))
  bspline <- c(bspline, mean(bspline_in))
  bols <- c(bols, mean(bols_in))

  N <- rep(140, 12)

  linear1 <- c(linear1, mean(linear_in[linear_in >= .1]))
  tsboost1 <- c(tsboost1, mean(tsboost_in[tsboost_in >= .1]))
  bspline1 <- c(bspline1, mean(bspline_in[bspline_in >= .1]))
  bols1 <- c(bols1, mean(bols_in[bols_in >= .1]))

  N1 <- c(N1, sum(tsboost_in >= .1))

  ### boosting is more accurate at least once

  linearb <- c(linearb, mean(linear_in[(linear_in <= tsboost_in | linear_in <= bspline_in) & linear_in >= 0]))
  tsboostb <- c(tsboostb, mean(tsboost_in[(tsboost_in >= linear_in | tsboost_in >= bols_in) & tsboost_in >= 0]))
  bsplineb <- c(bsplineb, mean(bspline_in[(bspline_in >= linear_in | bspline_in >= bols_in) & bspline_in >= 0]))
  bolsb <- c(bolsb, mean(bols_in[(bols_in <= tsboost_in | bols_in <= bspline_in) & bols_in >= 0]))

  N2 <- c(N2, sum((linear_in < tsboost_in | linear_in < bspline_in) & tsboost_in >= 0))

  linear1b <- c(linear1b, mean(linear_in[(linear_in < tsboost_in | linear_in < bspline_in) & linear_in >= .1]))
  tsboost1b <- c(tsboost1b, mean(tsboost_in[(tsboost_in >= linear_in | tsboost_in >= bols_in) & tsboost_in >= .1]))
  bspline1b <- c(bspline1b, mean(bspline_in[(bspline_in >= linear_in | bspline_in >= bols_in) & bspline_in >= .1]))
  bols1b <- c(bols1b, mean(bols_in[(bols_in < tsboost_in | bols_in < bspline_in) & bols_in >= .1]))

  N3 <- c(N3, sum((linear_in <= tsboost_in | linear_in <= bspline_in) & tsboost_in >= .1))

  ### linear is more accurate at least once

  linearl <- c(linearl, mean(linear_in[!(linear_in <= tsboost_in | linear_in <= bspline_in) & linear_in >= 0]))
  tsboostl <- c(tsboostl, mean(tsboost_in[!(tsboost_in >= linear_in | tsboost_in >= bols_in) & tsboost_in >= 0]))
  bsplinel <- c(bsplinel, mean(bspline_in[!(bspline_in >= linear_in | bspline_in >= bols_in) & bspline_in >= 0]))
  bolsl <- c(bolsl, mean(bols_in[!(bols_in <= tsboost_in | bols_in <= bspline_in) & bols_in >= 0]))

  N4 <- c(N4, sum(!(linear_in <= tsboost_in | linear_in <= bspline_in) & tsboost_in >= 0))

  linear1l <- c(linear1l, mean(linear_in[!(linear_in < tsboost_in | linear_in < bspline_in) & linear_in >= .1]))
  tsboost1l <- c(tsboost1l, mean(tsboost_in[!(tsboost_in >= linear_in | tsboost_in >= bols_in) & tsboost_in >= .1]))
  bspline1l <- c(bspline1l, mean(bspline_in[!(bspline_in >= linear_in | bspline_in >= bols_in) & bspline_in >= .1]))
  bols1l <- c(bols1l, mean(bols_in[!(bols_in < tsboost_in | bols_in < bspline_in) & bols_in >= .1]))

  N5 <- c(N5, sum(!(linear_in <= tsboost_in | linear_in <= bspline_in) & tsboost_in >= .1))
}


### TESTING FRED DATA

filenames <- list.files(paste("FREDTests", "/", sep = ""))
x <- c()

linear <- c()
tsboost <- c()
bspline <- c()

linear1 <- c()
tsboost1 <- c()
bspline1 <- c()
N1 <- c()

linearb <- c()
tsboostb <- c()
bsplineb <- c()
N2 <- c()

linear1b <- c()
tsboost1b <- c()
bspline1b <- c()
N3 <- c()

linearl <- c()
tsboostl <- c()
bsplinel <- c()
N4 <- c()

linear1l <- c()
tsboost1l <- c()
bspline1l <- c()
N5 <- c()

for (h in 1:12) {
  linear_in <- c()
  tsboost_in <- c()
  bspline_in <- c()

  for (name in filenames) {
    data_test <- readMat(paste("FREDTests", "/", name, sep = ""))


    start <- data_test$test.data.start[1]
    end <- tail(data_test$test.data.end, 1)


    xt <- data_test$bic.pred[h, ] # linear
    xt2 <- data_test$tsboost.pred[h, ] # tsboost
    xt3 <- data_test$bspline.pred[h, ] # bspline


    var <- data_test$const.pred[h, ] # for var
    yt <- data_test$true.Ytph[h, ] # true value

    tl <- start:(end - h)


    ind <- 1 - (mse(xt[tl], yt[tl]) / mse(var[tl], yt[tl])) # constant pred
    linear_in <- c(linear_in, ifelse(ind < 0, 0, ind))

    ind <- 1 - (mse(xt2[tl], yt[tl]) / mse(var[tl], yt[tl])) # const pred
    tsboost_in <- c(tsboost_in, ifelse(ind < 0, 0, ind))

    ind <- 1 - (mse(xt3[tl], yt[tl]) / mse(var[tl], yt[tl])) # constant pred
    bspline_in <- c(bspline_in, ifelse(ind < 0, 0, ind))
  }

  ### All cases
  linear <- c(linear, mean(linear_in))
  tsboost <- c(tsboost, mean(tsboost_in))
  bspline <- c(bspline, mean(bspline_in))

  N <- rep(140, 12)

  linear1 <- c(linear1, mean(linear_in[linear_in >= .1]))
  tsboost1 <- c(tsboost1, mean(tsboost_in[tsboost_in >= .1]))
  bspline1 <- c(bspline1, mean(bspline_in[bspline_in >= .1]))

  N1 <- c(N1, sum(tsboost_in >= .1))

  ### boosting is more accurate at least once

  linearb <- c(linearb, mean(linear_in[(linear_in <= tsboost_in & linear_in <= bspline_in) & linear_in >= 0]))
  tsboostb <- c(tsboostb, mean(tsboost_in[(tsboost_in >= linear_in & tsboost_in >= bols_in) & tsboost_in >= 0]))
  bsplineb <- c(bsplineb, mean(bspline_in[(bspline_in >= linear_in & bspline_in >= bols_in) & bspline_in >= 0]))

  N2 <- c(N2, sum((linear_in <= tsboost_in & linear_in <= bspline_in) & tsboost_in >= 0))

  linear1b <- c(linear1b, mean(linear_in[(linear_in < tsboost_in & linear_in < bspline_in) & linear_in >= .1]))
  tsboost1b <- c(tsboost1b, mean(tsboost_in[(tsboost_in >= linear_in) & tsboost_in >= .1]))
  bspline1b <- c(bspline1b, mean(bspline_in[(bspline_in >= linear_in) & bspline_in >= .1]))

  N3 <- c(N3, sum((linear_in <= tsboost_in & linear_in <= bspline_in) & linear_in >= .1))

  ### linear is more accurate at least once

  linearl <- c(linearl, mean(linear_in[!(linear_in <= tsboost_in & linear_in <= bspline_in) & linear_in >= 0]))
  tsboostl <- c(tsboostl, mean(tsboost_in[!(tsboost_in >= linear_in & tsboost_in >= bols_in) & tsboost_in >= 0]))
  bsplinel <- c(bsplinel, mean(bspline_in[!(bspline_in >= linear_in & bspline_in >= bols_in) & bspline_in >= 0]))

  N4 <- c(N4, sum(!(linear_in <= tsboost_in & linear_in <= bspline_in) & tsboost_in >= 0))


  # linearl <- c(linearl, mean(linear_in[(linear_in>tsboost_in & linear_in>bspline_in)& linear_in >0]))
  # tsboostl <- c(tsboostl, mean(tsboost_in[(tsboost_in<linear_in)& tsboost_in >0]))
  # bsplinel <- c(bsplinel, mean(bspline_in[(bspline_in<linear_in)& bspline_in >0]))
  # bolsl <- c(bolsl, mean(bols_in[(bols_in>tsboost_in & bols_in>bspline_in)& bols_in >0]))
  #
  # N4 <- c(N4, sum((linear_in>tsboost_in & linear_in>bspline_in) & linear_in > 0))

  linear1l <- c(linear1l, mean(linear_in[(linear_in > tsboost_in & linear_in > bspline_in) & linear_in >= .1]))
  tsboost1l <- c(tsboost1l, mean(tsboost_in[(tsboost_in < linear_in) & tsboost_in >= .1]))
  bspline1l <- c(bspline1l, mean(bspline_in[(bspline_in < linear_in) & bspline_in >= .1]))

  N5 <- c(N5, sum((linear_in > tsboost_in & linear_in > bspline_in) & linear_in >= .1))
}

### TESTING 2

M <- 300

filenames <- list.files(paste("IPEAtests2", M, "/", sep = ""))

linear <- rep(0, 140)
bols <- rep(0, 140)
bspline <- rep(0, 140)
tsboost <- rep(0, 140)

all_cases <- list()

for (i in 1:12) {
  all_cases[[i]] <- data.frame(linear, bols, bspline, tsboost)
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
  }
}

all <- data.frame(
  "linear" = rep(0, 12),
  "bols" = rep(0, 12),
  "bspline" = rep(0, 12),
  "tsboost" = rep(0, 12)
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
    
  N2 = c(N2, sum(cut))
  
  cut = (base['linear'] < base['tsboost'] | base['linear'] < base['bspline']) & base['linear'] >=0.1 & base['bols'] >=0.1 & base['bspline'] >=0.1 & base['tsboost'] >=0.1
  
  spl1[i,'linear'] = mean(base[cut,'linear'])
  spl1[i,'bols'] = mean(base[cut,'bols'])
  spl1[i,'bspline'] = mean(base[cut,'bspline'])
  spl1[i,'tsboost'] = mean(base[cut,'tsboost'])
  
  N3 = c(N3, sum(cut))
  
  ## Linear methods are better 
  
  cut = !(base['linear'] < base['tsboost'] | base['linear'] < base['bspline'])
  
  lin[i,'linear'] = mean(base[cut,'linear'])
  lin[i,'bols'] = mean(base[cut,'bols'])
  lin[i,'bspline'] = mean(base[cut,'bspline'])
  lin[i,'tsboost'] = mean(base[cut,'tsboost'])
  
  N4 = c(N4, sum(cut))
  
  cut = !(base['linear'] < base['tsboost'] | base['linear'] < base['bspline']) & base['linear'] >=0.1 & base['bols'] >=0.1 & base['bspline'] >=0.1 & base['tsboost'] >=0.1
  
  lin1[i,'linear'] = mean(base[cut,'linear'])
  lin1[i,'bols'] = mean(base[cut,'bols'])
  lin1[i,'bspline'] = mean(base[cut,'bspline'])
  lin1[i,'tsboost'] = mean(base[cut,'tsboost'])
  
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

