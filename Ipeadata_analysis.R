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
source("helper_functions_analysis.R")


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
bols_aic<- c()

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
    
    
    #ind <- 1 - (mse(xt[tl], yt[tl]) / mse(var[tl], yt[tl])) # constant pred
    
    #r_benchmark <- c(r_benchmark, ifelse(ind < 0, 0, ind))
    #ind <- 1 - (mse(xt2[tl], yt[tl]) / mse(var[tl], yt[tl])) # const pred
    
    #r_nonlinear <- c(r_nonlinear, ifelse(ind < 0, 0, ind))
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