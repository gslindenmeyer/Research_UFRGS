#############################################################
## Here we run the simulations using the simulation's functions and data
## gslindenmeyer@gmail.com - Oct 2021
#############################################################

##############################
# 0 - Load librairies
##############################
rm(list = ls())
require(mboost)
require(forecast)
require(R.matlab)
### for data visualization
require(ggplot2)
require(gridExtra)
require(ggthemes)
### fancy loading time
require(tcltk)
### adquiring our functions and objects

##############################
# 1 - Source file
##############################
dataPath <- "C:/Users/guili/Documents/GitHub/Research_UFRGS/simulation/"
source(paste0(dataPath, "functions_linear.R"))
source(paste0(dataPath, "functions_boosting.R"))
source(paste0(dataPath, "simulated_models.R"))
##############
## Settings ##
##############
save_directory <- "Mats/" # Add / to end or leave empty!
filename <- paste(dataPath, save_directory, "simulated_models.RData", sep = "")
filename <- paste(dataPath, save_directory, "simulated_models_boost.RData", sep = "")
load(filename)
###

##############################
# 2 - Start my code
##############################
modelo_a <- catch_r2(yta, name = "a", best = T)
modelo_b <- catch_r2(ytb, name = "b", best = T)
modelo_c <- catch_r2(ytc, name = "c", best = T)
modelo_d <- catch_r2(ytd, name = "d", best = T)
modelo_e <- catch_r2(yte, name = "e", best = T)
modelo_f <- catch_r2(ytf, name = "f", best = T)

models <- list("a" = modelo_a, "b" = modelo_b, "c" = modelo_c, "d" = modelo_d, "e" = modelo_e, "f" = modelo_f)

for (i in 1:length(models)) {
  df <- data.frame(
    "h" = models[[i]]$best_linear$h, "R2_best" = models[[i]]$best_linear$vetor_h,
    "R2" = models[[i]]$linear_estimate$vetor_h
  )
  models[[i]]$plot <- ggplot(data = df, mapping = aes(x = h)) +
    geom_line(aes(y = R2, color = "darkred")) +
    geom_point(aes(y = R2), shape = 21, colour = "black") +
    geom_line(aes(y = R2_best, color = "steelblue")) +
    geom_point(aes(y = R2_best), shape = 21, colour = "black") +
    ylab("R²") +
    ggtitle(paste("model_", substring(models[[i]]$best_linear$name, 1, 1), sep = "")) +
    ylim(0, 1) + # scale_y_continuous(labels=as.character(seq(0,1, by=0.1)),breaks=seq(0,1, by=0.1)) +
    scale_x_continuous(labels = as.character(models[[i]]$best_linear$h), breaks = models[[i]]$best_linear$h) +
    theme_calc() +
    scale_color_discrete(name = "Models", labels = c("Linear_Estimate", "Best_Linear"))
  #+theme(legend.position='bottom')
}

grid.arrange(models$a$plot, models$b$plot, models$c$plot, models$d$plot, models$e$plot, models$f$plot, ncol = 2)

hmax <- 12 # Def=12
hrange <- c(1:12) # Def=1:12
deg_freedom <- 4 # Def=4
cv_folds <- 5 # Def=5
boost_use_lags <- "FULL" # "FULL" for all 12 or "BIC" for same as with AR(BIC)
maxlearners <- 300 # Def 300
bols_maxlearners <- 100 # Def 300
spline_bctrl <- boost_control(mstop = maxlearners)
bols_bctrl <- boost_control(mstop = bols_maxlearners)

modelo_a_boost <- catch_r2_boosting(yta, name = "a", best = T, cv_folds = cv_folds, bctrl_out = spline_bctrl)
save(modelo_a_boost, file = paste(save_directory, "modelo_a_boost.RData", sep = ""))

modelo_b_boost <- catch_r2_boosting(ytb, name = "b", best = T, cv_folds = cv_folds, bctrl_out = spline_bctrl)
save(modelo_b_boost, file = paste(save_directory, "modelo_b_boost.RData", sep = ""))

modelo_c_boost <- catch_r2_boosting(ytc, name = "c", best = T, cv_folds = cv_folds, bctrl_out = spline_bctrl)
save(modelo_c_boost, file = paste(save_directory, "modelo_c_boost.RData", sep = ""))

modelo_d_boost <- catch_r2_boosting(ytd, name = "d", best = T, cv_folds = cv_folds, bctrl_out = spline_bctrl)
save(modelo_d_boost, file = paste(save_directory, "modelo_d_boost.RData", sep = ""))

modelo_e_boost <- catch_r2_boosting(yte, name = "e", best = T, cv_folds = cv_folds, bctrl_out = spline_bctrl)
save(modelo_e_boost, file = paste(save_directory, "modelo_e_boost.RData", sep = ""))

modelo_f_boost <- catch_r2_boosting(ytf, name = "f", best = T, cv_folds = cv_folds, bctrl_out = spline_bctrl)
save(modelo_f_boost, file = paste(save_directory, "modelo_f_boost.RData", sep = ""))


models_boost <- list(
  "a" = modelo_a_boost, "b" = modelo_b_boost,
  "c" = modelo_c_boost, "d" = modelo_d_boost,
  "e" = modelo_e_boost, "f" = modelo_f_boost
)

save(models, file = filename)

save(models_boost, file = filename)


dfa <- data.frame(
  "h" = models[[i]]$best_linear$h, "R2_best" = models[[i]]$best_linear$vetor_h,
  "R2" = models[[i]]$linear_estimate$vetor_h, "R2_best_nl" = modelo_a_boost$best_linear$vetor_h,
  "R2_nl" = modelo_a_boost$linear_estimate$vetor_h
)

superplot <- ggplot(data = dfa, mapping = aes(x = h)) +
  geom_line(aes(y = R2, color = "darkred")) +
  geom_point(aes(y = R2), shape = 21, colour = "black") +
  geom_line(aes(y = R2_best, color = "steelblue")) +
  geom_point(aes(y = R2_best), shape = 21, colour = "black") +
  geom_line(aes(y = R2_best_nl, color = "yellow")) +
  geom_line(aes(y = R2_nl, color = "green")) +
  ylab("R²") +
  ggtitle(paste("model_", substring(models[[i]]$best_linear$name, 1, 1), sep = "")) +
  ylim(0, 1) + # scale_y_continuous(labels=as.character(seq(0,1, by=0.1)),breaks=seq(0,1, by=0.1)) +
  scale_x_continuous(labels = as.character(models[[i]]$best_linear$h), breaks = models[[i]]$best_linear$h) +
  theme_calc() +
  scale_color_discrete(name = "Models", labels = c("Linear_Estimate", "Best_Linear", "Optimal", "Boosting Estimate"))
#+theme(legend.position='bottom')


## O que fazer? Eu tenho que fazer com que o best_linear capture o melhor p para cada h e salve essa informação
## Mas o que eu posso ir fazendo? Eu acho que eu já consigo aplicar esse algoritmo nos dados de verdade
## Posso ir fazendo a base de verdade.
