rm(list=ls())
require(mboost)
require(stats)
### for data visualization
require(ggplot2) 
require(gridExtra)
require(ggthemes)
### adquiring our functions and objects
source('functions.R')
source('simulated_models.R')

modelo_a = catch_r2(yta, name = 'modelo_a')
modelo_b = catch_r2(ytb, name = 'modelo_b')
modelo_c = catch_r2(ytc, name = 'modelo_c')
modelo_d = catch_r2(ytd, name = 'modelo_d')
modelo_e = catch_r2(yte, name = 'modelo_e')
modelo_f = catch_r2(ytf, name = 'modelo_f')

models = list('a' = modelo_a, 'b' = modelo_b, 'c' = modelo_c, 'd' = modelo_d, 'e'= modelo_e,'f'= modelo_f)

for(i in 1:6){
  models[[i]]$plot = ggplot(data = data.frame('h' = models[[i]]$h,'R2' = models[[i]]$vetor_h), mapping = aes(x = h, y = R2)) +
  geom_line() + geom_point(shape = 21, colour = "blue") + ylab("R²") + ggtitle(models[[i]]$name) + 
  ylim(0,1) + #scale_y_continuous(labels=as.character(seq(0,1, by=0.1)),breaks=seq(0,1, by=0.1)) +
  scale_x_continuous(labels=as.character(models[[i]]$h),breaks=models[[i]]$h) + theme_calc()
}

grid.arrange(models$a$plot, models$b$plot, models$c$plot, models$d$plot, models$e$plot,models$f$plot, ncol=2)

#df = data.frame(lag_1 = NA, lag_2=NA, lag_3 = NA, lag_4 = NA, lag_5 = NA, lag_6=NA, 
#                lag_7 = NA, lag_8 = NA,lag_9 = NA, lag_10=NA, lag_11 = NA, lag_12 = NA)
#for(j in 1:12){
#    for(i in 1:12){
#    modelo_1 = linear_model(yt=yt,h=i, p=j, n_out = 200)
#    df[i,j] <- mean(modelo_1$aic)
#}}
#which.min(lapply(df, mean))

