#set.seed(2506)
t = 500

yta = array(dim = t) ## Criação do vetor 1 linear
yta[c(1)] = rnorm(1, mean = 0, sd = 0.5)
for(i in 2:t){
  yta[i] = ifelse(yta[i-1]<=0.5,-0.5*yta[i-1],0) + ifelse(yta[i-1]>0.5,0.9*yta[i-1],0) + 0.5*rnorm(1, mean = 0, sd = 0.5)
  
}

ytb = array(dim = t) ## Criação do vetor 2 não linear
ytb[c(1)] = rnorm(1, mean = 0, sd = 0.25)
for(i in 2:t){
  ytb[i] = 0.4*(5-ytb[i-1]^2)/(1+ytb[i-1]^2)+rnorm(1, mean = 0, sd = 0.5)
  
}

ytc = array(dim = t) ## Criação do vetor 1 linear
ytc[c(1)] = rnorm(1, mean = 0, sd = 0.5)
for(i in 2:t){
  ytc[i] = (0.5+ 2*exp(-ytc[i-1]^2))*ytc[i-1] + rnorm(1, mean = 0, sd = 0.5)
}

ytd = array(dim = t) ## Criação do vetor 2 não linear
ytd[c(1,2,3,4,5,6,7,8,9,10)] = rnorm(10, mean = 0, sd = 0.1)
for(i in 11:t){
  ytd[i] = (0.4-2*exp(-50*ytd[i-6]^2))*ytd[i-6] + (0.5-0.5*exp(-50*ytd[i-10]^2))*ytd[i-10] +rnorm(1, mean = 0, sd = 0.1)
  
}

yte = array(dim = t) ## Criação do vetor 1 linear
yte[c(1,2)] = rnorm(2, mean = 0, sd = 0.1)
for(i in 3:t){
  yte[i] = -0.4*((3-yte[i-1]^2)/(1+yte[i-1]^2))+0.6*((3-(yte[i-2]-0.5)^3)/(1+(yte[i-2]-0.5)^4)) + rnorm(1, mean = 0, sd = 0.1)
  
}

ytf = array(dim = t) ## Criação do vetor 1 linear
ytf[c(1,2,3)] = rnorm(3, mean = 0, sd = 0.1)
for(i in 4:t){
  ytf[i] = 0.21*ytf[i-1] + 0.35*ytf[i-2] + 0.17*ytf[i-3] + rnorm(1, mean = 0, sd = 0.1)
  
}