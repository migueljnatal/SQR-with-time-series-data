## MC Simulation -  Location-shift model 1 Galvao (2013) - 50 quantile levels

#install.packages("knitr")
#install.packages("jtools")
#install.packages("conquer")
#install.packages('reshape2')
#install.packages('tidyr')
#install.packages("GoFKernel")
#install.packages('scales')
library(scales)
library(GoFKernel) 
library(reshape2)
library(tidyr)
library(ggplot2)
library(cowplot)
library(latex2exp) 
library(conquer)
library(jtools)
library(knitr)
library(quantreg)  
library(stargazer)
#source('sqr_function.r') 

#set.seed(1234)

# The variables

x = list()
nrep = 5000                                         #  number of replications
Tsize = 1100                                        #  sample size (before burn-in)
u = list()
Y = list()
ones = list()
l.x = list()
l.Y = list()
X = list()
X1 = list()
#taus <- c(0.1, 0.25, 0.5, 0.75, 0.9)                  # this 'taus' vector can be any of length we want; Galvao et al analyze only 5 quantiles
taus = seq(0.01,0.99,0.02)
A = array(0, dim = c(length(taus), 4, nrep))           # arrays for storing the simulated coef; A stands for rq, B & C for conquer and D for sqr functions
B = array(0, dim = c(length(taus), 4, nrep))
C = array(0, dim = c(length(taus), 4, nrep))
C2 = array(0, dim = c(length(taus), 4, nrep))
C3 = array(0, dim = c(length(taus), 4, nrep))
C_3 = array(0, dim = c(length(taus), 4, nrep))
C_2 = array(0, dim = c(length(taus), 4, nrep))
C8 = array(0, dim = c(length(taus), 4, nrep))
C16 = array(0, dim = c(length(taus), 4, nrep))
C_16 = array(0, dim = c(length(taus), 4, nrep))
C_8 = array(0, dim = c(length(taus), 4, nrep))
h_scale = 1
z = matrix(0, length(taus), nrep)
iqr_z =  matrix(0, length(taus), nrep) 
sigma_z =  matrix(0, length(taus), nrep) 
h =  matrix(0, length(taus), nrep) 
h2 =  matrix(0, length(taus), nrep) 
h3 =  matrix(0, length(taus), nrep) 
h_2 =  matrix(0, length(taus), nrep) 
h_3 =  matrix(0, length(taus), nrep) 
h8 =  matrix(0, length(taus), nrep) 
h16 =  matrix(0, length(taus), nrep) 
h_8 =  matrix(0, length(taus), nrep) 
h_16 =  matrix(0, length(taus), nrep) 


# Running the simulation

for(j in 1:nrep){ 
  x[[j]] = rnorm(n = Tsize, 0, 1)
  u[[j]] = rnorm(n = Tsize, 0, 1)
  Y[[j]] = 0 
  
  # The model to be estimated  
  
  for(t in 2:Tsize){
    
    Y[[j]][1] = 0
    
    Y[[j]][t] = 0.5*Y[[j]][t-1] + 0.5*x[[j]][t]  + 0.5*x[[j]][t-1] + u[[j]][t]
    
  }
  
  l.x[[j]] <- c(0, x[[j]][1:(Tsize-1)])   # c(0, x[[j]][1:(Tsize-1)]) ; c(0, x[1:(n-1)]) vs lag(x[[j]], n=1L)
  l.x[[j]][1] = 0                         # inserting a 0 for the 1st observation in the lagged x variable (instead of removing the NA value)
  
  l.Y[[j]] <- c(0, Y[[j]][1:(Tsize-1)])
  l.Y[[j]][1] = 0                         # same thing for Y  
  
  
  burn = 100                              # discarding the 1st 100 observations, as Galvao et al do in their paper
  Y[[j]] = Y[[j]][(burn+1):Tsize]
  x[[j]] = x[[j]][(burn+1):Tsize]
  u[[j]] = u[[j]][(burn+1):Tsize]
  l.Y[[j]] = l.Y[[j]][(burn+1):Tsize]
  l.x[[j]] = l.x[[j]][(burn+1):Tsize]
  
  ones[[j]] = rep(1, Tsize-burn)
  #ones[[j]] = rep(1, 1000)
  X[[j]] = cbind(l.Y[[j]], x[[j]], l.x[[j]]) 
  X1[[j]] = cbind(ones[[j]],l.Y[[j]], x[[j]], l.x[[j]]) 
  
  for (i in 1:length(taus)){
    qadl_rq <- rq(Y[[j]]~ l.Y[[j]] + x[[j]] + l.x[[j]] , tau=taus[i])  
    A[i,,j] = qadl_rq$coefficients # standard QR 
    
    #qadl_conquer <- conquer(X[[j]], Y[[j]], tau=taus[i], kernel = "Gaussian", h = 0)
    #B[i,,j] = conquer(as.matrix(X[[j]]), Y[[j]], tau=taus[i])$coef                 # 1st conquer estimation (NOT USED ANYMORE); we use A[i,,j] instead
    # Getting the optimal bandwidth
    
    z = Y[[j]] - cbind(rep(1,length(Y[[j]])),as.matrix(X[[j]]))%*%A[i,,j]#B[i,,j]
    iqr_z = quantile(z, .75) - quantile(z, .25)
    sigma_z = min(sd(z), ((iqr_z)/1.34898))
    h[i,j] = h_scale*(1.06*sigma_z)/(length(Y[[j]])^(1/5))
    
    # Demais niveis de $\zeta$ = h 'rule of thumb' 
    h_3[i,j]        = h[i,j]/4
    h_2[i,j]        = h[i,j]/2
    h2[i,j]         = h[i,j]*2
    h3[i,j]         = h[i,j]*4
    h_8[i,j]        = h[i,j]/8
    h8[i,j]         = h[i,j]*8
    h_16[i,j]        = h[i,j]/16
    h16[i,j]         = h[i,j]*16
    
    qadl_conquer_2 <- conquer(as.matrix(X[[j]]), Y[[j]], tau=taus[i], h = h[i,j])   # 2nd conquer estimation (optimal bandwidth = h)
    C[i,,j] = qadl_conquer_2$coeff
    
    qadl_conquer_2.1 <- conquer(as.matrix(X[[j]]), Y[[j]], tau=taus[i], h = h2[i,j])   # 2nd conquer estimation (optimal bandwidth = h*2)
    C2[i,,j] = qadl_conquer_2.1$coeff
    
    qadl_conquer_2.2 <- conquer(as.matrix(X[[j]]), Y[[j]], tau=taus[i], h = h3[i,j])   # 2nd conquer estimation (optimal bandwidth = h*4)
    C3[i,,j] = qadl_conquer_2.2$coeff
    
    qadl_conquer_2.3 <- conquer(as.matrix(X[[j]]), Y[[j]], tau=taus[i], h = h_3[i,j])   # 2nd conquer estimation (optimal bandwidth = h/4)
    C_3[i,,j] = qadl_conquer_2.3$coeff
    
    qadl_conquer_2.4 <- conquer(as.matrix(X[[j]]), Y[[j]], tau=taus[i], h = h_2[i,j])   # 2nd conquer estimation (optimal bandwidth = h/2)
    C_2[i,,j] = qadl_conquer_2.4$coeff
    
    qadl_conquer_2.5 <- conquer(as.matrix(X[[j]]), Y[[j]], tau=taus[i], h = h8[i,j])   # 2nd conquer estimation (optimal bandwidth = h*8)
    C8[i,,j] = qadl_conquer_2.5$coeff

    qadl_conquer_2.6 <- conquer(as.matrix(X[[j]]), Y[[j]], tau=taus[i], h = h_8[i,j])   # 2nd conquer estimation (optimal bandwidth = h*16)
    C_8[i,,j] = qadl_conquer_2.6$coeff

    qadl_conquer_2.7 <- conquer(as.matrix(X[[j]]), Y[[j]], tau=taus[i], h = h_16[i,j])   # 2nd conquer estimation (optimal bandwidth = h/16)
    C_16[i,,j] = qadl_conquer_2.7$coeff

    qadl_conquer_2.8 <- conquer(as.matrix(X[[j]]), Y[[j]], tau=taus[i], h = h16[i,j])   # 2nd conquer estimation (optimal bandwidth = h/8)
    C16[i,,j] = qadl_conquer_2.8$coeff
  }
}   

### Evaluating the estimators

# canonical estimator (Koenker)
intercept_rq = A[,1,]
alphas1_rq = A[,2,]
betas1_rq = A[,3,]
betas2_rq = A[,4,]

# smoothed estimator (FGH)
# h 
intercept_cq = C[,1,]
alphas1_cq = C[,2,]
betas1_cq = C[,3,]
betas2_cq = C[,4,] 
# h*2
intercept_cq2 = C2[,1,]
alphas1_cq2 = C2[,2,]
betas1_cq2 = C2[,3,]
betas2_cq2 = C2[,4,] 
# h*4
intercept_cq3 = C3[,1,]
alphas1_cq3 = C3[,2,]
betas1_cq3 = C3[,3,]
betas2_cq3 = C3[,4,]
# h/4
intercept_cq_3 = C_3[,1,]
alphas1_cq_3 = C_3[,2,]
betas1_cq_3 = C_3[,3,]
betas2_cq_3 = C_3[,4,] 
# h/2
intercept_cq_2 = C_2[,1,]
alphas1_cq_2 = C_2[,2,]
betas1_cq_2 = C_2[,3,]
betas2_cq_2 = C_2[,4,]
#h*8
intercept_cq8 = C8[,1,]
alphas1_cq8 = C8[,2,]
betas1_cq8 = C8[,3,]
betas2_cq8 = C8[,4,]
# h*16
intercept_cq16 = C16[,1,]
alphas1_cq16 = C16[,2,]
betas1_cq16 = C16[,3,]
betas2_cq16 = C16[,4,]
# h/16
intercept_cq_16 = C_16[,1,]
alphas1_cq_16 = C_16[,2,]
betas1_cq_16 = C_16[,3,]
betas2_cq_16 = C_16[,4,]
# h/8
intercept_cq_8 = C_8[,1,]
alphas1_cq_8 = C_8[,2,]
betas1_cq_8 = C_8[,3,]
betas2_cq_8 = C_8[,4,]

# Variables for storing Bias and MSE Values
intercept_rq_bias = c()
intercept_cq_bias = c()
alphas1_rq_bias = c()
alphas1_cq_bias = c()
betas1_rq_bias = c()
betas1_cq_bias = c()
betas2_rq_bias = c()
betas2_cq_bias = c()
# Other levels of bandwidth (h)
intercept_cq2_bias = c()      
alphas1_cq2_bias = c()
betas1_cq2_bias = c()
betas2_cq2_bias = c()

intercept_cq3_bias = c()      
alphas1_cq3_bias = c()
betas1_cq3_bias = c()
betas2_cq3_bias = c()

intercept_cq_3_bias = c()      
alphas1_cq_3_bias = c()
betas1_cq_3_bias = c()
betas2_cq_3_bias = c()

intercept_cq_2_bias = c()     
alphas1_cq_2_bias = c()
betas1_cq_2_bias = c()
betas2_cq_2_bias = c()

intercept_cq8_bias = c()
alphas1_cq8_bias = c()
betas1_cq8_bias = c()
betas2_cq8_bias = c()

intercept_cq16_bias = c()
alphas1_cq16_bias = c()
betas1_cq16_bias = c()
betas2_cq16_bias = c()

intercept_cq_16_bias = c()
alphas1_cq_16_bias = c()
betas1_cq_16_bias = c()
betas2_cq_16_bias = c()

intercept_cq_8_bias = c()
alphas1_cq_8_bias = c()
betas1_cq_8_bias = c()
betas2_cq_8_bias = c()

# MSE
intercept_rq_mse = c()
intercept_cq_mse = c()
alphas1_rq_mse = c()
alphas1_cq_mse = c()
betas1_rq_mse = c()
betas1_cq_mse = c()
betas2_rq_mse = c()
betas2_cq_mse = c()

# Other levels of bandwidth (h)
intercept_cq2_mse = c()
alphas1_cq2_mse = c()
betas1_cq2_mse = c()
betas2_cq2_mse = c()

intercept_cq3_mse = c()
alphas1_cq3_mse = c()
betas1_cq3_mse = c()
betas2_cq3_mse = c() 

intercept_cq_3_mse = c()
alphas1_cq_3_mse = c()
betas1_cq_3_mse = c()
betas2_cq_3_mse = c()

intercept_cq_2_mse = c()
alphas1_cq_2_mse = c()
betas1_cq_2_mse = c()
betas2_cq_2_mse = c()

intercept_cq8_mse = c()
alphas1_cq8_mse = c()
betas1_cq8_mse = c()
betas2_cq8_mse = c()

intercept_cq16_mse = c()
alphas1_cq16_mse = c()
betas1_cq16_mse = c()
betas2_cq16_mse = c()

intercept_cq_16_mse = c()
alphas1_cq_16_mse = c()
betas1_cq_16_mse = c()
betas2_cq_16_mse = c()

intercept_cq_8_mse = c()
alphas1_cq_8_mse = c()
betas1_cq_8_mse = c()
betas2_cq_8_mse = c()

### COMPUTING METRICS: Bias, RMSE and MSE

for(s in 1:length(taus)){
  
  # Bias 1st - 50th quantiles
  intercept_rq_bias[s] = mean(A[s,1,]) - qt(taus[s],df=3)
  intercept_cq_bias[s] = mean(C[s,1,]) - qt(taus[s],df=3)
  alphas1_rq_bias[s] = mean(A[s,2,]) - 0.5
  alphas1_cq_bias[s] = mean(C[s,2,]) - 0.5
  betas1_rq_bias[s] = mean(A[s,3,]) - 0.5
  betas1_cq_bias[s] = mean(C[s,3,]) - 0.5
  betas2_rq_bias[s] = mean(A[s,4,]) - 0.5
  betas2_cq_bias[s] = mean(C[s,4,]) - 0.5 
  # h*2
  intercept_cq2_bias[s] = mean(C2[s,1,]) - qt(taus[s],df=3)
  alphas1_cq2_bias[s] = mean(C2[s,2,]) - 0.5
  betas1_cq2_bias[s] = mean(C2[s,3,]) - 0.5
  betas2_cq2_bias[s] = mean(C2[s,4,]) - 0.5 
  # h*4
  intercept_cq3_bias[s] = mean(C3[s,1,]) - qt(taus[s],df=3)
  alphas1_cq3_bias[s] = mean(C3[s,2,]) - 0.5
  betas1_cq3_bias[s] = mean(C3[s,3,]) - 0.5
  betas2_cq3_bias[s] = mean(C3[s,4,]) - 0.5
  # h/4
  intercept_cq_3_bias[s] = mean(C_3[s,1,]) - qt(taus[s],df=3)
  alphas1_cq_3_bias[s] = mean(C_3[s,2,]) - 0.5
  betas1_cq_3_bias[s] = mean(C_3[s,3,]) - 0.5
  betas2_cq_3_bias[s] = mean(C_3[s,4,]) - 0.5
  # h/2
  intercept_cq_2_bias[s] = mean(C_2[s,1,]) - qt(taus[s],df=3)
  alphas1_cq_2_bias[s] = mean(C_2[s,2,]) - 0.5
  betas1_cq_2_bias[s] = mean(C_2[s,3,]) - 0.5
  betas2_cq_2_bias[s] = mean(C_2[s,4,]) - 0.5 
  # # h*8
  # intercept_cq8_bias[s] = mean(C8[s,1,]) - qt(taus[s],df=3)
  # alphas1_cq8_bias[s] = mean(C8[s,2,]) - 0.5
  # betas1_cq8_bias[s] = mean(C8[s,3,]) - 0.5
  # betas2_cq8_bias[s] = mean(C8[s,4,]) - 0.5 
  # # h*16
  # intercept_cq16_bias[s] = mean(C16[s,1,]) - qt(taus[s],df=3)
  # alphas1_cq16_bias[s] = mean(C16[s,2,]) - 0.5
  # betas1_cq16_bias[s] = mean(C16[s,3,]) - 0.5
  # betas2_cq16_bias[s] = mean(C16[s,4,]) - 0.5
  # # h/16
  # intercept_cq_16_bias[s] = mean(C_16[s,1,]) - qt(taus[s],df=3)
  # alphas1_cq_16_bias[s] = mean(C_16[s,2,]) - 0.5
  # betas1_cq_16_bias[s] = mean(C_16[s,3,]) - 0.5
  # betas2_cq_16_bias[s] = mean(C_16[s,4,]) - 0.5
  # # h/8
  # intercept_cq_8_bias[s] = mean(C_8[s,1,]) - qt(taus[s],df=3)
  # alphas1_cq_8_bias[s] = mean(C_8[s,2,]) - 0.5
  # betas1_cq_8_bias[s] = mean(C_8[s,3,]) - 0.5
  # betas2_cq_8_bias[s] = mean(C_8[s,4,]) - 0.5
  
  # MSE 1st - 50th quantiles
  intercept_rq_mse[s] = mean((A[s,1,] - qnorm(taus[s],0,1))^2)
  intercept_cq_mse[s] = mean((C[s,1,] - qnorm(taus[s],0,1))^2)
  alphas1_rq_mse[s] = mean((A[s,2,] - 0.5)^2)
  alphas1_cq_mse[s] = mean((C[s,2,] - 0.5)^2)
  betas1_rq_mse[s] = mean((A[s,3,] - 0.5)^2)
  betas1_cq_mse[s] = mean((C[s,3,] - 0.5)^2)
  betas2_rq_mse[s] = mean((A[s,4,] - 0.5)^2)
  betas2_cq_mse[s] = mean((C[s,4,] - 0.5)^2)
  # h*2
  intercept_cq2_mse[s] = mean((C2[s,1,] - qnorm(taus[s],0,1))^2)
  alphas1_cq2_mse[s] = mean((C2[s,2,] - 0.5)^2)
  betas1_cq2_mse[s] = mean((C2[s,3,] - 0.5)^2)
  betas2_cq2_mse[s] = mean((C2[s,4,] - 0.5)^2)
  # h*3
  intercept_cq3_mse[s] = mean((C3[s,1,] - qnorm(taus[s],0,1))^2)
  alphas1_cq3_mse[s] = mean((C3[s,2,] - 0.5)^2)
  betas1_cq3_mse[s] = mean((C3[s,3,] - 0.5)^2)
  betas2_cq3_mse[s] = mean((C3[s,4,] - 0.5)^2)
  # h*1/3
  intercept_cq_3_mse[s] = mean((C_3[s,1,] - qnorm(taus[s],0,1))^2)
  alphas1_cq_3_mse[s] = mean((C_3[s,2,] - 0.5)^2)
  betas1_cq_3_mse[s] = mean((C_3[s,3,] - 0.5)^2)
  betas2_cq_3_mse[s] = mean((C_3[s,4,] - 0.5)^2)
  # h*1/2
  intercept_cq_2_mse[s] = mean((C_2[s,1,] - qnorm(taus[s],0,1))^2)
  alphas1_cq_2_mse[s] = mean((C_2[s,2,] - 0.5)^2)
  betas1_cq_2_mse[s] = mean((C_2[s,3,] - 0.5)^2)
  betas2_cq_2_mse[s] = mean((C_2[s,4,] - 0.5)^2)
  # h*8
  intercept_cq8_mse[s] = mean((C8[s,1,] - qnorm(taus[s],0,1))^2)
  alphas1_cq8_mse[s] = mean((C8[s,2,] - 0.5)^2)
  betas1_cq8_mse[s] = mean((C8[s,3,] - 0.5)^2)
  betas2_cq8_mse[s] = mean((C8[s,4,] - 0.5)^2)
  # h*16
  intercept_cq16_mse[s] = mean((C16[s,1,] - qnorm(taus[s],0,1))^2)
  alphas1_cq16_mse[s] = mean((C16[s,2,] - 0.5)^2)
  betas1_cq16_mse[s] = mean((C16[s,3,] - 0.5)^2)
  betas2_cq16_mse[s] = mean((C16[s,4,] - 0.5)^2)
  # h/16
  intercept_cq_16_mse[s] = mean((C_16[s,1,] - qnorm(taus[s],0,1))^2)
  alphas1_cq_16_mse[s] = mean((C_16[s,2,] - 0.5)^2)
  betas1_cq_16_mse[s] = mean((C_16[s,3,] - 0.5)^2)
  betas2_cq_16_mse[s] = mean((C_16[s,4,] - 0.5)^2)
  # h/8
  intercept_cq_8_mse[s] = mean((C_8[s,1,] - qnorm(taus[s],0,1))^2)
  alphas1_cq_8_mse[s] = mean((C_8[s,2,] - 0.5)^2)
  betas1_cq_8_mse[s] = mean((C_8[s,3,] - 0.5)^2)
  betas2_cq_8_mse[s] = mean((C_8[s,4,] - 0.5)^2)
}

# Bias 
alphas0_rq_bias = intercept_rq_bias 
alphas0_cq_bias = intercept_cq_bias
# h*2
alphas0_cq2_bias = intercept_cq2_bias
# h*4
alphas0_cq3_bias = intercept_cq3_bias
# h/4
alphas0_cq_3_bias = intercept_cq_3_bias
# h/2
alphas0_cq_2_bias = intercept_cq_2_bias
# h*8
alphas0_cq8_bias = intercept_cq8_bias
# h*16
alphas0_cq16_bias = intercept_cq16_bias
# h/16
alphas0_cq_16_bias = intercept_cq_16_bias
# h/8
alphas0_cq_8_bias = intercept_cq_8_bias

# MSE
alphas0_rq_mse = intercept_rq_mse 
alphas0_cq_mse = intercept_cq_mse
# h*2
alphas0_cq2_mse = intercept_cq2_mse 
# h*4
alphas0_cq3_mse = intercept_cq3_mse 
# h/4
alphas0_cq_3_mse = intercept_cq_3_mse 
# h/2
alphas0_cq_2_mse = intercept_cq_2_mse 
# h*8
alphas0_cq8_mse = intercept_cq8_mse 
# h*16
alphas0_cq16_mse = intercept_cq16_mse 
# h/16
alphas0_cq_16_mse = intercept_cq_16_mse 
# h/8
alphas0_cq_8_mse = intercept_cq_8_mse 


### Stargazer tables 

### FULL TABLES -  Bias & MSE - h original

table_all_metrics  = cbind(intercept_rq_bias, alphas1_rq_bias, betas1_rq_bias, betas2_rq_bias, 
                           intercept_cq_bias, alphas1_cq_bias, betas1_cq_bias, betas2_cq_bias, 
                           intercept_rq_mse, alphas1_rq_mse, betas1_rq_mse, betas2_rq_mse,
                           intercept_cq_mse, alphas1_cq_mse, betas1_cq_mse, betas2_cq_mse)

table_all_metrics = round(table_all_metrics, 4)

colnames(table_all_metrics) <- c('alpha_0_rq','alpha_1_rq', 'beta_1_rq', 'beta_2_rq', 'alpha_0_cq','alpha_1_cq', 'beta_1_cq', 'beta_2_cq', 
                                 'alpha_0_rq','alpha_1_rq', 'beta_1_rq', 'beta_2_rq', 'alpha_0_cq','alpha_1_cq', 'beta_1_cq', 'beta_2_cq')

rownames(table_all_metrics) <- as.character(taus)

stargazer(table_all_metrics, title = 'Location-shift model x: bias and MSE of estimators', digits=4, digits.extra = 4, decimal.mark = '.')

# FULL TABLES -  Bias & mse - h*2
# 
# table_all_metrics_h2  = cbind(intercept_rq_bias, alphas1_rq_bias, betas1_rq_bias, betas2_rq_bias, 
#                               intercept_cq2_bias, alphas1_cq2_bias, betas1_cq2_bias, betas2_cq2_bias, 
#                               intercept_rq_mse, alphas1_rq_mse, betas1_rq_mse, betas2_rq_mse,
#                               intercept_cq2_mse, alphas1_cq2_mse, betas1_cq2_mse, betas2_cq2_mse)
# 
# table_all_metrics_h2 = round(table_all_metrics_h2, 4)
# 
# colnames(table_all_metrics_h2) <- c('alpha_0_rq','alpha_1_rq', 'beta_1_rq', 'beta_2_rq', 'alpha_0_cq','alpha_1_cq', 'beta_1_cq', 'beta_2_cq', 
#                                     'alpha_0_rq','alpha_1_rq', 'beta_1_rq', 'beta_2_rq', 'alpha_0_cq','alpha_1_cq', 'beta_1_cq', 'beta_2_cq')
# 
# rownames(table_all_metrics_h2) <- as.character(taus)
# 
# stargazer(table_all_metrics_h2, title = 'Location-shift model x: bias and mse of estimators - bandwidth = h*2', digits=4, digits.extra = 4, decimal.mark = '.')


### Plots - for each estimation, we plot the mean of the estimated values for each coefficient: mpq = mean per quantile (considering all replications)

# Storing coefficients for plotting
intercept_rq_mpq = c()
intercept_cq_mpq = c()
alphas1_rq_mpq = c()
alphas1_cq_mpq = c()
betas1_rq_mpq = c()
betas1_cq_mpq = c()
betas2_rq_mpq = c()
betas2_cq_mpq = c()
# h*2
intercept_cq2_mpq = c()
alphas1_cq2_mpq = c()
betas1_cq2_mpq = c()
betas2_cq2_mpq = c()
# h*4
intercept_cq3_mpq = c()
alphas1_cq3_mpq = c()
betas1_cq3_mpq = c()
betas2_cq3_mpq = c()
# h/4
intercept_cq_3_mpq = c()
alphas1_cq_3_mpq = c()
betas1_cq_3_mpq = c()
betas2_cq_3_mpq = c()
# h/2
intercept_cq_2_mpq = c()
alphas1_cq_2_mpq = c()
betas1_cq_2_mpq = c()
betas2_cq_2_mpq = c() 
# h*8
intercept_cq8_mpq = c()
alphas1_cq8_mpq = c()
betas1_cq8_mpq = c()
betas2_cq8_mpq = c()
# h*16
intercept_cq16_mpq = c()
alphas1_cq16_mpq = c()
betas1_cq16_mpq = c()
betas2_cq16_mpq = c()
# h/16
intercept_cq_16_mpq = c()
alphas1_cq_16_mpq = c()
betas1_cq_16_mpq = c()
betas2_cq_16_mpq = c()
# h/8
intercept_cq_8_mpq = c()
alphas1_cq_8_mpq = c()
betas1_cq_8_mpq = c()
betas2_cq_8_mpq = c() 

# Variance
intercept_rq_vpq = c()
intercept_cq_vpq = c()
alphas1_rq_vpq = c()
alphas1_cq_vpq = c()
betas1_rq_vpq = c()
betas1_cq_vpq = c()
betas2_rq_vpq = c()
betas2_cq_vpq = c()
# h*2
intercept_cq2_vpq = c()
alphas1_cq2_vpq = c()
betas1_cq2_vpq = c()
betas2_cq2_vpq = c()
# h*4
intercept_cq3_vpq = c()
alphas1_cq3_vpq = c()
betas1_cq3_vpq = c()
betas2_cq3_vpq = c()
# h/4
intercept_cq_3_vpq = c()
alphas1_cq_3_vpq = c()
betas1_cq_3_vpq = c()
betas2_cq_3_vpq = c()
# h/2
intercept_cq_2_vpq = c()
alphas1_cq_2_vpq = c()
betas1_cq_2_vpq = c()
betas2_cq_2_vpq = c()
# h*8
intercept_cq8_vpq = c()
alphas1_cq8_vpq = c()
betas1_cq8_vpq = c()
betas2_cq8_vpq = c()
# h*16
intercept_cq16_vpq = c()
alphas1_cq16_vpq = c()
betas1_cq16_vpq = c()
betas2_cq16_vpq = c()
# h/16
intercept_cq_16_vpq = c()
alphas1_cq_16_vpq = c()
betas1_cq_16_vpq = c()
betas2_cq_16_vpq = c()
# h/8
intercept_cq_8_vpq = c()
alphas1_cq_8_vpq = c()
betas1_cq_8_vpq = c()
betas2_cq_8_vpq = c()

for(g in 1:length(taus)){
  # Mean per quantiles
  intercept_rq_mpq[g] = mean(A[g,1,])
  intercept_cq_mpq[g] = mean(C[g,1,])
  alphas1_rq_mpq[g] = mean(A[g,2,])
  alphas1_cq_mpq[g] = mean(C[g,2,])
  betas1_rq_mpq[g] = mean(A[g,3,])
  betas1_cq_mpq[g] = mean(C[g,3,])
  betas2_rq_mpq[g] = mean(A[g,4,])
  betas2_cq_mpq[g] = mean(C[g,4,])
  # h*2
  intercept_cq2_mpq[g] = mean(C2[g,1,])
  alphas1_cq2_mpq[g] = mean(C2[g,2,])
  betas1_cq2_mpq[g] = mean(C2[g,3,])
  betas2_cq2_mpq[g] = mean(C2[g,4,])
  # h*3
  intercept_cq3_mpq[g] = mean(C3[g,1,])
  alphas1_cq3_mpq[g] = mean(C3[g,2,])
  betas1_cq3_mpq[g] = mean(C3[g,3,])
  betas2_cq3_mpq[g] = mean(C3[g,4,])
  # h/3
  intercept_cq_3_mpq[g] = mean(C_3[g,1,])
  alphas1_cq_3_mpq[g] = mean(C_3[g,2,])
  betas1_cq_3_mpq[g] = mean(C_3[g,3,])
  betas2_cq_3_mpq[g] = mean(C_3[g,4,])
  # h/2
  intercept_cq_2_mpq[g] = mean(C_2[g,1,])
  alphas1_cq_2_mpq[g] = mean(C_2[g,2,])
  betas1_cq_2_mpq[g] = mean(C_2[g,3,])
  betas2_cq_2_mpq[g] = mean(C_2[g,4,])
  # h*8
  intercept_cq8_mpq[g] = mean(C8[g,1,])
  alphas1_cq8_mpq[g] = mean(C8[g,2,])
  betas1_cq8_mpq[g] = mean(C8[g,3,])
  betas2_cq8_mpq[g] = mean(C8[g,4,])
  # h*16
  intercept_cq16_mpq[g] = mean(C16[g,1,])
  alphas1_cq16_mpq[g] = mean(C16[g,2,])
  betas1_cq16_mpq[g] = mean(C16[g,3,])
  betas2_cq16_mpq[g] = mean(C16[g,4,])
  # h/16
  intercept_cq_16_mpq[g] = mean(C_16[g,1,])
  alphas1_cq_16_mpq[g] = mean(C_16[g,2,])
  betas1_cq_16_mpq[g] = mean(C_16[g,3,])
  betas2_cq_16_mpq[g] = mean(C_16[g,4,])
  # h/8
  intercept_cq_8_mpq[g] = mean(C_8[g,1,])
  alphas1_cq_8_mpq[g] = mean(C_8[g,2,])
  betas1_cq_8_mpq[g] = mean(C_8[g,3,])
  betas2_cq_8_mpq[g] = mean(C_8[g,4,])
  
  # Variance per quantiles
  intercept_rq_vpq[g] = var(A[g,1,])
  intercept_cq_vpq[g] = var(C[g,1,])
  alphas1_rq_vpq[g] = var(A[g,2,])
  alphas1_cq_vpq[g] = var(C[g,2,])
  betas1_rq_vpq[g] = var(A[g,3,])
  betas1_cq_vpq[g] = var(C[g,3,])
  betas2_rq_vpq[g] = var(A[g,4,])
  betas2_cq_vpq[g] = var(C[g,4,])
  # h*2
  intercept_cq2_vpq[g] = var(C2[g,1,])
  alphas1_cq2_vpq[g] = var(C2[g,2,])
  betas1_cq2_vpq[g] = var(C2[g,3,])
  betas2_cq2_vpq[g] = var(C2[g,4,])
  # h*4
  intercept_cq3_vpq[g] = var(C3[g,1,])
  alphas1_cq3_vpq[g] = var(C3[g,2,])
  betas1_cq3_vpq[g] = var(C3[g,3,])
  betas2_cq3_vpq[g] = var(C3[g,4,])
  # h/4
  intercept_cq_3_vpq[g] = var(C_3[g,1,])
  alphas1_cq_3_vpq[g] = var(C_3[g,2,])
  betas1_cq_3_vpq[g] = var(C_3[g,3,])
  betas2_cq_3_vpq[g] = var(C_3[g,4,])
  # h/2
  intercept_cq_2_vpq[g] = var(C_2[g,1,])
  alphas1_cq_2_vpq[g] = var(C_2[g,2,])
  betas1_cq_2_vpq[g] = var(C_2[g,3,])
  betas2_cq_2_vpq[g] = var(C_2[g,4,])
  # h*8
  intercept_cq8_vpq[g] = var(C8[g,1,])
  alphas1_cq8_vpq[g] = var(C8[g,2,])
  betas1_cq8_vpq[g] = var(C8[g,3,])
  betas2_cq8_vpq[g] = var(C8[g,4,])
  # h*16
  intercept_cq16_vpq[g] = var(C16[g,1,])
  alphas1_cq16_vpq[g] = var(C16[g,2,])
  betas1_cq16_vpq[g] = var(C16[g,3,])
  betas2_cq16_vpq[g] = var(C16[g,4,])
  # h/16
  intercept_cq_16_vpq[g] = var(C_16[g,1,])
  alphas1_cq_16_vpq[g] = var(C_16[g,2,])
  betas1_cq_16_vpq[g] = var(C_16[g,3,])
  betas2_cq_16_vpq[g] = var(C_16[g,4,])
  # h/8
  intercept_cq_8_vpq[g] = var(C_8[g,1,])
  alphas1_cq_8_vpq[g] = var(C_8[g,2,])
  betas1_cq_8_vpq[g] = var(C_8[g,3,])
  betas2_cq_8_vpq[g] = var(C_8[g,4,])

}

# Final plots - intercept (alpha_0), alpha_1, beta_1 & beta_2  - MEAN per quantile - - h, h*2, h*4, h*1/4, h*1/2

coefs_intercept_mpq = data.frame(intercept_rq_mpq, intercept_cq_mpq, 
                                 intercept_cq2_mpq,intercept_cq3_mpq,intercept_cq_3_mpq,intercept_cq_2_mpq,
                                 intercept_cq8_mpq,intercept_cq16_mpq,intercept_cq_16_mpq,intercept_cq_8_mpq) 
coefs_alphas1_mpq = data.frame(alphas1_rq_mpq, alphas1_cq_mpq, 
                               alphas1_cq2_mpq,alphas1_cq3_mpq,alphas1_cq_3_mpq,alphas1_cq_2_mpq,
                               alphas1_cq8_mpq,alphas1_cq16_mpq,alphas1_cq_16_mpq,alphas1_cq_8_mpq) 
coefs_betas1_mpq = data.frame(betas1_rq_mpq, betas1_cq_mpq,
                              betas1_cq2_mpq,betas1_cq3_mpq,betas1_cq_3_mpq,betas1_cq_2_mpq,
                              betas1_cq8_mpq,betas1_cq16_mpq,betas1_cq_16_mpq,betas1_cq_8_mpq)
coefs_betas2_mpq = data.frame(betas2_rq_mpq, betas2_cq_mpq,
                              betas2_cq2_mpq,betas2_cq3_mpq,betas2_cq_3_mpq,betas2_cq_2_mpq,
                              betas2_cq8_mpq,betas2_cq16_mpq,betas2_cq_16_mpq,betas2_cq_8_mpq)

# xlab = TeX('$\\tau$')
# 
# ylab = TeX('$\\hat{\\alpha_0}(\\tau)$') 
# grafico_alpha0_mpq = ggplot(coefs_intercept_mpq, aes(x = taus, y = value)) + 
#   geom_line(aes(y = intercept_rq_mpq,linetype ='blank',color ='QR'),size=1.0) +
#   geom_line(aes(y = intercept_cq_mpq, linetype ='solid',color ='SQR'),size=1.0)+
#   geom_line(aes(y = intercept_cq2_mpq, linetype ='dotdash',color ='SQR_h*2'),size=1.0)+
#   geom_line(aes(y = intercept_cq3_mpq, linetype ='dotdash',color ='SQR_h*4'),size=1.0)+
#   geom_line(aes(y = intercept_cq_3_mpq, linetype ='dotdash',color ='SQR_h/4'),size=1.0)+
#   geom_line(aes(y = intercept_cq_2_mpq, linetype ='dotted',color ='SQR_h/2'),size=1.0)+
#   geom_line(aes(y = intercept_cq8_mpq, linetype ='dotdash',color ='SQR_h*8'),size=1.0)+
#   geom_line(aes(y = intercept_cq16_mpq, linetype ='dotdash',color ='SQR_h*16'),size=1.0)+
#   geom_line(aes(y = intercept_cq_16_mpq, linetype ='dotdash',color ='SQR_h/16'),size=1.0)+
#   geom_line(aes(y = intercept_cq_8_mpq, linetype ='dotted',color ='SQR_h/8'),size=1.0)+
#   scale_color_discrete(labels = unname(TeX(c("$QR$","$SQR_{\\zeta^*}$","$SQR_{16\\zeta^*}$","$SQR_{2\\zeta^*}$",
#                                              "$SQR_{4\\zeta^{*}}$","$SQR_{8\\zeta^{*}}$","$SQR_{\\zeta^*/16}$",
#                                              "$SQR_{\\zeta^*/2}$","$SQR_{\\zeta^{*}/4}$","$SQR_{\\zeta^{*}/8}$")))) +
#   geom_hline(yintercept = 0, color = 'yellow', size=1.3) +
#   labs(x =xlab,y = ylab) + #expand_limits(x=c(0, 1), y=c(-3,3)) +
#   scale_x_continuous(breaks=seq(0.01, 0.99, 0.07),limits = c(0.01, 0.99))  + 
#   theme_grey() + guides(linetype = FALSE) +
#   theme(plot.title = element_text(size = 15),
#         legend.title = element_blank(),
#         legend.text.align = 0,
#         legend.text = element_text(size = 12),
#         axis.title = element_text(size = 14))
# print(grafico_alpha0_mpq)
# 
# ylab = TeX('$\\hat{\\alpha_1}(\\tau)$') 
# grafico_alpha1_mpq = ggplot(coefs_alphas1_mpq, aes(x = taus, y = value)) + 
#   geom_line(aes(y = alphas1_rq_mpq,linetype ='blank',color ='QR'),size=1.0) +
#   geom_line(aes(y = alphas1_cq_mpq, linetype ='solid',color ='SQR'),size=1.0)+
#   geom_line(aes(y = alphas1_cq2_mpq, linetype ='dotdash',color ='SQR_h*2'),size=1.0)+
#   geom_line(aes(y = alphas1_cq3_mpq, linetype ='dotdash',color ='SQR_h*4'),size=1.0)+
#   geom_line(aes(y = alphas1_cq_3_mpq, linetype ='dotdash',color ='SQR_h/4'),size=1.0)+
#   geom_line(aes(y = alphas1_cq_2_mpq, linetype ='dotted',color ='SQR_h/2'),size=1.0)+
#   geom_line(aes(y = alphas1_cq8_mpq, linetype ='dotdash',color ='SQR_h*8'),size=1.0)+
#   geom_line(aes(y = alphas1_cq16_mpq, linetype ='dotdash',color ='SQR_h*16'),size=1.0)+
#   geom_line(aes(y = alphas1_cq_16_mpq, linetype ='dotdash',color ='SQR_h/16'),size=1.0)+
#   geom_line(aes(y = alphas1_cq_8_mpq, linetype ='dotted',color ='SQR_h/8'),size=1.0)+
#   scale_color_discrete(labels = unname(TeX(c("$QR$","$SQR_{\\zeta^*}$","$SQR_{16\\zeta^*}$","$SQR_{2\\zeta^*}$",
#                                              "$SQR_{4\\zeta^{*}}$","$SQR_{8\\zeta^{*}}$","$SQR_{\\zeta^*/16}$",
#                                              "$SQR_{\\zeta^*/2}$","$SQR_{\\zeta^{*}/4}$","$SQR_{\\zeta^{*}/8}$")))) +
#   geom_hline(yintercept = 0.5, color = 'yellow', size=1.3) +
#   labs(x =xlab,y = ylab) + #expand_limits(x=c(0, 1), y=c(0.4,0.6)) + 
#   scale_x_continuous(breaks=round(seq(0.01, 0.99, 0.07), 2),limits = c(0.01, 0.99))  + 
#   #scale_y_continuous(breaks=seq(0.4,0.6,0.02),limits = c(0.46, 0.54))  +
#   theme_grey() + guides(linetype = FALSE) +
#   theme(plot.title = element_text(size = 15),
#         legend.title = element_blank(),
#         legend.text.align = 0,
#         legend.text = element_text(size = 12),
#         axis.title = element_text(size = 14))
# print(grafico_alpha1_mpq)
# 
# ylab = TeX('$\\hat{\\theta_1}(\\tau)$')
# grafico_beta1_mpq = ggplot(coefs_betas1_mpq, aes(x = taus, y = value)) + 
#   geom_line(aes(y = betas1_rq_mpq,linetype ='blank',color ='QR'),size=1.0) +
#   geom_line(aes(y = betas1_cq_mpq, linetype ='solid',color ='SQR'),size=1.,)+
#   geom_line(aes(y = betas1_cq2_mpq, linetype ='dotdash',color ='SQR_h*2'),size=1.0)+
#   geom_line(aes(y = betas1_cq3_mpq, linetype ='dotdash',color ='SQR_h*4'),size=1.0)+
#   geom_line(aes(y = betas1_cq_3_mpq, linetype ='dotdash',color ='SQR_h/4'),size=1.0)+
#   geom_line(aes(y = betas1_cq_2_mpq, linetype ='dotted',color ='SQR_h/2'),size=1.0)+
#   geom_line(aes(y = betas1_cq8_mpq, linetype ='dotdash',color ='SQR_h*8'),size=1.0)+
#   geom_line(aes(y = betas1_cq16_mpq, linetype ='dotdash',color ='SQR_h*16'),size=1.0)+
#   geom_line(aes(y = betas1_cq_16_mpq, linetype ='dotdash',color ='SQR_h/16'),size=1.0)+
#   geom_line(aes(y = betas1_cq_8_mpq, linetype ='dotted',color ='SQR_h/8'),size=1.0)+
#   scale_color_discrete(labels = unname(TeX(c("$QR$","$SQR_{\\zeta^*}$","$SQR_{16\\zeta^*}$","$SQR_{2\\zeta^*}$",
#                                              "$SQR_{4\\zeta^{*}}$","$SQR_{8\\zeta^{*}}$","$SQR_{\\zeta^*/16}$",
#                                              "$SQR_{\\zeta^*/2}$","$SQR_{\\zeta^{*}/4}$","$SQR_{\\zeta^{*}/8}$")))) +
#   geom_hline(yintercept = 0.5, color = 'yellow', size=1.3) +
#   labs(x =xlab,y = ylab) + # ggtitle(bquote(beta[1][tau]^{tau})) + 
#   scale_x_continuous(breaks=seq(0.01, 0.99, 0.07),limits = c(0.01, 0.99))  + 
#   theme_grey() + guides(linetype = FALSE) +
#   theme(plot.title = element_text(size = 15),
#         legend.title = element_blank(),
#         legend.text.align = 0,
#         legend.text = element_text(size = 12),
#         axis.title = element_text(size = 14))
# print(grafico_beta1_mpq)
# 
# ylab = TeX('$\\hat{\\theta_2}(\\tau)$')
# grafico_beta2_mpq = ggplot(coefs_betas2_mpq, aes(x = taus, y = value)) + 
#   geom_line(aes(y = betas2_rq_mpq,linetype ='blank',color ='QR'),size=1.0) +
#   geom_line(aes(y = betas2_cq_mpq,linetype ='solid',color ='SQR'),size=1.0)+
#   geom_line(aes(y = betas2_cq2_mpq,linetype ='dotdash',color ='SQR_h*2'),size=1.0)+
#   geom_line(aes(y = betas2_cq3_mpq,linetype ='dotdash',color ='SQR_h*4'),size=1.0)+
#   geom_line(aes(y = betas2_cq_3_mpq,linetype ='dotdash',color ='SQR_h/4'),size=1.0)+
#   geom_line(aes(y = betas2_cq_2_mpq,linetype ='dotted',color ='SQR_h/2'),size=1.0)+
#   geom_line(aes(y = betas2_cq8_mpq,linetype ='dotdash',color ='SQR_h*8'),size=1.0)+
#   geom_line(aes(y = betas2_cq16_mpq,linetype ='dotdash',color ='SQR_h*16'),size=1.0)+
#   geom_line(aes(y = betas2_cq_16_mpq,linetype ='dotdash',color ='SQR_h/16'),size=1.0)+
#   geom_line(aes(y = betas2_cq_8_mpq,linetype ='dotted',color ='SQR_h/8'),size=1.0)+
#   scale_color_discrete(labels = unname(TeX(c("$QR$","$SQR_{\\zeta^*}$","$SQR_{16\\zeta^*}$","$SQR_{2\\zeta^*}$",
#                                              "$SQR_{4\\zeta^{*}}$","$SQR_{8\\zeta^{*}}$","$SQR_{\\zeta^*/16}$",
#                                              "$SQR_{\\zeta^*/2}$","$SQR_{\\zeta^{*}/4}$","$SQR_{\\zeta^{*}/8}$")))) +
#   geom_hline(yintercept = 0.5, color = 'yellow', size=1.3) +
#   labs(x =xlab,y = ylab)  + 
#   scale_x_continuous(breaks=seq(0.01, 0.99, 0.07),limits = c(0.01, 0.99))  + 
#   theme_grey() + guides(linetype = FALSE) +
#   theme(plot.title = element_text(size = 15),
#         legend.title = element_blank(),
#         legend.text.align = 0,
#         legend.text = element_text(size = 12),
#         axis.title = element_text(size = 14))
# print(grafico_beta2_mpq)

#plot_grid(grafico_alpha0_mpq,grafico_alpha1_mpq,grafico_beta1_mpq,grafico_beta2_mpq)


#### PLOTS MSE

mse_coefs = data.frame(cbind(alphas0_rq_mse,alphas0_cq_mse,
                             alphas0_cq2_mse,alphas0_cq3_mse,alphas0_cq_3_mse,alphas0_cq_2_mse,
                             alphas0_cq8_mse,alphas0_cq16_mse,alphas0_cq_16_mse,alphas0_cq_8_mse))
ylab = TeX('MSE $\\hat{\\alpha_0}(\\tau)$')
grafico_alphas0_mse = ggplot(mse_coefs, aes(x = taus, y = value)) + 
  geom_line(aes(y = alphas0_rq_mse,linetype ='blank',color ='QR'),size=1.0) +
  geom_line(aes(y = alphas0_cq_mse,linetype ='solid',color ='SQR'),size=1.0)+
  geom_line(aes(y = alphas0_cq2_mse,linetype ='dotdash',color ='SQR_h*2'),size=1.0)+
  geom_line(aes(y = alphas0_cq3_mse,linetype ='dotdash',color ='SQR_h*4'),size=1.0)+
  geom_line(aes(y = alphas0_cq_3_mse,linetype ='dotdash',color ='SQR_h/4'),size=1.0)+
  geom_line(aes(y = alphas0_cq_2_mse,linetype ='dotted',color ='SQR_h/2'),size=1.0)+
  # geom_line(aes(y = alphas0_cq8_mse,linetype ='dotdash',color ='SQR_h*8'),size=1.0)+
  # geom_line(aes(y = alphas0_cq16_mse,linetype ='dotdash',color ='SQR_h*16'),size=1.0)+
  # geom_line(aes(y = alphas0_cq_16_mse,linetype ='dotdash',color ='SQR_h/16'),size=1.0)+
  # geom_line(aes(y = alphas0_cq_8_mse,linetype ='dotted',color ='SQR_h/8'),size=1.0)+
  scale_color_discrete(labels = unname(TeX(c("$QR$","$SQR_{\\zeta^*}$","$SQR_{2\\zeta^*}$","$SQR_{4\\zeta^*}$","$SQR_{\\zeta^{*}/4}$","$SQR_{\\zeta^{*}/2}$")))) +
  geom_hline(yintercept = 0.0, color = 'yellow', size=1.3) +
  labs(x =xlab,y = ylab)  + 
  scale_x_continuous(breaks=seq(0.01, 0.99, 0.14),limits = c(0.01, 0.99))  + 
  #scale_y_log10(labels = trans_format("log10", math_format(10^.x))) + 
  scale_y_log10() + # Log-scale
  theme_grey() + guides(linetype = FALSE) +
  theme(plot.title = element_text(size = 15),
        legend.title = element_blank(),
        legend.text.align = 0,
        legend.text = element_text(size = 12),
        axis.title = element_text(size = 14))
print(grafico_alphas0_mse)

mse_coefs = data.frame(cbind(alphas1_rq_mse,alphas1_cq_mse,
                             alphas1_cq2_mse,alphas1_cq3_mse,alphas1_cq_3_mse,alphas1_cq_2_mse,
                             alphas1_cq8_mse,alphas1_cq16_mse,alphas1_cq_16_mse,alphas1_cq_8_mse))
ylab = TeX('MSE $\\hat{\\alpha_1}(\\tau)$')
grafico_alphas1_mse = ggplot(mse_coefs, aes(x = taus, y = value)) + 
  geom_line(aes(y = alphas1_rq_mse,linetype ='blank',color ='QR'),size=1.0) +
  geom_line(aes(y = alphas1_cq_mse,linetype ='solid',color ='SQR'),size=1.0)+
  geom_line(aes(y = alphas1_cq2_mse,linetype ='dotdash',color ='SQR_h*2'),size=1.0)+
  geom_line(aes(y = alphas1_cq3_mse,linetype ='dotdash',color ='SQR_h*4'),size=1.0)+
  geom_line(aes(y = alphas1_cq_3_mse,linetype ='dotdash',color ='SQR_h/4'),size=1.0)+
  geom_line(aes(y = alphas1_cq_2_mse,linetype ='dotted',color ='SQR_h/2'),size=1.0)+
  # geom_line(aes(y = alphas1_cq8_mse,linetype ='dotdash',color ='SQR_h*8'),size=1.0)+
  # geom_line(aes(y = alphas1_cq16_mse,linetype ='dotdash',color ='SQR_h*16'),size=1.0)+
  # geom_line(aes(y = alphas1_cq_16_mse,linetype ='dotdash',color ='SQR_h/16'),size=1.0)+
  # geom_line(aes(y = alphas1_cq_8_mse,linetype ='dotted',color ='SQR_h/8'),size=1.0)+
  scale_color_discrete(labels = unname(TeX(c("$QR$","$SQR_{\\zeta^*}$","$SQR_{2\\zeta^*}$","$SQR_{4\\zeta^*}$","$SQR_{\\zeta^{*}/4}$","$SQR_{\\zeta^{*}/2}$")))) +
  geom_hline(yintercept = 0.0, color = 'yellow', size=1.3) +
  labs(x =xlab,y = ylab)  + 
  scale_x_continuous(breaks=seq(0.01, 0.99, 0.14),limits = c(0.01, 0.99))  + 
  scale_y_log10() + # Log-scale
  #scale_y_log10(labels = trans_format("log10", math_format(10^.x))) + 
  theme_grey() + guides(linetype = FALSE) +
  theme(plot.title = element_text(size = 15),
        legend.title = element_blank(),
        legend.text.align = 0,
        legend.text = element_text(size = 12),
        axis.title = element_text(size = 14))
print(grafico_alphas1_mse)


mse_coefs = data.frame(cbind(betas1_rq_mse,betas1_cq_mse,
                             betas1_cq2_mse,betas1_cq3_mse,betas1_cq_3_mse,betas1_cq_2_mse,
                             betas1_cq8_mse,betas1_cq16_mse,betas1_cq_16_mse,betas1_cq_8_mse))
ylab = TeX('MSE $\\hat{\\theta_1}(\\tau)$')
grafico_betas1_mse = ggplot(mse_coefs, aes(x = taus, y = value)) + 
  geom_line(aes(y = betas1_rq_mse,linetype ='blank',color ='QR'),size=1.0) +
  geom_line(aes(y = betas1_cq_mse,linetype ='solid',color ='SQR'),size=1.0)+
  geom_line(aes(y = betas1_cq2_mse,linetype ='dotdash',color ='SQR_h*2'),size=1.0)+
  geom_line(aes(y = betas1_cq3_mse,linetype ='dotdash',color ='SQR_h*4'),size=1.0)+
  geom_line(aes(y = betas1_cq_3_mse,linetype ='dotdash',color ='SQR_h/4'),size=1.0)+
  geom_line(aes(y = betas1_cq_2_mse,linetype ='dotted',color ='SQR_h/2'),size=1.0)+
  # geom_line(aes(y = betas1_cq8_mse,linetype ='dotdash',color ='SQR_h*8'),size=1.0)+
  # geom_line(aes(y = betas1_cq16_mse,linetype ='dotdash',color ='SQR_h*16'),size=1.0)+
  # geom_line(aes(y = betas1_cq_16_mse,linetype ='dotdash',color ='SQR_h/16'),size=1.0)+
  # geom_line(aes(y = betas1_cq_8_mse,linetype ='dotted',color ='SQR_h/8'),size=1.0)+
  scale_color_discrete(labels = unname(TeX(c("$QR$","$SQR_{\\zeta^*}$","$SQR_{2\\zeta^*}$","$SQR_{4\\zeta^*}$","$SQR_{\\zeta^{*}/4}$","$SQR_{\\zeta^{*}/2}$")))) +
  geom_hline(yintercept = 0.0, color = 'yellow', size=1.3) +
  labs(x =xlab,y = ylab)  + 
  scale_x_continuous(breaks=seq(0.01, 0.99, 0.14),limits = c(0.01, 0.99))  + 
  scale_y_log10() + # Log-scale
  #scale_y_log10(labels = trans_format("log10", math_format(10^.x))) + 
  theme_grey() + guides(linetype = FALSE) +
  theme(plot.title = element_text(size = 15),
        legend.title = element_blank(),
        legend.text.align = 0,
        legend.text = element_text(size = 12),
        axis.title = element_text(size = 14))
print(grafico_betas1_mse)


mse_coefs = data.frame(cbind(betas2_rq_mse,betas2_cq_mse,
                             betas2_cq2_mse,betas2_cq3_mse,betas2_cq_3_mse,betas2_cq_2_mse,
                             betas2_cq8_mse,betas2_cq16_mse,betas2_cq_16_mse,betas2_cq_8_mse))
ylab = TeX('MSE $\\hat{\\theta_2}(\\tau)$')
grafico_betas2_mse = ggplot(mse_coefs, aes(x = taus, y = value)) + 
  geom_line(aes(y = betas2_rq_mse,linetype ='blank',color ='QR'),size=1.0) +
  geom_line(aes(y = betas2_cq_mse,linetype ='solid',color ='SQR'),size=1.0)+
  geom_line(aes(y = betas2_cq2_mse,linetype ='dotdash',color ='SQR_h*2'),size=1.0)+
  geom_line(aes(y = betas2_cq3_mse,linetype ='dotdash',color ='SQR_h*4'),size=1.0)+
  geom_line(aes(y = betas2_cq_3_mse,linetype ='dotdash',color ='SQR_h/4'),size=1.0)+
  geom_line(aes(y = betas2_cq_2_mse,linetype ='dotted',color ='SQR_h/2'),size=1.0)+
  # geom_line(aes(y = betas2_cq8_mse,linetype ='dotdash',color ='SQR_h*8'),size=1.0)+
  # geom_line(aes(y = betas2_cq16_mse,linetype ='dotdash',color ='SQR_h*16'),size=1.0)+
  # geom_line(aes(y = betas2_cq_16_mse,linetype ='dotdash',color ='SQR_h/16'),size=1.0)+
  # geom_line(aes(y = betas2_cq_8_mse,linetype ='dotted',color ='SQR_h/8'),size=1.0)+
  scale_color_discrete(labels = unname(TeX(c("$QR$","$SQR_{\\zeta^*}$","$SQR_{2\\zeta^*}$","$SQR_{4\\zeta^*}$","$SQR_{\\zeta^{*}/4}$","$SQR_{\\zeta^{*}/2}$")))) +
  geom_hline(yintercept = 0.0, color = 'yellow', size=1.3) +
  labs(x =xlab,y = ylab)  + 
  scale_x_continuous(breaks=seq(0.01, 0.99, 0.14),limits = c(0.01, 0.99))  + 
  #scale_y_log10(labels = trans_format("log10", math_format(10^.x))) + 
  scale_y_log10() + # Log-scale
  theme_grey() + guides(linetype = FALSE) +
  theme(plot.title = element_text(size = 15),
        legend.title = element_blank(),
        legend.text.align = 0,
        legend.text = element_text(size = 12),
        axis.title = element_text(size = 14))
print(grafico_betas2_mse)
#plot_grid(grafico_alphas0_mse,grafico_alphas1_mse,grafico_betas1_mse,grafico_betas2_mse)

#### PLOTS BIAS

bias_coefs = data.frame(cbind(alphas0_rq_bias,alphas0_cq_bias,
                              alphas0_cq2_bias,alphas0_cq3_bias,alphas0_cq_3_bias,alphas0_cq_2_bias,
                              alphas0_cq8_bias,alphas0_cq16_bias,alphas0_cq_16_bias,alphas0_cq_8_bias))
ylab = TeX('Bias $\\hat{\\alpha_0}(\\tau)$')
grafico_alphas0_bias = ggplot(bias_coefs, aes(x = taus, y = value)) + 
  geom_line(aes(y = alphas0_rq_bias,linetype ='blank',color ='QR'),size=1.0) +
  geom_line(aes(y = alphas0_cq_bias,linetype ='solid',color ='SQR'),size=1.0)+
  geom_line(aes(y = alphas0_cq2_bias,linetype ='dotdash',color ='SQR_h*2'),size=1.0)+
  geom_line(aes(y = alphas0_cq3_bias,linetype ='dotdash',color ='SQR_h*4'),size=1.0)+
  geom_line(aes(y = alphas0_cq_3_bias,linetype ='dotdash',color ='SQR_h/4'),size=1.0)+
  geom_line(aes(y = alphas0_cq_2_bias,linetype ='dotted',color ='SQR_h/2'),size=1.0)+
  # geom_line(aes(y = alphas0_cq8_bias,linetype ='dotdash',color ='SQR_h*8'),size=1.0)+
  # geom_line(aes(y = alphas0_cq16_bias,linetype ='dotdash',color ='SQR_h*16'),size=1.0)+
  # geom_line(aes(y = alphas0_cq_16_bias,linetype ='dotdash',color ='SQR_h/16'),size=1.0)+
  # geom_line(aes(y = alphas0_cq_8_bias,linetype ='dotted',color ='SQR_h/8'),size=1.0)+
  scale_color_discrete(labels = unname(TeX(c("$QR$","$SQR_{\\zeta^*}$","$SQR_{2\\zeta^*}$","$SQR_{4\\zeta^*}$","$SQR_{\\zeta^{*}/4}$","$SQR_{\\zeta^{*}/2}$")))) +
  geom_hline(yintercept = 0.0, color = 'yellow', size=1.3) +
  labs(x =xlab,y = ylab)  + 
  scale_x_continuous(breaks=seq(0.01, 0.99, 0.14),limits = c(0.01, 0.99))  + 
  theme_grey() + guides(linetype = FALSE) +
  theme(plot.title = element_text(size = 15),
        legend.title = element_blank(),
        legend.text.align = 0,
        legend.text = element_text(size = 12),
        axis.title = element_text(size = 14))
print(grafico_alphas0_bias)

bias_coefs = data.frame(cbind(alphas1_rq_bias,alphas1_cq_bias,
                              alphas1_cq2_bias,alphas1_cq3_bias,alphas1_cq_3_bias,alphas1_cq_2_bias,
                              alphas1_cq8_bias,alphas1_cq16_bias,alphas1_cq_16_bias,alphas1_cq_8_bias))
ylab = TeX('Bias $\\hat{\\alpha_1}(\\tau)$')
grafico_alphas1_bias = ggplot(bias_coefs, aes(x = taus, y = value)) + 
  geom_line(aes(y = alphas1_rq_bias,linetype ='blank',color ='QR'),size=1.0) +
  geom_line(aes(y = alphas1_cq_bias,linetype ='solid',color ='SQR'),size=1.0)+
  geom_line(aes(y = alphas1_cq2_bias,linetype ='dotdash',color ='SQR_h*2'),size=1.0)+
  geom_line(aes(y = alphas1_cq3_bias,linetype ='dotdash',color ='SQR_h*4'),size=1.0)+
  geom_line(aes(y = alphas1_cq_3_bias,linetype ='dotdash',color ='SQR_h/4'),size=1.0)+
  geom_line(aes(y = alphas1_cq_2_bias,linetype ='dotted',color ='SQR_h/2'),size=1.0)+
  # geom_line(aes(y = alphas1_cq8_bias,linetype ='dotdash',color ='SQR_h*8'),size=1.0)+
  # geom_line(aes(y = alphas1_cq16_bias,linetype ='dotdash',color ='SQR_h*16'),size=1.0)+
  # geom_line(aes(y = alphas1_cq_16_bias,linetype ='dotdash',color ='SQR_h/16'),size=1.0)+
  # geom_line(aes(y = alphas1_cq_8_bias,linetype ='dotted',color ='SQR_h/8'),size=1.0)+
  scale_color_discrete(labels = unname(TeX(c("$QR$","$SQR_{\\zeta^*}$","$SQR_{2\\zeta^*}$","$SQR_{4\\zeta^*}$","$SQR_{\\zeta^{*}/4}$","$SQR_{\\zeta^{*}/2}$")))) +
  geom_hline(yintercept = 0.0, color = 'yellow', size=1.3) +
  labs(x =xlab,y = ylab)  + 
  scale_x_continuous(breaks=seq(0.01, 0.99, 0.14),limits = c(0.01, 0.99))  + 
  theme_grey() + guides(linetype = FALSE) +
  theme(plot.title = element_text(size = 15),
        legend.title = element_blank(),
        legend.text.align = 0,
        legend.text = element_text(size = 12),
        axis.title = element_text(size = 14))
print(grafico_alphas1_bias)


bias_coefs = data.frame(cbind(betas1_rq_bias,betas1_cq_bias,
                              betas1_cq2_bias,betas1_cq3_bias,betas1_cq_3_bias,betas1_cq_2_bias,
                              betas1_cq8_bias,betas1_cq16_bias,betas1_cq_16_bias,betas1_cq_8_bias))
ylab = TeX('Bias $\\hat{\\theta_1}(\\tau)$')
grafico_betas1_bias = ggplot(bias_coefs, aes(x = taus, y = value)) + 
  geom_line(aes(y = betas1_rq_bias,linetype ='blank',color ='QR'),size=1.0) +
  geom_line(aes(y = betas1_cq_bias,linetype ='solid',color ='SQR'),size=1.0)+
  geom_line(aes(y = betas1_cq2_bias,linetype ='dotdash',color ='SQR_h*2'),size=1.0)+
  geom_line(aes(y = betas1_cq3_bias,linetype ='dotdash',color ='SQR_h*4'),size=1.0)+
  geom_line(aes(y = betas1_cq_3_bias,linetype ='dotdash',color ='SQR_h/4'),size=1.0)+
  geom_line(aes(y = betas1_cq_2_bias,linetype ='dotted',color ='SQR_h/2'),size=1.0)+
  # geom_line(aes(y = betas1_cq8_bias,linetype ='dotdash',color ='SQR_h*8'),size=1.0)+
  # geom_line(aes(y = betas1_cq16_bias,linetype ='dotdash',color ='SQR_h*16'),size=1.0)+
  # geom_line(aes(y = betas1_cq_16_bias,linetype ='dotdash',color ='SQR_h/16'),size=1.0)+
  # geom_line(aes(y = betas1_cq_8_bias,linetype ='dotted',color ='SQR_h/8'),size=1.0)+
  scale_color_discrete(labels = unname(TeX(c("$QR$","$SQR_{\\zeta^*}$","$SQR_{2\\zeta^*}$","$SQR_{4\\zeta^*}$","$SQR_{\\zeta^{*}/4}$","$SQR_{\\zeta^{*}/2}$")))) +
  geom_hline(yintercept = 0.0, color = 'yellow', size=1.3) +
  labs(x =xlab,y = ylab)  + 
  scale_x_continuous(breaks=seq(0.01, 0.99, 0.14),limits = c(0.01, 0.99))  + 
  theme_grey() + guides(linetype = FALSE) +
  theme(plot.title = element_text(size = 15),
        legend.title = element_blank(),
        legend.text.align = 0,
        legend.text = element_text(size = 12),
        axis.title = element_text(size = 14))
print(grafico_betas1_bias)


bias_coefs = data.frame(cbind(betas2_rq_bias,betas2_cq_bias,
                              betas2_cq2_bias,betas2_cq3_bias,betas2_cq_3_bias,betas2_cq_2_bias,
                              betas2_cq8_bias,betas2_cq16_bias,betas2_cq_16_bias,betas2_cq_8_bias))
ylab = TeX('Bias $\\hat{\\theta_2}(\\tau)$')
grafico_betas2_bias = ggplot(bias_coefs, aes(x = taus, y = value)) + 
  geom_line(aes(y = betas2_rq_bias,linetype ='blank',color ='QR'),size=1.0) +
  geom_line(aes(y = betas2_cq_bias,linetype ='solid',color ='SQR'),size=1.0)+
  geom_line(aes(y = betas2_cq2_bias,linetype ='dotdash',color ='SQR_h*2'),size=1.0)+
  geom_line(aes(y = betas2_cq3_bias,linetype ='dotdash',color ='SQR_h*4'),size=1.0)+
  geom_line(aes(y = betas2_cq_3_bias,linetype ='dotdash',color ='SQR_h/4'),size=1.0)+
  geom_line(aes(y = betas2_cq_2_bias,linetype ='dotted',color ='SQR_h/2'),size=1.0)+
  # geom_line(aes(y = betas2_cq8_bias,linetype ='dotdash',color ='SQR_h*8'),size=1.0)+
  # geom_line(aes(y = betas2_cq16_bias,linetype ='dotdash',color ='SQR_h*16'),size=1.0)+
  # geom_line(aes(y = betas2_cq_16_bias,linetype ='dotdash',color ='SQR_h/16'),size=1.0)+
  # geom_line(aes(y = betas2_cq_8_bias,linetype ='dotted',color ='SQR_h/8'),size=1.0)+
  scale_color_discrete(labels = unname(TeX(c("$QR$","$SQR_{\\zeta^*}$","$SQR_{2\\zeta^*}$","$SQR_{4\\zeta^*}$","$SQR_{\\zeta^{*}/4}$","$SQR_{\\zeta^{*}/2}$")))) +
  geom_hline(yintercept = 0.0, color = 'yellow', size=1.3) +
  labs(x =xlab,y = ylab)  + 
  scale_x_continuous(breaks=seq(0.01, 0.99, 0.14),limits = c(0.01, 0.99))  + 
  theme_grey() + guides(linetype = FALSE) +
  theme(plot.title = element_text(size = 15),
        legend.title = element_blank(),
        legend.text.align = 0,
        legend.text = element_text(size = 12),
        axis.title = element_text(size = 14))
print(grafico_betas2_bias)
#plot_grid(grafico_alphas0_bias,grafico_alphas1_bias,grafico_betas1_bias,grafico_betas2_bias)

### PLOTS vpq  - VARIANCE per quantile - - h, h*2, h*4, h*1/4, h*1/2

coefs_intercept_vpq = data.frame(intercept_rq_vpq, intercept_cq_vpq, 
                                 intercept_cq2_vpq,intercept_cq3_vpq,intercept_cq_3_vpq,intercept_cq_2_vpq,
                                 intercept_cq8_vpq,intercept_cq16_vpq,intercept_cq_16_vpq,intercept_cq_8_vpq) 
coefs_alphas1_vpq = data.frame(alphas1_rq_vpq, alphas1_cq_vpq, 
                               alphas1_cq2_vpq,alphas1_cq3_vpq,alphas1_cq_3_vpq,alphas1_cq_2_vpq,
                               alphas1_cq8_vpq,alphas1_cq16_vpq,alphas1_cq_16_vpq,alphas1_cq_8_vpq) 
coefs_betas1_vpq = data.frame(betas1_rq_vpq, betas1_cq_vpq,
                              betas1_cq2_vpq,betas1_cq3_vpq,betas1_cq_3_vpq,betas1_cq_2_vpq,
                              betas1_cq8_vpq,betas1_cq16_vpq,betas1_cq_16_vpq,betas1_cq_8_vpq)
coefs_betas2_vpq = data.frame(betas2_rq_vpq, betas2_cq_vpq,
                              betas2_cq2_vpq,betas2_cq3_vpq,betas2_cq_3_vpq,betas2_cq_2_vpq,
                              betas2_cq8_vpq,betas2_cq16_vpq,betas2_cq_16_vpq,betas2_cq_8_vpq)

xlab = TeX('$\\tau$')

ylab = TeX('Var $\\hat{\\alpha_0}(\\tau)$') 
grafico_alpha0_vpq = ggplot(coefs_intercept_vpq, aes(x = taus, y = value)) + 
  geom_line(aes(y = intercept_rq_vpq,linetype ='blank',color ='QR'),size=1.0) +
  geom_line(aes(y = intercept_cq_vpq, linetype ='solid',color ='SQR'),size=1.0)+
  geom_line(aes(y = intercept_cq2_vpq, linetype ='dotdash',color ='SQR_h*2'),size=1.0)+
  geom_line(aes(y = intercept_cq3_vpq, linetype ='dotdash',color ='SQR_h*4'),size=1.0)+
  geom_line(aes(y = intercept_cq_3_vpq, linetype ='dotdash',color ='SQR_h/4'),size=1.0)+
  geom_line(aes(y = intercept_cq_2_vpq, linetype ='dotted',color ='SQR_h/2'),size=1.0)+
  # geom_line(aes(y = intercept_cq8_vpq, linetype ='dotdash',color ='SQR_h*8'),size=1.0)+
  # geom_line(aes(y = intercept_cq16_vpq, linetype ='dotdash',color ='SQR_h*16'),size=1.0)+
  # geom_line(aes(y = intercept_cq_16_vpq, linetype ='dotdash',color ='SQR_h/16'),size=1.0)+
  # geom_line(aes(y = intercept_cq_8_vpq, linetype ='dotted',color ='SQR_h/8'),size=1.0)+
  scale_color_discrete(labels = unname(TeX(c("$QR$","$SQR_{\\zeta^*}$","$SQR_{2\\zeta^*}$","$SQR_{4\\zeta^*}$","$SQR_{\\zeta^{*}/4}$","$SQR_{\\zeta^{*}/2}$")))) +
  geom_hline(yintercept = 0, color = 'yellow', size=1.3) +
  labs(x =xlab,y = ylab) + #expand_limits(x=c(0, 1), y=c(-3,3)) +
  scale_x_continuous(breaks=seq(0.01, 0.99, 0.14),limits = c(0.01, 0.99))  + 
  scale_y_log10() + 
  theme_grey() + guides(linetype = FALSE) +
  theme(plot.title = element_text(size = 15),
        legend.title = element_blank(),
        legend.text.align = 0,
        legend.text = element_text(size = 12),
        axis.title = element_text(size = 14))
print(grafico_alpha0_vpq)

ylab = TeX('Var $\\hat{\\alpha_1}(\\tau)$') 
grafico_alpha1_vpq = ggplot(coefs_alphas1_vpq, aes(x = taus, y = value)) + 
  geom_line(aes(y = alphas1_rq_vpq,linetype ='blank',color ='QR'),size=1.0) +
  geom_line(aes(y = alphas1_cq_vpq, linetype ='solid',color ='SQR'),size=1.0)+
  geom_line(aes(y = alphas1_cq2_vpq, linetype ='dotdash',color ='SQR_h*2'),size=1.0)+
  geom_line(aes(y = alphas1_cq3_vpq, linetype ='dotdash',color ='SQR_h*4'),size=1.0)+
  geom_line(aes(y = alphas1_cq_3_vpq, linetype ='dotdash',color ='SQR_h/4'),size=1.0)+
  geom_line(aes(y = alphas1_cq_2_vpq, linetype ='dotted',color ='SQR_h/2'),size=1.0)+
  # geom_line(aes(y = alphas1_cq8_vpq, linetype ='dotdash',color ='SQR_h*8'),size=1.0)+
  # geom_line(aes(y = alphas1_cq16_vpq, linetype ='dotdash',color ='SQR_h*16'),size=1.0)+
  # geom_line(aes(y = alphas1_cq_16_vpq, linetype ='dotdash',color ='SQR_h/16'),size=1.0)+
  # geom_line(aes(y = alphas1_cq_8_vpq, linetype ='dotted',color ='SQR_h/8'),size=1.0)+
  scale_color_discrete(labels = unname(TeX(c("$QR$","$SQR_{\\zeta^*}$","$SQR_{2\\zeta^*}$","$SQR_{4\\zeta^*}$","$SQR_{\\zeta^{*}/4}$","$SQR_{\\zeta^{*}/2}$")))) +
  geom_hline(yintercept = 0, color = 'yellow', size=1.3) +
  labs(x =xlab,y = ylab) + #expand_limits(x=c(0, 1), y=c(0.4,0.6)) + 
  scale_x_continuous(breaks=round(seq(0.01, 0.99, 0.14), 2),limits = c(0.01, 0.99))  + 
  #scale_y_continuous(breaks=seq(0.4,0.6,0.02),limits = c(0.46, 0.54))  +
  scale_y_log10() + 
  theme_grey() + guides(linetype = FALSE) +
  theme(plot.title = element_text(size = 15),
        legend.title = element_blank(),
        legend.text.align = 0,
        legend.text = element_text(size = 12),
        axis.title = element_text(size = 14))
print(grafico_alpha1_vpq)

ylab = TeX('Var $\\hat{\\theta_1}(\\tau)$')
grafico_beta1_vpq = ggplot(coefs_betas1_vpq, aes(x = taus, y = value)) + 
  geom_line(aes(y = betas1_rq_vpq,linetype ='blank',color ='QR'),size=1.0) +
  geom_line(aes(y = betas1_cq_vpq, linetype ='solid',color ='SQR'),size=1.,)+
  geom_line(aes(y = betas1_cq2_vpq, linetype ='dotdash',color ='SQR_h*2'),size=1.0)+
  geom_line(aes(y = betas1_cq3_vpq, linetype ='dotdash',color ='SQR_h*4'),size=1.0)+
  geom_line(aes(y = betas1_cq_3_vpq, linetype ='dotdash',color ='SQR_h/4'),size=1.0)+
  geom_line(aes(y = betas1_cq_2_vpq, linetype ='dotted',color ='SQR_h/2'),size=1.0)+
  # geom_line(aes(y = betas1_cq8_vpq, linetype ='dotdash',color ='SQR_h*8'),size=1.0)+
  # geom_line(aes(y = betas1_cq16_vpq, linetype ='dotdash',color ='SQR_h*16'),size=1.0)+
  # geom_line(aes(y = betas1_cq_16_vpq, linetype ='dotdash',color ='SQR_h/16'),size=1.0)+
  # geom_line(aes(y = betas1_cq_8_vpq, linetype ='dotted',color ='SQR_h/8'),size=1.0)+
  scale_color_discrete(labels = unname(TeX(c("$QR$","$SQR_{\\zeta^*}$","$SQR_{2\\zeta^*}$","$SQR_{4\\zeta^*}$","$SQR_{\\zeta^{*}/4}$","$SQR_{\\zeta^{*}/2}$")))) +
  geom_hline(yintercept = 0, color = 'yellow', size=1.3) +
  labs(x =xlab,y = ylab) + # ggtitle(bquote(beta[1][tau]^{tau})) + 
  scale_x_continuous(breaks=seq(0.01, 0.99, 0.14),limits = c(0.01, 0.99))  + 
  scale_y_log10() + 
  theme_grey() + guides(linetype = FALSE) +
  theme(plot.title = element_text(size = 15),
        legend.title = element_blank(),
        legend.text.align = 0,
        legend.text = element_text(size = 12),
        axis.title = element_text(size = 14))
print(grafico_beta1_vpq)

ylab = TeX('Var $\\hat{\\theta_2}(\\tau)$')
grafico_beta2_vpq = ggplot(coefs_betas2_vpq, aes(x = taus, y = value)) + 
  geom_line(aes(y = betas2_rq_vpq,linetype ='blank',color ='QR'),size=1.0) +
  geom_line(aes(y = betas2_cq_vpq,linetype ='solid',color ='SQR'),size=1.0)+
  geom_line(aes(y = betas2_cq2_vpq,linetype ='dotdash',color ='SQR_h*2'),size=1.0)+
  geom_line(aes(y = betas2_cq3_vpq,linetype ='dotdash',color ='SQR_h*4'),size=1.0)+
  geom_line(aes(y = betas2_cq_3_vpq,linetype ='dotdash',color ='SQR_h/4'),size=1.0)+
  geom_line(aes(y = betas2_cq_2_vpq,linetype ='dotted',color ='SQR_h/2'),size=1.0)+
  # geom_line(aes(y = betas2_cq8_vpq,linetype ='dotdash',color ='SQR_h*8'),size=1.0)+
  # geom_line(aes(y = betas2_cq16_vpq,linetype ='dotdash',color ='SQR_h*16'),size=1.0)+
  # geom_line(aes(y = betas2_cq_16_vpq,linetype ='dotdash',color ='SQR_h/16'),size=1.0)+
  # geom_line(aes(y = betas2_cq_8_vpq,linetype ='dotted',color ='SQR_h/8'),size=1.0)+
  scale_color_discrete(labels = unname(TeX(c("$QR$","$SQR_{\\zeta^*}$","$SQR_{2\\zeta^*}$","$SQR_{4\\zeta^*}$","$SQR_{\\zeta^{*}/4}$","$SQR_{\\zeta^{*}/2}$")))) +
  geom_hline(yintercept = 0, color = 'yellow', size=1.3) +
  labs(x =xlab,y = ylab)  + 
  scale_x_continuous(breaks=seq(0.01, 0.99, 0.14),limits = c(0.01, 0.99))  + 
  scale_y_log10() + 
  theme_grey() + guides(linetype = FALSE) +
  theme(plot.title = element_text(size = 15),
        legend.title = element_blank(),
        legend.text.align = 0,
        legend.text = element_text(size = 12),
        axis.title = element_text(size = 14))
print(grafico_beta2_vpq)

#plot_grid(grafico_alpha0_vpq,grafico_alpha1_vpq,grafico_beta1_vpq,grafico_beta2_vpq)

options(scipen=999)
plot_grid(grafico_alphas0_mse,grafico_alphas0_bias,grafico_alpha0_vpq,
          grafico_alphas1_mse,grafico_alphas1_bias,grafico_alpha1_vpq,
          grafico_betas1_mse,grafico_betas1_bias,grafico_beta1_vpq,
          grafico_betas2_mse,grafico_betas2_bias, grafico_beta2_vpq,
          nrow = 4, ncol = 3, label_size = 8, align = 'v')

real_alpha0 = c()
#real_alpha1 = 0
real_alpha1 = c()
real_theta1 = c() 
real_theta2 = c() 
for(s in 1:length(taus)){
  
  real_alpha0[s] = qnorm(taus[s],0,1)
  real_alpha1[s] = 0.5
  real_theta1[s] = 0.5
  real_theta2[s] = 0.5
  
}


### SAVING MAIN RESULTS

# GLOBAL METRICS
alphas0_rq_global = c()
alphas1_rq_global = c()
thetas1_rq_global = c()
alphas0_cq_global = c()
alphas1_cq_global = c()
thetas1_cq_global = c()

alphas0_cq2_global = c()
alphas1_cq2_global = c()
thetas1_cq2_global = c()

alphas0_cq3_global = c()
alphas1_cq3_global = c()
thetas1_cq3_global = c()

alphas0_cq_3_global = c()
alphas1_cq_3_global = c()
thetas1_cq_3_global = c()

alphas0_cq_2_global = c()
alphas1_cq_2_global = c()
thetas1_cq_2_global = c()

alphas0_cq8_global = c()
alphas1_cq8_global = c()
thetas1_cq8_global = c()

alphas0_cq16_global = c()
alphas1_cq16_global = c()
thetas1_cq16_global = c()

alphas0_cq_16_global = c()
alphas1_cq_16_global = c()
thetas1_cq_16_global = c()

alphas0_cq_8_global = c()
alphas1_cq_8_global = c()
thetas1_cq_8_global = c()

thetas2_rq_global = c()
thetas2_cq_global = c()
thetas2_cq2_global = c()
thetas2_cq3_global = c()
thetas2_cq_3_global = c()
thetas2_cq_2_global = c()
thetas2_cq8_global = c()
thetas2_cq16_global = c()
thetas2_cq_8_global = c()
thetas2_cq_16_global = c()

for(j in 1:nrep){
  alphas0_rq_global[j] = sum((abs(A[,1,j] - real_alpha0))^2)*0.02
  alphas1_rq_global[j] = sum((abs(A[,2,j] - real_alpha1))^2)*0.02
  thetas1_rq_global[j] = sum((abs(A[,3,j] - real_theta1))^2)*0.02
  thetas2_rq_global[j] = sum((abs(A[,4,j] - real_theta2))^2)*0.02
  
  alphas0_cq_global[j] = sum((abs(C[,1,j] - real_alpha0))^2)*0.02
  alphas1_cq_global[j] = sum((abs(C[,2,j] - real_alpha1))^2)*0.02
  thetas1_cq_global[j] = sum((abs(C[,3,j] - real_theta1))^2)*0.02
  thetas2_cq_global[j] = sum((abs(C[,4,j] - real_theta2))^2)*0.02
  
  alphas0_cq2_global[j] = sum((abs(C2[,1,j] - real_alpha0))^2)*0.02
  alphas1_cq2_global[j] = sum((abs(C2[,2,j] - real_alpha1))^2)*0.02
  thetas1_cq2_global[j] = sum((abs(C2[,3,j] - real_theta1))^2)*0.02
  thetas2_cq2_global[j] = sum((abs(C2[,4,j] - real_theta2))^2)*0.02
  
  alphas0_cq3_global[j] = sum((abs(C3[,1,j] - real_alpha0))^2)*0.02
  alphas1_cq3_global[j] = sum((abs(C3[,2,j] - real_alpha1))^2)*0.02
  thetas1_cq3_global[j] = sum((abs(C3[,3,j] - real_theta1))^2)*0.02
  thetas2_cq3_global[j] = sum((abs(C3[,4,j] - real_theta2))^2)*0.02
  
  alphas0_cq_3_global[j] = sum((abs(C_3[,1,j] - real_alpha0))^2)*0.02
  alphas1_cq_3_global[j] = sum((abs(C_3[,2,j] - real_alpha1))^2)*0.02
  thetas1_cq_3_global[j] = sum((abs(C_3[,3,j] - real_theta1))^2)*0.02
  thetas2_cq_3_global[j] = sum((abs(C_3[,4,j] - real_theta2))^2)*0.02
  
  alphas0_cq_2_global[j] = sum((abs(C_2[,1,j] - real_alpha0))^2)*0.02
  alphas1_cq_2_global[j] = sum((abs(C_2[,2,j] - real_alpha1))^2)*0.02
  thetas1_cq_2_global[j] = sum((abs(C_2[,3,j] - real_theta1))^2)*0.02
  thetas2_cq_2_global[j] = sum((abs(C_2[,4,j] - real_theta2))^2)*0.02
  
  alphas0_cq8_global[j] = sum((abs(C8[,1,j] - real_alpha0))^2)*0.02
  alphas1_cq8_global[j] = sum((abs(C8[,2,j] - real_alpha1))^2)*0.02
  thetas1_cq8_global[j] = sum((abs(C8[,3,j] - real_theta1))^2)*0.02
  thetas2_cq8_global[j] = sum((abs(C8[,4,j] - real_theta2))^2)*0.02
  
  alphas0_cq16_global[j] = sum((abs(C16[,1,j] - real_alpha0))^2)*0.02
  alphas1_cq16_global[j] = sum((abs(C16[,2,j] - real_alpha1))^2)*0.02
  thetas1_cq16_global[j] = sum((abs(C16[,3,j] - real_theta1))^2)*0.02
  thetas2_cq16_global[j] = sum((abs(C16[,4,j] - real_theta2))^2)*0.02
  
  alphas0_cq_16_global[j] = sum((abs(C_16[,1,j] - real_alpha0))^2)*0.02
  alphas1_cq_16_global[j] = sum((abs(C_16[,2,j] - real_alpha1))^2)*0.02
  thetas1_cq_16_global[j] = sum((abs(C_16[,3,j] - real_theta1))^2)*0.02
  thetas2_cq_16_global[j] = sum((abs(C_16[,4,j] - real_theta2))^2)*0.02
  
  alphas0_cq_8_global[j] = sum((abs(C_8[,1,j] - real_alpha0))^2)*0.02
  alphas1_cq_8_global[j] = sum((abs(C_8[,2,j] - real_alpha1))^2)*0.02
  thetas1_cq_8_global[j] = sum((abs(C_8[,3,j] - real_theta1))^2)*0.02
  thetas2_cq_8_global[j] = sum((abs(C_8[,4,j] - real_theta2))^2)*0.02
}

global_rq = mean(alphas0_rq_global + alphas1_rq_global + thetas1_rq_global + thetas2_rq_global)
global_cq = mean(alphas0_cq_global + alphas1_cq_global + thetas1_cq_global + thetas2_cq_global)
global_cq2 = mean(alphas0_cq2_global + alphas1_cq2_global + thetas1_cq2_global + thetas2_cq2_global)
global_cq3 = mean(alphas0_cq3_global + alphas1_cq3_global + thetas1_cq3_global + thetas2_cq3_global)
global_cq_3 = mean(alphas0_cq_3_global + alphas1_cq_3_global + thetas1_cq_3_global + thetas2_cq_3_global)
global_cq_2 = mean(alphas0_cq_2_global + alphas1_cq_2_global + thetas1_cq_2_global + thetas2_cq_2_global)
global_cq8 = mean(alphas0_cq8_global + alphas1_cq8_global + thetas1_cq8_global + thetas2_cq8_global)
global_cq16 = mean(alphas0_cq16_global + alphas1_cq16_global + thetas1_cq16_global + thetas2_cq16_global)
global_cq_16 = mean(alphas0_cq_16_global + alphas1_cq_16_global + thetas1_cq_16_global + thetas2_cq_16_global)
global_cq_8 = mean(alphas0_cq_8_global + alphas1_cq_8_global + thetas1_cq_8_global + thetas2_cq_8_global)

globals_df = data.frame(cbind(global_rq,global_cq,
                              global_cq2,global_cq3,global_cq_3,global_cq_2,
                              global_cq8,global_cq16,global_cq_16,global_cq_8)) 

# FULL TABLES -  Bias & mse - h/2

table_all_metrics_h2  = cbind(intercept_rq_bias, alphas1_rq_bias, betas1_rq_bias, betas2_rq_bias, 
                              intercept_cq_2_bias, alphas1_cq_2_bias, betas1_cq_2_bias, betas2_cq_2_bias, 
                              intercept_rq_mse, alphas1_rq_mse, betas1_rq_mse, betas2_rq_mse,
                              intercept_cq_2_mse, alphas1_cq_2_mse, betas1_cq_2_mse, betas2_cq_2_mse)

table_all_metrics_h2 = round(table_all_metrics_h2, 4)

colnames(table_all_metrics_h2) <- c('alpha_0_rq','alpha_1_rq', 'beta_1_rq', 'beta_2_rq', 'alpha_0_cq','alpha_1_cq', 'beta_1_cq', 'beta_2_cq', 
                                    'alpha_0_rq','alpha_1_rq', 'beta_1_rq', 'beta_2_rq', 'alpha_0_cq','alpha_1_cq', 'beta_1_cq', 'beta_2_cq')

rownames(table_all_metrics_h2) <- as.character(taus)

stargazer(table_all_metrics_h2, title = 'Location-shift model x: bias and mse of estimators - bandwidth = h/2', digits=4, digits.extra = 4, decimal.mark = '.')