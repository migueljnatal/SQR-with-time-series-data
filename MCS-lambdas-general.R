setwd("C:/Users/MIGUEL/OneDrive - UFRGS/Mestrado PPGEst/DISSERTAÇÃO/Novas Simulações")

start_time = Sys.time()
#install.packages('Rfast')
#install.packages('urca')
#install.packages('tseries')
#install.packages("Metrics")
#library(repr) 
#library(extraDistr)
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
library(Rfast)
library(urca)
library(tseries)
library(Metrics)
source('choose_lambda.R')

# GLOBAL VARIABLES - USER MUST CHOOSE THESE PARAMETERS
nrep = 5000                                            # number of replications
T = 1101                                               # sample size
choice = 'hacovercos'                                  # choice of lambda function


# Quantile levels
tau.grid = seq(from=.01,to=.99, by=.02) 

# Parameters for the kumar quantile functions v10, v01 and v00
{a10 = 1; b10 = 3; a01 = 5; b01 = 1; a00 = 3; b00 = 1}

# Conditional quantile function of Y[t] given Y[t-1] = 1 and Z[t-1]=0
v10 = function (tau) qbeta(tau, a10,b10)

# Conditional quantile function of Y[t] given Y[t-1] = 0 and Z[t-1]=1
v01 = function(tau) qbeta(tau,a01,b01)

# Conditional quantile function of Y[t] given Y[t-1] = 0 and Z[t-1]=0 -- LAMBDAS <=> v00(\tau)

lambda_dv_qbeta = function(tau,a10,b10,a01,b01){
  
  lambda_tau = lambda(tau, choice=choice)
  d = D(lambda_tau, 'tau')
  
  Q_betas = (qbeta(tau,a10,b10) + qbeta(tau,a01,b01))
  
  dv_lambda = eval(d)
  
  return(dv_lambda*Q_betas)
  
}

v00 = function(tau) {
  
  lambda_tau = lambda(tau, choice=choice)
  
  Q_betas = (qbeta(tau,a10,b10) + qbeta(tau,a01,b01))
  
  itg = integrate(lambda_dv_qbeta, a10, b10, a01, b01, lower=0, upper=tau)$value
  
  full_itg = eval(lambda_tau)*Q_betas - itg
  
  return(full_itg)
} 

v00_init = v00

v00_u = function(u) v00_init(u)/v00_init(1)

v00_t = sapply(tau.grid,v00_u)


# Conditional quantile function of Y[t] given Y[t-1] = 1 and Z[t-1]=1
v11 = function(tau) v10(tau) + v01(tau) - v00_u(tau)

v11_t = sapply(tau.grid,v11)

# Functional parameters for the quantile regression equation
alpha0 = v00_u
alpha1 = function(tau) v10(tau) - v00_u(tau)
theta1 = function(tau) v01(tau) - v00_u(tau)

# Quantile function of Y given â„±[t-1]
Q = function(tau,Y.current,Z.current){
  alpha0(tau) + alpha1(tau)*Y.current + theta1(tau)*Z.current
}

# Simulating the sample paths:

A = array(0, dim = c(length(tau.grid), 3, nrep))           # arrays for storing the simulated coef; A stands for rq and C for conquer 
B = array(0, dim = c(length(tau.grid), 3, nrep))
C = array(0, dim = c(length(tau.grid), 3, nrep))
C2 = array(0, dim = c(length(tau.grid), 3, nrep))
C3 = array(0, dim = c(length(tau.grid), 3, nrep))
C_3 = array(0, dim = c(length(tau.grid), 3, nrep))
C_2 = array(0, dim = c(length(tau.grid), 3, nrep))
C8 = array(0, dim = c(length(tau.grid), 3, nrep))
C16 = array(0, dim = c(length(tau.grid), 3, nrep))
C_16 = array(0, dim = c(length(tau.grid), 3, nrep))
C_8 = array(0, dim = c(length(tau.grid), 3, nrep))
h_scale = 1
z = matrix(0, length(tau.grid), nrep)
iqr_z =  matrix(0, length(tau.grid), nrep) 
sigma_z =  matrix(0, length(tau.grid), nrep) 
h =  matrix(0, length(tau.grid), nrep) 
h2 =  matrix(0, length(tau.grid), nrep) 
h3 =  matrix(0, length(tau.grid), nrep) 
h_2 =  matrix(0, length(tau.grid), nrep) 
h_3 =  matrix(0, length(tau.grid), nrep) 
h8 =  matrix(0, length(tau.grid), nrep) 
# h16 =  matrix(0, length(tau.grid), nrep) 
# h_8 =  matrix(0, length(tau.grid), nrep) 
# h_16 =  matrix(0, length(tau.grid), nrep)

Y0 = list()
Z0 = list()
Y = Z = list()
unif_j = list()
unif_j2 = list()
Z.current = list()
Y.current = list()

for(j in 1:nrep){

# Arbitrary starting point (Y0,Z0)
 Y0[[j]] = runif(1)
 Z0[[j]] = runif(1)
 Y[[j]] = runif(1)
 Z[[j]] = runif(1)
# Y = Z = list()#numeric()
 Z.current[[j]] = Z0[[j]]
 Y.current[[j]] = Y0[[j]]

 
 unif_j[[j]] = runif(1)
 unif_j2[[j]] = runif(1)
 
 for(t in 1:T){

  Y[[j]][t] = Q(runif(1), Y.current[[j]], Z.current[[j]])
  
  Z[[j]][t] = Q(runif(1), Z.current[[j]], Y.current[[j]])
  
  Z.current[[j]] = Z[[j]][t]
  Y.current[[j]] = Y[[j]][t]

 }
  burn = 100                     # discarding the 1st 100 observations
  Y[[j]] = Y[[j]][(burn+1):T]
  Z[[j]] = Z[[j]][(burn+1):T]
  Y[[j]] = Y[[j]][-1]
  Z[[j]] = Z[[j]][-1] 
  
}


T = length(Y[[j]]) 

### Coefficients trial

Yvec = list()
Xmat = list()

for(j in 1:nrep){
  
 Yvec[[j]] = Y[[j]][2:T]
 Xmat[[j]] = cbind(Y[[j]][1:(T-1)], Z[[j]][1:(T-1)])
 
 for(i in 1:length(tau.grid)){
  
  qrfit = rq(Yvec[[j]]~Xmat[[j]], tau = tau.grid[i])
  A[i,,j] = qrfit$coefficients
  
  z = Yvec[[j]] - cbind(rep(1,length(Yvec[[j]])),as.matrix(Xmat[[j]]))%*%A[i,,j]#B[i,,j]
  iqr_z = quantile(z, .75) - quantile(z, .25)
  sigma_z = min(sd(z), ((iqr_z)/1.34898))
  h[i,j] = h_scale*(1.06*sigma_z)/(length(Yvec[[j]])^(1/5))
  
  # Demais niveis de $\zeta$ = h 'rule of thumb' como no paper do Artur
  h_3[i,j]        = h[i,j]/4
  h_2[i,j]        = h[i,j]/2
  h2[i,j]         = h[i,j]*2
  h3[i,j]         = h[i,j]*4
  h8[i,j]         = h[i,j]*8
  # h_8[i,j]        = h[i,j]/8
  # h_16[i,j]        = h[i,j]/16
  # h16[i,j]         = h[i,j]*16
  
  qadl_conquer_2 <- conquer(as.matrix(Xmat[[j]]), Yvec[[j]], tau=tau.grid[i], h = h[i,j])   # 2nd conquer estimation (optimal bandwidth = h)
  C[i,,j] = qadl_conquer_2$coeff
  
  qadl_conquer_2.1 <- conquer(as.matrix(Xmat[[j]]), Yvec[[j]], tau=tau.grid[i], h = h2[i,j])   # 2nd conquer estimation (optimal bandwidth = h*2)
  C2[i,,j] = qadl_conquer_2.1$coeff
  
  qadl_conquer_2.2 <- conquer(as.matrix(Xmat[[j]]), Yvec[[j]], tau=tau.grid[i], h = h3[i,j])   # 2nd conquer estimation (optimal bandwidth = h*4)
  C3[i,,j] = qadl_conquer_2.2$coeff
  
  qadl_conquer_2.3 <- conquer(as.matrix(Xmat[[j]]), Yvec[[j]], tau=tau.grid[i], h = h_3[i,j])   # 2nd conquer estimation (optimal bandwidth = h/4)
  C_3[i,,j] = qadl_conquer_2.3$coeff
  
  qadl_conquer_2.4 <- conquer(as.matrix(Xmat[[j]]), Yvec[[j]], tau=tau.grid[i], h = h_2[i,j])   # 2nd conquer estimation (optimal bandwidth = h/2)
  C_2[i,,j] = qadl_conquer_2.4$coeff
  
  qadl_conquer_2.5 <- conquer(as.matrix(Xmat[[j]]), Yvec[[j]], tau=tau.grid[i], h = h8[i,j])   # 2nd conquer estimation (optimal bandwidth = h*8)
  C8[i,,j] = qadl_conquer_2.5$coeff
  
  # qadl_conquer_2.6 <- conquer(as.matrix(Xmat[[j]]), Yvec[[j]], tau=tau.grid[i], h = h_8[i,j])   # 2nd conquer estimation (optimal bandwidth = h*16)
  # C_8[i,,j] = qadl_conquer_2.6$coeff
  # 
  # qadl_conquer_2.7 <- conquer(as.matrix(Xmat[[j]]), Yvec[[j]], tau=tau.grid[i], h = h_16[i,j])   # 2nd conquer estimation (optimal bandwidth = h/16)
  # C_16[i,,j] = qadl_conquer_2.7$coeff
  # 
  # qadl_conquer_2.8 <- conquer(as.matrix(Xmat[[j]]), Yvec[[j]], tau=tau.grid[i], h = h16[i,j])   # 2nd conquer estimation (optimal bandwidth = h/8)
  # C16[i,,j] = qadl_conquer_2.8$coeff
  
 }
}

# Functional parameters estimation 

# BIAS
alphas0_rq_bias = c()
alphas0_cq_bias = c()
alphas1_rq_bias = c()
alphas1_cq_bias = c()
thetas1_rq_bias = c()
thetas1_cq_bias = c()

# Other levels of bandwidth (h)
alphas0_cq2_bias = c()      
alphas1_cq2_bias = c()
thetas1_cq2_bias = c()

alphas0_cq3_bias = c()      
alphas1_cq3_bias = c()
thetas1_cq3_bias = c()

alphas0_cq_3_bias = c()      
alphas1_cq_3_bias = c()
thetas1_cq_3_bias = c()

alphas0_cq_2_bias = c()      
alphas1_cq_2_bias = c()
thetas1_cq_2_bias = c()

alphas0_cq8_bias = c()     
alphas1_cq8_bias = c()
thetas1_cq8_bias = c()

# alphas0_cq16_bias = c()      
# alphas1_cq16_bias = c()
# thetas1_cq16_bias = c()
# 
# alphas0_cq_16_bias = c()      
# alphas1_cq_16_bias = c()
# thetas1_cq_16_bias = c()
# 
# alphas0_cq_8_bias = c()      
# alphas1_cq_8_bias = c()
# thetas1_cq_8_bias = c()

# MSE
alphas0_rq_mse = c()
alphas0_cq_mse = c()
alphas1_rq_mse = c()
alphas1_cq_mse = c()
thetas1_rq_mse = c()
thetas1_cq_mse = c()

# Other levels of bandwidths (h)
alphas0_cq2_mse = c()
alphas1_cq2_mse = c()
thetas1_cq2_mse = c()

alphas0_cq3_mse = c()
alphas1_cq3_mse = c()
thetas1_cq3_mse = c()

alphas0_cq_3_mse = c()
alphas1_cq_3_mse = c()
thetas1_cq_3_mse = c()

alphas0_cq_2_mse = c()
alphas1_cq_2_mse = c()
thetas1_cq_2_mse = c()

alphas0_cq8_mse = c()
alphas1_cq8_mse = c()
thetas1_cq8_mse = c()

# alphas0_cq16_mse = c()
# alphas1_cq16_mse = c()
# thetas1_cq16_mse = c()
# 
# alphas0_cq_16_mse = c()
# alphas1_cq_16_mse = c()
# thetas1_cq_16_mse = c()
# 
# alphas0_cq_8_mse = c()
# alphas1_cq_8_mse = c()
# thetas1_cq_8_mse = c()

real_alpha0 = sapply(tau.grid, alpha0)
#real_alpha1 = 0
real_alpha1 = sapply(tau.grid, alpha1)
real_theta1 = sapply(tau.grid, theta1)


for(s in 1:length(tau.grid)){
  
  # Bias 1st - 50th quantiles
  alphas0_rq_bias[s] = mean(A[s,1,]) - real_alpha0[s]
  alphas0_cq_bias[s] = mean(C[s,1,]) - real_alpha0[s]
  alphas1_rq_bias[s] = mean(A[s,2,]) - real_alpha1[s]
  alphas1_cq_bias[s] = mean(C[s,2,]) - real_alpha1[s]
  thetas1_rq_bias[s] = mean(A[s,3,]) - real_theta1[s]
  thetas1_cq_bias[s] = mean(C[s,3,]) - real_theta1[s]
  # h*2
  alphas0_cq2_bias[s] = mean(C2[s,1,]) - real_alpha0[s]
  alphas1_cq2_bias[s] = mean(C2[s,2,]) - real_alpha1[s]
  thetas1_cq2_bias[s] = mean(C2[s,3,]) - real_theta1[s]
  # h*4
  alphas0_cq3_bias[s] = mean(C3[s,1,]) - real_alpha0[s]
  alphas1_cq3_bias[s] = mean(C3[s,2,]) - real_alpha1[s]
  thetas1_cq3_bias[s] = mean(C3[s,3,]) - real_theta1[s]
  # h/4
  alphas0_cq_3_bias[s] = mean(C_3[s,1,]) - real_alpha0[s]
  alphas1_cq_3_bias[s] = mean(C_3[s,2,]) - real_alpha1[s]
  thetas1_cq_3_bias[s] = mean(C_3[s,3,]) - real_theta1[s]
  # h/2
  alphas0_cq_2_bias[s] = mean(C_2[s,1,]) - real_alpha0[s]
  alphas1_cq_2_bias[s] = mean(C_2[s,2,]) - real_alpha1[s]
  thetas1_cq_2_bias[s] = mean(C_2[s,3,]) - real_theta1[s]
  # h*8
  alphas0_cq8_bias[s] = mean(C8[s,1,]) - real_alpha0[s]
  alphas1_cq8_bias[s] = mean(C8[s,2,]) - real_alpha1[s]
  thetas1_cq8_bias[s] = mean(C8[s,3,]) - real_theta1[s]
  # # h*16
  # alphas0_cq16_bias[s] = mean(C16[s,1,]) - real_alpha0[s]
  # alphas1_cq16_bias[s] = mean(C16[s,2,]) - real_alpha1[s]
  # thetas1_cq16_bias[s] = mean(C16[s,3,]) - real_theta1[s]
  # # h/16
  # alphas0_cq_16_bias[s] = mean(C_16[s,1,]) - real_alpha0[s]
  # alphas1_cq_16_bias[s] = mean(C_16[s,2,]) - real_alpha1[s]
  # thetas1_cq_16_bias[s] = mean(C_16[s,3,]) - real_theta1[s]
  # # h/8
  # alphas0_cq_8_bias[s] = mean(C_8[s,1,]) - real_alpha0[s]
  # alphas1_cq_8_bias[s] = mean(C_8[s,2,]) - real_alpha1[s]
  # thetas1_cq_8_bias[s] = mean(C_8[s,3,]) - real_theta1[s]
  
  # MSE 1st - 50th quantiles
  alphas0_rq_mse[s] = mean((A[s,1,] - real_alpha0[s])^2)
  alphas0_cq_mse[s] = mean((C[s,1,] - real_alpha0[s])^2)
  alphas1_rq_mse[s] = mean((A[s,2,] - real_alpha1[s])^2)
  alphas1_cq_mse[s] = mean((C[s,2,] - real_alpha1[s])^2)
  thetas1_rq_mse[s] = mean((A[s,3,] - real_theta1[s])^2)
  thetas1_cq_mse[s] = mean((C[s,3,] - real_theta1[s])^2)
  # h*2
  alphas0_cq2_mse[s] = mean((C2[s,1,] - real_alpha0[s])^2)
  alphas1_cq2_mse[s] = mean((C2[s,2,] - real_alpha1[s])^2)
  thetas1_cq2_mse[s] = mean((C2[s,3,] - real_theta1[s])^2)
  # h*4
  alphas0_cq3_mse[s] = mean((C3[s,1,] - real_alpha0[s])^2)
  alphas1_cq3_mse[s] = mean((C3[s,2,] - real_alpha1[s])^2)
  thetas1_cq3_mse[s] = mean((C3[s,3,] - real_theta1[s])^2)
  # h/4
  alphas0_cq_3_mse[s] = mean((C_3[s,1,] - real_alpha0[s])^2)
  alphas1_cq_3_mse[s] = mean((C_3[s,2,] - real_alpha1[s])^2)
  thetas1_cq_3_mse[s] = mean((C_3[s,3,] - real_theta1[s])^2)
  # h*1/2
  alphas0_cq_2_mse[s] = mean((C_2[s,1,] - real_alpha0[s])^2)
  alphas1_cq_2_mse[s] = mean((C_2[s,2,] - real_alpha1[s])^2)
  thetas1_cq_2_mse[s] = mean((C_2[s,3,] - real_theta1[s])^2)
  # h*8
  alphas0_cq8_mse[s] = mean((C8[s,1,] - real_alpha0[s])^2)
  alphas1_cq8_mse[s] = mean((C8[s,2,] - real_alpha1[s])^2)
  thetas1_cq8_mse[s] = mean((C8[s,3,] - real_theta1[s])^2)
  # # h*16
  # alphas0_cq16_mse[s] = mean((C16[s,1,] - real_alpha0[s])^2)
  # alphas1_cq16_mse[s] = mean((C16[s,2,] - real_alpha1[s])^2)
  # thetas1_cq16_mse[s] = mean((C16[s,3,] - real_theta1[s])^2)
  # # h/16
  # alphas0_cq_16_mse[s] = mean((C_16[s,1,] - real_alpha0[s])^2)
  # alphas1_cq_16_mse[s] = mean((C_16[s,2,] - real_alpha1[s])^2)
  # thetas1_cq_16_mse[s] = mean((C_16[s,3,] - real_theta1[s])^2)
  # # h/8
  # alphas0_cq_8_mse[s] = mean((C_8[s,1,] - real_alpha0[s])^2)
  # alphas1_cq_8_mse[s] = mean((C_8[s,2,] - real_alpha1[s])^2)
  # thetas1_cq_8_mse[s] = mean((C_8[s,3,] - real_theta1[s])^2)
  # 
  # alphas0_rq_global[s] = sum((abs(A[s,1,] - real_alpha0[s]))^2)*0.02
}

### Plots - for each estimation, we plot the mean of the estimated values for each coefficient: mpq = mean per quantile (considering all replications)

# Storing coefficients for plotting
alphas0_rq_mpq = c()
alphas0_cq_mpq = c()
alphas1_rq_mpq = c()
alphas1_cq_mpq = c()
thetas1_rq_mpq = c()
thetas1_cq_mpq = c()
# h*2
alphas0_cq2_mpq = c()
alphas1_cq2_mpq = c()
thetas1_cq2_mpq = c()
# h*4
alphas0_cq3_mpq = c()
alphas1_cq3_mpq = c()
thetas1_cq3_mpq = c()
# h/4
alphas0_cq_3_mpq = c()
alphas1_cq_3_mpq = c()
thetas1_cq_3_mpq = c()
# h/2
alphas0_cq_2_mpq = c()
alphas1_cq_2_mpq = c()
thetas1_cq_2_mpq = c()
# h*8
alphas0_cq8_mpq = c()
alphas1_cq8_mpq = c()
thetas1_cq8_mpq = c()
# # h*16
# alphas0_cq16_mpq = c()
# alphas1_cq16_mpq = c()
# thetas1_cq16_mpq = c()
# # h/16
# alphas0_cq_16_mpq = c()
# alphas1_cq_16_mpq = c()
# thetas1_cq_16_mpq = c()
# # h/8
# alphas0_cq_8_mpq = c()
# alphas1_cq_8_mpq = c()
# thetas1_cq_8_mpq = c()

# Variance
alphas0_rq_vpq = c()
alphas0_cq_vpq = c()
alphas1_rq_vpq = c()
alphas1_cq_vpq = c()
thetas1_rq_vpq = c()
thetas1_cq_vpq = c()
# h*2
alphas0_cq2_vpq = c()
alphas1_cq2_vpq = c()
thetas1_cq2_vpq = c()
# h*4
alphas0_cq3_vpq = c()
alphas1_cq3_vpq = c()
thetas1_cq3_vpq = c()
# h/4
alphas0_cq_3_vpq = c()
alphas1_cq_3_vpq = c()
thetas1_cq_3_vpq = c()
# h/2
alphas0_cq_2_vpq = c()
alphas1_cq_2_vpq = c()
thetas1_cq_2_vpq = c()
# h*8
alphas0_cq8_vpq = c()
alphas1_cq8_vpq = c()
thetas1_cq8_vpq = c()
# # h*16
# alphas0_cq16_vpq = c()
# alphas1_cq16_vpq = c()
# thetas1_cq16_vpq = c()
# # h/16
# alphas0_cq_16_vpq = c()
# alphas1_cq_16_vpq = c()
# thetas1_cq_16_vpq = c()
# # h/8
# alphas0_cq_8_vpq = c()
# alphas1_cq_8_vpq = c()
# thetas1_cq_8_vpq = c()

for(g in 1:length(tau.grid)){
  # Mean per quantiles
  alphas0_rq_mpq[g] = mean(A[g,1,])
  alphas0_cq_mpq[g] = mean(C[g,1,])
  alphas1_rq_mpq[g] = mean(A[g,2,])
  alphas1_cq_mpq[g] = mean(C[g,2,])
  thetas1_rq_mpq[g] = mean(A[g,3,])
  thetas1_cq_mpq[g] = mean(C[g,3,])
  # h*2
  alphas0_cq2_mpq[g] = mean(C2[g,1,])
  alphas1_cq2_mpq[g] = mean(C2[g,2,])
  thetas1_cq2_mpq[g] = mean(C2[g,3,])
  # h*3
  alphas0_cq3_mpq[g] = mean(C3[g,1,])
  alphas1_cq3_mpq[g] = mean(C3[g,2,])
  thetas1_cq3_mpq[g] = mean(C3[g,3,])
  # h/3
  alphas0_cq_3_mpq[g] = mean(C_3[g,1,])
  alphas1_cq_3_mpq[g] = mean(C_3[g,2,])
  thetas1_cq_3_mpq[g] = mean(C_3[g,3,])
  # h/2
  alphas0_cq_2_mpq[g] = mean(C_2[g,1,])
  alphas1_cq_2_mpq[g] = mean(C_2[g,2,])
  thetas1_cq_2_mpq[g] = mean(C_2[g,3,])
  # h*8
  alphas0_cq8_mpq[g] = mean(C8[g,1,])
  alphas1_cq8_mpq[g] = mean(C8[g,2,])
  thetas1_cq8_mpq[g] = mean(C8[g,3,])
  # # h*16
  # alphas0_cq16_mpq[g] = mean(C16[g,1,])
  # alphas1_cq16_mpq[g] = mean(C16[g,2,])
  # thetas1_cq16_mpq[g] = mean(C16[g,3,])
  # # h/16
  # alphas0_cq_16_mpq[g] = mean(C_16[g,1,])
  # alphas1_cq_16_mpq[g] = mean(C_16[g,2,])
  # thetas1_cq_16_mpq[g] = mean(C_16[g,3,])
  # # h/8
  # alphas0_cq_8_mpq[g] = mean(C_8[g,1,])
  # alphas1_cq_8_mpq[g] = mean(C_8[g,2,])
  # thetas1_cq_8_mpq[g] = mean(C_8[g,3,])
  
  # Variance per quantiles
  alphas0_rq_vpq[g] = var(A[g,1,])
  alphas0_cq_vpq[g] = var(C[g,1,])
  alphas1_rq_vpq[g] = var(A[g,2,])
  alphas1_cq_vpq[g] = var(C[g,2,])
  thetas1_rq_vpq[g] = var(A[g,3,])
  thetas1_cq_vpq[g] = var(C[g,3,])
  # h*2
  alphas0_cq2_vpq[g] = var(C2[g,1,])
  alphas1_cq2_vpq[g] = var(C2[g,2,])
  thetas1_cq2_vpq[g] = var(C2[g,3,])
  # h*4
  alphas0_cq3_vpq[g] = var(C3[g,1,])
  alphas1_cq3_vpq[g] = var(C3[g,2,])
  thetas1_cq3_vpq[g] = var(C3[g,3,])
  # h/4
  alphas0_cq_3_vpq[g] = var(C_3[g,1,])
  alphas1_cq_3_vpq[g] = var(C_3[g,2,])
  thetas1_cq_3_vpq[g] = var(C_3[g,3,])
  # h/2
  alphas0_cq_2_vpq[g] = var(C_2[g,1,])
  alphas1_cq_2_vpq[g] = var(C_2[g,2,])
  thetas1_cq_2_vpq[g] = var(C_2[g,3,])
  # h*8
  alphas0_cq8_vpq[g] = var(C8[g,1,])
  alphas1_cq8_vpq[g] = var(C8[g,2,])
  thetas1_cq8_vpq[g] = var(C8[g,3,])
  # # h*16
  # alphas0_cq16_vpq[g] = var(C16[g,1,])
  # alphas1_cq16_vpq[g] = var(C16[g,2,])
  # thetas1_cq16_vpq[g] = var(C16[g,3,])
  # # h/16
  # alphas0_cq_16_vpq[g] = var(C_16[g,1,])
  # alphas1_cq_16_vpq[g] = var(C_16[g,2,])
  # thetas1_cq_16_vpq[g] = var(C_16[g,3,])
  # # h/8
  # alphas0_cq_8_vpq[g] = var(C_8[g,1,])
  # alphas1_cq_8_vpq[g] = var(C_8[g,2,])
  # thetas1_cq_8_vpq[g] = var(C_8[g,3,])
}

# Final plots - alphas0 (alpha_0), alpha_1, beta_1 & beta_2  - MEAN per quantile - - h, h*2, h*4, h*1/4, h*1/2

# coefs_alphas0_mpq = data.frame(alphas0_rq_mpq, alphas0_cq_mpq, 
#                                alphas0_cq2_mpq,alphas0_cq3_mpq,alphas0_cq_3_mpq,alphas0_cq_2_mpq,
#                                alphas0_cq8_mpq,alphas0_cq16_mpq,alphas0_cq_16_mpq,alphas0_cq_8_mpq) 
# coefs_alphas1_mpq = data.frame(alphas1_rq_mpq, alphas1_cq_mpq, 
#                                alphas1_cq2_mpq,alphas1_cq3_mpq,alphas1_cq_3_mpq,alphas1_cq_2_mpq,
#                                alphas1_cq8_mpq,alphas1_cq16_mpq,alphas1_cq_16_mpq,alphas1_cq_8_mpq) 
# coefs_thetas1_mpq = data.frame(thetas1_rq_mpq, thetas1_cq_mpq,
#                                thetas1_cq2_mpq,thetas1_cq3_mpq,thetas1_cq_3_mpq,thetas1_cq_2_mpq,
#                                thetas1_cq8_mpq,thetas1_cq16_mpq,thetas1_cq_16_mpq,thetas1_cq_8_mpq)

# xlab = TeX('$\\tau$')
# 
# ylab = TeX('$\\hat{\\alpha_0}(\\tau)$') 
# grafico_alpha0_mpq = ggplot(coefs_alphas0_mpq, aes(x = tau.grid, y = value)) + 
#   geom_line(aes(y = alphas0_rq_mpq,linetype ='blank',color ='QR'),size=1.0) +
#   geom_line(aes(y = alphas0_cq_mpq, linetype ='solid',color ='SQR'),size=1.0)+
#   geom_line(aes(y = alphas0_cq2_mpq, linetype ='dotdash',color ='SQR_h*2'),size=1.0)+
#   geom_line(aes(y = alphas0_cq3_mpq, linetype ='dotdash',color ='SQR_h*4'),size=1.0)+
#   geom_line(aes(y = alphas0_cq_3_mpq, linetype ='dotdash',color ='SQR_h/4'),size=1.0)+
#   geom_line(aes(y = alphas0_cq_2_mpq, linetype ='dotted',color ='SQR_h/2'),size=1.0)+
#   geom_line(aes(y = alphas0_cq8_mpq, linetype ='dotdash',color ='SQR_h*8'),size=1.0)+
#   geom_line(aes(y = alphas0_cq16_mpq, linetype ='dotdash',color ='SQR_h*16'),size=1.0)+
#   geom_line(aes(y = alphas0_cq_16_mpq, linetype ='dotdash',color ='SQR_h/16'),size=1.0)+
#   geom_line(aes(y = alphas0_cq_8_mpq, linetype ='dotted',color ='SQR_h/8'),size=1.0)+
#   geom_hline(yintercept = 0, color = 'yellow', size=1.3) +
#   labs(x =xlab,y = ylab) + #expand_limits(x=c(0, 1), y=c(-3,3)) +
#   scale_x_continuous(breaks=seq(0.01, 0.99, 0.07),limits = c(0.01, 0.99))  + 
#   theme_grey() + guides(linetype = FALSE) +
#   theme(plot.title = element_text(size = 15),
#         legend.title = element_blank(),
#         legend.text = element_text(size = 12),
#         axis.title = element_text(size = 14))
# print(grafico_alpha0_mpq)
# 
# ylab = TeX('$\\hat{\\alpha_1}(\\tau)$') 
# grafico_alpha1_mpq = ggplot(coefs_alphas1_mpq, aes(x = tau.grid, y = value)) + 
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
#   geom_hline(yintercept = 0.5, color = 'yellow', size=1.3) +
#   labs(x =xlab,y = ylab) + #expand_limits(x=c(0, 1), y=c(0.4,0.6)) + 
#   scale_x_continuous(breaks=round(seq(0.01, 0.99, 0.07), 2),limits = c(0.01, 0.99))  + 
#   #scale_y_continuous(breaks=seq(0.4,0.6,0.02),limits = c(0.46, 0.54))  +
#   theme_grey() + guides(linetype = FALSE) +
#   theme(plot.title = element_text(size = 15),
#         legend.title = element_blank(),
#         legend.text = element_text(size = 12),
#         axis.title = element_text(size = 14))
# print(grafico_alpha1_mpq)
# 
# ylab = TeX('$\\hat{\\theta_1}(\\tau)$')
# grafico_beta1_mpq = ggplot(coefs_thetas1_mpq, aes(x = tau.grid, y = value)) + 
#   geom_line(aes(y = thetas1_rq_mpq,linetype ='blank',color ='QR'),size=1.0) +
#   geom_line(aes(y = thetas1_cq_mpq, linetype ='solid',color ='SQR'),size=1.,)+
#   geom_line(aes(y = thetas1_cq2_mpq, linetype ='dotdash',color ='SQR_h*2'),size=1.0)+
#   geom_line(aes(y = thetas1_cq3_mpq, linetype ='dotdash',color ='SQR_h*4'),size=1.0)+
#   geom_line(aes(y = thetas1_cq_3_mpq, linetype ='dotdash',color ='SQR_h/4'),size=1.0)+
#   geom_line(aes(y = thetas1_cq_2_mpq, linetype ='dotted',color ='SQR_h/2'),size=1.0)+
#   geom_line(aes(y = thetas1_cq8_mpq, linetype ='dotdash',color ='SQR_h*8'),size=1.0)+
#   geom_line(aes(y = thetas1_cq16_mpq, linetype ='dotdash',color ='SQR_h*16'),size=1.0)+
#   geom_line(aes(y = thetas1_cq_16_mpq, linetype ='dotdash',color ='SQR_h/16'),size=1.0)+
#   geom_line(aes(y = thetas1_cq_8_mpq, linetype ='dotted',color ='SQR_h/8'),size=1.0)+
#   geom_hline(yintercept = 0.5, color = 'yellow', size=1.3) +
#   labs(x =xlab,y = ylab) + # ggtitle(bquote(beta[1][tau]^{tau})) + 
#   scale_x_continuous(breaks=seq(0.01, 0.99, 0.07),limits = c(0.01, 0.99))  + 
#   theme_grey() + guides(linetype = FALSE) +
#   theme(plot.title = element_text(size = 15),
#         legend.title = element_blank(),
#         legend.text = element_text(size = 12),
#         axis.title = element_text(size = 14))
# print(grafico_beta1_mpq)

#plot_grid(grafico_alpha0_mpq,grafico_alpha1_mpq,grafico_beta1_mpq,grafico_beta2_mpq)

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

# alphas0_cq16_global = c()
# alphas1_cq16_global = c()
# thetas1_cq16_global = c()
# 
# alphas0_cq_16_global = c()
# alphas1_cq_16_global = c()
# thetas1_cq_16_global = c()
# 
# alphas0_cq_8_global = c()
# alphas1_cq_8_global = c()
# thetas1_cq_8_global = c()

for(j in 1:nrep){
  alphas0_rq_global[j] = sum((abs(A[,1,j] - real_alpha0))^2)*0.02
  alphas1_rq_global[j] = sum((abs(A[,2,j] - real_alpha1))^2)*0.02
  thetas1_rq_global[j] = sum((abs(A[,3,j] - real_theta1))^2)*0.02
  
  alphas0_cq_global[j] = sum((abs(C[,1,j] - real_alpha0))^2)*0.02
  alphas1_cq_global[j] = sum((abs(C[,2,j] - real_alpha1))^2)*0.02
  thetas1_cq_global[j] = sum((abs(C[,3,j] - real_theta1))^2)*0.02
  
  alphas0_cq2_global[j] = sum((abs(C2[,1,j] - real_alpha0))^2)*0.02
  alphas1_cq2_global[j] = sum((abs(C2[,2,j] - real_alpha1))^2)*0.02
  thetas1_cq2_global[j] = sum((abs(C2[,3,j] - real_theta1))^2)*0.02
  
  alphas0_cq3_global[j] = sum((abs(C3[,1,j] - real_alpha0))^2)*0.02
  alphas1_cq3_global[j] = sum((abs(C3[,2,j] - real_alpha1))^2)*0.02
  thetas1_cq3_global[j] = sum((abs(C3[,3,j] - real_theta1))^2)*0.02
  
  alphas0_cq_3_global[j] = sum((abs(C_3[,1,j] - real_alpha0))^2)*0.02
  alphas1_cq_3_global[j] = sum((abs(C_3[,2,j] - real_alpha1))^2)*0.02
  thetas1_cq_3_global[j] = sum((abs(C_3[,3,j] - real_theta1))^2)*0.02
  
  alphas0_cq_2_global[j] = sum((abs(C_2[,1,j] - real_alpha0))^2)*0.02
  alphas1_cq_2_global[j] = sum((abs(C_2[,2,j] - real_alpha1))^2)*0.02
  thetas1_cq_2_global[j] = sum((abs(C_2[,3,j] - real_theta1))^2)*0.02
  
  alphas0_cq8_global[j] = sum((abs(C8[,1,j] - real_alpha0))^2)*0.02
  alphas1_cq8_global[j] = sum((abs(C8[,2,j] - real_alpha1))^2)*0.02
  thetas1_cq8_global[j] = sum((abs(C8[,3,j] - real_theta1))^2)*0.02
  
  # alphas0_cq16_global[j] = sum((abs(C16[,1,j] - real_alpha0))^2)*0.02
  # alphas1_cq16_global[j] = sum((abs(C16[,2,j] - real_alpha1))^2)*0.02
  # thetas1_cq16_global[j] = sum((abs(C16[,3,j] - real_theta1))^2)*0.02
  # 
  # alphas0_cq_16_global[j] = sum((abs(C_16[,1,j] - real_alpha0))^2)*0.02
  # alphas1_cq_16_global[j] = sum((abs(C_16[,2,j] - real_alpha1))^2)*0.02
  # thetas1_cq_16_global[j] = sum((abs(C_16[,3,j] - real_theta1))^2)*0.02
  # 
  # alphas0_cq_8_global[j] = sum((abs(C_8[,1,j] - real_alpha0))^2)*0.02
  # alphas1_cq_8_global[j] = sum((abs(C_8[,2,j] - real_alpha1))^2)*0.02
  # thetas1_cq_8_global[j] = sum((abs(C_8[,3,j] - real_theta1))^2)*0.02
}

global_rq = mean(alphas0_rq_global + alphas1_rq_global + thetas1_rq_global)
global_cq = mean(alphas0_cq_global + alphas1_cq_global + thetas1_cq_global)
global_cq2 = mean(alphas0_cq2_global + alphas1_cq2_global + thetas1_cq2_global)
global_cq3 = mean(alphas0_cq3_global + alphas1_cq3_global + thetas1_cq3_global)
global_cq_3 = mean(alphas0_cq_3_global + alphas1_cq_3_global + thetas1_cq_3_global)
global_cq_2 = mean(alphas0_cq_2_global + alphas1_cq_2_global + thetas1_cq_2_global)
global_cq8 = mean(alphas0_cq8_global + alphas1_cq8_global + thetas1_cq8_global)
# global_cq16 = mean(alphas0_cq16_global + alphas1_cq16_global + thetas1_cq16_global)
# global_cq_16 = mean(alphas0_cq_16_global + alphas1_cq_16_global + thetas1_cq_16_global)
# global_cq_8 = mean(alphas0_cq_8_global + alphas1_cq_8_global + thetas1_cq_8_global)

globals_df = data.frame(cbind(global_rq,global_cq,
                              global_cq2,global_cq3,global_cq_3,global_cq_2,
                              global_cq8))#,global_cq16,global_cq_16,global_cq_8))

# WRITING CHECKPOINTS FOR THE OUTPUTS

coefs_alphas0_mse = data.frame(cbind(alphas0_rq_mse,alphas0_cq_mse,
                             alphas0_cq2_mse,alphas0_cq3_mse,alphas0_cq_3_mse,alphas0_cq_2_mse,
                             alphas0_cq8_mse))#,alphas0_cq16_mse,alphas0_cq_16_mse,alphas0_cq_8_mse)) 

coefs_alphas1_mse = data.frame(cbind(alphas1_rq_mse,alphas1_cq_mse,
                             alphas1_cq2_mse,alphas1_cq3_mse,alphas1_cq_3_mse,alphas1_cq_2_mse,
                             alphas1_cq8_mse))#,alphas1_cq16_mse,alphas1_cq_16_mse,alphas1_cq_8_mse)) 

coefs_thetas1_mse = data.frame(cbind(thetas1_rq_mse,thetas1_cq_mse,
                             thetas1_cq2_mse,thetas1_cq3_mse,thetas1_cq_3_mse,thetas1_cq_2_mse,
                             thetas1_cq8_mse))#,thetas1_cq16_mse,thetas1_cq_16_mse,thetas1_cq_8_mse))

coefs_alphas0_bias = data.frame(cbind(alphas0_rq_bias,alphas0_cq_bias,
                              alphas0_cq2_bias,alphas0_cq3_bias,alphas0_cq_3_bias,alphas0_cq_2_bias,
                              alphas0_cq8_bias))#,alphas0_cq16_bias,alphas0_cq_16_bias,alphas0_cq_8_bias))

coefs_alphas1_bias = data.frame(cbind(alphas1_rq_bias,alphas1_cq_bias,
                              alphas1_cq2_bias,alphas1_cq3_bias,alphas1_cq_3_bias,alphas1_cq_2_bias,
                              alphas1_cq8_bias))#,alphas1_cq16_bias,alphas1_cq_16_bias,alphas1_cq_8_bias))

coefs_thetas1_bias = data.frame(cbind(thetas1_rq_bias,thetas1_cq_bias,
                              thetas1_cq2_bias,thetas1_cq3_bias,thetas1_cq_3_bias,thetas1_cq_2_bias,
                              thetas1_cq8_bias))#,thetas1_cq16_bias,thetas1_cq_16_bias,thetas1_cq_8_bias)) 

coefs_alphas0_vpq = data.frame(cbind(alphas0_rq_vpq, alphas0_cq_vpq,
                               alphas0_cq2_vpq,alphas0_cq3_vpq,alphas0_cq_3_vpq,alphas0_cq_2_vpq,
                               alphas0_cq8_vpq))#,alphas0_cq16_vpq,alphas0_cq_16_vpq,alphas0_cq_8_vpq)

coefs_alphas1_vpq = data.frame(cbind(alphas1_rq_vpq, alphas1_cq_vpq,
                               alphas1_cq2_vpq,alphas1_cq3_vpq,alphas1_cq_3_vpq,alphas1_cq_2_vpq,
                               alphas1_cq8_vpq))
                               #,alphas1_cq16_vpq,alphas1_cq_16_vpq,alphas1_cq_8_vpq)
                               
coefs_thetas1_vpq = data.frame(cbind(thetas1_rq_vpq, thetas1_cq_vpq,
                               thetas1_cq2_vpq,thetas1_cq3_vpq,thetas1_cq_3_vpq,thetas1_cq_2_vpq,
                               thetas1_cq8_vpq))#,thetas1_cq16_vpq,thetas1_cq_16_vpq,thetas1_cq_8_vpq)


path_outputs = "C:/Users/MIGUEL/OneDrive - UFRGS/Mestrado PPGEst/DISSERTAÇÃO/Novas Simulações/outputs_checkpoints/"
path_lambda = choice
path_sample = as.character(T)
path_write = paste(path_outputs, choice, sep="")
setwd(path_write)
# Coefficientes - mpq
write.csv(coefs_alphas0_mpq, file= paste('coefs_alphas0_mpq-n',path_sample,'.csv',sep=""), row.names=FALSE)
write.csv(coefs_alphas1_mpq, file= paste('coefs_alphas1_mpq-n',path_sample,'.csv',sep=""), row.names=FALSE)
write.csv(coefs_thetas1_mpq, file= paste('coefs_thetas1_mpq-n',path_sample,'.csv',sep=""), row.names=FALSE)
# Coefficientes - vpq
write.csv(coefs_alphas0_vpq, file= paste('coefs_alphas0_vpq-n',path_sample,'.csv',sep=""), row.names=FALSE)
write.csv(coefs_alphas1_vpq, file= paste('coefs_alphas1_vpq-n',path_sample,'.csv',sep=""), row.names=FALSE)
write.csv(coefs_thetas1_vpq, file= paste('coefs_thetas_vpq-n',path_sample,'.csv',sep=""), row.names=FALSE)
# Coefficientes - mse
write.csv(coefs_alphas0_mse, file= paste('coefs_alphas0_mse-n',path_sample,'.csv',sep=""), row.names=FALSE)
write.csv(coefs_alphas1_mse, file= paste('coefs_alphas1_mse-n',path_sample,'.csv',sep=""), row.names=FALSE)
write.csv(coefs_thetas1_mse, file= paste('coefs_thetas1_mse-n',path_sample,'.csv',sep=""), row.names=FALSE)
# Coefficientes - bias
write.csv(coefs_alphas0_bias, file= paste('coefs_alphas0_bias-n',path_sample,'.csv',sep=""), row.names=FALSE)
write.csv(coefs_alphas1_bias, file= paste('coefs_alphas1_bias-n',path_sample,'.csv',sep=""), row.names=FALSE)
write.csv(coefs_thetas1_bias, file= paste('coefs_thetas1_bias-n',path_sample,'.csv',sep=""), row.names=FALSE)
# Global metrics
write.csv(globals_df, file= paste('globals_df-n',path_sample,'.csv',sep=""))

#Saving arrays for the estimated coefficients
saveRDS(A, "A.rds")
saveRDS(C, "C.rds")
saveRDS(C2, "C2.rds")
saveRDS(C3, "C4.rds")
saveRDS(C_3, "C_4.rds")
saveRDS(C_2, "C_2.rds")
saveRDS(C8, "C8.rds")
# saveRDS(C16, "C16.rds")
# saveRDS(C_16, "C_16.rds")
# saveRDS(C_8, "C_8.rds")
#Atest <- readRDS("A.rds")

end_time = Sys.time()

### FINAL PLOTS - WITHOUT *8,*16,/8,/16 bandwidths

#### PLOTS MSE

mse_coefs = data.frame(cbind(alphas0_rq_mse,alphas0_cq_mse,
                             alphas0_cq2_mse,alphas0_cq3_mse,alphas0_cq_3_mse,alphas0_cq_2_mse,
                             alphas0_cq8_mse))#,alphas0_cq16_mse,alphas0_cq_16_mse,alphas0_cq_8_mse))
ylab = TeX('MSE $\\hat{\\alpha_0}(\\tau)$')
grafico_alphas0_mse = ggplot(mse_coefs, aes(x = tau.grid, y = value)) + 
  geom_line(aes(y = alphas0_rq_mse,linetype ='blank',color ='QR'),size=1.0) +
  geom_line(aes(y = alphas0_cq_mse,linetype ='solid',color ='SQR'),size=1.0)+
  geom_line(aes(y = alphas0_cq2_mse,linetype ='dotdash',color ='SQR_h*2'),size=1.0)+
  geom_line(aes(y = alphas0_cq3_mse,linetype ='dotdash',color ='SQR_h*4'),size=1.0)+
  geom_line(aes(y = alphas0_cq_3_mse,linetype ='dotdash',color ='SQR_h/4'),size=1.0)+
  geom_line(aes(y = alphas0_cq_2_mse,linetype ='dotted',color ='SQR_h/2'),size=1.0)+
  scale_color_discrete(labels = unname(TeX(c("$QR$","$SQR_{\\zeta^*}$","$SQR_{2\\zeta^*}$","$SQR_{4\\zeta^*}$","$SQR_{\\zeta^{*}/4}$","$SQR_{\\zeta^{*}/2}$")))) + 
  # geom_line(aes(y = alphas0_cq8_mse,linetype ='dotdash',color ='SQR_h*8'),size=1.0)+
  # geom_line(aes(y = alphas0_cq16_mse,linetype ='dotdash',color ='SQR_h*16'),size=1.0)+
  # geom_line(aes(y = alphas0_cq_16_mse,linetype ='dotdash',color ='SQR_h/16'),size=1.0)+
  # geom_line(aes(y = alphas0_cq_8_mse,linetype ='dotted',color ='SQR_h/8'),size=1.0)+
  geom_hline(yintercept = 0.0, color = 'yellow', size=1.3) +
  labs(x =xlab,y = ylab)  + 
  scale_x_continuous(breaks=seq(0.01, 0.99, 0.07),limits = c(0.01, 0.99))  + 
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
                             alphas1_cq8_mse))#,alphas1_cq16_mse,alphas1_cq_16_mse,alphas1_cq_8_mse))
ylab = TeX('MSE $\\hat{\\alpha_1}(\\tau)$')
grafico_alphas1_mse = ggplot(mse_coefs, aes(x = tau.grid, y = value)) + 
  geom_line(aes(y = alphas1_rq_mse,linetype ='blank',color ='QR'),size=1.0) +
  geom_line(aes(y = alphas1_cq_mse,linetype ='solid',color ='SQR'),size=1.0)+
  geom_line(aes(y = alphas1_cq2_mse,linetype ='dotdash',color ='SQR_h*2'),size=1.0)+
  geom_line(aes(y = alphas1_cq3_mse,linetype ='dotdash',color ='SQR_h*4'),size=1.0)+
  geom_line(aes(y = alphas1_cq_3_mse,linetype ='dotdash',color ='SQR_h/4'),size=1.0)+
  geom_line(aes(y = alphas1_cq_2_mse,linetype ='dotted',color ='SQR_h/2'),size=1.0)+
  scale_color_discrete(labels = unname(TeX(c("$QR$","$SQR_{\\zeta^*}$","$SQR_{2\\zeta^*}$","$SQR_{4\\zeta^*}$","$SQR_{\\zeta^{*}/4}$","$SQR_{\\zeta^{*}/2}$")))) + 
  # geom_line(aes(y = alphas1_cq8_mse,linetype ='dotdash',color ='SQR_h*8'),size=1.0)+
  # geom_line(aes(y = alphas1_cq16_mse,linetype ='dotdash',color ='SQR_h*16'),size=1.0)+
  # geom_line(aes(y = alphas1_cq_16_mse,linetype ='dotdash',color ='SQR_h/16'),size=1.0)+
  # geom_line(aes(y = alphas1_cq_8_mse,linetype ='dotted',color ='SQR_h/8'),size=1.0)+
  geom_hline(yintercept = 0.0, color = 'yellow', size=1.3) +
  labs(x =xlab,y = ylab)  + 
  scale_x_continuous(breaks=seq(0.01, 0.99, 0.07),limits = c(0.01, 0.99))  + 
  scale_y_log10() + # Log-scale
  #scale_y_log10(labels = trans_format("log10", math_format(10^.x))) + 
  theme_grey() + guides(linetype = FALSE) +
  theme(plot.title = element_text(size = 15),
        legend.title = element_blank(),
        legend.text.align = 0,
        legend.text = element_text(size = 12),
        axis.title = element_text(size = 14))
print(grafico_alphas1_mse)

mse_coefs = data.frame(cbind(thetas1_rq_mse,thetas1_cq_mse,
                             thetas1_cq2_mse,thetas1_cq3_mse,thetas1_cq_3_mse,thetas1_cq_2_mse,
                             thetas1_cq8_mse))#,thetas1_cq16_mse,thetas1_cq_16_mse,thetas1_cq_8_mse))
ylab = TeX('MSE $\\hat{\\theta_1}(\\tau)$')
grafico_thetas1_mse = ggplot(mse_coefs, aes(x = tau.grid, y = value)) + 
  geom_line(aes(y = thetas1_rq_mse,linetype ='blank',color ='QR'),size=1.0) +
  geom_line(aes(y = thetas1_cq_mse,linetype ='solid',color ='SQR'),size=1.0)+
  geom_line(aes(y = thetas1_cq2_mse,linetype ='dotdash',color ='SQR_h*2'),size=1.0)+
  geom_line(aes(y = thetas1_cq3_mse,linetype ='dotdash',color ='SQR_h*4'),size=1.0)+
  geom_line(aes(y = thetas1_cq_3_mse,linetype ='dotdash',color ='SQR_h/4'),size=1.0)+
  geom_line(aes(y = thetas1_cq_2_mse,linetype ='dotted',color ='SQR_h/2'),size=1.0)+
  scale_color_discrete(labels = unname(TeX(c("$QR$","$SQR_{\\zeta^*}$","$SQR_{2\\zeta^*}$","$SQR_{4\\zeta^*}$","$SQR_{\\zeta^{*}/4}$","$SQR_{\\zeta^{*}/2}$")))) + 
  # geom_line(aes(y = thetas1_cq8_mse,linetype ='dotdash',color ='SQR_h*8'),size=1.0)+
  # geom_line(aes(y = thetas1_cq16_mse,linetype ='dotdash',color ='SQR_h*16'),size=1.0)+
  # geom_line(aes(y = thetas1_cq_16_mse,linetype ='dotdash',color ='SQR_h/16'),size=1.0)+
  # geom_line(aes(y = thetas1_cq_8_mse,linetype ='dotted',color ='SQR_h/8'),size=1.0)+
  geom_hline(yintercept = 0.0, color = 'yellow', size=1.3) +
  labs(x =xlab,y = ylab)  + 
  scale_x_continuous(breaks=seq(0.01, 0.99, 0.07),limits = c(0.01, 0.99))  + 
  scale_y_log10() + # Log-scale
  #scale_y_log10(labels = trans_format("log10", math_format(10^.x))) + 
  theme_grey() + guides(linetype = FALSE) +
  theme(plot.title = element_text(size = 15),
        legend.title = element_blank(),
        legend.text.align = 0,
        legend.text = element_text(size = 12),
        axis.title = element_text(size = 14))
print(grafico_thetas1_mse)
#plot_grid(grafico_alphas0_mse,grafico_alphas1_mse,grafico_thetas1_mse,grafico_betas2_mse)

#### PLOTS BIAS

bias_coefs = data.frame(cbind(alphas0_rq_bias,alphas0_cq_bias,
                              alphas0_cq2_bias,alphas0_cq3_bias,alphas0_cq_3_bias,alphas0_cq_2_bias,
                              alphas0_cq8_bias))#,alphas0_cq16_bias,alphas0_cq_16_bias,alphas0_cq_8_bias))
ylab = TeX('Bias $\\hat{\\alpha_0}(\\tau)$')
grafico_alphas0_bias = ggplot(bias_coefs, aes(x = tau.grid, y = value)) + 
  geom_line(aes(y = alphas0_rq_bias,linetype ='blank',color ='QR'),size=1.0) +
  geom_line(aes(y = alphas0_cq_bias,linetype ='solid',color ='SQR'),size=1.0)+
  geom_line(aes(y = alphas0_cq2_bias,linetype ='dotdash',color ='SQR_h*2'),size=1.0)+
  geom_line(aes(y = alphas0_cq3_bias,linetype ='dotdash',color ='SQR_h*4'),size=1.0)+
  geom_line(aes(y = alphas0_cq_3_bias,linetype ='dotdash',color ='SQR_h/4'),size=1.0)+
  geom_line(aes(y = alphas0_cq_2_bias,linetype ='dotted',color ='SQR_h/2'),size=1.0)+
  scale_color_discrete(labels = unname(TeX(c("$QR$","$SQR_{\\zeta^*}$","$SQR_{2\\zeta^*}$","$SQR_{4\\zeta^*}$","$SQR_{\\zeta^{*}/4}$","$SQR_{\\zeta^{*}/2}$")))) + 
  # geom_line(aes(y = alphas0_cq8_bias,linetype ='dotdash',color ='SQR_h*8'),size=1.0)+
  # geom_line(aes(y = alphas0_cq16_bias,linetype ='dotdash',color ='SQR_h*16'),size=1.0)+
  # geom_line(aes(y = alphas0_cq_16_bias,linetype ='dotdash',color ='SQR_h/16'),size=1.0)+
  # geom_line(aes(y = alphas0_cq_8_bias,linetype ='dotted',color ='SQR_h/8'),size=1.0)+
  geom_hline(yintercept = 0.0, color = 'yellow', size=1.3) +
  labs(x =xlab,y = ylab)  + 
  scale_x_continuous(breaks=seq(0.01, 0.99, 0.07),limits = c(0.01, 0.99))  + 
  theme_grey() + guides(linetype = FALSE) +
  theme(plot.title = element_text(size = 15),
        legend.title = element_blank(),
        legend.text.align = 0,
        legend.text = element_text(size = 12),
        axis.title = element_text(size = 14))
print(grafico_alphas0_bias)

bias_coefs = data.frame(cbind(alphas1_rq_bias,alphas1_cq_bias,
                              alphas1_cq2_bias,alphas1_cq3_bias,alphas1_cq_3_bias,alphas1_cq_2_bias,
                              alphas1_cq8_bias))#,alphas1_cq16_bias,alphas1_cq_16_bias,alphas1_cq_8_bias))
ylab = TeX('Bias $\\hat{\\alpha_1}(\\tau)$')
grafico_alphas1_bias = ggplot(bias_coefs, aes(x = tau.grid, y = value)) + 
  geom_line(aes(y = alphas1_rq_bias,linetype ='blank',color ='QR'),size=1.0) +
  geom_line(aes(y = alphas1_cq_bias,linetype ='solid',color ='SQR'),size=1.0)+
  geom_line(aes(y = alphas1_cq2_bias,linetype ='dotdash',color ='SQR_h*2'),size=1.0)+
  geom_line(aes(y = alphas1_cq3_bias,linetype ='dotdash',color ='SQR_h*4'),size=1.0)+
  geom_line(aes(y = alphas1_cq_3_bias,linetype ='dotdash',color ='SQR_h/4'),size=1.0)+
  geom_line(aes(y = alphas1_cq_2_bias,linetype ='dotted',color ='SQR_h/2'),size=1.0)+
  scale_color_discrete(labels = unname(TeX(c("$QR$","$SQR_{\\zeta^*}$","$SQR_{2\\zeta^*}$","$SQR_{4\\zeta^*}$","$SQR_{\\zeta^{*}/4}$","$SQR_{\\zeta^{*}/2}$")))) + 
  # geom_line(aes(y = alphas1_cq8_bias,linetype ='dotdash',color ='SQR_h*8'),size=1.0)+
  # geom_line(aes(y = alphas1_cq16_bias,linetype ='dotdash',color ='SQR_h*16'),size=1.0)+
  # geom_line(aes(y = alphas1_cq_16_bias,linetype ='dotdash',color ='SQR_h/16'),size=1.0)+
  # geom_line(aes(y = alphas1_cq_8_bias,linetype ='dotted',color ='SQR_h/8'),size=1.0)+
  geom_hline(yintercept = 0.0, color = 'yellow', size=1.3) +
  labs(x =xlab,y = ylab)  + 
  scale_x_continuous(breaks=seq(0.01, 0.99, 0.07),limits = c(0.01, 0.99))  + 
  theme_grey() + guides(linetype = FALSE) +
  theme(plot.title = element_text(size = 15),
        legend.title = element_blank(),
        legend.text.align = 0,
        legend.text = element_text(size = 12),
        axis.title = element_text(size = 14))
print(grafico_alphas1_bias)


bias_coefs = data.frame(cbind(thetas1_rq_bias,thetas1_cq_bias,
                              thetas1_cq2_bias,thetas1_cq3_bias,thetas1_cq_3_bias,thetas1_cq_2_bias,
                              thetas1_cq8_bias))#,thetas1_cq16_bias,thetas1_cq_16_bias,thetas1_cq_8_bias))
ylab = TeX('Bias $\\hat{\\theta_1}(\\tau)$')
grafico_thetas1_bias = ggplot(bias_coefs, aes(x = tau.grid, y = value)) + 
  geom_line(aes(y = thetas1_rq_bias,linetype ='blank',color ='QR'),size=1.0) +
  geom_line(aes(y = thetas1_cq_bias,linetype ='solid',color ='SQR'),size=1.0)+
  geom_line(aes(y = thetas1_cq2_bias,linetype ='dotdash',color ='SQR_h*2'),size=1.0)+
  geom_line(aes(y = thetas1_cq3_bias,linetype ='dotdash',color ='SQR_h*4'),size=1.0)+
  geom_line(aes(y = thetas1_cq_3_bias,linetype ='dotdash',color ='SQR_h/4'),size=1.0)+
  geom_line(aes(y = thetas1_cq_2_bias,linetype ='dotted',color ='SQR_h/2'),size=1.0)+
  scale_color_discrete(labels = unname(TeX(c("$QR$","$SQR_{\\zeta^*}$","$SQR_{2\\zeta^*}$","$SQR_{4\\zeta^*}$","$SQR_{\\zeta^{*}/4}$","$SQR_{\\zeta^{*}/2}$")))) + 
  # geom_line(aes(y = thetas1_cq8_bias,linetype ='dotdash',color ='SQR_h*8'),size=1.0)+
  # geom_line(aes(y = thetas1_cq16_bias,linetype ='dotdash',color ='SQR_h*16'),size=1.0)+
  # geom_line(aes(y = thetas1_cq_16_bias,linetype ='dotdash',color ='SQR_h/16'),size=1.0)+
  # geom_line(aes(y = thetas1_cq_8_bias,linetype ='dotted',color ='SQR_h/8'),size=1.0)+
  geom_hline(yintercept = 0.0, color = 'yellow', size=1.3) +
  labs(x =xlab,y = ylab)  + 
  scale_x_continuous(breaks=seq(0.01, 0.99, 0.07),limits = c(0.01, 0.99))  + 
  theme_grey() + guides(linetype = FALSE) +
  theme(plot.title = element_text(size = 15),
        legend.title = element_blank(),
        legend.text.align = 0,
        legend.text = element_text(size = 12),
        axis.title = element_text(size = 14))
print(grafico_thetas1_bias)
#plot_grid(grafico_alphas0_bias,grafico_alphas1_bias,grafico_thetas1_bias,grafico_betas2_bias)

### PLOTS vpq  - VARIANCE per quantile - - h, h*2, h*4, h*1/4, h*1/2

xlab = TeX('$\\tau$')

ylab = TeX('Var $\\hat{\\alpha_0}(\\tau)$') 
grafico_alpha0_vpq = ggplot(coefs_alphas0_vpq, aes(x = tau.grid, y = value)) + 
  geom_line(aes(y = alphas0_rq_vpq,linetype ='blank',color ='QR'),size=1.0) +
  geom_line(aes(y = alphas0_cq_vpq, linetype ='solid',color ='SQR'),size=1.0)+
  geom_line(aes(y = alphas0_cq2_vpq, linetype ='dotdash',color ='SQR_h*2'),size=1.0)+
  geom_line(aes(y = alphas0_cq3_vpq, linetype ='dotdash',color ='SQR_h*4'),size=1.0)+
  geom_line(aes(y = alphas0_cq_3_vpq, linetype ='dotdash',color ='SQR_h/4'),size=1.0)+
  geom_line(aes(y = alphas0_cq_2_vpq, linetype ='dotted',color ='SQR_h/2'),size=1.0)+
  scale_color_discrete(labels = unname(TeX(c("$QR$","$SQR_{\\zeta^*}$","$SQR_{2\\zeta^*}$","$SQR_{4\\zeta^*}$","$SQR_{\\zeta^{*}/4}$","$SQR_{\\zeta^{*}/2}$")))) + 
  # geom_line(aes(y = alphas0_cq8_vpq, linetype ='twodash',color ='SQR_h*8'),size=1.0)+
  # geom_line(aes(y = alphas0_cq16_vpq, linetype ='twodash',color ='SQR_h*16'),size=1.0)+
  # geom_line(aes(y = alphas0_cq_16_vpq, linetype ='twodash',color ='SQR_h/16'),size=1.0)+
  # geom_line(aes(y = alphas0_cq_8_vpq, linetype ='dashed',color ='SQR_h/8'),size=1.0)+
  geom_hline(yintercept = 0, color = 'yellow', size=1.3) +
  labs(x =xlab,y = ylab) + #expand_limits(x=c(0, 1), y=c(-3,3)) +
  scale_x_continuous(breaks=seq(0.01, 0.99, 0.07),limits = c(0.01, 0.99))  + 
  #scale_y_log10() + 
  theme_grey() + guides(linetype = FALSE) +
  theme(plot.title = element_text(size = 15),
        legend.title = element_blank(),
        legend.text.align = 0,
        legend.text = element_text(size = 12),
        axis.title = element_text(size = 14))
print(grafico_alpha0_vpq)

ylab = TeX('Var $\\hat{\\alpha_1}(\\tau)$') 
grafico_alpha1_vpq = ggplot(coefs_alphas1_vpq, aes(x = tau.grid, y = value)) + 
  geom_line(aes(y = alphas1_rq_vpq,linetype ='blank',color ='QR'),size=1.0) +
  geom_line(aes(y = alphas1_cq_vpq, linetype ='solid',color ='SQR'),size=1.0)+
  geom_line(aes(y = alphas1_cq2_vpq, linetype ='dotdash',color ='SQR_h*2'),size=1.0)+
  geom_line(aes(y = alphas1_cq3_vpq, linetype ='dotdash',color ='SQR_h*4'),size=1.0)+
  geom_line(aes(y = alphas1_cq_3_vpq, linetype ='dotdash',color ='SQR_h/4'),size=1.0)+
  geom_line(aes(y = alphas1_cq_2_vpq, linetype ='dotted',color ='SQR_h/2'),size=1.0)+
  scale_color_discrete(labels = unname(TeX(c("$QR$","$SQR_{\\zeta^*}$","$SQR_{2\\zeta^*}$","$SQR_{4\\zeta^*}$","$SQR_{\\zeta^{*}/4}$","$SQR_{\\zeta^{*}/2}$")))) + 
  # geom_line(aes(y = alphas1_cq8_vpq, linetype ='dotdash',color ='SQR_h*8'),size=1.0)+
  # geom_line(aes(y = alphas1_cq16_vpq, linetype ='dotdash',color ='SQR_h*16'),size=1.0)+
  # geom_line(aes(y = alphas1_cq_16_vpq, linetype ='dotdash',color ='SQR_h/16'),size=1.0)+
  # geom_line(aes(y = alphas1_cq_8_vpq, linetype ='dotted',color ='SQR_h/8'),size=1.0)+
  geom_hline(yintercept = 0.0, color = 'yellow', size=1.3) +
  labs(x =xlab,y = ylab) + #expand_limits(x=c(0, 1), y=c(0.4,0.6)) + 
  scale_x_continuous(breaks=round(seq(0.01, 0.99, 0.07), 2),limits = c(0.01, 0.99))  + 
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
grafico_beta1_vpq = ggplot(coefs_thetas1_vpq, aes(x = tau.grid, y = value)) + 
  geom_line(aes(y = thetas1_rq_vpq,linetype ='blank',color ='QR'),size=1.0) +
  geom_line(aes(y = thetas1_cq_vpq, linetype ='solid',color ='SQR'),size=1.,)+
  geom_line(aes(y = thetas1_cq2_vpq, linetype ='dotdash',color ='SQR_h*2'),size=1.0)+
  geom_line(aes(y = thetas1_cq3_vpq, linetype ='dotdash',color ='SQR_h*4'),size=1.0)+
  geom_line(aes(y = thetas1_cq_3_vpq, linetype ='dotdash',color ='SQR_h/4'),size=1.0)+
  geom_line(aes(y = thetas1_cq_2_vpq, linetype ='dotted',color ='SQR_h/2'),size=1.0)+
  scale_color_discrete(labels = unname(TeX(c("$QR$","$SQR_{\\zeta^*}$","$SQR_{2\\zeta^*}$","$SQR_{4\\zeta^*}$","$SQR_{\\zeta^{*}/4}$","$SQR_{\\zeta^{*}/2}$")))) + 
  # geom_line(aes(y = thetas1_cq8_vpq, linetype ='dotdash',color ='SQR_h*8'),size=1.0)+
  # geom_line(aes(y = thetas1_cq16_vpq, linetype ='dotdash',color ='SQR_h*16'),size=1.0)+
  # geom_line(aes(y = thetas1_cq_16_vpq, linetype ='dotdash',color ='SQR_h/16'),size=1.0)+
  # geom_line(aes(y = thetas1_cq_8_vpq, linetype ='dotted',color ='SQR_h/8'),size=1.0)+
  geom_hline(yintercept = 0.0, color = 'yellow', size=1.3) +
  labs(x =xlab,y = ylab) + # ggtitle(bquote(beta[1][tau]^{tau})) + 
  scale_x_continuous(breaks=seq(0.01, 0.99, 0.07),limits = c(0.01, 0.99))  + 
  scale_y_log10() + 
  theme_grey() + guides(linetype = FALSE) +
  theme(plot.title = element_text(size = 15),
        legend.title = element_blank(),
        legend.text.align = 0,
        legend.text = element_text(size = 12),
        axis.title = element_text(size = 14))
print(grafico_beta1_vpq)
#plot_grid(grafico_alpha0_vpq,grafico_alpha1_vpq,grafico_beta1_vpq,grafico_beta2_vpq)

options(scipen=999)
plot_grid(grafico_alphas0_mse,grafico_alphas0_bias,grafico_alpha0_vpq,
          grafico_alphas1_mse,grafico_alphas1_bias,grafico_alpha1_vpq,
          grafico_thetas1_mse,grafico_thetas1_bias,grafico_beta1_vpq,
          nrow = 3, ncol = 3, label_size = 8) 



end_time - start_time
