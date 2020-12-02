# Decision Theory and Bayesian Inference 
# R Project 2
# Marie Ternes, Patrick Finke and Xiaotong Meng

graphics.off()    # Close graphs
rm(list=ls())     # Clean memory
cat("\014")       # Clear Console


# IMPORTANT!!! The working directory needs to be adjusted 
#setwd("~/Documents/Year 5_Master/Decision Theory and Bayesian Inference/Project")

#####################################
############# Task (a) ##############
#####################################

# set the seed and default to the PID
seed_self <<- Sys.getpid()
set.seed_self = function(seed) {
  if (missing(seed)) {
    seed = Sys.getpid()
  }
  seed_self <<- seed
}

# calculate n (pseudo) random samples from U(0,1) using LCG
runif_self = function(n) {
  a <- 171
  b <- 0
  m <- 30269
  
  random = numeric(n)
  for (j in 1:n) {
    seed_self <<- (a * seed_self + b) %% m
    random[j] <- seed_self
  }
  return(random / m)
}

# approximate PI using n U(0,1) samples
# Note: we don't generate in the unit circle but only in the first quadrant of the unit circle.
# However, that does affect the ratio 
approx_pi = function(n) {
  x = runif_self(n)
  y = runif_self(n)
  m = sum(x**2 + y**2 < 1) 
  return(4 * m / n)
}

graphics.off()
#library(plotrix)

temp_func = function(x){
    y=sqrt(1-x^2)
   return(y)
} 

x = seq(0,1, 0.00001) 
plot(x,temp_func(x), xlab='x', ylab='y', type = 'l')
n = 10
x1 = runif_self(n)
y1 = runif_self(n)
points(x1[x1^2+y1^2<=1],y1[x1^2+y1^2<=1],col="green",cex =0.5)
points(x1[x1^2+y1^2>1],y1[x1^2+y1^2>1],col="red", cex=0.5)

#draw.arc(x=0,y=0,angle1=0*pi/180,angle2=90*pi/180 )





task.a = function() {
  set.seed_self(1234) # fix seed
  
  print("Generating 5 random numbers")
  print(runif_self(5))
  print("Testing performance of runif")
  print(system.time(for (i in 1:1000) runif(10000)))
  print("Testing performance of runif_self")
  print(system.time(for (i in 1:1000) runif_self(10000)))
  print("Approximating PI")
  approx = mapply(approx_pi, c(10, 100, 1000, 10000, 100000, 1000000))
  approx.error = abs(approx - pi)
  print(round(approx, 4))
  print(round(approx.error, 4))
  
  # visualization of approximation
  x = seq(100, 100000, 200)
  y = mapply(approx_pi, x)
  plot(x, y, type="l", ylab=bquote("approximation of "*pi), xlab="sample size", ylim =c(3,3.2))
  abline(h=pi, col="red")
}
  
task.a()

#####################################
############# Task (b) ##############
#####################################

#######################################
################ Prior (i) #############
##### All relevant functions for Task b) under Prior (i)######

# PDF of Gamma distribution
my_dgamma <- function(x, alpha, beta){
  fx = (beta^alpha/gamma(alpha))*x^(alpha-1)*exp(-beta*x)
  return(fx)
}

# CDF of Laplace distribution
my_plaplace <- function(p, mu, lambda){
  if(p < mu) Fp = 0.5*exp(lambda*(p-mu))
  else Fp = 1-0.5*exp(-lambda*(p-mu))
  
  return(Fp)
}

# PDF of truncated Laplace distribution
my_dlaplace_trunc <- function(x, mu, lambda){
  p = 1-my_plaplace(0, mu, lambda)       # =Probability(x > 0), where x is Laplace distribution
  fx = (lambda/2)*exp(-lambda*abs(x-mu))/p
  return(fx)
}

# CDF of truncated Laplace distribution
my_plaplace_trunc <- function(p, mu, lambda){
  prob = 1-my_plaplace(0, mu, lambda)   # =Probability(x > 0), where x is Laplace distribution  
  if(p < mu) Fp = (0.5*exp(lambda*(p-mu))-0.5*exp(-lambda*mu))/prob
  else Fp = (1-0.5*exp(-lambda*(p-mu))-0.5*exp(-lambda*mu))/prob
  
  return(Fp)
}

# Randon number generator for truncated Laplace distribution (using inverse density method)
my_rlaplace_trunc <- function(n, mu, lambda) {
  U=runif_self(n)
  prob = 1-my_plaplace(0, mu, lambda)         # =Probability(x > 0), where x is Laplace distribution
  cutoff = (1/2-1/2*exp(-lambda*mu))/prob     # cutoff point for CDF calculations 
  for (j in 1:n) {
    if (U[j]< cutoff) U[j]=mu + log(2*prob*U[j]+exp(-lambda*mu))/lambda
    else U[j]=mu-log(2-2*prob*U[j]-exp(-lambda*mu))/lambda
  }
  return(U)
}


# Random number generator for Gamma density (using rejection method)
rgamma_reject <- function(n, alpha, beta) {
  X=numeric(n)
  M=1.4
  lambda = 120
  Counter = 0             # not necessary since normalizing constant known
  for (j in 1:n) {
    repeat{
      Counter = Counter +1   
      Y = my_rlaplace_trunc(1, mu = (alpha-1)/beta, lambda=lambda)
      U = runif_self(1)
      Mg = M*my_dlaplace_trunc(Y, mu = (alpha-1)/beta, lambda=lambda)
      UMg = U*Mg
      if (UMg<my_dgamma(Y, alpha, beta)) break
    }
    X[j]=Y
  }
  rel_accept = n/Counter
  C = 1/(M*rel_accept)    
  return(list(X=X,C=C,rel_accept=rel_accept))
}



#############################################
##### Prior (i) - actual calculations #######

# Read in cars data file
data = read.csv('cars.csv', sep = "", dec = ",")

# Define variables
Y = data$speed 
n = length(Y)
alpha = 0.55
sumY = sum(Y)

# How to find right M and g
x = seq(0,0.1,0.0001)

# With the help of graphs
plot(x,my_dgamma(x, n*alpha+2, sumY+0.2),col="black", type ='l', ylim=c(0,85) , ylab = 'density',  xlab = bquote(beta*"|y"))
lines(x, my_dlaplace_trunc(x, mu = (n*alpha+2-1)/(sumY+0.2), lambda = 120), col = 'blue')
lines(x, 1.4*my_dlaplace_trunc(x, mu = (n*alpha+2-1)/(sumY+0.2), lambda = 120), col = 'red')
legend('topright', c('f', 'g', 'M*g'), lty=c(1), lwd=c(1,1,1),col=c('black', 'blue', 'red'), cex=0.9)

plot(x,my_dgamma(x, n*alpha+2, sumY+0.2),col="black", type ='l', ylim=c(32,42), xlim = c(0.044,0.046),  ylab = 'density',  xlab = bquote(beta*"|y")) #main = 'Plot of f (black), g (blue), and M*g (red)',
lines(x, 1.4*my_dlaplace_trunc(x, mu = (n*alpha+2-1)/(sumY+0.2), lambda = 120), col = 'red')
legend('topright', c('f', 'M*g'), lty=c(1), lwd=c(1,1,1),col=c('black', 'red'), cex=0.9)

# With the help of calculations
f = my_dgamma(x, n*alpha+2, sumY+0.2)
Mg = 1.4*my_dlaplace_trunc(x, mu = (n*alpha+2-1)/(sumY+0.2), lambda = 120)
sum(f-Mg>0) 


# Applying rejection algorithm to sample from gamma distribution
set.seed_self(2)
nSim = 100000
posterior_reject = rgamma_reject(nSim,n*alpha+2, sumY+0.2)
head(posterior_reject$X)
posterior_reject$C # approx equal to 1 (because integration constant already known)
posterior_reject$rel_accept

# Plots of the Gamma distribution simulated by rejection algorithm
hist(posterior_reject$X, nclass=100, probability = TRUE, xlab = bquote(beta*"|y"), main ='') # main=paste("Histogram of Rejection sample of size n =",nSim)
lines(density(posterior_reject$X), col = 'blue')
lines(x, my_dgamma(x, n*alpha+2, sumY+0.2), col="red")
legend('topright', c('True G(29.5, 755.15) density ', 'Estimated G(29.5, 755.15) density'), lty=c(1), lwd=c(1,1),col=c('red', 'blue'), cex=0.76)



#########################################
################ Prior (ii) #############
##### All relevant functions for Task b) under Prior (ii)######

# Unnormalized PDF of posterior density under prior (ii) given parameters n, alpha, sumY
posterior_ii <- function(x){
  n = 50
  alpha = 0.55
  sumY = 754.95
  val = x^(n*alpha)*exp(-x*sumY)*exp(-0.2*x)*(4*sin(2*x+1)^2+1+(x*exp(-15*(x-15)))/(2*(1+exp(-15*(x-15)))^2))
  return(val)
}

# Random number generator for unknown density (with rejection method)
rdist_reject <- function(nSim, posterior_ii) {
  X=numeric(nSim)

  M = 3.1e-53
  mu = 0.0368
  lambda = 125
  
  Counter = 0            
  for (j in 1:nSim) {
    repeat{
      Counter = Counter +1   
      
      Y = my_rlaplace_trunc(1, mu = mu, lambda=lambda)
      U = runif_self(1)
      Mg = M*my_dlaplace_trunc(Y, mu = mu, lambda=lambda)
      UMg = U*Mg
      if (UMg<posterior_ii(Y)) break
    }
    X[j]=Y
  }
  rel_accept = nSim/Counter
  C = 1/(M*rel_accept)    
  return(list(X=X,C=C,rel_accept=rel_accept))
}

#############################################
##### Prior (i) - actual calculations #######

# How to find right M and g

# With the help of Graphs
M = 3.1e-53
mu = 0.0368
lambda = 125
x=seq(0,0.1,0.0001)
plot(x, posterior_ii(x),type='l', ylim=c(0,2e-51),  ylab = 'density',  xlab = bquote(beta*"|y") )
lines(x,M*my_dlaplace_trunc(x, mu = mu, lambda = lambda),col='red')
legend('topright', c('f', 'M*g'), lty=c(1), lwd=c(1,1,1),col=c('black', 'red'), cex=0.9)

# With the help of Calculations
f = posterior_ii(x)
Mg = M*my_dlaplace_trunc(x, mu = mu, lambda = lambda)
sum(f-Mg >0)


# Applying rejection algorithm to sample from unknown distribution
set.seed_self(2)
nSim = 100000
posterior_ii_reject = rdist_reject(nSim, posterior_ii)
head(posterior_ii_reject$X)
posterior_ii_reject$rel_accept
C = posterior_ii_reject$C
C

# Plots of the posterior distribution simulated by rejection algorithm
hist(posterior_ii_reject$X,breaks=100,probability=TRUE,xlab = bquote(beta*"|y"), main ='', ylim = c(0,60))
lines(density(posterior_ii_reject$X), col ='blue')
lines(x,C*posterior_ii(x),col='red')
legend('topright', c('True density ', 'Estimated density'), lty=c(1), lwd=c(1,1),col=c('red', 'blue'), cex=0.8)



########################################################################
##### Comparison of influence of the different prior distributions #####
hist(posterior_reject$X, nclass=100, probability = TRUE, col=rgb(1,0,0,1/2), xlab = bquote(beta*"|y"), border = F, main ='')
hist(posterior_ii_reject$X,nclass=100,probability=TRUE, col=rgb(0,0,1,1/3), add =TRUE, xlab = bquote(beta*"|y"), border = F)
lines(x, my_dgamma(x, n*alpha+2, sumY+0.2), col="red")
lines(x,C*posterior_ii(x),col='blue')
legend('topright', c('Prior (1)', 'Prior (2)'), lty=c(1), lwd=c(1,1),col=c('red', 'blue'), cex=0.7)



##########################################
############### Task (c) #################
##########################################
Y = data$speed 
n = length(Y)
alpha = 0.55
sumY = sum(Y)

##########################################
############ Prior (i) ###################

# Exact estimation of Bayes Estimator 
bayes_estimator_gamma_linex_loss<- function(c, alpha, beta){
  E = (1-c/beta)^(-alpha) #MGF of a Gamma distribution
  
  d = (1/c)*log(E)  
  return(d)
}

c <- c(-15,-5,5,15) 
# Exact Bayes estimator 
BE_true = bayes_estimator_gamma_linex_loss(c, n*alpha+2, sumY+0.2)
round(BE_true,6)


# Estimation by Monte Carlo Simulation
bayes_estimator_MCgamma_linex_loss <- function(MC.N, c, alpha, beta){
  x = rgamma_reject(MC.N, alpha, beta)
  hx = exp(c*(x$X))
  #MGF_estim = mean(hx)
  MGF_estim = cumsum(hx)/(1:MC.N) # sample mean of e^(c*beta)
  
  d  = 1/c * log(MGF_estim)
  return(d) 
}

set.seed_self(2)
MC.N = 10000

# c = -15
BE1 = bayes_estimator_MCgamma_linex_loss(MC.N, -15, n*alpha+2, sumY+0.2)
round(BE1[MC.N],6)

plot(x = c(1:MC.N), BE1, type ='l', xlab = 'number of MC simulations', ylab = bquote("Bayes estimator of "*beta), main = '')
abline(h = BE_true[1], col ='red')
legend('bottomright', c('MC simulations', 'True value'), lty=c(1,1),col=c('black','red'), cex=0.8)

plot(x = c(1000:MC.N), BE1[c(1000:MC.N)], xlim = c(1000,MC.N), type ='l', xlab = 'number of MC simulations', ylab = bquote("Bayes estimator of "*beta), main = '')
abline(h = BE_true[1], col ='red')

# c = -5
BE2 = bayes_estimator_MCgamma_linex_loss(MC.N, -5, n*alpha+2, sumY+0.2 )
round(BE2[MC.N],6)

plot(x = c(1:MC.N), BE2, type ='l', xlab = 'number of MC simulations', ylab = bquote("Bayes estimator of "*beta), main = '')
abline(h = BE_true[2], col ='red')
legend('bottomright', c('MC simulations', 'True value'), lty=c(1,1),col=c('black','red'), cex=0.8)

plot(x = c(1000:MC.N), BE2[c(1000:MC.N)], xlim = c(1000,MC.N), type ='l', xlab = 'number of MC simulations', ylab = bquote("Bayes estimator of "*beta), main = '')
abline(h = BE_true[2], col ='red')

# c = 5
BE3 = bayes_estimator_MCgamma_linex_loss(MC.N,  5, n*alpha+2, sumY+0.2 )
round(BE3[MC.N],6)

plot(x = c(1:MC.N), BE3, type ='l', xlab = 'number of MC simulations', ylab = bquote("Bayes estimator of "*beta), main = ' ')
abline(h = BE_true[3], col ='red')
legend('bottomright', c('MC simulations', 'True value'), lty=c(1,1),col=c('black','red'), cex=0.8)

plot(x = c(1000:MC.N), BE3[c(1000:MC.N)], xlim = c(1000,MC.N), type ='l', xlab = 'number of MC simulations', ylab = bquote("Bayes estimator of "*beta), main = '')
abline(h = BE_true[3], col ='red')

# c = 15
BE4 = bayes_estimator_MCgamma_linex_loss(MC.N, 15, n*alpha+2, sumY+0.2 )
round(BE4[MC.N],6)

plot(x = c(1:MC.N), BE4, type ='l', xlab = 'number of MC simulations', ylab = bquote("Bayes estimator of "*beta), main = ' ')
abline(h = BE_true[4], col ='red')
legend('bottomright', c('MC simulations', 'True value'), lty=c(1,1),col=c('black','red'), cex=0.8)

plot(x = c(1000:MC.N), BE4[c(1000:MC.N)], xlim = c(1000,MC.N), type ='l', xlab = 'number of MC simulations', ylab = bquote("Bayes estimator of "*beta), main = '')
abline(h = BE_true[4], col ='red')



#####################################
######## Prior (ii) ################

# Estimation by Monte Carlo Simulation
bayes_estimator_MCdist_linex_loss <- function(MC.N, c, posterior_ii){
  x = rdist_reject(MC.N, posterior_ii)
  hx = exp(c*(x$X))
  #MGF_estim = mean(hx)
  #S = sum((hx-MGF_estim)^2)/(MC.N-1)
  MGF_estim = cumsum(hx)/(1:MC.N)                             # sample mean of e^(c*beta)
  var_estim = cumsum((hx-cumsum(hx)/(1:MC.N))^2)/(1:MC.N)^2   # sample variance of e^(c*beta)
  
  d  = 1/c * log(MGF_estim)

  return(list(d = d, mean_estim = MGF_estim , var_estim = var_estim)) 
}


# Estimation by Monte Carlo Simulation
MC.N = 10000
set.seed_self(2)
# c = -15
BE1_ii = bayes_estimator_MCdist_linex_loss(MC.N, -15, posterior_ii)
round(BE1_ii$d[MC.N],6)

# c = -5
BE2_ii = bayes_estimator_MCdist_linex_loss(MC.N, -5, posterior_ii)
round(BE2_ii$d[MC.N],6)

# c = 5
BE3_ii = bayes_estimator_MCdist_linex_loss(MC.N, 5, posterior_ii)
round(BE3_ii$d[MC.N],6)

# c = 15
BE4_ii = bayes_estimator_MCdist_linex_loss(MC.N, 15, posterior_ii)
round(BE4_ii$d[MC.N],6)

BE1_ii$var_estim[MC.N]
BE2_ii$var_estim[MC.N]
BE3_ii$var_estim[MC.N]
BE4_ii$var_estim[MC.N]


# Variance Reduction Method  -> Importance sampling with Gamma Distribution
Importance_Sampling <- function(MC.N, c, IntC, posterior_ii, alpha, beta){
  
  y = rgamma_reject(MC.N, alpha, beta)$X                             # importance distribution Gamma
  weights = IntC*exp(c*y)*posterior_ii(y)/my_dgamma(y, alpha, beta)  # importance weights (h=1)
  
  #mean_weights = mean(weights)
  #S_weights = sum((weights-mean_weights)^2)/(MC.N-1)
  mean_weights = cumsum(weights)/(1:MC.N)                                # sample mean of e^(c*beta)
  var_estim = cumsum((weights-cumsum(weights)/(1:MC.N))^2)/(1:MC.N)^2    # sample variance of e^(c*beta)
  
  d = 1/c * log(mean_weights)                                         # Bayes estimator 
  effective_ss = ((sum(weights))^2)/(sum(weights^2))                  # effective sample size
  
  return(list(d = d, mean = mean_weights, var = var_estim, ess = effective_ss))
}

######################################
######################################
#####################################
#####################################


MC.N = 10000
x=seq(0,0.1,0.0001)
C = posterior_ii_reject$C
set.seed_self(2)

##### Plot of the importance weights for all c
c1 = -15
c2 = -5
c3 = 5
c4 = 15

par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plot(x,(C*exp(c1*x)*posterior_ii(x))/ (my_dgamma(x,28.6, 771)), ylim = c(0,3), type ='l', ylab = 'weights = f/g', main = '', xlab = bquote(beta*"|y"))
lines(x,(C*exp(c2*x)*posterior_ii(x))/ (my_dgamma(x,28.6, 761)), col = 'red')
lines(x,(C*exp(c3*x)*posterior_ii(x))/ (my_dgamma(x,28.6, 751)), col = 'blue')
lines(x,(C*exp(c4*x)*posterior_ii(x))/ (my_dgamma(x,28.6, 741)), col ='darkgreen')
legend("topright", inset=c(-0.32,0), legend = c('c = -15', 'c = - 5', 'c = 5', 'c = 15'), lty=c(1,1,1,1),col=c('black','red', 'blue', 'darkgreen'), cex=0.9)
par(mar=c(5, 4, 4, 2) + 0.1, xpd =FALSE)

# Plots of the target distribution and importance distribution for different c's and applying the importance sampling algorithm

##### c = -15
# Plots
c1 = -15
plot(x, C*exp(c1*x)*posterior_ii(x),type='l',ylab="density",main="", ylim = c(0,60), xlab = bquote(beta*"|y"))
lines(x, my_dgamma(x,28.6, 771), col = 'green')
legend('topright', legend = c('f', 'g ~ G(28.6,771)'), lty=c(1,1),col=c('black','green'), cex=0.9)

# Applying Importance Sampling
BE_importSamp1 = Importance_Sampling(MC.N, c = -15, C, posterior_ii, alpha = 28.6, beta = 771)
round(BE_importSamp1$d[MC.N],6)
BE_importSamp1$var[MC.N]
BE_importSamp1$ess


####### c = -5
# Plots
c1 = -5
plot(x, C*exp(c1*x)*posterior_ii(x),type='l',ylab="density",main="", ylim = c(0,60), xlab = bquote(beta*"|y"))
lines(x, my_dgamma(x,28.6, 761), col = 'green')
legend('topright', legend = c('f', 'g ~ G(28.6,761)'), lty=c(1,1),col=c('black','green'), cex=0.9)


# Applying Importance Sampling
BE_importSamp2 = Importance_Sampling(MC.N, c = -5, C, posterior_ii, alpha = 28.6, beta = 761)
round(BE_importSamp2$d[MC.N],6)
BE_importSamp2$var[MC.N]
BE_importSamp2$ess


##### c = 5
# Plots
c1 = 5
plot(x, C*exp(c1*x)*posterior_ii(x),type='l',ylab="density",main="", ylim = c(0,70), xlab = bquote(beta*"|y"))
lines(x, my_dgamma(x,28.6, 751), col = 'green')
legend('topright', legend = c('f', 'g ~ G(28.6,751)'), lty=c(1,1),col=c('black','green'), cex=0.9)


# Applying Importance Sampling
BE_importSamp3 = Importance_Sampling(MC.N, c = 5, C, posterior_ii, alpha = 28.6, beta = 751)
round(BE_importSamp3$d[MC.N],6)
BE_importSamp3$var[MC.N]
BE_importSamp3$ess

##### c = 15
# Plots
c1 = 15
plot(x, C*exp(c1*x)*posterior_ii(x),type='l',ylab="density",main="", xlab = bquote(beta*"|y"), )
lines(x, my_dgamma(x,28.6, 741), col = 'green')
legend('topright', legend = c('f', 'g ~ G(28.6,741)'), lty=c(1,1),col=c('black','green'), cex=0.9)


# Applying Importance Sampling
BE_importSamp4 = Importance_Sampling(MC.N, c = 15, C, posterior_ii, alpha = 28.6, beta = 741)
round(BE_importSamp4$d[MC.N],6)
BE_importSamp4$var[MC.N]
BE_importSamp4$ess

#####################
# Variance + Variance Reduction
BE_importSamp1$var[MC.N]
BE_importSamp2$var[MC.N]
BE_importSamp3$var[MC.N]
BE_importSamp4$var[MC.N]

round(100*((BE1_ii$var_estim[MC.N]-BE_importSamp1$var[MC.N])/BE1_ii$var_estim[MC.N]),4)
round(100*((BE2_ii$var_estim[MC.N]-BE_importSamp2$var[MC.N])/BE2_ii$var_estim[MC.N]),4)
round(100*((BE3_ii$var_estim[MC.N]-BE_importSamp3$var[MC.N])/BE3_ii$var_estim[MC.N]),4)
round(100*((BE4_ii$var_estim[MC.N]-BE_importSamp4$var[MC.N])/BE4_ii$var_estim[MC.N]),4)

# Difference in Bayes estimator MC and Importance Sampling
round(BE1_ii$d[MC.N]-BE_importSamp1$d[MC.N],6)
round(BE2_ii$d[MC.N]-BE_importSamp2$d[MC.N],6)
round(BE3_ii$d[MC.N]-BE_importSamp3$d[MC.N],6)
round(BE4_ii$d[MC.N]-BE_importSamp4$d[MC.N],6)



##### Comparison Plots between MC integration and Importance Sampling

## c = -15
#### Plot of running bayes estimator:
plot(x = c(1:MC.N), y = BE1_ii$d, type='l', ylim = c(0.029,0.039), xlab="number of simulations",ylab=bquote("Bayes estimator of "*beta), 
     main="")
lines(x = c(1:MC.N), y = BE_importSamp1$d, type='l',ylab="",col="red")
legend('bottomright', legend = c('MC integration', 'Importance Sampling'), lty=c(1,1),col=c('black','red'), cex=0.9)


#### Plot of estimated variances: 
plot(x = c(1:MC.N), y = BE1_ii$var_estim, type='l',xlab="number of simulations",ylab="Estimated Variance",
     main="")#, xlim = c(2000,MC.N), ylim = c(0,0.0001))
lines(x = c(1:MC.N), y = BE_importSamp1$var,type='l',col="red")
legend('topright', legend = c('MC integration', 'Importance Sampling'), lty=c(1,1),col=c('black','red'), cex=0.9)


plot(x = c(1000:MC.N), y = BE1_ii$var_estim[1000:MC.N], type='l',xlab="number of simulations",ylab="Estimated Variance",
     main="", ylim = c(0,3.5e-06))#, xlim = c(2000,MC.N), ylim = c(0,0.0001))
lines(x = c(1000:MC.N), y = BE_importSamp1$var[1000:MC.N],type='l',col="red")


############################
## c = -5
plot(x = c(1:MC.N), y = BE2_ii$d, type='l',xlab="number of simulations",ylab=bquote("Bayes estimator of "*beta) )
lines(x = c(1:MC.N), y = BE_importSamp2$d, type='l',xlab="Iterations",ylab="",col="red")
legend('topright', legend = c('MC integration', 'Importance Sampling'), lty=c(1,1),col=c('black','red'), cex=0.9)


#### Plot of estimated variances: 
plot(x = c(1:MC.N), y = BE2_ii$var_estim, type='l',xlab="number of simulations",ylab="Estimated Variance")#, xlim = c(2000,MC.N), ylim = c(0,0.0001))
lines(x = c(1:MC.N), y = BE_importSamp2$var,type='l',col="red")
legend('topright', legend = c('MC integration', 'Importance Sampling'), lty=c(1,1),col=c('black','red'), cex=0.9)

plot(x = c(1000:MC.N), y = BE2_ii$var_estim[1000:MC.N], type='l',xlab="number of simulations",ylab="Estimated Variance",
     main="", ylim = c(0,1.1e-06))#, xlim = c(2000,MC.N), ylim = c(0,0.0001))
lines(x = c(1000:MC.N), y = BE_importSamp2$var[1000:MC.N],type='l',col="red")


######################
## c = 5
plot(x = c(1:MC.N), y = BE3_ii$d, type='l',xlab="number of simulations",ylab=bquote("Bayes estimator of "*beta))
lines(x = c(1:MC.N), y = BE_importSamp3$d, type='l',ylab="",col="red")
legend('topright', legend = c('MC integration', 'Importance Sampling'), lty=c(1,1),col=c('black','red'), cex=0.9)

#### Plot of estimated variances: 
plot(x = c(1:MC.N), y = BE3_ii$var_estim, type='l',xlab="number of simulations",ylab="Estimated Variance")#, xlim = c(2000,MC.N), ylim = c(0,0.0001))
lines(x = c(1:MC.N), y = BE_importSamp3$var,type='l',col="red")
legend('topright', legend = c('MC integration', 'Importance Sampling'), lty=c(1,1),col=c('black','red'), cex=0.9)

plot(x = c(1000:MC.N), y = BE3_ii$var_estim[1000:MC.N], type='l',xlab="number of simulations",ylab="Estimated Variance",
     main="", ylim = c(0,2e-06))#, xlim = c(2000,MC.N), ylim = c(0,0.0001))
lines(x = c(1000:MC.N), y = BE_importSamp3$var[1000:MC.N],type='l',col="red")

######################
## c = 15
plot(x = c(1:MC.N), y = BE4_ii$d, type='l',xlab="number of simulations",ylab=bquote("Bayes estimator of "*beta))
lines(x = c(1:MC.N), y = BE_importSamp4$d, type='l',ylab="",col="red")
legend('bottomright', legend = c('MC integration', 'Importance Sampling'), lty=c(1,1),col=c('black','red'), cex=0.9)

#### Plot of estimated variances: 
plot(x = c(1:MC.N), y = BE4_ii$var_estim, type='l',xlab="number of simulations",ylab="Estimated Variance")#, xlim = c(2000,MC.N), ylim = c(0,0.0001))
lines(x = c(1:MC.N), y = BE_importSamp4$var,type='l',col="red")
legend('topright', legend = c('MC integration', 'Importance Sampling'), lty=c(1,1),col=c('black','red'), cex=0.9)

plot(x = c(1000:MC.N), y = BE4_ii$var_estim[1000:MC.N], type='l',xlab="number of simulations",ylab="Estimated Variance",
     ylim = c(0,4e-05))#, xlim = c(2000,MC.N), ylim = c(0,0.0001))
lines(x = c(1000:MC.N), y = BE_importSamp4$var[1000:MC.N],type='l',col="red")



##########################################
############### Task (d) #################
##########################################

# Calculate (1-alpha)100% equal tail credibility interval based on samples
credibility.interval = function(samples, alpha) {
  samples = sort(samples)
  m = alpha/2 * length(samples)
  return(c(samples[m+1], samples[length(samples)-m+1])) #samples[m+1], samples[length(samples)-m+1]
}


#### Prior (i) #####
Y = data$speed
n = length(Y)
alpha = 0.55
sumY = sum(Y)
x = seq(0,0.1,0.0001)

set.seed_self(2)
nSim = 10000
y = rgamma_reject(nSim, n*alpha+2, sumY+0.2)$X
credInt_gamma = credibility.interval(y, 0.1); credInt_gamma
plot(x, my_dgamma(x, n*alpha+2, sumY+0.2), type = 'l', ylab = '', xlab = bquote(beta*"|y")) 
abline(v = credInt_gamma[1], col = 'blue')
abline(v = credInt_gamma[2], col = 'blue')
abline(v = BE_true[1], col = 'orange', lty = 1)
abline(v = BE_true[2], col = 'red', lty = 2)
abline(v = BE_true[3], col = 'darkgreen', lty = 3)
abline(v = BE_true[4], col = 'green', lty = 4)
legend('topright', c('90% Credible Int', 'c = -15', 'c = -5', 'c = 5', 'c = 15'), lty=c(1,1,2,3,4),col=c('blue','orange','red', 'darkgreen', 'green'), cex=0.7)


#### Prior (ii) ####
C = posterior_ii_reject$C
y1 = rdist_reject(nSim, posterior_ii)$X
credInt_dist = credibility.interval(y1, 0.1)
plot(x, C*posterior_ii(x), type = 'l', ylab = '', xlab = bquote(beta*"|y")) 
abline(v = credInt_dist[1], col = 'blue')
abline(v = credInt_dist[2], col = 'blue')
abline(v = BE_importSamp1$d[MC.N], col = 'orange', lty = 1)
abline(v = BE_importSamp2$d[MC.N], col = 'red', lty = 2)
abline(v = BE_importSamp3$d[MC.N], col = 'darkgreen', lty = 3)
abline(v = BE_importSamp4$d[MC.N], col = 'green', lty = 4)
legend('topright', c('90% Credible Int', 'c = -15', 'c = -5', 'c = 5', 'c = 15'), lty=c(1,1,2,3,4),col=c('blue','orange','red', 'darkgreen', 'green'), cex=0.7)







