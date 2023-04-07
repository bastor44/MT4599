# MT4599 Chapter 5: EM (simulations)

setwd('~/MT4599')
library(missMethods)
library(Metrics)
library(rms)
library(em)
library(mixtools)

wine <- read.table('wine.data', dec='.', sep=',')
colnames(wine) <- c('class', 'alcohol', 'malic acid', 'ash', 'alkalinity of ash',
                    'magnesium', 'total phenols', 'flavanoids', 'nonflavanoid phenols', 
                    'proanthocyanins', 'color intensity', 'hue', 'OD2800/OD315', 
                    'proline')

############################ COMPLETE DATA ##################################
x <- wine$alcohol
mean(x) #1 3.00
var(x) # 0.66

plot(density(x))
hist(x, breaks=30)
# initial parameter estimates for normal mixture model: N(12.5, 0.5), N(13.7, 0.5)
EM_comp <- normalmixEM(x, k=2, lambda=c(0.5, 0.5), mu=c(12.5, 13.7), 
                       sigma=c(0.5, 0.5))
# maximised parameter estimates
EM_comp$mu
EM_comp$sigma
EM_comp$lambda

# calculate new overall mean and variance 
EM_comp_mu <- EM_comp$lambda[1] * EM_comp$mu[1] + EM_comp$lambda[2]*EM_comp$mu[2]
EM_comp_mu # 13.00

EM_comp_var <- EM_comp$lambda[1]*(EM_comp$sigma[1]^2 + EM_comp$mu[1]^2) + 
               EM_comp$lambda[2]*(EM_comp$sigma[2]^2 + EM_comp$mu[2]^2) -
               (EM_comp$lambda[1]*EM_comp$mu[1] + EM_comp$lambda[2]*EM_comp$mu[2])^2
EM_comp_var # 0.66

# function to calculate the overall mean of a normal mixture model 
emmean <- function(emobject, k) {
  mean <- 0
  for (i in 1:k) {
    mean <- mean + emobject$lambda[i]*emobject$mu[i]
  }
  return(mean)
}
emmean(EM_comp, 2) # 13.00

# function to calculate the overall variance of a normal mixture model 
emvar <- function(emobject, k) {
  var <- 0
  minus <- 0 
  for (i in 1:k) {
    var <- var + emobject$lambda[i]*(emobject$sigma[i]^2 + emobject$mu[i]^2)
    minus <- minus + emobject$lambda[i]*emobject$mu[i]
  }
  var <- var - minus^2
  return(var)
}
emvar(EM_comp, 2) # 0.66

# function to run a t-test by hand (without data)
pttest <- function(mu1, mu2, sigma1, sigma2, d1, d2) {
  n1 <- length(d1)
  n2 <- length(d2)
  df <- n1 + n2 - 2
  se <- sqrt((sigma1^2 / n1) + (sigma2^2 / n2))
  t <- (mu1-mu2)/se
  pval <- pt(t, df)
  return(pval)
}

pFtest <- function(d1, d2, sigma1, sigma2) {
  df1 <- length(d1) - 1
  df2 <- length(d2) - 1
  F <- sigma1^2 / sigma2^2
  pval <- pf(F, df1, df2)
  return(pval)
}


################################ MCAR ########################################
# one iteration
# MCAR, p=0.1
# plot one simulation to construct an initial model 
mcar_1 <- delete_MCAR(ds=wine, p=0.1, cols_mis='alcohol')
mean(mcar_1$alcohol, na.rm=TRUE) # 12.99
plot(density(mcar_1$alcohol, na.rm=TRUE)) # looks like a normal mixture 
hist(mcar_1$alcohol, breaks=30)

x_mcar_1 <- mcar_1$alcohol[!is.na(mcar_1$alcohol)]

# initial parameters: N(12.4, 0.5) & N(13.7, 0.5)
EM_mcar_1 <- normalmixEM(x_mcar_1, k=2, lambda=c(0.5, 0.5), mu=c(12.4, 13.7), 
                         sigma=c(0.5, 0.5)) # 121 iterations 
# new parameter estimates for normal mixture model:
EM_mcar_1$mu
EM_mcar_1$sigma
EM_mcar_1$lambda

# overall mean and variance
emmean(EM_mcar_1, 2) # 13.02
emvar(EM_mcar_1, 2) # 0.66


# 100 simulations of MCAR, p=0.1
nsim <- 100
means <- rep(0, nsim)
vars <- rep(0, nsim)
ts <- rep(0, nsim)
Fs <- rep(0, nsim)
s <- 1646 

for (i in 1:nsim) {
  mcar_1 <- delete_MCAR(ds=wine, p=0.1, cols_mis='alcohol')
  x_mcar_1 <- mcar_1$alcohol[!is.na(mcar_1$alcohol)]
  
  # initial parameter estimates: N(12.4, 0.5), N(13.7, 0.5), equally weighted
  EM_mcar_1 <- normalmixEM(x_mcar_1, k=2, lambda=c(0.5, 0.5), mu=c(12.4, 13.7),
                           sigma=c(0.5, 0.5)) 
  
  # overall means and variances
  means[i] <- emmean(EM_mcar_1, 2)
  vars[i] <- emvar(EM_mcar_1, 2)
  
  # t-tests and F-tests
  ts[i] <- pttest(mean(wine$alcohol), means[i], sqrt(var(wine$alcohol)), sqrt(vars[i]), wine$alcohol, x_mcar_1)
  Fs[i] <- pFtest(wine$alcohol, x_mcar_1, sqrt(var(wine$alcohol)), sqrt(vars[i]))
  s <- s+1
}
length(ts[ts<=0.05]) # 0
length(Fs[Fs<=0.05]) # 0 


# MCAR, p=0.3
# plot once to determine initial parameter estimates
mcar_3 <- delete_MCAR(ds=wine, p=0.3, cols_mis='alcohol')
plot(density(mcar_3$alcohol, na.rm=TRUE))
# N(12.2, 0.5), N(13.1, 1)

nsim <- 100
means <- rep(0, nsim)
vars <- rep(0, nsim)
ts <- rep(0, nsim)
Fs <- rep(0, nsim)
s <- 1646 

for (i in 1:nsim) {
  mcar_3 <- delete_MCAR(ds=wine, p=0.3, cols_mis='alcohol')
  x_mcar_3 <- mcar_3$alcohol[!is.na(mcar_3$alcohol)]
  
  # initial parameter estimates: N(12.4, 0.5), N(13.7, 0.5), equally weighted
  EM_mcar_3 <- normalmixEM(x_mcar_3, k=2, lambda=c(0.5, 0.5), mu=c(12.2, 13.1),
                           sigma=c(0.5, 1)) 
  
  # overall means and variances
  means[i] <- emmean(EM_mcar_3, 2)
  vars[i] <- emvar(EM_mcar_3, 2)
  
  # t-tests and F-tests
  ts[i] <- pttest(mean(wine$alcohol), means[i], sqrt(var(wine$alcohol)), sqrt(vars[i]), wine$alcohol, x_mcar_3)
  Fs[i] <- pFtest(wine$alcohol, x_mcar_3, sqrt(var(wine$alcohol)), sqrt(vars[i]))
  s <- s+1
}
length(ts[ts<=0.05]) # 0
length(Fs[Fs<=0.05]) # 0


# MCAR, p=0.5
# plot once to determine initial parameter estimates
mcar_5 <- delete_MCAR(ds=wine, p=0.5, cols_mis='alcohol')
plot(density(mcar_5$alcohol, na.rm=TRUE))
# N(12.5, 0.5), N(13.7, 0.5)

nsim <- 100
means <- rep(0, nsim)
vars <- rep(0, nsim)
ts <- rep(0, nsim)
Fs <- rep(0, nsim)
s <- 1646 

for (i in 1:nsim) {
  mcar_5 <- delete_MCAR(ds=wine, p=0.5, cols_mis='alcohol')
  x_mcar_5 <- mcar_5$alcohol[!is.na(mcar_5$alcohol)]
  
  # initial parameter estimates: N(12.4, 0.5), N(13.7, 0.5), equally weighted
  EM_mcar_5 <- normalmixEM(x_mcar_5, k=2, lambda=c(0.5, 0.5), mu=c(12.5, 13.7),
                           sigma=c(0.5, 0.5)) 
  
  # overall means and variances
  means[i] <- emmean(EM_mcar_5, 2)
  vars[i] <- emvar(EM_mcar_5, 2)
  
  # t-tests and F-tests
  ts[i] <- pttest(mean(wine$alcohol), means[i], sqrt(var(wine$alcohol)), sqrt(vars[i]), wine$alcohol, x_mcar_3)
  Fs[i] <- pFtest(wine$alcohol, x_mcar_3, sqrt(var(wine$alcohol)), sqrt(vars[i]))
  s <- s+1
}
length(ts[ts<=0.05]) # 0
length(Fs[Fs<=0.05]) # 1


# MCAR, p=0.7
# plot once to determine initial parameter estimates
mcar_7 <- delete_MCAR(ds=wine, p=0.7, cols_mis='alcohol')
plot(density(mcar_7$alcohol, na.rm=TRUE))
# N(12.5, 0.5), N(13.5, 0.5)

nsim <- 100
means <- rep(0, nsim)
vars <- rep(0, nsim)
ts <- rep(0, nsim)
Fs <- rep(0, nsim)
s <- 1646 

for (i in 1:nsim) {
  mcar_7 <- delete_MCAR(ds=wine, p=0.7, cols_mis='alcohol')
  x_mcar_7 <- mcar_7$alcohol[!is.na(mcar_7$alcohol)]
  
  # initial parameter estimates: N(12.4, 0.5), N(13.7, 0.5), equally weighted
  EM_mcar_7 <- normalmixEM(x_mcar_7, k=2, lambda=c(0.5, 0.5), mu=c(12.5, 13.5),
                           sigma=c(0.5, 0.5)) 
  
  # overall means and variances
  means[i] <- emmean(EM_mcar_7, 2)
  vars[i] <- emvar(EM_mcar_7, 2)
  
  # t-tests and F-tests
  ts[i] <- pttest(mean(wine$alcohol), means[i], sqrt(var(wine$alcohol)), sqrt(vars[i]), wine$alcohol, x_mcar_3)
  Fs[i] <- pFtest(wine$alcohol, x_mcar_3, sqrt(var(wine$alcohol)), sqrt(vars[i]))
  s <- s+1
}
length(ts[ts<=0.05]) # 6
length(Fs[Fs<=0.05]) # 0



################################### MAR ######################################
# MAR, p=0.1
mar_1 <- delete_MAR_censoring(ds=wine, p=0.1, cols_mis='alcohol', cols_ctrl='ash')
plot(density(mar_1$alcohol, na.rm=TRUE))
# N(12.3, 0.5), N(13.7, 0.5)

nsim <- 100
means <- rep(0, nsim)
vars <- rep(0, nsim)
ts <- rep(0, nsim)
Fs <- rep(0, nsim)
s <- 1646 

for (i in 1:nsim) {
  mar_1 <- delete_MAR_censoring(ds=wine, p=0.1, cols_mis='alcohol', cols_ctrl='ash')
  x_mar_1 <- mar_1$alcohol[!is.na(mar_1$alcohol)]
  
  # initial parameter estimates: N(12.4, 0.5), N(13.7, 0.5), equally weighted
  EM_mar_1 <- normalmixEM(x_mar_1, k=2, lambda=c(0.5, 0.5), mu=c(12.3, 13.7),
                           sigma=c(0.5, 0.5)) 
  
  # overall means and variances
  means[i] <- emmean(EM_mar_1, 2)
  vars[i] <- emvar(EM_mar_1, 2)
  
  # t-tests and F-tests
  ts[i] <- pttest(mean(wine$alcohol), means[i], sqrt(var(wine$alcohol)), sqrt(vars[i]), wine$alcohol, x_mar_1)
  Fs[i] <- pFtest(wine$alcohol, x_mar_1, sqrt(var(wine$alcohol)), sqrt(vars[i]))
  s <- s+1
}
length(ts[ts<=0.05]) # 0
length(Fs[Fs<=0.05]) # 0 


# MAR, p=0.3
mar_3 <- delete_MAR_censoring(ds=wine, p=0.3, cols_mis='alcohol', cols_ctrl='ash')
plot(density(mar_3$alcohol, na.rm=TRUE))
# N(12.0, 1), N(13.9, 0.7)

nsim <- 100
means <- rep(0, nsim)
vars <- rep(0, nsim)
ts <- rep(0, nsim)
Fs <- rep(0, nsim)
s <- 1646 

for (i in 1:nsim) {
  mar_3 <- delete_MAR_censoring(ds=wine, p=0.3, cols_mis='alcohol', cols_ctrl='ash')
  x_mar_3 <- mar_3$alcohol[!is.na(mar_3$alcohol)]
  
  # initial parameter estimates
  EM_mar_3 <- normalmixEM(x_mar_3, k=2, lambda=c(0.5, 0.5), mu=c(12.0, 13.9),
                          sigma=c(1, 0.7)) 
  
  # overall means and variances
  means[i] <- emmean(EM_mar_3, 2)
  vars[i] <- emvar(EM_mar_3, 2)
  
  # t-tests and F-tests
  ts[i] <- pttest(mean(wine$alcohol), means[i], sqrt(var(wine$alcohol)), sqrt(vars[i]), wine$alcohol, x_mar_3)
  Fs[i] <- pFtest(wine$alcohol, x_mar_3, sqrt(var(wine$alcohol)), sqrt(vars[i]))
  s <- s+1
}
length(ts[ts<=0.05]) # 0
length(Fs[Fs<=0.05]) # 0 


# MAR, p=0.5
mar_5 <- delete_MAR_censoring(ds=wine, p=0.5, cols_mis='alcohol', cols_ctrl='ash')
plot(density(mar_5$alcohol, na.rm=TRUE))
# N(11.9, 1), N(13.9, 1), weighted 0.3, 0.7

nsim <- 100
means <- rep(0, nsim)
vars <- rep(0, nsim)
ts <- rep(0, nsim)
Fs <- rep(0, nsim)
s <- 1646 

for (i in 1:nsim) {
  mar_5 <- delete_MAR_censoring(ds=wine, p=0.5, cols_mis='alcohol', cols_ctrl='ash')
  x_mar_5 <- mar_5$alcohol[!is.na(mar_5$alcohol)]
  
  # initial parameter estimates
  EM_mar_5 <- normalmixEM(x_mar_5, k=2, lambda=c(0.3, 0.7), mu=c(11.9, 13.9),
                          sigma=c(1, 1)) 
  
  # overall means and variances
  means[i] <- emmean(EM_mar_5, 2)
  vars[i] <- emvar(EM_mar_5, 2)
  
  # t-tests and F-tests
  ts[i] <- pttest(mean(wine$alcohol), means[i], sqrt(var(wine$alcohol)), sqrt(vars[i]), wine$alcohol, x_mar_5)
  Fs[i] <- pFtest(wine$alcohol, x_mar_5, sqrt(var(wine$alcohol)), sqrt(vars[i]))
  s <- s+1
}
length(ts[ts<=0.05]) # 0
length(Fs[Fs<=0.05]) # 0 


# MAR, p=0.7
mar_7 <- delete_MAR_censoring(ds=wine, p=0.7, cols_mis='alcohol', cols_ctrl='ash')
plot(density(mar_7$alcohol, na.rm=TRUE))
# N(12.0, 1), N(13.7, 0.5), weighted 0.3, 0.7

nsim <- 100
means <- rep(0, nsim)
vars <- rep(0, nsim)
ts <- rep(0, nsim)
Fs <- rep(0, nsim)
s <- 1646 

for (i in 1:nsim) {
  mar_7 <- delete_MAR_censoring(ds=wine, p=0.7, cols_mis='alcohol', cols_ctrl='ash')
  x_mar_7 <- mar_5$alcohol[!is.na(mar_7$alcohol)]
  
  # initial parameter estimates
  EM_mar_7 <- normalmixEM(x_mar_7, k=2, lambda=c(0.3, 0.7), mu=c(12.0, 13.7),
                          sigma=c(1, 0.5)) 
  
  # overall means and variances
  means[i] <- emmean(EM_mar_7, 2)
  vars[i] <- emvar(EM_mar_7, 2)
  
  # t-tests and F-tests
  ts[i] <- pttest(mean(wine$alcohol), means[i], sqrt(var(wine$alcohol)), sqrt(vars[i]), wine$alcohol, x_mar_5)
  Fs[i] <- pFtest(wine$alcohol, x_mar_7, sqrt(var(wine$alcohol)), sqrt(vars[i]))
  s <- s+1
}
length(ts[ts<=0.05]) # 0
length(Fs[Fs<=0.05]) # 0 



################################## MNAR ######################################
# MNAR, p=0.1
mnar_1 <- delete_MNAR_censoring(ds=wine, p=0.1, cols_mis='alcohol')
plot(density(mnar_1$alcohol, na.rm=TRUE))
# N(12.3, 0.3), N(13.7, 0.5)

nsim <- 100
means <- rep(0, nsim)
vars <- rep(0, nsim)
ts <- rep(0, nsim)
Fs <- rep(0, nsim)
s <- 1646 

for (i in 1:nsim) {
  mnar_1 <- delete_MNAR_censoring(ds=wine, p=0.1, cols_mis='alcohol')
  x_mnar_1 <- mnar_1$alcohol[!is.na(mnar_1$alcohol)]
  
  # initial parameter estimates: N(12.4, 0.5), N(13.7, 0.5), equally weighted
  EM_mnar_1 <- normalmixEM(x_mnar_1, k=2, lambda=c(0.5, 0.5), mu=c(12.3, 13.7),
                           sigma=c(0.3, 0.5)) 
  
  # overall means and variances
  means[i] <- emmean(EM_mnar_1, 2)
  vars[i] <- emvar(EM_mnar_1, 2)
  
  # t-tests and F-tests
  ts[i] <- pttest(mean(wine$alcohol), means[i], sqrt(var(wine$alcohol)), sqrt(vars[i]), wine$alcohol, x_mnar_1)
  Fs[i] <- pFtest(wine$alcohol, x_mnar_1, sqrt(var(wine$alcohol)), sqrt(vars[i]))
  s <- s+1
}
length(ts[ts<=0.05]) # 100
length(Fs[Fs<=0.05]) # 0 


# MNAR, p=0.3
mnar_3 <- delete_MNAR_censoring(ds=wine, p=0.3, cols_mis='alcohol')
plot(density(mnar_3$alcohol, na.rm=TRUE))
# N(13, 2), N(13.7, 0.3)

nsim <- 100
means <- rep(0, nsim)
vars <- rep(0, nsim)
ts <- rep(0, nsim)
Fs <- rep(0, nsim)
s <- 1646 

for (i in 1:nsim) {
  mnar_3 <- delete_MNAR_censoring(ds=wine, p=0.3, cols_mis='alcohol')
  x_mnar_3 <- mnar_3$alcohol[!is.na(mnar_3$alcohol)]
  
  # initial parameter estimates: N(12.4, 0.5), N(13.7, 0.5), equally weighted
  EM_mnar_3 <- normalmixEM(x_mnar_3, k=2, lambda=c(0.5, 0.5), mu=c(13, 13.7),
                           sigma=c(2, 0.3)) 
  
  # overall means and variances
  means[i] <- emmean(EM_mnar_3, 2)
  vars[i] <- emvar(EM_mnar_3, 2)
  
  # t-tests and F-tests
  ts[i] <- pttest(mean(wine$alcohol), means[i], sqrt(var(wine$alcohol)), sqrt(vars[i]), wine$alcohol, x_mnar_3)
  Fs[i] <- pFtest(wine$alcohol, x_mnar_3, sqrt(var(wine$alcohol)), sqrt(vars[i]))
  s <- s+1
}
length(ts[ts<=0.05]) # 100
length(Fs[Fs<=0.05]) # 0 


# MNAR, p=0.5
mnar_5 <- delete_MNAR_censoring(ds=wine, p=0.5, cols_mis='alcohol')
plot(density(mnar_5$alcohol, na.rm=TRUE))
# N(13.3, 0.5), N(13.7, 0.3)

nsim <- 100
means <- rep(0, nsim)
vars <- rep(0, nsim)
ts <- rep(0, nsim)
Fs <- rep(0, nsim)
s <- 1646 

for (i in 1:nsim) {
  mnar_5 <- delete_MNAR_censoring(ds=wine, p=0.5, cols_mis='alcohol')
  x_mnar_5 <- mnar_5$alcohol[!is.na(mnar_5$alcohol)]
  
  # initial parameter estimates: N(12.4, 0.5), N(13.7, 0.5), equally weighted
  EM_mnar_5 <- normalmixEM(x_mnar_5, k=2, lambda=c(0.5, 0.5), mu=c(13.3, 13.7),
                           sigma=c(0.5, 0.3)) 
  
  # overall means and variances
  means[i] <- emmean(EM_mnar_5, 2)
  vars[i] <- emvar(EM_mnar_5, 2)
  
  # t-tests and F-tests
  ts[i] <- pttest(mean(wine$alcohol), means[i], sqrt(var(wine$alcohol)), sqrt(vars[i]), wine$alcohol, x_mnar_5)
  Fs[i] <- pFtest(wine$alcohol, x_mnar_5, sqrt(var(wine$alcohol)), sqrt(vars[i]))
  s <- s+1
}
length(ts[ts<=0.05]) # 100
length(Fs[Fs<=0.05]) # 0 


# MNAR, p=0.7
mnar_7 <- delete_MNAR_censoring(ds=wine, p=0.7, cols_mis='alcohol')
plot(density(mnar_7$alcohol, na.rm=TRUE))
# N(13.7, 0.5), N(14.3, 0.5)

nsim <- 100
means <- rep(0, nsim)
vars <- rep(0, nsim)
ts <- rep(0, nsim)
Fs <- rep(0, nsim)
s <- 1646 

for (i in 1:nsim) {
  mnar_7 <- delete_MNAR_censoring(ds=wine, p=0.7, cols_mis='alcohol')
  x_mnar_7 <- mnar_7$alcohol[!is.na(mnar_7$alcohol)]
  
  # initial parameter estimates: N(12.4, 0.5), N(13.7, 0.5), equally weighted
  EM_mnar_7 <- normalmixEM(x_mnar_7, k=2, lambda=c(0.5, 0.5), mu=c(13.5, 14.3),
                           sigma=c(2, 0.3)) 
  
  # overall means and variances
  means[i] <- emmean(EM_mnar_7, 2)
  vars[i] <- emvar(EM_mnar_7, 2)
  
  # t-tests and F-tests
  ts[i] <- pttest(mean(wine$alcohol), means[i], sqrt(var(wine$alcohol)), sqrt(vars[i]), wine$alcohol, x_mnar_7)
  Fs[i] <- pFtest(wine$alcohol, x_mnar_7, sqrt(var(wine$alcohol)), sqrt(vars[i]))
  s <- s+1
}
length(ts[ts<=0.05]) # 100
length(Fs[Fs<=0.05]) # 0 



######################### illustrations #######################################
mcar_3
set.seed(1157)
em_ill <- delete_MCAR(wine, p=0.3, cols_mis='alcohol')
hist(em_ill$alcohol, breaks=30, xlim=c(11,15), freq=FALSE, main='Histogram of alcohol under MCAR with p=0.3',
     xlab='alcohol content', cex.lab=1.3)
x <- seq(11,15,length=1000)

wnorm <- function(x, lambda, mu, sigma, k) {
  d <- 0
  for (i in 1:k) {
    d <- d + lambda*((1/sqrt(2*pi*sigma[i]^2)) * exp(-0.5*(((x-mu[i])/sigma[i])^2)))
  }
  return(d)
}

ymix <- wnorm(x, lambda=c(0.5,0.5), mu=c(12.3, 13.5), sigma=c(0.7, 0.7), k=2)
lines(x, ymix, lwd=2, col='blue')


EMx <- em_ill$alcohol[!is.na(em_ill$alcohol)]
EMm <- normalmixEM(EMx, k=2, lambda=c(0.5, 0.5), mu=c(12.3, 13.5), sigma=c(0.7, 0.7))
EMm$mu
EMm$sigma
EMm$lambda

em <- wnorm(x, lambda=EMm$lambda, mu=EMm$mu, sigma=EMm$sigma, k=2)
lines(x, em, lwd=2, col='firebrick')

legend('topleft', legend=c('Normal mix model with initial parameters: N(12.3, 0.7) and N(13.5, 0.7)', 
                           'Normal mix model after EM: N(12.33, 0.52) and N(13.55, 0.46)'), 
       col=c('blue', 'firebrick'), lwd=2, cex=1.2)

