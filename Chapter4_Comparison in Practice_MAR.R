# Chapter 4: Comparison in Practice (MAR)

setwd('~/MT4599')
library(missMethods)
library(Metrics)
library(rms)

# dataset to run simulation studies on:
wine <- read.table('wine.data', dec='.', sep=',') # read in data
# name columns according to documentation for dataset
colnames(wine) <- c('class', 'alcohol', 'malic acid', 'ash', 'alkalinity of ash',
                    'magnesium', 'total phenols', 'flavanoids', 'nonflavanoid phenols', 
                    'proanthocyanins', 'color intensity', 'hue', 'OD2800/OD315', 
                    'proline')

truemu <- mean(wine$alcohol)
truevar <- var(wine$alcohol)


##################### COMPLETE CASE ANALYSIS (CCA) #############################
# univariate, MAR, p=0.1, CCA
nsim <- 100
ts <- rep(0, nsim)
Fs <- rep(0, nsim)
s <- 1157

for (i in 1:nsim) {
  set.seed(s)
  MAR_1 <- delete_MAR_censoring(ds=wine, p=0.1, cols_mis='alcohol', cols_ctrl='ash')
  mar_1_cca <- na.omit(MAR_1)
  ts[i] <- t.test(mar_1_cca$alcohol, wine$alcohol)$p.value
  Fs[i] <- var.test(mar_1_cca$alcohol, wine$alcohol)$p.value
  # cannot calculate RMSE for CCA
  s <- s+1
}
length(ts[ts<=0.05]) # 0
length(Fs[Fs<=0.05]) # 0


# univariate, MAR, p=0.3, CCA
nsim <- 100
ts <- rep(0, nsim)
Fs <- rep(0, nsim)
s <- 1157

for (i in 1:nsim) {
  set.seed(s)
  MAR_3 <- delete_MAR_censoring(ds=wine, p=0.3, cols_mis='alcohol', cols_ctrl='ash')
  mar_3_cca <- na.omit(MAR_3)
  ts[i] <- t.test(mar_1_cca$alcohol, wine$alcohol)$p.value
  Fs[i] <- var.test(mar_1_cca$alcohol, wine$alcohol)$p.value
  # cannot calculate RMSE for CCA
  s <- s+1
}
length(ts[ts<=0.05]) # 0
length(Fs[Fs<=0.05]) # 0


# univariate, MAR, p=0.5, CCA
nsim <- 100
ts <- rep(0, nsim)
Fs <- rep(0, nsim)
s <- 1157

for (i in 1:nsim) {
  set.seed(s)
  MAR_5 <- delete_MAR_censoring(ds=wine, p=0.5, cols_mis='alcohol', cols_ctrl='ash')
  mar_5_cca <- na.omit(MAR_5)
  ts[i] <- t.test(mar_5_cca$alcohol, wine$alcohol)$p.value
  Fs[i] <- var.test(mar_5_cca$alcohol, wine$alcohol)$p.value
  # cannot calculate RMSE for CCA
  s <- s+1
}
length(ts[ts<=0.05]) # 0
length(Fs[Fs<=0.05]) # 0


# univariate, MAR, p=0.7, CCA
nsim <- 100
ts <- rep(0, nsim)
Fs <- rep(0, nsim)
s <- 1157

for (i in 1:nsim) {
  set.seed(s)
  MAR_7 <- delete_MAR_censoring(ds=wine, p=0.7, cols_mis='alcohol', cols_ctrl='ash')
  mar_7_cca <- na.omit(MAR_7)
  ts[i] <- t.test(mar_7_cca$alcohol, wine$alcohol)$p.value
  Fs[i] <- var.test(mar_7_cca$alcohol, wine$alcohol)$p.value
  # cannot calculate RMSE for CCA
  s <- s+1
}
length(ts[ts<=0.05]) # 0
length(Fs[Fs<=0.05]) # 0



########################## MEAN IMPUTATION ###################################
# univariate, MAR, p=0.1, Mean imputation
nsim <- 100
ts <- rep(0, nsim)
Fs <- rep(0, nsim)
RMSE <- rep(0, nsim)
s <- 1157

for (i in 1:nsim) {
  set.seed(s)
  MAR_1 <- delete_MAR_censoring(ds=wine, p=0.1, cols_mis='alcohol', cols_ctrl='ash')
  mar_1_mi <- MAR_1
  mar_1_mi$alcohol <- replace(mar_1_mi$alcohol, is.na(mar_1_mi$alcohol),
                               round(mean(mar_1_mi$alcohol, na.rm=TRUE), 2))
  ts[i] <- t.test(mar_1_mi$alcohol, wine$alcohol)$p.value
  Fs[i] <- var.test(mar_1_mi$alcohol, wine$alcohol)$p.value
  RMSE[i] <- rmse(wine$alcohol, mar_1_mi$alcohol)
  s <- s+1
}
length(ts[ts<=0.05]) # 0
length(Fs[Fs<=0.05]) # 0
RMSEmu <- round(mean(RMSE), 2) 
RMSEmu               # 0.28
RMSEci <- round(c(quantile(RMSE, 0.025), quantile(RMSE, 0.975)), 2)
RMSEci               # 0.28, 0.28
# note: acceptable when between 0.2 and 0.5 (rule of thumb, best used comparatively)

# univariate, MAR, p=0.3, mean imputation 
nsim <- 100
ts <- rep(0, nsim)
Fs <- rep(0, nsim)
RMSE <- rep(0, nsim)
s <- 1157

for (i in 1:nsim) {
  set.seed(s)
  MAR_3 <- delete_MAR_censoring(ds=wine, p=0.3, cols_mis='alcohol', cols_ctrl='ash')
  mar_3_mi <- MAR_3
  mar_3_mi$alcohol <- replace(mar_3_mi$alcohol, is.na(mar_3_mi$alcohol),
                               round(mean(mar_3_mi$alcohol, na.rm=TRUE), 2))
  ts[i] <- t.test(mar_3_mi$alcohol, wine$alcohol)$p.value
  Fs[i] <- var.test(mar_3_mi$alcohol, wine$alcohol)$p.value
  RMSE[i] <- rmse(wine$alcohol, mar_3_mi$alcohol)
  s <- s+1
}
length(ts[ts<=0.05]) # 0
length(Fs[Fs<=0.05]) # 100
RMSEmu <- round(mean(RMSE), 2)
RMSEmu               # 0.46
RMSEci <- round(c(quantile(RMSE, 0.025), quantile(RMSE, 0.975)), 2) 
RMSEci               # 0.46, 0.46

# univariate, MAR, p=0.5, mean imputation 
nsim <- 100
ts <- rep(0, nsim)
Fs <- rep(0, nsim)
RMSE <- rep(0, nsim)
s <- 1157

for (i in 1:nsim) {
  set.seed(s)
  MAR_5 <- delete_MAR_censoring(ds=wine, p=0.5, cols_mis='alcohol', cols_ctrl='ash')
  mar_5_mi <- MAR_5
  mar_5_mi$alcohol <- replace(mar_5_mi$alcohol, is.na(mar_5_mi$alcohol),
                               round(mean(mar_5_mi$alcohol, na.rm=TRUE), 2))
  ts[i] <- t.test(mar_5_mi$alcohol, wine$alcohol)$p.value
  Fs[i] <- var.test(mar_5_mi$alcohol, wine$alcohol)$p.value
  RMSE[i] <- rmse(wine$alcohol, mar_5_mi$alcohol)
  s <- s+1
}
length(ts[ts<=0.05]) # 100
length(Fs[Fs<=0.05]) # 100
RMSEmu <- round(mean(RMSE), 2)
RMSEmu               # 0.58
RMSEci <- round(c(quantile(RMSE, 0.025), quantile(RMSE, 0.975)), 2) 
RMSEci               # 0.53, 0.62

# univariate MAR, p=0.7, mean imputation 
nsim <- 100
ts <- rep(0, nsim)
Fs <- rep(0, nsim)
RMSE <- rep(0, nsim)
# s <- 1157 

for (i in 1:nsim) {
  # set.seed(s)
  MAR_7 <- delete_MAR_censoring(ds=wine, p=0.7, cols_mis='alcohol', cols_ctrl='ash')
  mar_7_mi <- MAR_7
  mar_7_mi$alcohol <- replace(mar_7_mi$alcohol, is.na(mar_7_mi$alcohol),
                               round(mean(mar_7_mi$alcohol, na.rm=TRUE), 2))
  ts[i] <- t.test(mar_7_mi$alcohol, wine$alcohol)$p.value
  Fs[i] <- var.test(mar_7_mi$alcohol, wine$alcohol)$p.value
  RMSE[i] <- rmse(wine$alcohol, mar_7_mi$alcohol)
  # s <- s+1
}
length(ts[ts<=0.05]) # 17
length(Fs[Fs<=0.05]) # 100
RMSEmu <- round(mean(RMSE), 2)
RMSEmu               # 0.69
RMSEci <- round(c(quantile(RMSE, 0.025), quantile(RMSE, 0.975)), 2) 
RMSEci               # 0.66, 0.72


########################## REGRESSION-BASED IMPUTATION ########################
# univariate, MAR, p=0.1, regression-based imputation 
nsim <- 100
ts <- rep(0, nsim)
Fs <- rep(0, nsim)
RMSE <- rep(0, nsim)
s <- 1157

for (i in 1:nsim) {
  set.seed(s)
  MAR_1 <- delete_MAR_censoring(ds=wine, p=0.1, cols_mis='alcohol', cols_ctrl='ash')
  mar_1_reg <- MAR_1
  mar_reg <- lm(alcohol ~ ash, data=mar_1_reg)
  newdat <- mar_1_reg[is.na(mar_1_reg$alcohol), ]
  predmar1 <- predict(mar_reg, newdata=newdat)
  mar_1_reg$alcohol <- replace(mar_1_reg$alcohol, is.na(mar_1_reg$alcohol), 
                                round(predmar1, 2))
  ts[i] <- t.test(mar_1_reg$alcohol, wine$alcohol)$p.value
  Fs[i] <- var.test(mar_1_reg$alcohol, wine$alcohol)$p.value
  RMSE[i] <- rmse(wine$alcohol, mar_1_reg$alcohol)
  s <- s+1
}
length(ts[ts<=0.05]) # 0
length(Fs[Fs<=0.05]) # 0
RMSEmu <- round(mean(RMSE), 2)
RMSEmu               # 0.26
RMSEci <- round(c(quantile(RMSE, 0.025), quantile(RMSE, 0.975)), 2) 
RMSEci               # 0.20, 0.32


# univariate, MAR, p=0.3, regression-based imputation
nsim <- 100
ts <- rep(0, nsim)
Fs <- rep(0, nsim)
RMSE <- rep(0, nsim)
s <- 1157

for (i in 1:nsim) {
  set.seed(s)
  MAR_3 <- delete_MAR_censoring(ds=wine, p=0.3, cols_mis='alcohol', cols_ctrl='ash')
  mar_3_reg <- MAR_3
  mar_reg <- lm(alcohol ~ ash, data=mar_3_reg)
  newdat <- mar_3_reg[is.na(mar_3_reg$alcohol), ]
  predmar3 <- predict(mar_reg, newdata=newdat)
  mar_3_reg$alcohol <- replace(mar_3_reg$alcohol, is.na(mar_3_reg$alcohol), 
                                round(predmar3, 2))
  ts[i] <- t.test(mar_3_reg$alcohol, wine$alcohol)$p.value
  Fs[i] <- var.test(mar_3_reg$alcohol, wine$alcohol)$p.value
  RMSE[i] <- rmse(wine$alcohol, mar_3_reg$alcohol)
  s <- s+1
}
length(ts[ts<=0.05]) # 0
length(Fs[Fs<=0.05]) # 84
RMSEmu <- round(mean(RMSE), 2)
RMSEmu               # 0.44
RMSEci <- round(c(quantile(RMSE, 0.025), quantile(RMSE, 0.975)), 2) 
RMSEci               # 0.38, 0.49


# univariate, MAR, p=0.5, regression-based imputation 
nsim <- 100
ts <- rep(0, nsim)
Fs <- rep(0, nsim)
RMSE <- rep(0, nsim)
s <- 1157

for (i in 1:nsim) {
  set.seed(s)
  MAR_5 <- delete_MAR_censoring(ds=wine, p=0.5, cols_mis='alcohol', cols_ctrl='ash')
  mar_5_reg <- MAR_5
  mar_reg <- lm(alcohol ~ ash, data=mar_5_reg)
  newdat <- mar_5_reg[is.na(mar_5_reg$alcohol), ]
  predmar5 <- predict(mar_reg, newdata=newdat)
  mar_5_reg$alcohol <- replace(mar_5_reg$alcohol, is.na(mar_5_reg$alcohol), 
                                round(predmar5, 2))
  ts[i] <- t.test(mar_5_reg$alcohol, wine$alcohol)$p.value
  Fs[i] <- var.test(mar_5_reg$alcohol, wine$alcohol)$p.value
  RMSE[i] <- rmse(wine$alcohol, mar_5_reg$alcohol)
  s <- s+1
}
length(ts[ts<=0.05]) # 100
length(Fs[Fs<=0.05]) # 100
RMSEmu <- round(mean(RMSE), 2)
RMSEmu               # 0.68
RMSEci <- round(c(quantile(RMSE, 0.025), quantile(RMSE, 0.975)), 2) 
RMSEci               # 


# univariate, MAR, p=0.7, regression-based imputation 
nsim <- 100
ts <- rep(0, nsim)
Fs <- rep(0, nsim)
RMSE <- rep(0, nsim)
s <- 1157

for (i in 1:nsim) {
  set.seed(s)
  MAR_7 <- delete_MAR_censoring(ds=wine, p=0.7, cols_mis='alcohol', cols_ctrl='ash')
  mar_7_reg <- MAR_7
  mar_reg <- lm(alcohol ~ ash, data=mar_7_reg)
  newdat <- mar_7_reg[is.na(mar_7_reg$alcohol), ]
  predmar7 <- predict(mar_reg, newdata=newdat)
  mar_7_reg$alcohol <- replace(mar_7_reg$alcohol, is.na(mar_7_reg$alcohol), 
                                round(predmar7, 2))
  ts[i] <- t.test(mar_7_reg$alcohol, wine$alcohol)$p.value
  Fs[i] <- var.test(mar_7_reg$alcohol, wine$alcohol)$p.value
  RMSE[i] <- rmse(wine$alcohol, mar_7_reg$alcohol)
  s <- s+1
}
length(ts[ts<=0.05]) # 100
length(Fs[Fs<=0.05]) # 100
RMSEmu <- round(mean(RMSE), 2)
RMSEmu               # 0.81
RMSEci <- round(c(quantile(RMSE, 0.025), quantile(RMSE, 0.975)), 2) 
RMSEci               # 


#################### HOT DECK IMPUTATION #####################################
# univariate, MAR, p=0.1, hot deck imputation
nsim <- 100
ts <- rep(0, nsim)
Fs <- rep(0, nsim)
RMSE <- rep(0, nsim)
s <- 1157

for (i in 1:nsim) {
  set.seed(s)
  MAR_1 <- delete_MAR_censoring(ds=wine, p=0.1, cols_mis='alcohol', cols_ctrl='ash')
  mar_1_hd <- MAR_1
  mar_1_hd$alcohol <- replace(mar_1_hd$alcohol, is.na(mar_1_hd$alcohol),
                               sample(mar_1_hd$alcohol[!is.na(mar_1_hd$alcohol)],
                                      length(mar_1_hd$alcohol[is.na(mar_1_hd$alcohol)]),
                                      replace=TRUE))
  ts[i] <- t.test(mar_1_hd$alcohol, wine$alcohol)$p.value
  Fs[i] <- var.test(mar_1_hd$alcohol, wine$alcohol)$p.value
  RMSE[i] <- rmse(wine$alcohol, mar_1_hd$alcohol)
  s <- s+1
}
length(ts[ts<=0.05]) # 0
length(Fs[Fs<=0.05]) # 0
RMSEmu <- round(mean(RMSE), 2)
RMSEmu               # 0.39
RMSEci <- round(c(quantile(RMSE, 0.025), quantile(RMSE, 0.975)), 2) 
RMSEci               # 0.29, 0.48


# univariate, MAR, p=0.3, hot deck imputation
nsim <- 100
ts <- rep(0, nsim)
Fs <- rep(0, nsim)
RMSE <- rep(0, nsim)
s <- 1157

for (i in 1:nsim) {
  set.seed(s)
  MAR_3 <- delete_MAR_censoring(ds=wine, p=0.3, cols_mis='alcohol', cols_ctrl='ash')
  mar_3_hd <- MAR_3
  mar_3_hd$alcohol <- replace(mar_3_hd$alcohol, is.na(mar_3_hd$alcohol),
                               sample(mar_3_hd$alcohol[!is.na(mar_3_hd$alcohol)],
                                      length(mar_3_hd$alcohol[is.na(mar_3_hd$alcohol)]),
                                      replace=TRUE))
  ts[i] <- t.test(mar_3_hd$alcohol, wine$alcohol)$p.value
  Fs[i] <- var.test(mar_3_hd$alcohol, wine$alcohol)$p.value
  RMSE[i] <- rmse(wine$alcohol, mar_3_hd$alcohol)
  s <- s+1
}
length(ts[ts<=0.05]) # 26
length(Fs[Fs<=0.05]) # 0
RMSEmu <- round(mean(RMSE), 2)
RMSEmu               # 0.64
RMSEci <- round(c(quantile(RMSE, 0.025), quantile(RMSE, 0.975)), 2) 
RMSEci               # 0.56, 0.72


# univariate, MAR, p=0.5, hot deck imputation 
nsim <- 100
ts <- rep(0, nsim)
Fs <- rep(0, nsim)
RMSE <- rep(0, nsim)
s <- 1157

for (i in 1:nsim) {
  set.seed(s)
  MAR_5 <- delete_MAR_censoring(ds=wine, p=0.5, cols_mis='alcohol', cols_ctrl='ash')
  mar_5_hd <- MAR_5
  mar_5_hd$alcohol <- replace(mar_5_hd$alcohol, is.na(mar_5_hd$alcohol),
                               sample(mar_5_hd$alcohol[!is.na(mar_5_hd$alcohol)],
                                      length(mar_5_hd$alcohol[is.na(mar_5_hd$alcohol)]),
                                      replace=TRUE))
  ts[i] <- t.test(mar_5_hd$alcohol, wine$alcohol)$p.value
  Fs[i] <- var.test(mar_5_hd$alcohol, wine$alcohol)$p.value
  RMSE[i] <- rmse(wine$alcohol, mar_5_hd$alcohol)
  s <- s+1
}
length(ts[ts<=0.05]) # 28
length(Fs[Fs<=0.05]) # 0
RMSEmu <- round(mean(RMSE), 2)
RMSEmu               # 0.82
RMSEci <- round(c(quantile(RMSE, 0.025), quantile(RMSE, 0.975)), 2) 
RMSEci               # 0.74, 0.91


# univariate, MAR, p=0.7, hot deck imputation
nsim <- 100
ts <- rep(0, nsim)
Fs <- rep(0, nsim)
RMSE <- rep(0, nsim)
s <- 1157

for (i in 1:nsim) {
  set.seed(s)
  MAR_7 <- delete_MAR_censoring(ds=wine, p=0.7, cols_mis='alcohol', cols_ctrl='ash')
  mar_7_hd <- MAR_7
  mar_7_hd$alcohol <- replace(mar_7_hd$alcohol, is.na(mar_7_hd$alcohol),
                               sample(mar_7_hd$alcohol[!is.na(mar_7_hd$alcohol)],
                                      length(mar_7_hd$alcohol[is.na(mar_7_hd$alcohol)]),
                                      replace=TRUE))
  ts[i] <- t.test(mar_7_hd$alcohol, wine$alcohol)$p.value
  Fs[i] <- var.test(mar_7_hd$alcohol, wine$alcohol)$p.value
  RMSE[i] <- rmse(wine$alcohol, mar_7_hd$alcohol)
  s <- s+1
}
length(ts[ts<=0.05]) # 7
length(Fs[Fs<=0.05]) # 0
RMSEmu <- round(mean(RMSE), 2)
RMSEmu               # 0.95
RMSEci <- round(c(quantile(RMSE, 0.025), quantile(RMSE, 0.975)), 2) 
RMSEci               # 0.87, 1.04


################# IMPUTATION FROM A CONDITIONAL DISTRIBUTION ##################
# univariate, MAR, p=0.1, imputation from a conditional distribution 
# Gibbs sampler to find predictive distribution (to then sample from)
# alcohol ~ N(mu, sigma^2), mu and sigma unknown
# priors: mu ~ N(phi, tau^2), sigma ~ invGamma(alpha, beta)
nsim <- 100
ts <- rep(0, nsim)
Fs <- rep(0, nsim)
RMSE <- rep(0, nsim)
s <- 1157

for (i in 1:nsim) {
  set.seed(s)
  MAR_1 <- delete_MAR_censoring(ds=wine, p=0.1, cols_mis='alcohol', cols_ctrl='ash')
  # Priors
  phi <- round(mean(MAR_1$alcohol, na.rm=TRUE), 2)
  tau2 <- 1
  alpha <- 0.1
  beta <- 0.1
  
  # simulate data
  n <- length(MAR_1$alcohol[is.na(MAR_1$alcohol)])
  x <- rnorm(n, phi, tau2)
  xbar <- mean(x, na.rm=TRUE)
  
  # samples to take in the chain 
  T <- 1000
  mu <- rep(0, 2)
  sigma2 <- rep(0, 2)
  
  # set starting values
  mu[1] <- 13
  sigma2[1] <- 0.1
  
  # run Gibbs sampler
  for (t in 1:(T-1)) {
    mu[t+1]<-rnorm(1,(tau2*n*xbar + sigma2[t]*phi)/(tau2*n + sigma2[t]),
                   sqrt(sigma2[t]*tau2/(tau2*n+sigma2[t])))
    sigma2[t+1]<-1/rgamma(1,shape=(n/2 + alpha),
                          rate=(1/2*sum((x-mu[t+1])^2) + beta))
  }
  
  # posteriors
  post_mu <- mean(mu)
  post_sig <- mean(sigma2)
  
  # sample from the posterior predictive distribution 
  mar_1_cd <- MAR_1
  mar_1_cd$alcohol <- replace(mar_1_cd$alcohol, is.na(mar_1_cd$alcohol),
                               round(rnorm(length(mar_1_cd$alcohol[is.na(mar_1_cd$alcohol)]), 
                                           post_mu, post_sig), 2))
  
  # tests
  ts[i] <- t.test(wine$alcohol, mar_1_cd$alcohol)
  Fs[i] <- var.test(wine$alcohol, mar_1_cd$alcohol)
  RMSE[i] <- rmse(wine$alcohol, mar_1_cd$alcohol)
  s <- s+1
}

length(ts[ts<=0.05]) # 97
length(Fs[Fs<=0.05]) # 0
RMSEmu <- round(mean(RMSE), 2)
RMSEmu               # 0.46
RMSEci <- round(c(quantile(RMSE, 0.025), quantile(RMSE, 0.975)), 2) 
RMSEci               # 0.26  0.71


# univariate, MAR, p=0.3, imputation from a conditonal distribution 
nsim <- 100
ts <- rep(0, nsim)
Fs <- rep(0, nsim)
RMSE <- rep(0, nsim)
s <- 1157

for (i in 1:nsim) {
  set.seed(s)
  MAR_3 <- delete_MAR_censoring(ds=wine, p=0.3, cols_mis='alcohol', cols_ctrl='ash')
  # Priors
  phi <- round(mean(MAR_3$alcohol, na.rm=TRUE), 2)
  tau2 <- 1
  alpha <- 0.1
  beta <- 0.1
  
  # simulate data
  n <- length(MAR_3$alcohol[is.na(MAR_3$alcohol)])
  x <- rnorm(n, phi, tau2)
  xbar <- mean(x, na.rm=TRUE)
  
  # samples to take in the chain 
  T <- 1000
  mu <- rep(0, 2)
  sigma2 <- rep(0, 2)
  
  # set starting values
  mu[1] <- 13
  sigma2[1] <- 0.1
  
  # run Gibbs sampler
  for (t in 1:(T-1)) {
    mu[t+1]<-rnorm(1,(tau2*n*xbar + sigma2[t]*phi)/(tau2*n + sigma2[t]),
                   sqrt(sigma2[t]*tau2/(tau2*n+sigma2[t])))
    sigma2[t+1]<-1/rgamma(1,shape=(n/2 + alpha),
                          rate=(1/2*sum((x-mu[t+1])^2) + beta))
  }
  
  # posteriors
  post_mu <- mean(mu)
  post_sig <- mean(sigma2)
  
  # sample from the posterior predictive distribution 
  mar_3_cd <- MAR_3
  mar_3_cd$alcohol <- replace(mar_3_cd$alcohol, is.na(mar_3_cd$alcohol),
                               round(rnorm(length(mar_3_cd$alcohol[is.na(mar_3_cd$alcohol)]), 
                                           post_mu, post_sig), 2))
  
  # tests
  ts[i] <- t.test(wine$alcohol, mar_3_cd$alcohol)
  Fs[i] <- var.test(wine$alcohol, mar_3_cd$alcohol)
  RMSE[i] <- rmse(wine$alcohol, mar_3_cd$alcohol)
  s <- s+1
}

length(ts[ts<=0.05]) # 98
length(Fs[Fs<=0.05]) # 0
RMSEmu <- round(mean(RMSE), 2)
RMSEmu               # 0.75
RMSEci <- round(c(quantile(RMSE, 0.025), quantile(RMSE, 0.975)), 2) 
RMSEci               # 0.57, 0.95


# univariate, MAR, p=0.5, imputation from a conditonal distribution 
nsim <- 100
ts <- rep(0, nsim)
Fs <- rep(0, nsim)
RMSE <- rep(0, nsim)
s <- 1157

for (i in 1:nsim) {
  set.seed(s)
  MAR_5 <- delete_MAR_censoring(ds=wine, p=0.5, cols_mis='alcohol', cols_ctrl='ash')
  # Priors
  phi <- round(mean(MAR_5$alcohol, na.rm=TRUE), 2)
  tau2 <- 1
  alpha <- 0.1
  beta <- 0.1
  
  # simulate data
  n <- length(MAR_5$alcohol[is.na(MAR_5$alcohol)])
  x <- rnorm(n, phi, tau2)
  xbar <- mean(x, na.rm=TRUE)
  
  # samples to take in the chain 
  T <- 1000
  mu <- rep(0, 2)
  sigma2 <- rep(0, 2)
  
  # set starting values
  mu[1] <- 13
  sigma2[1] <- 0.1
  
  # run Gibbs sampler
  for (t in 1:(T-1)) {
    mu[t+1]<-rnorm(1,(tau2*n*xbar + sigma2[t]*phi)/(tau2*n + sigma2[t]),
                   sqrt(sigma2[t]*tau2/(tau2*n+sigma2[t])))
    sigma2[t+1]<-1/rgamma(1,shape=(n/2 + alpha),
                          rate=(1/2*sum((x-mu[t+1])^2) + beta))
  }
  
  # posteriors
  post_mu <- mean(mu)
  post_sig <- mean(sigma2)
  
  # sample from the posterior predictive distribution 
  mar_5_cd <- MAR_5
  mar_5_cd$alcohol <- replace(mar_5_cd$alcohol, is.na(mar_5_cd$alcohol),
                               round(rnorm(length(mar_5_cd$alcohol[is.na(mar_5_cd$alcohol)]), 
                                           post_mu, post_sig), 2))
  
  # tests
  ts[i] <- t.test(wine$alcohol, mar_5_cd$alcohol)
  Fs[i] <- var.test(wine$alcohol, mar_5_cd$alcohol)
  RMSE[i] <- rmse(wine$alcohol, mar_5_cd$alcohol)
  s <- s+1
}

length(ts[ts<=0.05]) # 100
length(Fs[Fs<=0.05]) # 0
RMSEmu <- round(mean(RMSE), 2)
RMSEmu               # 0.93
RMSEci <- round(c(quantile(RMSE, 0.025), quantile(RMSE, 0.975)), 2) 
RMSEci               # 0.74, 1.16


# univariate, MAR, p=0.7, imputation from a conditional distribution 
nsim <- 100
ts <- rep(0, nsim)
Fs <- rep(0, nsim)
RMSE <- rep(0, nsim)
s <- 1157

for (i in 1:nsim) {
  set.seed(s)
  MAR_7 <- delete_MAR_censoring(ds=wine, p=0.7, cols_mis='alcohol', cols_ctrl='ash')
  # Priors
  phi <- round(mean(MAR_7$alcohol, na.rm=TRUE), 2)
  tau2 <- 1
  alpha <- 0.1
  beta <- 0.1
  
  # simulate data
  n <- length(MAR_7$alcohol[is.na(MAR_7$alcohol)])
  x <- rnorm(n, phi, tau2)
  xbar <- mean(x, na.rm=TRUE)
  
  # samples to take in the chain 
  T <- 1000
  mu <- rep(0, 2)
  sigma2 <- rep(0, 2)
  
  # set starting values
  mu[1] <- 13
  sigma2[1] <- 0.1
  
  # run Gibbs sampler
  for (t in 1:(T-1)) {
    mu[t+1]<-rnorm(1,(tau2*n*xbar + sigma2[t]*phi)/(tau2*n + sigma2[t]),
                   sqrt(sigma2[t]*tau2/(tau2*n+sigma2[t])))
    sigma2[t+1]<-1/rgamma(1,shape=(n/2 + alpha),
                          rate=(1/2*sum((x-mu[t+1])^2) + beta))
  }
  
  # posteriors
  post_mu <- mean(mu)
  post_sig <- mean(sigma2)
  
  # sample from the posterior predictive distribution 
  mar_7_cd <- MAR_7
  mar_7_cd$alcohol <- replace(mar_7_cd$alcohol, is.na(mar_7_cd$alcohol),
                               round(rnorm(length(mar_7_cd$alcohol[is.na(mar_7_cd$alcohol)]), 
                                           post_mu, post_sig), 2))
  
  # tests
  ts[i] <- t.test(wine$alcohol, mar_7_cd$alcohol)
  Fs[i] <- var.test(wine$alcohol, mar_7_cd$alcohol)
  RMSE[i] <- rmse(wine$alcohol, mar_7_cd$alcohol)
  s <- s+1
}

length(ts[ts<=0.05]) # 85
length(Fs[Fs<=0.05]) # 0
RMSEmu <- round(mean(RMSE), 2)
RMSEmu               # 1.1
RMSEci <- round(c(quantile(RMSE, 0.025), quantile(RMSE, 0.975)), 2) 
RMSEci               # 0.92, 1.29


################################ IPW ##########################################
# univariate, MAR, p=0.1, IPW
nsim <- 100
ts <- rep(0, nsim)
Fs <- rep(0, nsim)
RMSE <- rep(0, nsim)
s <- 1157

for (i in 1:nsim) {
  set.seed(s)
  MAR_1 <- delete_MAR_censoring(ds=wine, p=0.1, cols_mis='alcohol', cols_ctrl='ash')
  mar_1_ipw <- MAR_1
  # add two columns to the data: 1) 0/1 missingness, 2) weight
  mar_1_ipw$R <- rep(0, 178)
  mar_1_ipw$weight <- rep(0, 178)
  for (j in 1:178) {
    if (is.na(mar_1_ipw$alcohol[j])) {
      mar_1_ipw$R[j] <- 1
    }
  }
  
  # regression on missingness to determine weights 
  # ipw <- lm(R ~ ash, data=mar_1_ipw)
  ipw <- lm(R ~ ash, data=mar_1_ipw)
  ipw_newdat <- mar_1_ipw
  
  # calculate weights based on logistic regression model
  weight <- predict(ipw, ipw_newdat)
  mar_1_ipw$weight <- weight
  
  # calculate weighted mean
  wmean <- round(weighted.mean(mar_1_ipw$alcohol[!is.na(mar_1_ipw$alcohol)],
                               mar_1_ipw$weight[!is.na(mar_1_ipw$alcohol)]), 2)
  
  # replace missing values with weighted mean
  mar_1_ipw$alcohol <- replace(mar_1_ipw$alcohol, is.na(mar_1_ipw$alcohol), wmean)
  
  ts[i] <- t.test(mar_1_ipw$alcohol, wine$alcohol)$p.value
  Fs[i] <- var.test(mar_1_ipw$alcohol, wine$alcohol)$p.value
  RMSE[i] <- rmse(wine$alcohol, mar_1_ipw$alcohol)
  s <- s+1
}

length(ts[ts<=0.05]) # 100
length(Fs[Fs<=0.05]) # 100
RMSEmu <- round(mean(RMSE), 2)
RMSEmu               # 0.25
RMSEci <- round(c(quantile(RMSE, 0.025), quantile(RMSE, 0.975)), 2) 
RMSEci               # 0.25, 0.25


# univariate, MAR, p=0.3, IPW
nsim <- 100
ts <- rep(0, nsim)
Fs <- rep(0, nsim)
RMSE <- rep(0, nsim)
s <- 1157

for (i in 1:nsim) {
  set.seed(s)
  MAR_3 <- delete_MAR_censoring(ds=wine, p=0.3, cols_mis='alcohol', cols_ctrl='ash')
  mar_3_ipw <- MAR_3
  # add two columns to the data: 1) 0/1 missingness, 2) weight
  mar_3_ipw$R <- rep(0, 178)
  mar_3_ipw$weight <- rep(0, 178)
  for (j in 1:178) {
    if (is.na(mar_3_ipw$alcohol[j])) {
      mar_3_ipw$R[j] <- 1
    }
  }
  
  # logistic regression on missingness to determine weights 
  # ipw <- lm(R ~ ash, data=mar_3_ipw)
  ipw <- lm(R ~ ash, data=mar_3_ipw)
  ipw_newdat <- mar_3_ipw
  
  # calculate weights based on logistic regression model
  weight <- predict(ipw, ipw_newdat)
  mar_3_ipw$weight <- weight
  
  # calculate weighted mean
  wmean <- round(weighted.mean(mar_3_ipw$alcohol[!is.na(mar_3_ipw$alcohol)],
                               mar_3_ipw$weight[!is.na(mar_3_ipw$alcohol)]), 2)
  
  # replace missing values with weighted mean
  mar_3_ipw$alcohol <- replace(mar_3_ipw$alcohol, is.na(mar_3_ipw$alcohol), wmean)
  
  ts[i] <- t.test(mar_3_ipw$alcohol, wine$alcohol)$p.value
  Fs[i] <- var.test(mar_3_ipw$alcohol, wine$alcohol)$p.value
  RMSE[i] <- rmse(wine$alcohol, mar_3_ipw$alcohol)
  s <- s+1
}

length(ts[ts<=0.05]) # 100
length(Fs[Fs<=0.05]) # 100
RMSEmu <- round(mean(RMSE), 2)
RMSEmu               # 0.49
RMSEci <- round(c(quantile(RMSE, 0.025), quantile(RMSE, 0.975)), 2) 
RMSEci               # 0.49, 0.49


# univariate, MAR, p=0.5, IPW
nsim <- 100
ts <- rep(0, nsim)
Fs <- rep(0, nsim)
RMSE <- rep(0, nsim)
s <- 1157

for (i in 1:nsim) {
  set.seed(s)
  MAR_5 <- delete_MAR_censoring(ds=wine, p=0.5, cols_mis='alcohol', cols_ctrl='ash')
  mar_5_ipw <- MAR_5
  # add two columns to the data: 1) 0/1 missingness, 2) weight
  mar_5_ipw$R <- rep(0, 178)
  mar_5_ipw$weight <- rep(0, 178)
  for (j in 1:178) {
    if (is.na(mar_5_ipw$alcohol[j])) {
      mar_5_ipw$R[j] <- 1
    }
  }
  
  # logistic regression on missingness to determine weights 
  # ipw <- lm(R ~ ash, data=mar_5_ipw)
  ipw <- lm(R ~ ash, data=mar_5_ipw)
  ipw_newdat <- mar_5_ipw
  
  # calculate weights based on logistic regression model
  weight <- predict(ipw, ipw_newdat)
  mar_5_ipw$weight <- weight
  
  # calculate weighted mean
  wmean <- round(weighted.mean(mar_5_ipw$alcohol[!is.na(mar_5_ipw$alcohol)],
                               mar_5_ipw$weight[!is.na(mar_5_ipw$alcohol)]), 2)
  
  # replace missing values with weighted mean
  mar_5_ipw$alcohol <- replace(mar_5_ipw$alcohol, is.na(mar_5_ipw$alcohol), wmean)
  
  ts[i] <- t.test(mar_5_ipw$alcohol, wine$alcohol)$p.value
  Fs[i] <- var.test(mar_5_ipw$alcohol, wine$alcohol)$p.value
  RMSE[i] <- rmse(wine$alcohol, mar_5_ipw$alcohol)
  s <- s+1
}

length(ts[ts<=0.05]) # 100
length(Fs[Fs<=0.05]) # 100
RMSEmu <- round(mean(RMSE), 2)
RMSEmu               # 0.61
RMSEci <- round(c(quantile(RMSE, 0.025), quantile(RMSE, 0.975)), 2) 
RMSEci               # 0.61, 0.61


# univariate, MAR, p=0.7, IPW
nsim <- 100
ts <- rep(0, nsim)
Fs <- rep(0, nsim)
RMSE <- rep(0, nsim)
s <- 1157

for (i in 1:nsim) {
  set.seed(s)
  MAR_7 <- delete_MAR_censoring(ds=wine, p=0.7, cols_mis='alcohol', cols_ctrl='ash')
  mar_7_ipw <- MAR_7
  # add two columns to the data: 1) 0/1 missingness, 2) weight
  mar_7_ipw$R <- rep(0, 178)
  mar_7_ipw$weight <- rep(0, 178)
  for (j in 1:178) {
    if (is.na(mar_7_ipw$alcohol[j])) {
      mar_7_ipw$R[j] <- 1
    }
  }
  
  # logistic regression on missingness to determine weights 
  # ipw <- lm(R ~ ash + `color intensity` + magnesium + `total phenols` + `flavanoids`, data=mar_7_ipw)
  ipw <- lm(R ~ ash, data=mar_7_ipw)
  ipw_newdat <- mar_7_ipw
  
  # calculate weights based on logistic regression model
  weight <- predict(ipw, ipw_newdat)
  mar_7_ipw$weight <- weight
  
  # calculate weighted mean
  wmean <- round(weighted.mean(mar_7_ipw$alcohol[!is.na(mar_7_ipw$alcohol)],
                               mar_7_ipw$weight[!is.na(mar_7_ipw$alcohol)]), 2)
  
  # replace missing values with weighted mean
  mar_7_ipw$alcohol <- replace(mar_7_ipw$alcohol, is.na(mar_7_ipw$alcohol), wmean)
  
  ts[i] <- t.test(mar_7_ipw$alcohol, wine$alcohol)$p.value
  Fs[i] <- var.test(mar_7_ipw$alcohol, wine$alcohol)$p.value
  RMSE[i] <- rmse(wine$alcohol, mar_7_ipw$alcohol)
  s <- s+1
}

length(ts[ts<=0.05]) # 0
length(Fs[Fs<=0.05]) # 100
RMSEmu <- round(mean(RMSE), 2)
RMSEmu               # 0.70
RMSEci <- round(c(quantile(RMSE, 0.025), quantile(RMSE, 0.975)), 2) 
RMSEci               # 0.70, 0.70