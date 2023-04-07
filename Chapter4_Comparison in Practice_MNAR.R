# Chapter 4: Comparison in Practice (MNAR)

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
# univariate, MNAR, p=0.1, CCA
nsim <- 100
ts <- rep(0, nsim)
Fs <- rep(0, nsim)
s <- 1157

for (i in 1:nsim) {
  set.seed(s)
  MNAR_1 <- delete_MNAR_censoring(ds=wine, p=0.1, cols_mis='alcohol')
  mnar_1_cca <- na.omit(MNAR_1)
  ts[i] <- t.test(mnar_1_cca$alcohol, wine$alcohol)$p.value
  Fs[i] <- var.test(mnar_1_cca$alcohol, wine$alcohol)$p.value
  # cannot calculate RMSE for CCA
  s <- s+1
}
length(ts[ts<=0.05]) # 0
length(Fs[Fs<=0.05]) # 0


# univariate, MNAR, p=0.3, CCA
nsim <- 100
ts <- rep(0, nsim)
Fs <- rep(0, nsim)
s <- 1157

for (i in 1:nsim) {
  set.seed(s)
  MNAR_3 <- delete_MNAR_censoring(ds=wine, p=0.3, cols_mis='alcohol')
  mnar_3_cca <- na.omit(MNAR_3)
  ts[i] <- t.test(mnar_1_cca$alcohol, wine$alcohol)$p.value
  Fs[i] <- var.test(mnar_1_cca$alcohol, wine$alcohol)$p.value
  # cannot calculate RMSE for CCA
  s <- s+1
}
length(ts[ts<=0.05]) # 0
length(Fs[Fs<=0.05]) # 0


# univariate, MNAR, p=0.5, CCA
nsim <- 100
ts <- rep(0, nsim)
Fs <- rep(0, nsim)
s <- 1157

for (i in 1:nsim) {
  set.seed(s)
  MNAR_5 <- delete_MNAR_censoring(ds=wine, p=0.5, cols_mis='alcohol')
  mnar_5_cca <- na.omit(MNAR_5)
  ts[i] <- t.test(mnar_5_cca$alcohol, wine$alcohol)$p.value
  Fs[i] <- var.test(mnar_5_cca$alcohol, wine$alcohol)$p.value
  # cannot calculate RMSE for CCA
  s <- s+1
}
length(ts[ts<=0.05]) # 100
length(Fs[Fs<=0.05]) # 100


# univariate, MNAR, p=0.7, CCA
nsim <- 100
ts <- rep(0, nsim)
Fs <- rep(0, nsim)
s <- 1157

for (i in 1:nsim) {
  set.seed(s)
  MNAR_7 <- delete_MNAR_censoring(ds=wine, p=0.7, cols_mis='alcohol')
  mnar_7_cca <- na.omit(MNAR_7)
  ts[i] <- t.test(mnar_7_cca$alcohol, wine$alcohol)$p.value
  Fs[i] <- var.test(mnar_7_cca$alcohol, wine$alcohol)$p.value
  # cannot calculate RMSE for CCA
  s <- s+1
}
length(ts[ts<=0.05]) # 100
length(Fs[Fs<=0.05]) # 100



########################## MEAN IMPUTATION ###################################
# univariate, MNAR, p=0.1, Mean imputation
nsim <- 100
ts <- rep(0, nsim)
Fs <- rep(0, nsim)
RMSE <- rep(0, nsim)
s <- 1157

for (i in 1:nsim) {
  set.seed(s)
  MNAR_1 <- delete_MNAR_censoring(ds=wine, p=0.1, cols_mis='alcohol')
  mnar_1_mi <- MNAR_1
  mnar_1_mi$alcohol <- replace(mnar_1_mi$alcohol, is.na(mnar_1_mi$alcohol),
                               round(mean(mnar_1_mi$alcohol, na.rm=TRUE), 2))
  ts[i] <- t.test(mnar_1_mi$alcohol, wine$alcohol)$p.value
  Fs[i] <- var.test(mnar_1_mi$alcohol, wine$alcohol)$p.value
  RMSE[i] <- rmse(wine$alcohol, mnar_1_mi$alcohol)
  s <- s+1
}
length(ts[ts<=0.05]) # 0
length(Fs[Fs<=0.05]) # 0
RMSEmu <- round(mean(RMSE), 2) 
RMSEmu               # 0.26
RMSEci <- round(c(quantile(RMSE, 0.025), quantile(RMSE, 0.975)), 2)
RMSEci               # 0.21, 0.32
# note: acceptable when between 0.2 and 0.5 (rule of thumb, best used comparatively)

# univariate, MNAR, p=0.3, mean imputation 
nsim <- 100
ts <- rep(0, nsim)
Fs <- rep(0, nsim)
RMSE <- rep(0, nsim)
s <- 1157

for (i in 1:nsim) {
  set.seed(s)
  MNAR_3 <- delete_MNAR_censoring(ds=wine, p=0.3, cols_mis='alcohol')
  mnar_3_mi <- MNAR_3
  mnar_3_mi$alcohol <- replace(mnar_3_mi$alcohol, is.na(mnar_3_mi$alcohol),
                               round(mean(mnar_3_mi$alcohol, na.rm=TRUE), 2))
  ts[i] <- t.test(mnar_3_mi$alcohol, wine$alcohol)$p.value
  Fs[i] <- var.test(mnar_3_mi$alcohol, wine$alcohol)$p.value
  RMSE[i] <- rmse(wine$alcohol, mnar_3_mi$alcohol)
  s <- s+1
}
length(ts[ts<=0.05]) # 0
length(Fs[Fs<=0.05]) # 92
RMSEmu <- round(mean(RMSE), 2)
RMSEmu               # 0.45
RMSEci <- round(c(quantile(RMSE, 0.025), quantile(RMSE, 0.975)), 2) 
RMSEci               # 0.39, 0.49

# univariate, MNAR, p=0.5, mean imputation 
nsim <- 100
ts <- rep(0, nsim)
Fs <- rep(0, nsim)
RMSE <- rep(0, nsim)
s <- 1157

for (i in 1:nsim) {
  set.seed(s)
  MNAR_5 <- delete_MNAR_censoring(ds=wine, p=0.5, cols_mis='alcohol')
  mnar_5_mi <- MNAR_5
  mnar_5_mi$alcohol <- replace(mnar_5_mi$alcohol, is.na(mnar_5_mi$alcohol),
                               round(mean(mnar_5_mi$alcohol, na.rm=TRUE), 2))
  ts[i] <- t.test(mnar_5_mi$alcohol, wine$alcohol)$p.value
  Fs[i] <- var.test(mnar_5_mi$alcohol, wine$alcohol)$p.value
  RMSE[i] <- rmse(wine$alcohol, mnar_5_mi$alcohol)
  s <- s+1
}
length(ts[ts<=0.05]) # 1
length(Fs[Fs<=0.05]) # 100
RMSEmu <- round(mean(RMSE), 2)
RMSEmu               # 0.58
RMSEci <- round(c(quantile(RMSE, 0.025), quantile(RMSE, 0.975)), 2) 
RMSEci               # 0.53, 0.62

# univariate MNAR, p=0.7, mean imputation 
nsim <- 100
ts <- rep(0, nsim)
Fs <- rep(0, nsim)
RMSE <- rep(0, nsim)
s <- 1157 

for (i in 1:nsim) {
  set.seed(s)
  MNAR_7 <- delete_MNAR_censoring(ds=wine, p=0.7, cols_mis='alcohol')
  mnar_7_mi <- MNAR_7
  mnar_7_mi$alcohol <- replace(mnar_7_mi$alcohol, is.na(mnar_7_mi$alcohol),
                               round(mean(mnar_7_mi$alcohol, na.rm=TRUE), 2))
  ts[i] <- t.test(mnar_7_mi$alcohol, wine$alcohol)$p.value
  Fs[i] <- var.test(mnar_7_mi$alcohol, wine$alcohol)$p.value
  RMSE[i] <- rmse(wine$alcohol, mnar_7_mi$alcohol)
  s <- s+1
}
length(ts[ts<=0.05]) # 17
length(Fs[Fs<=0.05]) # 100
RMSEmu <- round(mean(RMSE), 2)
RMSEmu               # 0.69
RMSEci <- round(c(quantile(RMSE, 0.025), quantile(RMSE, 0.975)), 2) 
RMSEci               # 0.66, 0.72


########################## REGRESSION-BASED IMPUTATION ########################
# univariate, MNAR, p=0.1, regression-based imputation 
nsim <- 100
ts <- rep(0, nsim)
Fs <- rep(0, nsim)
RMSE <- rep(0, nsim)
s <- 1157

for (i in 1:nsim) {
  set.seed(s)
  MNAR_1 <- delete_MNAR_censoring(ds=wine, p=0.1, cols_mis='alcohol')
  mnar_1_reg <- MNAR_1
  mnar_reg <- lm(alcohol ~ ash, data=mnar_1_reg)
  newdat <- mnar_1_reg[is.na(mnar_1_reg$alcohol), ]
  predmnar1 <- predict(mnar_reg, newdata=newdat)
  mnar_1_reg$alcohol <- replace(mnar_1_reg$alcohol, is.na(mnar_1_reg$alcohol), 
                                round(predmnar1, 2))
  ts[i] <- t.test(mnar_1_reg$alcohol, wine$alcohol)$p.value
  Fs[i] <- var.test(mnar_1_reg$alcohol, wine$alcohol)$p.value
  RMSE[i] <- rmse(wine$alcohol, mnar_1_reg$alcohol)
  s <- s+1
}
length(ts[ts<=0.05]) # 0
length(Fs[Fs<=0.05]) # 0
RMSEmu <- round(mean(RMSE), 2)
RMSEmu               # 0.26
RMSEci <- round(c(quantile(RMSE, 0.025), quantile(RMSE, 0.975)), 2) 
RMSEci               # 0.20, 0.32


# univariate, MNAR, p=0.3, regression-based imputation
nsim <- 100
ts <- rep(0, nsim)
Fs <- rep(0, nsim)
RMSE <- rep(0, nsim)
s <- 1157

for (i in 1:nsim) {
  set.seed(s)
  MNAR_3 <- delete_MNAR_censoring(ds=wine, p=0.3, cols_mis='alcohol')
  mnar_3_reg <- MNAR_3
  mnar_reg <- lm(alcohol ~ ash, data=mnar_3_reg)
  newdat <- mnar_3_reg[is.na(mnar_3_reg$alcohol), ]
  predmnar3 <- predict(mnar_reg, newdata=newdat)
  mnar_3_reg$alcohol <- replace(mnar_3_reg$alcohol, is.na(mnar_3_reg$alcohol), 
                                round(predmnar3, 2))
  ts[i] <- t.test(mnar_3_reg$alcohol, wine$alcohol)$p.value
  Fs[i] <- var.test(mnar_3_reg$alcohol, wine$alcohol)$p.value
  RMSE[i] <- rmse(wine$alcohol, mnar_3_reg$alcohol)
  s <- s+1
}
length(ts[ts<=0.05]) # 0
length(Fs[Fs<=0.05]) # 84
RMSEmu <- round(mean(RMSE), 2)
RMSEmu               # 0.44
RMSEci <- round(c(quantile(RMSE, 0.025), quantile(RMSE, 0.975)), 2) 
RMSEci               # 0.38, 0.49


# univariate, MNAR, p=0.5, regression-based imputation 
nsim <- 100
ts <- rep(0, nsim)
Fs <- rep(0, nsim)
RMSE <- rep(0, nsim)
s <- 1157

for (i in 1:nsim) {
  set.seed(s)
  MNAR_5 <- delete_MNAR_censoring(ds=wine, p=0.5, cols_mis='alcohol')
  mnar_5_reg <- MNAR_5
  mnar_reg <- lm(alcohol ~ ash, data=mnar_5_reg)
  newdat <- mnar_5_reg[is.na(mnar_5_reg$alcohol), ]
  predmnar5 <- predict(mnar_reg, newdata=newdat)
  mnar_5_reg$alcohol <- replace(mnar_5_reg$alcohol, is.na(mnar_5_reg$alcohol), 
                                round(predmnar5, 2))
  ts[i] <- t.test(mnar_5_reg$alcohol, wine$alcohol)$p.value
  Fs[i] <- var.test(mnar_5_reg$alcohol, wine$alcohol)$p.value
  RMSE[i] <- rmse(wine$alcohol, mnar_5_reg$alcohol)
  s <- s+1
}
length(ts[ts<=0.05]) # 1
length(Fs[Fs<=0.05]) # 100
RMSEmu <- round(mean(RMSE), 2)
RMSEmu               # 0.57
RMSEci <- round(c(quantile(RMSE, 0.025), quantile(RMSE, 0.975)), 2) 
RMSEci               # 0.52, 0.62


# univariate, MNAR, p=0.7, regression-based imputation 
nsim <- 100
ts <- rep(0, nsim)
Fs <- rep(0, nsim)
RMSE <- rep(0, nsim)
s <- 1157

for (i in 1:nsim) {
  set.seed(s)
  MNAR_7 <- delete_MNAR_censoring(ds=wine, p=0.7, cols_mis='alcohol')
  mnar_7_reg <- MNAR_7
  mnar_reg <- lm(alcohol ~ ash, data=mnar_7_reg)
  newdat <- mnar_7_reg[is.na(mnar_7_reg$alcohol), ]
  predmnar7 <- predict(mnar_reg, newdata=newdat)
  mnar_7_reg$alcohol <- replace(mnar_7_reg$alcohol, is.na(mnar_7_reg$alcohol), 
                                round(predmnar7, 2))
  ts[i] <- t.test(mnar_7_reg$alcohol, wine$alcohol)$p.value
  Fs[i] <- var.test(mnar_7_reg$alcohol, wine$alcohol)$p.value
  RMSE[i] <- rmse(wine$alcohol, mnar_7_reg$alcohol)
  s <- s+1
}
length(ts[ts<=0.05]) # 15
length(Fs[Fs<=0.05]) # 100
RMSEmu <- round(mean(RMSE), 2)
RMSEmu               # 0.68
RMSEci <- round(c(quantile(RMSE, 0.025), quantile(RMSE, 0.975)), 2) 
RMSEci               # 0.65, 0.75


#################### HOT DECK IMPUTATION #####################################
# univariate, MNAR, p=0.1, hot deck imputation
nsim <- 100
ts <- rep(0, nsim)
Fs <- rep(0, nsim)
RMSE <- rep(0, nsim)
s <- 1157

for (i in 1:nsim) {
  set.seed(s)
  MNAR_1 <- delete_MNAR_censoring(ds=wine, p=0.1, cols_mis='alcohol')
  mnar_1_hd <- MNAR_1
  mnar_1_hd$alcohol <- replace(mnar_1_hd$alcohol, is.na(mnar_1_hd$alcohol),
                               sample(mnar_1_hd$alcohol[!is.na(mnar_1_hd$alcohol)],
                                      length(mnar_1_hd$alcohol[is.na(mnar_1_hd$alcohol)]),
                                      replace=TRUE))
  ts[i] <- t.test(mnar_1_hd$alcohol, wine$alcohol)$p.value
  Fs[i] <- var.test(mnar_1_hd$alcohol, wine$alcohol)$p.value
  RMSE[i] <- rmse(wine$alcohol, mnar_1_hd$alcohol)
  s <- s+1
}
length(ts[ts<=0.05]) # 0
length(Fs[Fs<=0.05]) # 0
RMSEmu <- round(mean(RMSE), 2)
RMSEmu               # 0.37
RMSEci <- round(c(quantile(RMSE, 0.025), quantile(RMSE, 0.975)), 2) 
RMSEci               # 0.28, 0.47


# univariate, MNAR, p=0.3, hot deck imputation
nsim <- 100
ts <- rep(0, nsim)
Fs <- rep(0, nsim)
RMSE <- rep(0, nsim)
s <- 1157

for (i in 1:nsim) {
  set.seed(s)
  MNAR_3 <- delete_MNAR_censoring(ds=wine, p=0.3, cols_mis='alcohol')
  mnar_3_hd <- MNAR_3
  mnar_3_hd$alcohol <- replace(mnar_3_hd$alcohol, is.na(mnar_3_hd$alcohol),
                               sample(mnar_3_hd$alcohol[!is.na(mnar_3_hd$alcohol)],
                                      length(mnar_3_hd$alcohol[is.na(mnar_3_hd$alcohol)]),
                                      replace=TRUE))
  ts[i] <- t.test(mnar_3_hd$alcohol, wine$alcohol)$p.value
  Fs[i] <- var.test(mnar_3_hd$alcohol, wine$alcohol)$p.value
  RMSE[i] <- rmse(wine$alcohol, mnar_3_hd$alcohol)
  s <- s+1
}
length(ts[ts<=0.05]) # 0
length(Fs[Fs<=0.05]) # 0
RMSEmu <- round(mean(RMSE), 2)
RMSEmu               # 0.62
RMSEci <- round(c(quantile(RMSE, 0.025), quantile(RMSE, 0.975)), 2) 
RMSEci               # 0.53, 0.70


# univariate, MNAR, p=0.5, hot deck imputation 
nsim <- 100
ts <- rep(0, nsim)
Fs <- rep(0, nsim)
RMSE <- rep(0, nsim)
s <- 1157

for (i in 1:nsim) {
  set.seed(s)
  MNAR_5 <- delete_MNAR_censoring(ds=wine, p=0.5, cols_mis='alcohol')
  mnar_5_hd <- MNAR_5
  mnar_5_hd$alcohol <- replace(mnar_5_hd$alcohol, is.na(mnar_5_hd$alcohol),
                               sample(mnar_5_hd$alcohol[!is.na(mnar_5_hd$alcohol)],
                                      length(mnar_5_hd$alcohol[is.na(mnar_5_hd$alcohol)]),
                                      replace=TRUE))
  ts[i] <- t.test(mnar_5_hd$alcohol, wine$alcohol)$p.value
  Fs[i] <- var.test(mnar_5_hd$alcohol, wine$alcohol)$p.value
  RMSE[i] <- rmse(wine$alcohol, mnar_5_hd$alcohol)
  s <- s+1
}
length(ts[ts<=0.05]) # 3
length(Fs[Fs<=0.05]) # 1
RMSEmu <- round(mean(RMSE), 2)
RMSEmu               # 0.81
RMSEci <- round(c(quantile(RMSE, 0.025), quantile(RMSE, 0.975)), 2) 
RMSEci               # 0.71, 0.91


# univariate, MNAR, p=0.7, hot deck imputation
nsim <- 100
ts <- rep(0, nsim)
Fs <- rep(0, nsim)
RMSE <- rep(0, nsim)
s <- 1157

for (i in 1:nsim) {
  set.seed(s)
  MNAR_7 <- delete_MNAR_censoring(ds=wine, p=0.7, cols_mis='alcohol')
  mnar_7_hd <- MNAR_7
  mnar_7_hd$alcohol <- replace(mnar_7_hd$alcohol, is.na(mnar_7_hd$alcohol),
                               sample(mnar_7_hd$alcohol[!is.na(mnar_7_hd$alcohol)],
                                      length(mnar_7_hd$alcohol[is.na(mnar_7_hd$alcohol)]),
                                      replace=TRUE))
  ts[i] <- t.test(mnar_7_hd$alcohol, wine$alcohol)$p.value
  Fs[i] <- var.test(mnar_7_hd$alcohol, wine$alcohol)$p.value
  RMSE[i] <- rmse(wine$alcohol, mnar_7_hd$alcohol)
  s <- s+1
}
length(ts[ts<=0.05]) # 14
length(Fs[Fs<=0.05]) # 4
RMSEmu <- round(mean(RMSE), 2)
RMSEmu               # 0.95
RMSEci <- round(c(quantile(RMSE, 0.025), quantile(RMSE, 0.975)), 2) 
RMSEci               # 0.85, 1.05


################# IMPUTATION FROM A CONDITIONAL DISTRIBUTION ##################
# univariate, MNAR, p=0.1, imputation from a conditional distribution 
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
  MNAR_1 <- delete_MNAR_censoring(ds=wine, p=0.1, cols_mis='alcohol')
  # Priors
  phi <- round(mean(MNAR_1$alcohol, na.rm=TRUE), 2)
  tau2 <- 1
  alpha <- 0.1
  beta <- 0.1
  
  # simulate data
  n <- length(MNAR_1$alcohol[is.na(MNAR_1$alcohol)])
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
  mnar_1_cd <- MNAR_1
  mnar_1_cd$alcohol <- replace(mnar_1_cd$alcohol, is.na(mnar_1_cd$alcohol),
                               round(rnorm(length(mnar_1_cd$alcohol[is.na(mnar_1_cd$alcohol)]), 
                                           post_mu, post_sig), 2))
  
  # tests
  ts[i] <- t.test(wine$alcohol, mnar_1_cd$alcohol)
  Fs[i] <- var.test(wine$alcohol, mnar_1_cd$alcohol)
  RMSE[i] <- rmse(wine$alcohol, mnar_1_cd$alcohol)
  s <- s+1
}

length(ts[ts<=0.05]) # 51
length(Fs[Fs<=0.05]) # 0
RMSEmu <- round(mean(RMSE), 2)
RMSEmu               # 0.47
RMSEci <- round(c(quantile(RMSE, 0.025), quantile(RMSE, 0.975)), 2) 
RMSEci               # 0.28, 0.76


# univariate, MNAR, p=0.3, imputation from a conditonal distribution 
nsim <- 100
ts <- rep(0, nsim)
Fs <- rep(0, nsim)
RMSE <- rep(0, nsim)
s <- 1157

for (i in 1:nsim) {
  set.seed(s)
  MNAR_3 <- delete_MNAR_censoring(ds=wine, p=0.3, cols_mis='alcohol')
  # Priors
  phi <- round(mean(MNAR_3$alcohol, na.rm=TRUE), 2)
  tau2 <- 1
  alpha <- 0.1
  beta <- 0.1
  
  # simulate data
  n <- length(MNAR_3$alcohol[is.na(MNAR_3$alcohol)])
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
  mnar_3_cd <- MNAR_3
  mnar_3_cd$alcohol <- replace(mnar_3_cd$alcohol, is.na(mnar_3_cd$alcohol),
                               round(rnorm(length(mnar_3_cd$alcohol[is.na(mnar_3_cd$alcohol)]), 
                                           post_mu, post_sig), 2))
  
  # tests
  ts[i] <- t.test(wine$alcohol, mnar_3_cd$alcohol)
  Fs[i] <- var.test(wine$alcohol, mnar_3_cd$alcohol)
  RMSE[i] <- rmse(wine$alcohol, mnar_3_cd$alcohol)
  s <- s+1
}

length(ts[ts<=0.05]) # 61
length(Fs[Fs<=0.05]) # 0
RMSEmu <- round(mean(RMSE), 2)
RMSEmu               # 0.72
RMSEci <- round(c(quantile(RMSE, 0.025), quantile(RMSE, 0.975)), 2) 
RMSEci               # 0.56, 1.01


# univariate, MNAR, p=0.5, imputation from a conditonal distribution 
nsim <- 100
ts <- rep(0, nsim)
Fs <- rep(0, nsim)
RMSE <- rep(0, nsim)
s <- 1157

for (i in 1:nsim) {
  set.seed(s)
  MNAR_5 <- delete_MNAR_censoring(ds=wine, p=0.5, cols_mis='alcohol')
  # Priors
  phi <- round(mean(MNAR_5$alcohol, na.rm=TRUE), 2)
  tau2 <- 1
  alpha <- 0.1
  beta <- 0.1
  
  # simulate data
  n <- length(MNAR_5$alcohol[is.na(MNAR_5$alcohol)])
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
  mnar_5_cd <- MNAR_5
  mnar_5_cd$alcohol <- replace(mnar_5_cd$alcohol, is.na(mnar_5_cd$alcohol),
                               round(rnorm(length(mnar_5_cd$alcohol[is.na(mnar_5_cd$alcohol)]), 
                                           post_mu, post_sig), 2))
  
  # tests
  ts[i] <- t.test(wine$alcohol, mnar_5_cd$alcohol)
  Fs[i] <- var.test(wine$alcohol, mnar_5_cd$alcohol)
  RMSE[i] <- rmse(wine$alcohol, mnar_5_cd$alcohol)
  s <- s+1
}

length(ts[ts<=0.05]) # 45
length(Fs[Fs<=0.05]) # 0
RMSEmu <- round(mean(RMSE), 2)
RMSEmu               # 0.93
RMSEci <- round(c(quantile(RMSE, 0.025), quantile(RMSE, 0.975)), 2) 
RMSEci               # 0.75, 1.16


# univariate, MNAR, p=0.7, imputation from a conditional distribution 
nsim <- 100
ts <- rep(0, nsim)
Fs <- rep(0, nsim)
RMSE <- rep(0, nsim)
s <- 1157

for (i in 1:nsim) {
  set.seed(s)
  MNAR_7 <- delete_MNAR_censoring(ds=wine, p=0.7, cols_mis='alcohol')
  # Priors
  phi <- round(mean(MNAR_7$alcohol, na.rm=TRUE), 2)
  tau2 <- 1
  alpha <- 0.1
  beta <- 0.1
  
  # simulate data
  n <- length(MNAR_7$alcohol[is.na(MNAR_7$alcohol)])
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
  mnar_7_cd <- MNAR_7
  mnar_7_cd$alcohol <- replace(mnar_7_cd$alcohol, is.na(mnar_7_cd$alcohol),
                               round(rnorm(length(mnar_7_cd$alcohol[is.na(mnar_7_cd$alcohol)]), 
                                           post_mu, post_sig), 2))
  
  # tests
  ts[i] <- t.test(wine$alcohol, mnar_7_cd$alcohol)
  Fs[i] <- var.test(wine$alcohol, mnar_7_cd$alcohol)
  RMSE[i] <- rmse(wine$alcohol, mnar_7_cd$alcohol)
  s <- s+1
}

length(ts[ts<=0.05]) # 51
length(Fs[Fs<=0.05]) # 0
RMSEmu <- round(mean(RMSE), 2)
RMSEmu               # 1.07
RMSEci <- round(c(quantile(RMSE, 0.025), quantile(RMSE, 0.975)), 2) 
RMSEci               # 0.91, 1.29


################################ IPW ##########################################
# univariate, MNAR, p=0.1, IPW
nsim <- 100
ts <- rep(0, nsim)
Fs <- rep(0, nsim)
RMSE <- rep(0, nsim)
s <- 1157

for (i in 1:nsim) {
  set.seed(s)
  MNAR_1 <- delete_MNAR_censoring(ds=wine, p=0.1, cols_mis='alcohol')
  mnar_1_ipw <- MNAR_1
  # add two columns to the data: 1) 0/1 missingness, 2) weight
  mnar_1_ipw$R <- rep(0, 178)
  mnar_1_ipw$weight <- rep(0, 178)
  for (j in 1:178) {
    if (is.na(mnar_1_ipw$alcohol[j])) {
      mnar_1_ipw$R[j] <- 1
    }
  }
  
  # logistic regression on missingness to determine weights 
  # ipw <- lrm(R ~ ash + `color intensity` + magnesium + `total phenols` + `flavanoids`, data=mnar_1_ipw)
  ipw <- lrm(R ~ ash, data=mnar_1_ipw)
  ipw_newdat <- mnar_1_ipw
  
  # calculate weights based on logistic regression model
  weight <- predict(ipw, ipw_newdat)
  mnar_1_ipw$weight <- weight
  
  # calculate weighted mean
  wmean <- round(weighted.mean(mnar_1_ipw$alcohol[!is.na(mnar_1_ipw$alcohol)],
                               mnar_1_ipw$weight[!is.na(mnar_1_ipw$alcohol)]), 2)
  
  # replace missing values with weighted mean
  mnar_1_ipw$alcohol <- replace(mnar_1_ipw$alcohol, is.na(mnar_1_ipw$alcohol), wmean)
  
  ts[i] <- t.test(mnar_1_ipw$alcohol, wine$alcohol)$p.value
  Fs[i] <- var.test(mnar_1_ipw$alcohol, wine$alcohol)$p.value
  RMSE[i] <- rmse(wine$alcohol, mnar_1_ipw$alcohol)
  s <- s+1
}

length(ts[ts<=0.05]) # 0
length(Fs[Fs<=0.05]) # 0 
RMSEmu <- round(mean(RMSE), 2)
RMSEmu               # 0.26
RMSEci <- round(c(quantile(RMSE, 0.025), quantile(RMSE, 0.975)), 2) 
RMSEci               # 0.21, 0.32


# univariate, MNAR, p=0.3, IPW
nsim <- 100
ts <- rep(0, nsim)
Fs <- rep(0, nsim)
RMSE <- rep(0, nsim)
s <- 1157

for (i in 1:nsim) {
  set.seed(s)
  MNAR_3 <- delete_MNAR_censoring(ds=wine, p=0.3, cols_mis='alcohol')
  mnar_3_ipw <- MNAR_3
  # add two columns to the data: 1) 0/1 missingness, 2) weight
  mnar_3_ipw$R <- rep(0, 178)
  mnar_3_ipw$weight <- rep(0, 178)
  for (j in 1:178) {
    if (is.na(mnar_3_ipw$alcohol[j])) {
      mnar_3_ipw$R[j] <- 1
    }
  }
  
  # logistic regression on missingness to determine weights 
  # ipw <- lrm(R ~ ash + `color intensity` + magnesium + `total phenols` + `flavanoids`, data=mnar_3_ipw)
  ipw <- lrm(R ~ ash, data=mnar_3_ipw)
  ipw_newdat <- mnar_3_ipw
  
  # calculate weights based on logistic regression model
  weight <- predict(ipw, ipw_newdat)
  mnar_3_ipw$weight <- weight
  
  # calculate weighted mean
  wmean <- round(weighted.mean(mnar_3_ipw$alcohol[!is.na(mnar_3_ipw$alcohol)],
                               mnar_3_ipw$weight[!is.na(mnar_3_ipw$alcohol)]), 2)
  
  # replace missing values with weighted mean
  mnar_3_ipw$alcohol <- replace(mnar_3_ipw$alcohol, is.na(mnar_3_ipw$alcohol), wmean)
  
  ts[i] <- t.test(mnar_3_ipw$alcohol, wine$alcohol)$p.value
  Fs[i] <- var.test(mnar_3_ipw$alcohol, wine$alcohol)$p.value
  RMSE[i] <- rmse(wine$alcohol, mnar_3_ipw$alcohol)
  s <- s+1
}

length(ts[ts<=0.05]) # 0
length(Fs[Fs<=0.05]) # 92
RMSEmu <- round(mean(RMSE), 2)
RMSEmu               # 0.45
RMSEci <- round(c(quantile(RMSE, 0.025), quantile(RMSE, 0.975)), 2) 
RMSEci               # 0.39, 0.49


# univariate, MNAR, p=0.5, IPW
nsim <- 100
ts <- rep(0, nsim)
Fs <- rep(0, nsim)
RMSE <- rep(0, nsim)
s <- 1157

for (i in 1:nsim) {
  set.seed(s)
  MNAR_5 <- delete_MNAR_censoring(ds=wine, p=0.5, cols_mis='alcohol')
  mnar_5_ipw <- MNAR_5
  # add two columns to the data: 1) 0/1 missingness, 2) weight
  mnar_5_ipw$R <- rep(0, 178)
  mnar_5_ipw$weight <- rep(0, 178)
  for (j in 1:178) {
    if (is.na(mnar_5_ipw$alcohol[j])) {
      mnar_5_ipw$R[j] <- 1
    }
  }
  
  # logistic regression on missingness to determine weights 
  # ipw <- lrm(R ~ ash, data=mnar_5_ipw)
  ipw <- lrm(R ~ ash, data=mnar_5_ipw)
  ipw_newdat <- mnar_5_ipw
  
  # calculate weights based on logistic regression model
  weight <- predict(ipw, ipw_newdat)
  mnar_5_ipw$weight <- weight
  
  # calculate weighted mean
  wmean <- round(weighted.mean(mnar_5_ipw$alcohol[!is.na(mnar_5_ipw$alcohol)],
                               mnar_5_ipw$weight[!is.na(mnar_5_ipw$alcohol)]), 2)
  
  # replace missing values with weighted mean
  mnar_5_ipw$alcohol <- replace(mnar_5_ipw$alcohol, is.na(mnar_5_ipw$alcohol), wmean)
  
  ts[i] <- t.test(mnar_5_ipw$alcohol, wine$alcohol)$p.value
  Fs[i] <- var.test(mnar_5_ipw$alcohol, wine$alcohol)$p.value
  RMSE[i] <- rmse(wine$alcohol, mnar_5_ipw$alcohol)
  s <- s+1
}

length(ts[ts<=0.05]) # 98
length(Fs[Fs<=0.05]) # 85
RMSEmu <- round(mean(RMSE), 2)
RMSEmu               # 7.45
RMSEci <- round(c(quantile(RMSE, 0.025), quantile(RMSE, 0.975)), 2) 
RMSEci               # 0.87, 51.01


# univariate, MNAR, p=0.7, IPW
nsim <- 100
ts <- rep(0, nsim)
Fs <- rep(0, nsim)
RMSE <- rep(0, nsim)
s <- 1157

for (i in 1:nsim) {
  set.seed(s)
  MNAR_7 <- delete_MNAR_censoring(ds=wine, p=0.7, cols_mis='alcohol')
  mnar_7_ipw <- MNAR_7
  # add two columns to the data: 1) 0/1 missingness, 2) weight
  mnar_7_ipw$R <- rep(0, 178)
  mnar_7_ipw$weight <- rep(0, 178)
  for (j in 1:178) {
    if (is.na(mnar_7_ipw$alcohol[j])) {
      mnar_7_ipw$R[j] <- 1
    }
  }
  
  # logistic regression on missingness to determine weights 
  # ipw <- lrm(R ~ ash, data=mnar_7_ipw)
  ipw <- lrm(R ~ ash, data=mnar_7_ipw)
  ipw_newdat <- mnar_7_ipw
  
  # calculate weights based on logistic regression model
  weight <- predict(ipw, ipw_newdat)
  mnar_7_ipw$weight <- weight
  
  # calculate weighted mean
  wmean <- round(weighted.mean(mnar_7_ipw$alcohol[!is.na(mnar_7_ipw$alcohol)],
                               mnar_7_ipw$weight[!is.na(mnar_7_ipw$alcohol)]), 2)
  
  # replace missing values with weighted mean
  mnar_7_ipw$alcohol <- replace(mnar_7_ipw$alcohol, is.na(mnar_7_ipw$alcohol), wmean)
  
  ts[i] <- t.test(mnar_7_ipw$alcohol, wine$alcohol)$p.value
  Fs[i] <- var.test(mnar_7_ipw$alcohol, wine$alcohol)$p.value
  RMSE[i] <- rmse(wine$alcohol, mnar_7_ipw$alcohol)
  s <- s+1
}

length(ts[ts<=0.05]) # 100
length(Fs[Fs<=0.05]) # 100
RMSEmu <- round(mean(RMSE), 2)
RMSEmu               # 1.25
RMSEci <- round(c(quantile(RMSE, 0.025), quantile(RMSE, 0.975)), 2) 
RMSEci               # 1.25, 1.25