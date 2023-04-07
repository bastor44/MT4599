# MT4599 - Chapter 4: Comparison in Practice (MCAR)

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
# univariate, MCAR, p=0.1, CCA
nsim <- 100
ts <- rep(0, nsim)
Fs <- rep(0, nsim)
s <- 1157

for (i in 1:nsim) {
  set.seed(s)
  MCAR_1 <- delete_MCAR(ds=wine, p=0.1, cols_mis='alcohol')
  mcar_1_cca <- na.omit(MCAR_1)
  ts[i] <- t.test(mcar_1_cca$alcohol, wine$alcohol)$p.value
  Fs[i] <- var.test(mcar_1_cca$alcohol, wine$alcohol)$p.value
  # cannot calculate RMSE for CCA
  s <- s+1
}
length(ts[ts<=0.05]) # 0
length(Fs[Fs<=0.05]) # 0


# univariate, MCAR, p=0.3, CCA
nsim <- 100
ts <- rep(0, nsim)
Fs <- rep(0, nsim)
s <- 1157

for (i in 1:nsim) {
  set.seed(s)
  MCAR_3 <- delete_MCAR(ds=wine, p=0.3, cols_mis='alcohol')
  mcar_3_cca <- na.omit(MCAR_3)
  ts[i] <- t.test(mcar_1_cca$alcohol, wine$alcohol)$p.value
  Fs[i] <- var.test(mcar_1_cca$alcohol, wine$alcohol)$p.value
  # cannot calculate RMSE for CCA
  s <- s+1
}
length(ts[ts<=0.05]) # 0
length(Fs[Fs<=0.05]) # 0


# univariate, MCAR, p=0.5, CCA
nsim <- 100
ts <- rep(0, nsim)
Fs <- rep(0, nsim)
s <- 1157

for (i in 1:nsim) {
  set.seed(s)
  MCAR_5 <- delete_MCAR(ds=wine, p=0.5, cols_mis='alcohol')
  mcar_5_cca <- na.omit(MCAR_5)
  ts[i] <- t.test(mcar_5_cca$alcohol, wine$alcohol)$p.value
  Fs[i] <- var.test(mcar_5_cca$alcohol, wine$alcohol)$p.value
  # cannot calculate RMSE for CCA
  s <- s+1
}
length(ts[ts<=0.05]) # 0
length(Fs[Fs<=0.05]) # 0


# univariate, MCAR, p=0.7, CCA
nsim <- 100
ts <- rep(0, nsim)
Fs <- rep(0, nsim)
s <- 1157

for (i in 1:nsim) {
  set.seed(s)
  MCAR_7 <- delete_MCAR(ds=wine, p=0.7, cols_mis='alcohol')
  mcar_7_cca <- na.omit(MCAR_7)
  ts[i] <- t.test(mcar_7_cca$alcohol, wine$alcohol)$p.value
  Fs[i] <- var.test(mcar_7_cca$alcohol, wine$alcohol)$p.value
  # cannot calculate RMSE for CCA
  s <- s+1
}
length(ts[ts<=0.05]) # 1
length(Fs[Fs<=0.05]) # 0



########################## MEAN IMPUTATION ###################################
# univariate, MCAR, p=0.1, Mean imputation
nsim <- 100
ts <- rep(0, nsim)
Fs <- rep(0, nsim)
RMSE <- rep(0, nsim)
s <- 1157

for (i in 1:nsim) {
  set.seed(s)
  MCAR_1 <- delete_MCAR(ds=wine, p=0.1, cols_mis='alcohol')
  mcar_1_mi <- MCAR_1
  mcar_1_mi$alcohol <- replace(mcar_1_mi$alcohol, is.na(mcar_1_mi$alcohol),
                               round(mean(mcar_1_mi$alcohol, na.rm=TRUE), 2))
  ts[i] <- t.test(mcar_1_mi$alcohol, wine$alcohol)$p.value
  Fs[i] <- var.test(mcar_1_mi$alcohol, wine$alcohol)$p.value
  RMSE[i] <- rmse(wine$alcohol, mcar_1_mi$alcohol)
  s <- s+1
}
length(ts[ts<=0.05]) # 0
length(Fs[Fs<=0.05]) # 0
RMSEmu <- round(mean(RMSE), 2) 
RMSEmu               # 0.26
RMSEci <- round(c(quantile(RMSE, 0.025), quantile(RMSE, 0.975)), 2)
RMSEci               # 0.21, 0.32
# note: acceptable when between 0.2 and 0.5 (rule of thumb, best used comparatively)

# univariate, MCAR, p=0.3, mean imputation 
nsim <- 100
ts <- rep(0, nsim)
Fs <- rep(0, nsim)
RMSE <- rep(0, nsim)
s <- 1157

for (i in 1:nsim) {
  set.seed(s)
  MCAR_3 <- delete_MCAR(ds=wine, p=0.3, cols_mis='alcohol')
  mcar_3_mi <- MCAR_3
  mcar_3_mi$alcohol <- replace(mcar_3_mi$alcohol, is.na(mcar_3_mi$alcohol),
                               round(mean(mcar_3_mi$alcohol, na.rm=TRUE), 2))
  ts[i] <- t.test(mcar_3_mi$alcohol, wine$alcohol)$p.value
  Fs[i] <- var.test(mcar_3_mi$alcohol, wine$alcohol)$p.value
  RMSE[i] <- rmse(wine$alcohol, mcar_3_mi$alcohol)
  s <- s+1
}
length(ts[ts<=0.05]) # 0
length(Fs[Fs<=0.05]) # 92
RMSEmu <- round(mean(RMSE), 2)
RMSEmu               # 0.45
RMSEci <- round(c(quantile(RMSE, 0.025), quantile(RMSE, 0.975)), 2) 
RMSEci               # 0.39, 0.49

# univariate, MCAR, p=0.5, mean imputation 
nsim <- 100
ts <- rep(0, nsim)
Fs <- rep(0, nsim)
RMSE <- rep(0, nsim)
s <- 1157

for (i in 1:nsim) {
  set.seed(s)
  MCAR_5 <- delete_MCAR(ds=wine, p=0.5, cols_mis='alcohol')
  mcar_5_mi <- MCAR_5
  mcar_5_mi$alcohol <- replace(mcar_5_mi$alcohol, is.na(mcar_5_mi$alcohol),
                               round(mean(mcar_5_mi$alcohol, na.rm=TRUE), 2))
  ts[i] <- t.test(mcar_5_mi$alcohol, wine$alcohol)$p.value
  Fs[i] <- var.test(mcar_5_mi$alcohol, wine$alcohol)$p.value
  RMSE[i] <- rmse(wine$alcohol, mcar_5_mi$alcohol)
  s <- s+1
}
length(ts[ts<=0.05]) # 1
length(Fs[Fs<=0.05]) # 100
RMSEmu <- round(mean(RMSE), 2)
RMSEmu               # 0.58
RMSEci <- round(c(quantile(RMSE, 0.025), quantile(RMSE, 0.975)), 2) 
RMSEci               # 0.53, 0.62

# univariate MCAR, p=0.7, mean imputation 
nsim <- 100
ts <- rep(0, nsim)
Fs <- rep(0, nsim)
RMSE <- rep(0, nsim)
s <- 1157 

for (i in 1:nsim) {
  set.seed(s)
  MCAR_7 <- delete_MCAR(ds=wine, p=0.7, cols_mis='alcohol')
  mcar_7_mi <- MCAR_7
  mcar_7_mi$alcohol <- replace(mcar_7_mi$alcohol, is.na(mcar_7_mi$alcohol),
                               round(mean(mcar_7_mi$alcohol, na.rm=TRUE), 2))
  ts[i] <- t.test(mcar_7_mi$alcohol, wine$alcohol)$p.value
  Fs[i] <- var.test(mcar_7_mi$alcohol, wine$alcohol)$p.value
  RMSE[i] <- rmse(wine$alcohol, mcar_7_mi$alcohol)
  s <- s+1
}
length(ts[ts<=0.05]) # 17
length(Fs[Fs<=0.05]) # 100
RMSEmu <- round(mean(RMSE), 2)
RMSEmu               # 0.69
RMSEci <- round(c(quantile(RMSE, 0.025), quantile(RMSE, 0.975)), 2) 
RMSEci               # 0.66, 0.72


########################## REGRESSION-BASED IMPUTATION ########################
# univariate, MCAR, p=0.1, regression-based imputation 
nsim <- 100
ts <- rep(0, nsim)
Fs <- rep(0, nsim)
RMSE <- rep(0, nsim)
s <- 1157

for (i in 1:nsim) {
  set.seed(s)
  MCAR_1 <- delete_MCAR(ds=wine, p=0.1, cols_mis='alcohol')
  mcar_1_reg <- MCAR_1
  mcar_reg <- lm(alcohol ~ ash, data=mcar_1_reg)
  newdat <- mcar_1_reg[is.na(mcar_1_reg$alcohol), ]
  predmcar1 <- predict(mcar_reg, newdata=newdat)
  mcar_1_reg$alcohol <- replace(mcar_1_reg$alcohol, is.na(mcar_1_reg$alcohol), 
                                round(predmcar1, 2))
  ts[i] <- t.test(mcar_1_reg$alcohol, wine$alcohol)$p.value
  Fs[i] <- var.test(mcar_1_reg$alcohol, wine$alcohol)$p.value
  RMSE[i] <- rmse(wine$alcohol, mcar_1_reg$alcohol)
  s <- s+1
}
length(ts[ts<=0.05]) # 0
length(Fs[Fs<=0.05]) # 0
RMSEmu <- round(mean(RMSE), 2)
RMSEmu               # 0.26
RMSEci <- round(c(quantile(RMSE, 0.025), quantile(RMSE, 0.975)), 2) 
RMSEci               # 0.20, 0.32


# univariate, MCAR, p=0.3, regression-based imputation
nsim <- 100
ts <- rep(0, nsim)
Fs <- rep(0, nsim)
RMSE <- rep(0, nsim)
s <- 1157

for (i in 1:nsim) {
  set.seed(s)
  MCAR_3 <- delete_MCAR(ds=wine, p=0.3, cols_mis='alcohol')
  mcar_3_reg <- MCAR_3
  mcar_reg <- lm(alcohol ~ ash, data=mcar_3_reg)
  newdat <- mcar_3_reg[is.na(mcar_3_reg$alcohol), ]
  predmcar3 <- predict(mcar_reg, newdata=newdat)
  mcar_3_reg$alcohol <- replace(mcar_3_reg$alcohol, is.na(mcar_3_reg$alcohol), 
                                round(predmcar3, 2))
  ts[i] <- t.test(mcar_3_reg$alcohol, wine$alcohol)$p.value
  Fs[i] <- var.test(mcar_3_reg$alcohol, wine$alcohol)$p.value
  RMSE[i] <- rmse(wine$alcohol, mcar_3_reg$alcohol)
  s <- s+1
}
length(ts[ts<=0.05]) # 0
length(Fs[Fs<=0.05]) # 84
RMSEmu <- round(mean(RMSE), 2)
RMSEmu               # 0.44
RMSEci <- round(c(quantile(RMSE, 0.025), quantile(RMSE, 0.975)), 2) 
RMSEci               # 0.38, 0.49


# univariate, MCAR, p=0.5, regression-based imputation 
nsim <- 100
ts <- rep(0, nsim)
Fs <- rep(0, nsim)
RMSE <- rep(0, nsim)
s <- 1157

for (i in 1:nsim) {
  set.seed(s)
  MCAR_5 <- delete_MCAR(ds=wine, p=0.5, cols_mis='alcohol')
  mcar_5_reg <- MCAR_5
  mcar_reg <- lm(alcohol ~ ash, data=mcar_5_reg)
  newdat <- mcar_5_reg[is.na(mcar_5_reg$alcohol), ]
  predmcar5 <- predict(mcar_reg, newdata=newdat)
  mcar_5_reg$alcohol <- replace(mcar_5_reg$alcohol, is.na(mcar_5_reg$alcohol), 
                                round(predmcar5, 2))
  ts[i] <- t.test(mcar_5_reg$alcohol, wine$alcohol)$p.value
  Fs[i] <- var.test(mcar_5_reg$alcohol, wine$alcohol)$p.value
  RMSE[i] <- rmse(wine$alcohol, mcar_5_reg$alcohol)
  s <- s+1
}
length(ts[ts<=0.05]) # 1
length(Fs[Fs<=0.05]) # 100
RMSEmu <- round(mean(RMSE), 2)
RMSEmu               # 0.57
RMSEci <- round(c(quantile(RMSE, 0.025), quantile(RMSE, 0.975)), 2) 
RMSEci               # 0.52, 0.62


# univariate, MCAR, p=0.7, regression-based imputation 
nsim <- 100
ts <- rep(0, nsim)
Fs <- rep(0, nsim)
RMSE <- rep(0, nsim)
s <- 1157

for (i in 1:nsim) {
  set.seed(s)
  MCAR_7 <- delete_MCAR(ds=wine, p=0.7, cols_mis='alcohol')
  mcar_7_reg <- MCAR_7
  mcar_reg <- lm(alcohol ~ ash, data=mcar_7_reg)
  newdat <- mcar_7_reg[is.na(mcar_7_reg$alcohol), ]
  predmcar7 <- predict(mcar_reg, newdata=newdat)
  mcar_7_reg$alcohol <- replace(mcar_7_reg$alcohol, is.na(mcar_7_reg$alcohol), 
                                round(predmcar7, 2))
  ts[i] <- t.test(mcar_7_reg$alcohol, wine$alcohol)$p.value
  Fs[i] <- var.test(mcar_7_reg$alcohol, wine$alcohol)$p.value
  RMSE[i] <- rmse(wine$alcohol, mcar_7_reg$alcohol)
  s <- s+1
}
length(ts[ts<=0.05]) # 15
length(Fs[Fs<=0.05]) # 100
RMSEmu <- round(mean(RMSE), 2)
RMSEmu               # 0.68
RMSEci <- round(c(quantile(RMSE, 0.025), quantile(RMSE, 0.975)), 2) 
RMSEci               # 0.65, 0.75


#################### HOT DECK IMPUTATION #####################################
# univariate, MCAR, p=0.1, hot deck imputation
nsim <- 100
ts <- rep(0, nsim)
Fs <- rep(0, nsim)
RMSE <- rep(0, nsim)
s <- 1157

for (i in 1:nsim) {
  set.seed(s)
  MCAR_1 <- delete_MCAR(ds=wine, p=0.1, cols_mis='alcohol')
  mcar_1_hd <- MCAR_1
  mcar_1_hd$alcohol <- replace(mcar_1_hd$alcohol, is.na(mcar_1_hd$alcohol),
                               sample(mcar_1_hd$alcohol[!is.na(mcar_1_hd$alcohol)],
                                      length(mcar_1_hd$alcohol[is.na(mcar_1_hd$alcohol)]),
                                      replace=TRUE))
  ts[i] <- t.test(mcar_1_hd$alcohol, wine$alcohol)$p.value
  Fs[i] <- var.test(mcar_1_hd$alcohol, wine$alcohol)$p.value
  RMSE[i] <- rmse(wine$alcohol, mcar_1_hd$alcohol)
  s <- s+1
}
length(ts[ts<=0.05]) # 0
length(Fs[Fs<=0.05]) # 0
RMSEmu <- round(mean(RMSE), 2)
RMSEmu               # 0.37
RMSEci <- round(c(quantile(RMSE, 0.025), quantile(RMSE, 0.975)), 2) 
RMSEci               # 0.28, 0.47


# univariate, MCAR, p=0.3, hot deck imputation
nsim <- 100
ts <- rep(0, nsim)
Fs <- rep(0, nsim)
RMSE <- rep(0, nsim)
s <- 1157

for (i in 1:nsim) {
  set.seed(s)
  MCAR_3 <- delete_MCAR(ds=wine, p=0.3, cols_mis='alcohol')
  mcar_3_hd <- MCAR_3
  mcar_3_hd$alcohol <- replace(mcar_3_hd$alcohol, is.na(mcar_3_hd$alcohol),
                               sample(mcar_3_hd$alcohol[!is.na(mcar_3_hd$alcohol)],
                                      length(mcar_3_hd$alcohol[is.na(mcar_3_hd$alcohol)]),
                                      replace=TRUE))
  ts[i] <- t.test(mcar_3_hd$alcohol, wine$alcohol)$p.value
  Fs[i] <- var.test(mcar_3_hd$alcohol, wine$alcohol)$p.value
  RMSE[i] <- rmse(wine$alcohol, mcar_3_hd$alcohol)
  s <- s+1
}
length(ts[ts<=0.05]) # 0
length(Fs[Fs<=0.05]) # 0
RMSEmu <- round(mean(RMSE), 2)
RMSEmu               # 0.62
RMSEci <- round(c(quantile(RMSE, 0.025), quantile(RMSE, 0.975)), 2) 
RMSEci               # 0.53, 0.70


# univariate, MCAR, p=0.5, hot deck imputation 
nsim <- 100
ts <- rep(0, nsim)
Fs <- rep(0, nsim)
RMSE <- rep(0, nsim)
s <- 1157

for (i in 1:nsim) {
  set.seed(s)
  MCAR_5 <- delete_MCAR(ds=wine, p=0.5, cols_mis='alcohol')
  mcar_5_hd <- MCAR_5
  mcar_5_hd$alcohol <- replace(mcar_5_hd$alcohol, is.na(mcar_5_hd$alcohol),
                               sample(mcar_5_hd$alcohol[!is.na(mcar_5_hd$alcohol)],
                                      length(mcar_5_hd$alcohol[is.na(mcar_5_hd$alcohol)]),
                                      replace=TRUE))
  ts[i] <- t.test(mcar_5_hd$alcohol, wine$alcohol)$p.value
  Fs[i] <- var.test(mcar_5_hd$alcohol, wine$alcohol)$p.value
  RMSE[i] <- rmse(wine$alcohol, mcar_5_hd$alcohol)
  s <- s+1
}
length(ts[ts<=0.05]) # 3
length(Fs[Fs<=0.05]) # 1
RMSEmu <- round(mean(RMSE), 2)
RMSEmu               # 0.81
RMSEci <- round(c(quantile(RMSE, 0.025), quantile(RMSE, 0.975)), 2) 
RMSEci               # 0.71, 0.91


# univariate, MCAR, p=0.7, hot deck imputation
nsim <- 100
ts <- rep(0, nsim)
Fs <- rep(0, nsim)
RMSE <- rep(0, nsim)
s <- 1157

for (i in 1:nsim) {
  set.seed(s)
  MCAR_7 <- delete_MCAR(ds=wine, p=0.7, cols_mis='alcohol')
  mcar_7_hd <- MCAR_7
  mcar_7_hd$alcohol <- replace(mcar_7_hd$alcohol, is.na(mcar_7_hd$alcohol),
                               sample(mcar_7_hd$alcohol[!is.na(mcar_7_hd$alcohol)],
                                      length(mcar_7_hd$alcohol[is.na(mcar_7_hd$alcohol)]),
                                      replace=TRUE))
  ts[i] <- t.test(mcar_7_hd$alcohol, wine$alcohol)$p.value
  Fs[i] <- var.test(mcar_7_hd$alcohol, wine$alcohol)$p.value
  RMSE[i] <- rmse(wine$alcohol, mcar_7_hd$alcohol)
  s <- s+1
}
length(ts[ts<=0.05]) # 14
length(Fs[Fs<=0.05]) # 4
RMSEmu <- round(mean(RMSE), 2)
RMSEmu               # 0.95
RMSEci <- round(c(quantile(RMSE, 0.025), quantile(RMSE, 0.975)), 2) 
RMSEci               # 0.85, 1.05


################# IMPUTATION FROM A CONDITIONAL DISTRIBUTION ##################
# univariate, MCAR, p=0.1, imputation from a conditional distribution 
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
  MCAR_1 <- delete_MCAR(ds=wine, p=0.1, cols_mis='alcohol')
  # Priors
  phi <- round(mean(MCAR_1$alcohol, na.rm=TRUE), 2)
  tau2 <- 1
  alpha <- 0.1
  beta <- 0.1
  
  # simulate data
  n <- length(MCAR_1$alcohol[is.na(MCAR_1$alcohol)])
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
  mcar_1_cd <- MCAR_1
  mcar_1_cd$alcohol <- replace(mcar_1_cd$alcohol, is.na(mcar_1_cd$alcohol),
                               round(rnorm(length(mcar_1_cd$alcohol[is.na(mcar_1_cd$alcohol)]), 
                                           post_mu, post_sig), 2))
  
  # tests
  ts[i] <- t.test(wine$alcohol, mcar_1_cd$alcohol)
  Fs[i] <- var.test(wine$alcohol, mcar_1_cd$alcohol)
  RMSE[i] <- rmse(wine$alcohol, mcar_1_cd$alcohol)
  s <- s+1
}

length(ts[ts<=0.05]) # 51
length(Fs[Fs<=0.05]) # 0
RMSEmu <- round(mean(RMSE), 2)
RMSEmu               # 0.47
RMSEci <- round(c(quantile(RMSE, 0.025), quantile(RMSE, 0.975)), 2) 
RMSEci               # 0.28, 0.76


# univariate, MCAR, p=0.3, imputation from a conditonal distribution 
nsim <- 100
ts <- rep(0, nsim)
Fs <- rep(0, nsim)
RMSE <- rep(0, nsim)
s <- 1157

for (i in 1:nsim) {
  set.seed(s)
  MCAR_3 <- delete_MCAR(ds=wine, p=0.3, cols_mis='alcohol')
  # Priors
  phi <- round(mean(MCAR_3$alcohol, na.rm=TRUE), 2)
  tau2 <- 1
  alpha <- 0.1
  beta <- 0.1
  
  # simulate data
  n <- length(MCAR_3$alcohol[is.na(MCAR_3$alcohol)])
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
  mcar_3_cd <- MCAR_3
  mcar_3_cd$alcohol <- replace(mcar_3_cd$alcohol, is.na(mcar_3_cd$alcohol),
                               round(rnorm(length(mcar_3_cd$alcohol[is.na(mcar_3_cd$alcohol)]), 
                                           post_mu, post_sig), 2))
  
  # tests
  ts[i] <- t.test(wine$alcohol, mcar_3_cd$alcohol)
  Fs[i] <- var.test(wine$alcohol, mcar_3_cd$alcohol)
  RMSE[i] <- rmse(wine$alcohol, mcar_3_cd$alcohol)
  s <- s+1
}

length(ts[ts<=0.05]) # 61
length(Fs[Fs<=0.05]) # 0
RMSEmu <- round(mean(RMSE), 2)
RMSEmu               # 0.72
RMSEci <- round(c(quantile(RMSE, 0.025), quantile(RMSE, 0.975)), 2) 
RMSEci               # 0.56, 1.01


# univariate, MCAR, p=0.5, imputation from a conditonal distribution 
nsim <- 100
ts <- rep(0, nsim)
Fs <- rep(0, nsim)
RMSE <- rep(0, nsim)
s <- 1157

for (i in 1:nsim) {
  set.seed(s)
  MCAR_5 <- delete_MCAR(ds=wine, p=0.5, cols_mis='alcohol')
  # Priors
  phi <- round(mean(MCAR_5$alcohol, na.rm=TRUE), 2)
  tau2 <- 1
  alpha <- 0.1
  beta <- 0.1
  
  # simulate data
  n <- length(MCAR_5$alcohol[is.na(MCAR_5$alcohol)])
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
  mcar_5_cd <- MCAR_5
  mcar_5_cd$alcohol <- replace(mcar_5_cd$alcohol, is.na(mcar_5_cd$alcohol),
                               round(rnorm(length(mcar_5_cd$alcohol[is.na(mcar_5_cd$alcohol)]), 
                                           post_mu, post_sig), 2))
  
  # tests
  ts[i] <- t.test(wine$alcohol, mcar_5_cd$alcohol)
  Fs[i] <- var.test(wine$alcohol, mcar_5_cd$alcohol)
  RMSE[i] <- rmse(wine$alcohol, mcar_5_cd$alcohol)
  s <- s+1
}

length(ts[ts<=0.05]) # 45
length(Fs[Fs<=0.05]) # 0
RMSEmu <- round(mean(RMSE), 2)
RMSEmu               # 0.93
RMSEci <- round(c(quantile(RMSE, 0.025), quantile(RMSE, 0.975)), 2) 
RMSEci               # 0.75, 1.16


# univariate, MCAR, p=0.7, imputation from a conditional distribution 
nsim <- 100
ts <- rep(0, nsim)
Fs <- rep(0, nsim)
RMSE <- rep(0, nsim)
s <- 1157

for (i in 1:nsim) {
  set.seed(s)
  MCAR_7 <- delete_MCAR(ds=wine, p=0.7, cols_mis='alcohol')
  # Priors
  phi <- round(mean(MCAR_7$alcohol, na.rm=TRUE), 2)
  tau2 <- 1
  alpha <- 0.1
  beta <- 0.1
  
  # simulate data
  n <- length(MCAR_7$alcohol[is.na(MCAR_7$alcohol)])
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
  mcar_7_cd <- MCAR_7
  mcar_7_cd$alcohol <- replace(mcar_7_cd$alcohol, is.na(mcar_7_cd$alcohol),
                               round(rnorm(length(mcar_7_cd$alcohol[is.na(mcar_7_cd$alcohol)]), 
                                           post_mu, post_sig), 2))
  
  # tests
  ts[i] <- t.test(wine$alcohol, mcar_7_cd$alcohol)
  Fs[i] <- var.test(wine$alcohol, mcar_7_cd$alcohol)
  RMSE[i] <- rmse(wine$alcohol, mcar_7_cd$alcohol)
  s <- s+1
}

length(ts[ts<=0.05]) # 51
length(Fs[Fs<=0.05]) # 0
RMSEmu <- round(mean(RMSE), 2)
RMSEmu               # 1.07
RMSEci <- round(c(quantile(RMSE, 0.025), quantile(RMSE, 0.975)), 2) 
RMSEci               # 0.91, 1.29


################################ IPW ##########################################
# univariate, MCAR, p=0.1, IPW
nsim <- 100
ts <- rep(0, nsim)
Fs <- rep(0, nsim)
RMSE <- rep(0, nsim)
s <- 1157

for (i in 1:nsim) {
  set.seed(s)
  MCAR_1 <- delete_MCAR(ds=wine, p=0.1, cols_mis='alcohol')
  mcar_1_ipw <- MCAR_1
  # add two columns to the data: 1) 0/1 missingness, 2) weight
  mcar_1_ipw$R <- rep(0, 178)
  mcar_1_ipw$weight <- rep(0, 178)
  for (j in 1:178) {
    if (is.na(mcar_1_ipw$alcohol[j])) {
      mcar_1_ipw$R[j] <- 1
    }
  }
  
  # logistic regression on missingness to determine weights 
  # ipw <- lrm(R ~ ash + `color intensity` + magnesium + `total phenols` + `flavanoids`, data=mcar_1_ipw)
  ipw <- lrm(R ~ ash, data=mcar_1_ipw)
  ipw_newdat <- mcar_1_ipw
  
  # calculate weights based on logistic regression model
  weight <- predict(ipw, ipw_newdat)
  mcar_1_ipw$weight <- weight

  # calculate weighted mean
  wmean <- round(weighted.mean(mcar_1_ipw$alcohol[!is.na(mcar_1_ipw$alcohol)],
                               mcar_1_ipw$weight[!is.na(mcar_1_ipw$alcohol)]), 2)
  
  # replace missing values with weighted mean
  mcar_1_ipw$alcohol <- replace(mcar_1_ipw$alcohol, is.na(mcar_1_ipw$alcohol), wmean)
  
  ts[i] <- t.test(mcar_1_ipw$alcohol, wine$alcohol)$p.value
  Fs[i] <- var.test(mcar_1_ipw$alcohol, wine$alcohol)$p.value
  RMSE[i] <- rmse(wine$alcohol, mcar_1_ipw$alcohol)
  s <- s+1
}

length(ts[ts<=0.05]) # 0
length(Fs[Fs<=0.05]) # 0 
RMSEmu <- round(mean(RMSE), 2)
RMSEmu               # 0.26
RMSEci <- round(c(quantile(RMSE, 0.025), quantile(RMSE, 0.975)), 2) 
RMSEci               # 0.21, 0.32


# univariate, MCAR, p=0.3, IPW
nsim <- 100
ts <- rep(0, nsim)
Fs <- rep(0, nsim)
RMSE <- rep(0, nsim)
s <- 1157

for (i in 1:nsim) {
  set.seed(s)
  MCAR_3 <- delete_MCAR(ds=wine, p=0.3, cols_mis='alcohol')
  mcar_3_ipw <- MCAR_3
  # add two columns to the data: 1) 0/1 missingness, 2) weight
  mcar_3_ipw$R <- rep(0, 178)
  mcar_3_ipw$weight <- rep(0, 178)
  for (j in 1:178) {
    if (is.na(mcar_3_ipw$alcohol[j])) {
      mcar_3_ipw$R[j] <- 1
    }
  }
  
  # logistic regression on missingness to determine weights 
  # ipw <- lrm(R ~ ash + `color intensity` + magnesium + `total phenols` + `flavanoids`, data=mcar_3_ipw)
  ipw <- lrm(R ~ ash, data=mcar_3_ipw)
  ipw_newdat <- mcar_3_ipw
  
  # calculate weights based on logistic regression model
  weight <- predict(ipw, ipw_newdat)
  mcar_3_ipw$weight <- weight
  
  # calculate weighted mean
  wmean <- round(weighted.mean(mcar_3_ipw$alcohol[!is.na(mcar_3_ipw$alcohol)],
                               mcar_3_ipw$weight[!is.na(mcar_3_ipw$alcohol)]), 2)
  
  # replace missing values with weighted mean
  mcar_3_ipw$alcohol <- replace(mcar_3_ipw$alcohol, is.na(mcar_3_ipw$alcohol), wmean)
  
  ts[i] <- t.test(mcar_3_ipw$alcohol, wine$alcohol)$p.value
  Fs[i] <- var.test(mcar_3_ipw$alcohol, wine$alcohol)$p.value
  RMSE[i] <- rmse(wine$alcohol, mcar_3_ipw$alcohol)
  s <- s+1
}

length(ts[ts<=0.05]) # 0
length(Fs[Fs<=0.05]) # 92
RMSEmu <- round(mean(RMSE), 2)
RMSEmu               # 0.45
RMSEci <- round(c(quantile(RMSE, 0.025), quantile(RMSE, 0.975)), 2) 
RMSEci               # 0.39, 0.49


# univariate, MCAR, p=0.5, IPW
nsim <- 100
ts <- rep(0, nsim)
Fs <- rep(0, nsim)
RMSE <- rep(0, nsim)
s <- 1157

for (i in 1:nsim) {
  set.seed(s)
  MCAR_5 <- delete_MCAR(ds=wine, p=0.5, cols_mis='alcohol')
  mcar_5_ipw <- MCAR_5
  # add two columns to the data: 1) 0/1 missingness, 2) weight
  mcar_5_ipw$R <- rep(0, 178)
  mcar_5_ipw$weight <- rep(0, 178)
  for (j in 1:178) {
    if (is.na(mcar_5_ipw$alcohol[j])) {
      mcar_5_ipw$R[j] <- 1
    }
  }
  
  # logistic regression on missingness to determine weights 
  # ipw <- lrm(R ~ ash + `color intensity` + magnesium + `total phenols` + `flavanoids`, data=mcar_5_ipw)
  ipw <- lrm(R ~ ash, data=mcar_5_ipw)
  ipw_newdat <- mcar_5_ipw
  
  # calculate weights based on logistic regression model
  weight <- predict(ipw, ipw_newdat)
  mcar_5_ipw$weight <- weight
  
  # calculate weighted mean
  wmean <- round(weighted.mean(mcar_5_ipw$alcohol[!is.na(mcar_5_ipw$alcohol)],
                               mcar_5_ipw$weight[!is.na(mcar_5_ipw$alcohol)]), 2)
  
  # replace missing values with weighted mean
  mcar_5_ipw$alcohol <- replace(mcar_5_ipw$alcohol, is.na(mcar_5_ipw$alcohol), wmean)
  
  ts[i] <- t.test(mcar_5_ipw$alcohol, wine$alcohol)$p.value
  Fs[i] <- var.test(mcar_5_ipw$alcohol, wine$alcohol)$p.value
  RMSE[i] <- rmse(wine$alcohol, mcar_5_ipw$alcohol)
  s <- s+1
}

length(ts[ts<=0.05]) # 98
length(Fs[Fs<=0.05]) # 85
RMSEmu <- round(mean(RMSE), 2)
RMSEmu               # 7.45
RMSEci <- round(c(quantile(RMSE, 0.025), quantile(RMSE, 0.975)), 2) 
RMSEci               # 0.87, 51.01


# univariate, MCAR, p=0.7, IPW
nsim <- 100
ts <- rep(0, nsim)
Fs <- rep(0, nsim)
RMSE <- rep(0, nsim)
s <- 1157

for (i in 1:nsim) {
  set.seed(s)
  MCAR_7 <- delete_MCAR(ds=wine, p=0.7, cols_mis='alcohol')
  mcar_7_ipw <- MCAR_7
  # add two columns to the data: 1) 0/1 missingness, 2) weight
  mcar_7_ipw$R <- rep(0, 178)
  mcar_7_ipw$weight <- rep(0, 178)
  for (j in 1:178) {
    if (is.na(mcar_7_ipw$alcohol[j])) {
      mcar_7_ipw$R[j] <- 1
    }
  }
  
  # logistic regression on missingness to determine weights 
  # ipw <- lrm(R ~ ash + `color intensity` + magnesium + `total phenols` + `flavanoids`, data=mcar_7_ipw)
  ipw <- lrm(R ~ ash, data=mcar_7_ipw)
  ipw_newdat <- mcar_7_ipw
  
  # calculate weights based on logistic regression model
  weight <- predict(ipw, ipw_newdat)
  mcar_7_ipw$weight <- weight
  
  # calculate weighted mean
  wmean <- round(weighted.mean(mcar_7_ipw$alcohol[!is.na(mcar_7_ipw$alcohol)],
                               mcar_7_ipw$weight[!is.na(mcar_7_ipw$alcohol)]), 2)
  
  # replace missing values with weighted mean
  mcar_7_ipw$alcohol <- replace(mcar_7_ipw$alcohol, is.na(mcar_7_ipw$alcohol), wmean)
  
  ts[i] <- t.test(mcar_7_ipw$alcohol, wine$alcohol)$p.value
  Fs[i] <- var.test(mcar_7_ipw$alcohol, wine$alcohol)$p.value
  RMSE[i] <- rmse(wine$alcohol, mcar_7_ipw$alcohol)
  s <- s+1
}

length(ts[ts<=0.05]) # 16
length(Fs[Fs<=0.05]) # 100
RMSEmu <- round(mean(RMSE), 2)
RMSEmu               # 0.69
RMSEci <- round(c(quantile(RMSE, 0.025), quantile(RMSE, 0.975)), 2) 
RMSEci               # 0.66, 0.72