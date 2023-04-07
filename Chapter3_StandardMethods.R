# Chapter 3: Standard Methods (Illustrations)

library(missMethods)
library(rms)
setwd('~/MT4599')
wine <- read.table('wine.data', dec='.', sep=',') # read in data
# name columns according to documentation for dataset
colnames(wine) <- c('class', 'alcohol', 'malic acid', 'ash', 'alkalinity of ash',
                    'magnesium', 'total phenols', 'flavanoids', 'nonflavanoid phenols', 
                    'proanthocyanins', 'color intensity', 'hue', 'OD2800/OD315', 
                    'proline')
wine
miniwine <- wine[1:10, c('ash', 'alcohol', 'hue', 'color intensity')]

#### 3.2.1 Mean Imputation
# MCAR - univariate
# simulate MCAR missingness with p=0.3, cols_mis='alcohol'
set.seed(1326)
miniwine_MCAR <- delete_MCAR(ds=miniwine, p=0.3, cols_mis='alcohol')
miniwine_MCAR

miniwine_mi <- miniwine_MCAR
mu_alc <- mean(miniwine_mi$alcohol, na.rm=TRUE)
miniwine_mi$alcohol <- replace(miniwine_mi$alcohol, is.na(miniwine_mi$alcohol), mu_alc)


# MCAR - multivariate 
set.seed(1327)
miniwine_multiMCAR <- delete_MCAR(ds=miniwine, p=0.3, cols_mis=c('alcohol', 'ash'))
miniwine_multiMCAR
mu_alc_mult <- mean(miniwine_multiMCAR$alcohol, na.rm=TRUE)
mu_ash_mult <- mean(miniwine_multiMCAR$ash, na.rm=TRUE)
mu_alc_mult
mu_ash_mult


#### 3.2.2 Regression-based Single Imputation 
minireg <- lm(alcohol ~ ash, miniwine_MCAR)
newdat <- miniwine_MCAR[is.na(miniwine_MCAR$alcohol), ]
pred <- predict(minireg, newdat)
pred <- round(pred,2)
miniwine_reg <- miniwine_MCAR
miniwine_reg$alcohol <- replace(miniwine_reg$alcohol, is.na(miniwine_reg$alcohol), 
                                pred)
miniwine_reg
minireg



#### 3.2.3 Hot Deck Imputation
miniwine_hd <- miniwine_MCAR
n <- length(miniwine_hd$alcohol[is.na(miniwine_hd$alcohol)])
set.seed(1320)
hd <- sample(miniwine_hd$alcohol[!is.na(miniwine_hd$alcohol)], n, replace=TRUE)
hd


#### 3.2.4 Imputation from a Conditional Distribution
set.seed(1157)
mcar <- delete_MCAR(ds=wine, p=0.3, cols_mis='alcohol')
hist(mcar$alcohol, breaks=30, freq=FALSE, xlim=c(11, 15), xlab='Alcohol content', 
     main='Histogram of alcohol content under MCAR with p=0.3', cex.axis=1.5, cex.lab=1.5)
x <- seq(11, 15, length=1000)
y <- dnorm(x, 13, .75)
lines(x,y, lwd=2, col='blue')
legend('topleft', legend='N(13, 0.75) distribution', col='blue', lwd=2, cex=1.3)

miniwine_cd <- miniwine_MCAR
# Priors
phi <- round(mean(miniwine_cd$alcohol, na.rm=TRUE), 2)
tau2 <- 1
alpha <- 0.1
beta <- 0.1

# simulate data
n <- length(miniwine_cd$alcohol[is.na(miniwine_cd$alcohol)])
x <- rnorm(n, phi, tau2)
xbar <- mean(x, na.rm=TRUE)

# samples to take in the chain 
T <- 1000
mu <- rep(0, 2)
sigma2 <- rep(0, 2)

# set starting values
mu[1] <- 13
sigma2[1] <- 0.75

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
miniwine_cd$alcohol <- replace(miniwine_cd$alcohol, is.na(miniwine_cd$alcohol),
                             round(rnorm(length(miniwine_cd$alcohol[is.na(miniwine_cd$alcohol)]), 
                                         post_mu, post_sig), 2))
miniwine_cd
miniwine_MCAR$alcohol
miniwine_cd$alcohol


#### 3.2.5 Inverse Probability Weighting
miniwine
miniwine_MCAR
miniwine_ipw <- miniwine_MCAR
miniwine_ipw$R <- rep(0, 10)
miniwine_ipw$weight <- rep(0, 10)
for (i in 1:10) {
  if (is.na(miniwine_ipw$alcohol[i])) {
    miniwine_ipw$R[i] <- 1
  }
}
ipw <- lrm(R ~ ash, miniwine_ipw)
newdat <- miniwine_ipw

w <- predict(ipw, newdat)
miniwine_ipw$weight <- w

wmean <- round(weighted.mean(miniwine_ipw$alcohol[!is.na(miniwine_ipw$alcohol)], 
                             miniwine_ipw$weight[!is.na(miniwine_ipw$alcohol)]),2)
wmean
