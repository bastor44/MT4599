# Chapter 2: Missing Data Mechanisms

library(missMethods)
setwd('~/MT4599')
wine <- read.table('wine.data', dec='.', sep=',') # read in data
# name columns according to documentation for dataset
colnames(wine) <- c('class', 'alcohol', 'malic acid', 'ash', 'alkalinity of ash',
                    'magnesium', 'total phenols', 'flavanoids', 'nonflavanoid phenols', 
                    'proanthocyanins', 'color intensity', 'hue', 'OD2800/OD315', 
                    'proline')
wine
miniwine <- wine[1:10, c('alcohol', 'ash', 'color intensity', 'hue')]

# 2.1 MCAR - univariate
# simulate MCAR missingness with p=0.3, cols_mis='alcohol'
set.seed(1326)
miniwine_MCAR <- delete_MCAR(ds=miniwine, p=0.3, cols_mis='alcohol')
miniwine_MCAR

# 2.1 MCAR - multivariate
# simulate MCAR missingness with p=0.3 in each column, cols_mis=c('alcohol', 'ash')
set.seed(1327)
miniwine_multiMCAR <- delete_MCAR(ds=miniwine, p=0.3, cols_mis=c('alcohol', 'ash'))
miniwine_multiMCAR


# 2.2 MAR - univariate
miniwine_MAR <- delete_MAR_censoring(ds=miniwine, p=0.3, cols_mis='alcohol', cols_ctrl='ash')
miniwine_MAR

# 2.2 MAR - multivariate
miniwine_multiMAR <- delete_MAR_censoring(ds=miniwine, p=0.3, cols_mis=c('alcohol', 'color intensity'), 
                                          cols_ctrl=c('ash', 'hue'))
miniwine_multiMAR


# 2.3 MNAR - univariate
miniwine_MNAR <- delete_MNAR_censoring(ds=miniwine, p=0.3, cols_mis='alcohol')
miniwine_MNAR

# 2.3 MNAR - multivariate
miniwine_multiMNAR <- delete_MNAR_censoring(ds=miniwine, p=0.3, cols_mis=c('alcohol', 'ash'))
miniwine_multiMNAR

