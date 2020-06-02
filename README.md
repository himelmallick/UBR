# UBR (Unified Bayesian Regularization via Scale Mixture of Uniform Distribution)

## Introduction

This repository houses the `R` package `UBR`, which provides Bayesian regularization methods for high-dimensional linear regression. The methods implemented in this package leverage a novel data augmentation scheme based on the scale mixture of uniform (SMU) distribution ([Mallick and Yi, 2014](http://intlpress.com/site/pub/files/_fulltext/journals/sii/2014/0007/0004/SII-2014-0007-0004-a012.pdf); [Mallick, 2015](https://www.researchgate.net/publication/290002016_Some_Contributions_to_Bayesian_Regularization_Methods_with_Applications_to_Genetics_and_Clinical_Trials); [Mallick and Yi, 2017](https://www.ncbi.nlm.nih.gov/pubmed/28943688); [Mallick and Yi, 2018](https://www.researchgate.net/publication/315695016_Bayesian_Bridge_Regression)) that leads to a set of efficient Gibbs samplers with tractable full conditional posterior distributions and superior performance over existing methods. 

## Dependencies

`UBR` requires the following `R` package: `devtools` (for installation only). Please install it before installing `UBR`, which can be done as follows (execute from within a fresh R session):
```r
install.packages("devtools")
library(devtools)
```

## Installation

Once the dependencies are installed, `UBR` can be loaded using the following command:
```r
devtools::install_github("himelmallick/UBR")
library(UBR)
```
 
## Basic Usage
```r
BayesBridgeSMU(x, y, alpha = 0.5, max.steps = 20000, n.burn = 10000, n.thin = 1, tau = 0, r = 1, delta = 0.1, posterior.summary.beta = "mean", posterior.summary.lambda = "mean", posterior.summary.alpha = "mean", beta.ci.level = 0.95, lambda.ci.level = 0.95, alpha.ci.level = 0.95, seed = 1234)
```

```r
BayesGroupBridgeSMU(x, y, alpha = 0.5, group = group, max.steps = 20000, n.burn = 10000, n.thin = 1, r = 1, delta = 0.1, posterior.summary.beta = "mean", posterior.summary.lambda = "mean", posterior.summary.alpha = "mean", beta.ci.level = 0.95, lambda.ci.level = 0.95, alpha.ci.level = 0.95, seed = 1234, alpha.prior = "none", a = 10, b = 10, update.sigma2 = T)
```

## Input

- **x**:  A data frame or matrix with standardized predictors in columns and samples in rows.
- **y**: A mean-centered continuous response variable with matched sample rows with `x`.
- **alpha**: Concavity parameter of the bridge or group bridge penalty (the exponent to which the L1 norm of the coefficients are raised). Default is 0.5, the square root. Must be between 0 and 1 (0 exclusive). It can also be estimated from the data by specifying an arbitrary non-positive value such as 0 (currently only implemented for group bridge).
- **max.steps**: Number of MCMC iterations. Default is `20000`.
- **n.burn**: Number of burn-in iterations. Default is `10000`.
- **n.thin**: Lag at which thinning should be done. Default is `1` (no thinning).
- **tau**: Tolerance parameter to be used when inverting the covariance matrix. Default is `0` (bridge).
- **r**: Shape hyperparameter for the Gamma prior on lambda. Default is `1`.
- **delta**: Rate hyperparameter for the Gamma prior on lambda. Default is `0.1`.
- **posterior.summary.beta**: Posterior summary measure for beta (mean or median). Default is `'mean'`.
- **posterior.summary.lambda**:	Posterior summary measure for lambda (mean or median). Default is `'mean'`.
- **posterior.summary.alpha**: Posterior summary measure for alpha (mean or median). Default is `'mean'` (group bridge).
- **beta.ci.level**: Credible interval level for beta. Default is `0.95` (95\%).
- **lambda.ci.level**: Credible interval level for lambda. Default is `0.95` (95\%).
- **alpha.ci.level**: Credible interval level for alpha when alpha is not pre-fixed. Default is `0.95` (95\%) (group bridge).
- **seed**: Seed value for reproducibility. Default is `1234`.
- **alpha.prior**: Prior distribution on alpha when alpha is unknown. Must be one of 'uniform' and 'beta' (group bridge).
- **a**: Beta prior distribution shape parameter 1 when `alpha.prior` is 'beta'. Default is 10 (group bridge).
- **B**: Beta prior distribution shape parameter 2 when `alpha.prior` is 'beta'. Default is 10 (group bridge).
- **update.sigma2**: Whether sigma2 should be updated. Default is TRUE (group bridge).



## Output

A list containing the following components is returned:

- **time**: Computational time in minutes.
- **beta**: Posterior mean or median estimates of beta.
- **lowerbeta**: Lower limit credible interval of beta.
- **upperbeta**: Upper limit credible interval of beta.
- **alpha**: Posterior estimate of alpha if alpha is unknown (group bridge).
- **alphaci**: Posterior credible interval of alpha if alpha is unknown (group bridge).
- **lambda**: Posterior estimate of lambda.
- **lowerlambda**: Lower limit credible interval of lambda (group bridge).
- **upperlambda**: Upper limit credible interval of lambda (group bridge).
- **lambdaci**: Posterior credible interval of lambda.
- **beta.post**: Post-burn-in posterior samples of beta.
- **sigma2.post**: Post-burn-in posterior samples of sigma2.
- **lambda.post**: Post-burn-in posterior samples of lambda
- **alpha.post**: Post-burn-in posterior samples of alpha if alpha is unknown (group bridge).


## Example 1 (Bayesian Bridge Regression)

Let's consider the `prostate` dataset from the `ElemStatLearn` package.

```r

#########################
# Load Prostate dataset #
#########################

library(ElemStatLearn)
prost<-prostate

```

First, we standardize the variables and center the response variable. We also separate out training and test samples.

```r

###########################################
# Scale data and prepare train/test split #
###########################################

prost.std <- data.frame(cbind(scale(prost[,1:8]),prost$lpsa))
names(prost.std)[9] <- 'lpsa'
data.train <- prost.std[prost$train,]
data.test <- prost.std[!prost$train,]

##################################
# Extract standardized variables #
##################################

y.train   = data.train$lpsa - mean(data.train$lpsa)
y.test <- data.test$lpsa - mean(data.test$lpsa)
x.train = scale(as.matrix(data.train[,1:8], ncol=8))
x.test = scale(as.matrix(data.test[,1:8], ncol=8))

```

Let's apply Bayesian Bridge with a concavity parameter of `0.5` and calculate predictions in the test set.

```r

##############################################
# Bayesian Bridge with alpha = 0.5 using SMU #
##############################################

library(UBR)
BayesBridge<- BayesBridgeSMU(as.data.frame(x.train), y.train, alpha = 0.5)
y.pred.BayesBridge <-x.test%*%BayesBridge$beta

```

Note that, Bayesian Bridge with a concavity parameter of `1` is equivalent to a Bayesian LASSO regression. Let's apply Bayesian LASSO using SMU and calculate predictions in the test set.

```r

#############################################################
# Bayesian Bridge with alpha = 1 using SMU (Bayesian LASSO) #
#############################################################

BayesLASSO<- BayesBridgeSMU(as.data.frame(x.train), y.train, alpha = 1)
y.pred.BayesLASSO <-x.test%*%BayesLASSO$beta

```

Let's compare with various frequentist methods such as OLS, LASSO, adaptive LASSO, classical bridge, and elastic net.

```r 

####################
# Compare with OLS #
####################

ols.lm<-lm(y.train ~ x.train)
y.pred.ols<-x.test%*%ols.lm$coefficients[-1]

##################################
# Compare with Frequentist LASSO #
##################################

set.seed(1234)
library(glmnet)
lasso.cv=cv.glmnet(x.train,y.train,alpha = 1)  
lambda.cv.lasso=lasso.cv$lambda.min          
lasso.sol=glmnet(x.train,y.train,alpha = 1)    
lasso.coeff=predict(lasso.sol,type="coef",s=lambda.cv.lasso)
y.pred.lasso=x.test%*%lasso.coeff[-1]

###########################################
# Compare with Frequentist Adaptive LASSO #
###########################################

alasso.cv=cv.glmnet(x.train,y.train,alpha = 1,penalty.factor=1/abs(ols.lm$coefficients))  
lambda.cv.alasso=alasso.cv$lambda.min          
alasso.sol=glmnet(x.train,y.train,alpha = 1,penalty.factor=1/abs(ols.lm$coefficients))    
alasso.coeff=predict(alasso.sol,type="coef",s=lambda.cv.alasso)
y.pred.alasso=x.test%*%alasso.coeff[-1]

#################################
# Compare with Classical Bridge #
#################################

library(grpreg)
group<-c(rep(1,8))
fit.bridge <- gBridge(x.train, y.train,group=group)
coef.bridge<-grpreg::select(fit.bridge,'AIC')$beta
y.pred.bridge<-x.test%*%coef.bridge[-1]

########################################
# Compare with Frequentist Elastic Net #
########################################

library(glmnet)
enet.cv=cv.glmnet(x.train,y.train,alpha = 0.5)  
lambda.cv.enet=enet.cv$lambda.min          
enet.sol=glmnet(x.train,y.train,alpha = 0.5)    
enet.coeff=predict(enet.sol,type="coef",s=lambda.cv.enet)
y.pred.enet=x.test%*%enet.coeff[-1]

```

Let's compare the performance in test data.

```r 

####################
# MSE on Test Data #
####################

mean((y.pred.alasso-y.test)^2) # Frequentist Adaptive LASSO: 0.7419881
mean((y.pred.ols-y.test)^2) # OLS: 0.5421042
mean((y.pred.lasso-y.test)^2) # Frequentist LASSO: 0.5219545
mean((y.pred.bridge-y.test)^2) # Classical Bridge: 0.5168394
mean((y.pred.enet-y.test)^2) # Frequentist Elastic Net: 0.5135257
mean((y.pred.BayesLASSO - y.test)^2) # Bayesian LASSO with SMU: 0.4835839
mean((y.pred.BayesBridge - y.test)^2) # Bayesian Bridge with SMU: 0.4725654

```

Let's check MCMC convergence of the Bayesian Bridge estimator based on SMU through two visualizations: trace plots and histograms.


```r 

######################################
# Visualization of Posterior Samples #
######################################

##############
# Trace Plot #
##############

library(coda)
plot(mcmc(BayesBridge$beta.post),density=FALSE,smooth=TRUE)

```

![plot of chunk traceplot](https://github.com/himelmallick/UBR/blob/master/misc/traceplot_prostate.png)

```r 

#############
# Histogram #
#############

library(psych)
multi.hist(BayesBridge$beta.post,density=TRUE,main="")

```

![plot of chunk histogram](https://github.com/himelmallick/UBR/blob/master/misc/histogram_prostate.png)

## Example 2 (Bayesian Group Bridge)

```r
################################
# Load the Birthweight Dataset #
################################

library(grpreg)
data(birthwt.grpreg) # Hosmer and lemeshow (1989)
x <- scale(as.matrix(birthwt.grpreg[,-c(1,2)]))
y <- birthwt.grpreg$bwt - mean(birthwt.grpreg$bwt)
group <- c(1,1,1,2,2,2,3,3,4,5,5,6,7,8,8,8)

###########################
# Classical Group Bridge  #
###########################

fit.gbridge <- gBridge(x, y, group=group)
coef.gbridge<-select(fit.gbridge,'BIC')$beta
coef.gbridge[-1][coef.gbridge[-1]!=0]

#       age2        age3        lwt1        lwt3       white       black       smoke        ptl1 
# 0.06509616  0.02053118  0.07431276  0.05355660  0.13286843 -0.01387433 -0.10794526 -0.05735128 
# ht          ui 
# -0.07009045 -0.13688881 

###############
# Group LASSO #
###############
fit.glasso <- grpreg(x, y, group=group,penalty='grLasso')
coef.glasso<-select(fit.glasso,'BIC')$beta
coef.glasso[-1][coef.glasso[-1]!=0]

# age1         age2         age3         lwt1         lwt2         lwt3        white 
# 0.009441572  0.037498245  0.022409360  0.045851715 -0.011028536  0.036018025  0.084093467 
# black        smoke         ptl1        ptl2m           ht           ui 
# -0.018460475 -0.085362146 -0.052697766  0.007808026 -0.065206615 -0.131695837 


##########################
# Bayesian Group Bridge  #
##########################

library(UBR)
fit.BayesGroupBridge<-BayesGroupBridgeSMU(x, y, group = group) 
fit.BayesGroupBridge$beta[which(sign(fit.BayesGroupBridge$lowerbeta) == sign(fit.BayesGroupBridge$upperbeta))]

# white      smoke         ht         ui 
# 0.1388594 -0.1314631 -0.1122126 -0.1656339 

```



## Citation

If you use `UBR` in your work, please cite the following papers:

Mallick, H., Yi, N. (2014). [A New Bayesian LASSO](http://intlpress.com/site/pub/files/_fulltext/journals/sii/2014/0007/0004/SII-2014-0007-0004-a012.pdf). Statistics and Its Interface, 7(4): 571-582.

Mallick, H. (2015). [Some Contributions to Bayesian Regularization Methods with Applications to Genetics and Clinical Trials](https://www.researchgate.net/publication/290002016_Some_Contributions_to_Bayesian_Regularization_Methods_with_Applications_to_Genetics_and_Clinical_Trials). Doctoral Dissertation, University of Alabama at Birmingham.

Mallick, H., Yi, N. (2017). [Bayesian Group Bridge for Bi-level Variable Selection](https://www.ncbi.nlm.nih.gov/pubmed/28943688). Computational Statistics and Data Analysis, 110(6): 115-133.

Mallick, H., Yi, N. (2018). [Bayesian Bridge Regression](https://www.researchgate.net/publication/315695016_Bayesian_Bridge_Regression). Journal of Applied Statistics, 45(6): 988-1008.
