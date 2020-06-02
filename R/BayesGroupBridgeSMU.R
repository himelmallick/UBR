#' Unified Bayesian Regularization via Scale Mixture of Uniform Distribution
#'
#' Bayesian Group Bridge regression
#' 
#' @param x A numeric matrix with standardized predictors in columns and samples in rows.
#' @param y A mean-centered continuous response variable with matched sample rows with x.
#' @param alpha Concavity parameter of the group bridge penalty (the exponent to which the L1 norm of the coefficients are raised). 
#' Default is 0.5, the square root. Must be between 0 and 1 (0 exclusive). 
#' It can also be estimated from the data by specifying arbitrary non-postive value such as 0. 
#' @param group Index of grouping information among the predictors.
#' @param max.steps Number of MCMC iterations. Default is 20000.
#' @param n.burn Number of burn-in iterations. Default is 10000.
#' @param n.thin Lag at which thinning should be done. Default is 1 (no thinning).
#' @param r Shape hyperparameter for the Gamma prior on lambda. Default is 1.
#' @param delta Rate hyperparameter for the Gamma prior on lambda. Default is 0.1.
#' @param posterior.summary.beta Posterior summary measure for beta (mean, median, or mode). Default is 'mean'. 
#' @param posterior.summary.lambda Posterior summary measure for lambda (mean or median). Default is 'mean'. 
#' @param posterior.summary.alpha Posterior summary measure for alpha (mean or median) when alpha is unknown. Default is 'mean'. 
#' @param beta.ci.level Credible interval level for beta. Default is 0.95 (95\%).
#' @param lambda.ci.level Credible interval level for lambda. Default is 0.95 (95\%).
#' @param alpha.ci.level Credible interval level for alpha when alpha is unknown. 
#' Default is 0.95 (95\%).
#' @param seed Seed value for reproducibility. Default is 1234.
#' @param alpha.prior Prior distribution on alpha when alpha is unknown. 
#' Must be one of 'uniform' and 'beta'.
#' @param a Beta prior distribution shape parameter 1 when `alpha.prior` is 'beta'. Default is 10.
#' @param b Beta prior distribution shape parameter 2 when `alpha.prior` is 'beta'. Default is 10.
#' @param update.sigma2 Whether sigma2 should be updated. Default is TRUE.
#' 
#' @return A list containing the following components is returned:
#' \item{time}{Computational time in minutes.}
#' \item{beta}{Posterior mean or median estimates of beta.}
#' \item{lowerbeta}{Lower limit credible interval of beta.}
#' \item{upperbeta}{Upper limit credible interval of beta.}
#' \item{lambda}{Posterior mean or median estimates of lambda.}
#' \item{lowerlambda}{Lower limit credible interval of lambda}
#' \item{upperlambda}{Upper limit credible interval of lambda}
#' \item{alpha}{Posterior estimate of alpha if alpha is unknown.}
#' \item{alphaci}{Posterior credible interval of alpha if alpha is unknown.}
#' \item{beta.post}{Post-burn-in posterior samples of beta.}
#' \item{sigma2.post}{Post-burn-in posterior samples of sigma2.}
#' \item{lambda.post}{Post-burn-in posterior samples of lambda}
#' \item{alpha.post}{Post-burn-in posterior samples of alpha if alpha is unknown.}
#'
#' @examples
#' \dontrun{
#' 
#' ################################
#' # Load the Birthweight Dataset #
#' ################################
#' 
#' library(grpreg)
#' data(birthwt.grpreg) # Hosmer and lemeshow (1989)
#' x <- scale(as.matrix(birthwt.grpreg[,-c(1,2)]))
#' y <- birthwt.grpreg$bwt - mean(birthwt.grpreg$bwt)
#' group <- c(1,1,1,2,2,2,3,3,4,5,5,6,7,8,8,8)
#' 
#' ###########################
#' # Classical Group Bridge  #
#' ###########################
#' 
#' fit.gbridge <- gBridge(x, y, group=group)
#' coef.gbridge<-select(fit.gbridge,'BIC')$beta
#' coef.gbridge[-1][coef.gbridge[-1]!=0]
#' 
#' ###############
#' # Group LASSO #
#' ###############
#' 
#' fit.glasso <- grpreg(x, y, group=group,penalty='grLasso')
#' coef.glasso<-select(fit.glasso,'BIC')$beta
#' coef.glasso[-1][coef.glasso[-1]!=0]
#' 
#' ##########################
#' # Bayesian Group Bridge  #
#' ##########################
#' 
#' library(UBR)
#' fit.BayesGroupBridge<-BayesGroupBridgeSMU(x, y, group = group) 
#' fit.BayesGroupBridge$beta[which(sign(fit.BayesGroupBridge$lowerbeta) == sign(fit.BayesGroupBridge$upperbeta))]
#' 
#' }
#'
#'
#' @export
#' @keywords SMU, MCMC, Bayesian regularization, Bridge, Group Bridge
BayesGroupBridgeSMU = function(x, 
                y,
                alpha = 0.5, 
                group = group, 
                max.steps = 20000,
                n.burn = 10000,
                n.thin = 1, 
                r = 1,
                delta = 0.1,
                posterior.summary.beta = 'mean',
                posterior.summary.lambda = 'mean',
                posterior.summary.alpha = 'mean',
                beta.ci.level = 0.95,
                lambda.ci.level = 0.95,
                alpha.ci.level = 0.95,
                seed = 1234,
                alpha.prior = 'none', 
                a = 10, 
                b = 10, 
                update.sigma2 = T){
  
  # Set random seed
  set.seed(seed)
  
  # Create grouping structure from the provided group information
  group<-create.group(group)
  G<-length(group$gsize)
  index<-create.index(group)
  
  # Set parameter values and extract relevant quantities
  known.alpha = alpha > 0
  x <- as.matrix(x)
  n <- nrow(x)
  m <- ncol(x)

  # Calculate OLS quantities
  XtX <- t(x) %*% x	
  xy <- t(x) %*% y
  ixx <- chol2inv(chol(XtX))
  xy <- t(x) %*% y
  bhat <- ixx%*%xy;

  # Initialize the Coefficient Vectors and Matrices
  betaSamples <- matrix(0, max.steps, m)
  sigma2Samples <- rep(0, max.steps)
  uSamples <- matrix(0, max.steps, G)
  lambdaSamples <- matrix(0, max.steps, G)
  alphaSamples<-rep(0,max.steps)
  
  # Initial Values 
  sigma2<-(1/n)*sum((y-x%*%bhat)^2)
  u<-lambda<-rep(1,G)
  beta<-rep(0,m);
  if (!known.alpha) alpha<-runif(1,0,1);
  residue <-drop(y-x%*%beta)
  
  # Counter for the MCMC sampling iterations
  k <- 0
  
  # Track time
  cat(c("Job started at:",date()),fill=TRUE)
  start.time <- Sys.time()
  
  # Gibbs Sampler
  while (k < max.steps) {
    k <- k + 1
    if (k %% 1000 == 0) {
      cat('Iteration:', k, "\r")
    }

    # Sample beta
    bound<-u^(1/alpha)
    Mu<-bhat
    Sigma<-sigma2*ixx
    betag<-c()
    for (g in 1:G){
      case<-unlist(group$gindex[g]);
      given<-seq(1,m,1)[-case];
      if(length(case)==1){
        obj<-condNormal(beta[given],mu=Mu, sigma=Sigma, req=case, given=given)
        meang<-as.vector(obj$condMean)
        covg<-as.vector(obj$condVar)
        betag<-tmvmixnorm::rtuvn(1,mean=meang, sd=sqrt(covg),lower=-bound[g],upper=bound[g])
        beta[case]<-betag
      }
      else{
        obj<-condNormal(beta[given],mu=Mu, sigma=Sigma, req=case, given=given)
        meang<-obj$condMean
        covg<-obj$condVar
        betag<-as.vector(rtmvn.simp.1(1,meang,covg,bound[g],beta[case]))
        beta[case]<-betag
      }
    }
    betaSamples[k,] <- beta;
    
    
    # Sample sigma2
    if (update.sigma2==T)
    {
      shape.sig <- n/2
      residue <- drop(y - x%*%beta)
      rate.sig <- (t(residue) %*% (residue))/2
      sigma2 <-1/(rgamma(1,shape.sig,rate.sig))
    }
    sigma2Samples[k] <- sigma2
    
    # Sample u
    lambdaPrime <- lambda
    for (i in 1:G) {
      norm<-sum(abs(beta[group$gindex[[i]]]))
      T<- norm^alpha
      u[i] <- as.vector(rexp(1,lambdaPrime[i])) + T
    }
    uSamples[k, ] <- u
    
    # Update lambda
    for (i in 1:G) {
      shape.lamb =  r + (group$gsize[i]/alpha)
      rate.lamb = delta + u[i]
      lambda[i] <- rgamma(1, shape=shape.lamb, rate=rate.lamb)
    }
    lambdaSamples[k,] <- lambda
    
    # Update alpha	
    if (!known.alpha) {
      if (alpha.prior=='beta') alpha<-draw.alpha(alpha,beta,lambda, group=group,pr.a=a, pr.b=b,ep=0.1)
      if (alpha.prior=='uniform') alpha<-draw.alpha.2(alpha,beta,lambda, group=group,ep=0.1)
    }
    alphaSamples[k] <- alpha
  }
  
  # Collect all quantities and prepare outputs
  beta.post=betaSamples[seq(n.burn+1, max.steps, n.thin),]
  sigma2.post=sigma2Samples[seq(n.burn+1, max.steps, n.thin)]
  lambda.post=lambdaSamples[seq(n.burn+1, max.steps, n.thin),]
  alpha.post<-alphaSamples[seq(n.burn+1, max.steps, n.thin)]
  
  # Posterior estimates of beta
  if (posterior.summary.beta == 'mean') {
    beta <-apply(beta.post, 2, mean)
  } else if (posterior.summary.beta == 'median') {
    beta <-apply(beta.post, 2, median)
  } else if (posterior.summary.beta == 'mode') {
    beta <-posterior.mode(mcmc(beta.post))
  } else{
    stop('posterior.summary.beta must be either mean, median, or mode.')
  }
  
  # Posterior estimates of lambda
  if (posterior.summary.lambda == 'mean') {
    lambda <-apply(lambda.post, 2, mean)
  } else if (posterior.summary.lambda == 'median') {
    lambda <-apply(lambda.post, 2, median)
  } else{
    stop('posterior.summary.lambda must be either mean or median.')
  }
  

  # Credible intervals for beta
  lowerbeta<-colQuantiles(beta.post,prob=(1 - beta.ci.level)/2)
  upperbeta<-colQuantiles(beta.post,prob=(1 + beta.ci.level)/2)
  
  # Assign names
  names(beta)<-names(lowerbeta)<-names(upperbeta)<-colnames(beta.post)<-colnames(x)
  
  # Credible intervals for lambda
  lowerlambda<-colQuantiles(lambda.post,prob=(1 - lambda.ci.level)/2)
  upperlambda<-colQuantiles(lambda.post,prob=(1 + lambda.ci.level)/2)
  
  # Assign names
  names(lambda)<-names(lowerlambda)<-names(upperlambda)<-colnames(lambda.post)<-paste('Group', 1:G, sep='')
  
  # Credible intervals for alpha
  if (!known.alpha){
    alphaci<-credInt(alpha.post, cdf = NULL,conf = alpha.ci.level, type="twosided")      
  } else{
    alphaci<-NULL
    alpha.post<-NULL
  }
  

  # Track stopping time
  stop.time <- Sys.time()
  time<-round(difftime(stop.time, start.time, units="min"),3)
  cat(c("Job finished at:",date()),fill=TRUE)
  
  
  # Return output
  return(list(time = format(time),
              beta = beta,
              lowerbeta = lowerbeta,
              upperbeta = upperbeta,
              lambda = lambda,
              lowerlambda = lowerlambda,
              upperlambda = upperlambda,
              alpha = alpha,
              alphaci = alphaci,
              beta.post = beta.post,
              sigma2.post = sigma2.post,
              lambda.post = lambda.post,
              alpha.post = alpha.post))
}



