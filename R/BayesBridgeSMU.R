#' Unified Bayesian Regularization via Scale Mixture of Uniform Distribution
#' 
#' Bayesian Bridge Regression
#' 
#' @param x A numeric matrix with standardized predictors in columns and samples in rows.
#' @param y A mean-centered continuous response variable with matched sample rows with x.
#' @param alpha Concavity parameter of the bridge penalty (the exponent to which the L1 norm of the coefficients are raised). 
#' Default is 0.5, the square root. Must be in (0,1]. 
#' @param max.steps Number of MCMC iterations. Default is 20000.
#' @param n.burn Number of burn-in iterations. Default is 10000.
#' @param n.thin Lag at which thinning should be done. Default is 1 (no thinning).
#' @param eps Small tolerance parameter to be used when inverting X,
#' if it is less than full rank. Default is 1e-06.
#' @param r Shape hyperparameter for the Gamma prior on lambda. Default is 1.
#' @param delta Rate hyperparameter for the Gamma prior on lambda. Default is 0.1.
#' @param posterior.summary.beta Posterior summary measure for beta (mean, median, or mode). Default is 'mean'. 
#' @param posterior.summary.lambda Posterior summary measure for lambda (mean or median). Default is 'mean'. 
#' @param beta.ci.level Credible interval level for beta. Default is 0.95 (95\%).
#' @param lambda.ci.level Credible interval level for lambda. Default is 0.95 (95\%).
#' @param seed Seed value for reproducibility. Default is 1234.
#' 
#' @return A list containing the following components is returned:
#' \item{time}{Computational time in minutes.}
#' \item{beta}{Posterior mean or median estimates of beta.}
#' \item{lowerbeta}{Lower limit credible interval of beta.}
#' \item{upperbeta}{Upper limit credible interval of beta.}
#' \item{lambda}{Posterior estimate of lambda.}
#' \item{lambdaci}{Posterior credible interval of lambda.}
#' \item{beta.post}{Post-burn-in posterior samples of beta.}
#' \item{sigma2.post}{Post-burn-in posterior samples of sigma2.}
#' \item{lambda.post}{Post-burn-in posterior samples of lambda}
#' 
#' @examples
#' \dontrun{
#' 
#' #########################
#' # Load Diabetes dataset #
#' #########################
#' 
#' library(ElemStatLearn)
#' prost<-prostate
#' 
#' ###########################################
#' # Scale data and prepare train/test split #
#' ###########################################
#' 
#' prost.std <- data.frame(cbind(scale(prost[,1:8]),prost$lpsa))
#' names(prost.std)[9] <- 'lpsa'
#' data.train <- prost.std[prost$train,]
#' data.test <- prost.std[!prost$train,]
#'
#' ##################################
#' # Extract standardized variables #
#' ##################################
#' 
#' y.train   = data.train$lpsa - mean(data.train$lpsa)
#' y.test <- data.test$lpsa - mean(data.test$lpsa)
#' x.train = scale(as.matrix(data.train[,1:8], ncol=8))
#' x.test = scale(as.matrix(data.test[,1:8], ncol=8))
#' 
#' ##############################################
#' # Bayesian Bridge with alpha = 0.5 using SMU #
#' ##############################################
#' 
#' library(UBR)
#' BayesBridge<- BayesBridgeSMU(as.data.frame(x.train), y.train, alpha = 0.5)
#' y.pred.BayesBridge<-x.test%*%BayesBridge$beta
#' mean((y.pred.BayesBridge - y.test)^2) # Performance on test data
#'
#' #############################################################
#' # Bayesian Bridge with alpha = 1 using SMU (Bayesian LASSO) #
#' #############################################################
#' 
#' BayesLASSO<- BayesBridgeSMU(as.data.frame(x.train), y.train, alpha = 1)
#' y.pred.BayesLASSO<-x.test%*%BayesLASSO$beta
#' mean((y.pred.BayesLASSO - y.test)^2) # Performance on test data
#'
#' ########################################################
#' # Visualization of Posterior Samples (Bayesian Bridge) #
#' ########################################################
#' 
#' ##############
#' # Trace Plot #
#' ##############
#' 
#' library(coda)
#' plot(mcmc(BayesBridge$beta.post),density=FALSE,smooth=TRUE)
#' 
#' #############
#' # Histogram #
#' #############
#' 
#' library(psych)
#' multi.hist(BayesBridge$beta.post,density=TRUE,main="")
#' 
#' }
#'
#' @export
#' @keywords SMU, MCMC, Bayesian regularization

BayesBridgeSMU = function(x,
                          y, 
                          alpha = 0.5, 
                          max.steps=20000,
                          n.burn=10000,
                          n.thin=1, 
                          eps=1e-06,
                          r = 1,
                          delta = 0.1,
                          posterior.summary.beta = 'mean',
                          posterior.summary.lambda = 'mean',
                          posterior.summary.alpha = 'mean',
                          beta.ci.level = 0.95,
                          lambda.ci.level = 0.95,
                          alpha.ci.level = 0.95,
                          seed = 1234) {
  
  # Set random seed
  set.seed(seed)
  
  # Set parameter values and extract relevant quantities
	x <- as.matrix(x)
 	n <- nrow(x)
	m <- ncol(x)

	# Calculate OLS quantities
	XtX <- t(x) %*% x	
	xy <- t(x) %*% y

	# Initialize the Coefficient Vectors and Matrices
	betaSamples <- matrix(0, max.steps, m)
	sigma2Samples <- rep(0, max.steps)
	uSamples <- matrix(0, max.steps, m)
	lambdaSamples <- rep(0, max.steps)

	# Initial Values 
 	lambda <- rgamma(1,shape=r, rate=delta)
 	sigma2<- runif(1,0.1,10)
	beta   <- rep(0,m)
	u<-rep(0,m)
	for (i in 1:m){
		u[i]     <- rgamma(1,shape=(1/alpha+1), rate=lambda) 
			}
    for (i in 1:m){
		beta[i]     <- runif(1,-sqrt(sigma2)*u[i]^(1/alpha), sqrt(sigma2)*u[i]^(1/alpha))
			}
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
	  if (qr(x)$rank == m) eps<-0
		invA <- solve(XtX+diag(eps,nrow(XtX),ncol(XtX)))
		mean <- drop(invA%*%xy)
		varcov <- sigma2*invA
		z<-sqrt(sigma2)*u^(1/alpha)
		identity<-diag(1,m)
		beta<- tmvmixnorm::rtmvn(n = 1,
		                   Mean = mean,
		                   Sigma = varcov,
		                   D = identity,
		                   lower = -z,
		                   upper = z, 
		                   int = beta) 
		
		betaSamples[k,] <- beta

		# Sample sigma2
		shape.sig <- (n+m-1)/2
		residue <- drop(y -x%*%beta)
		rate.sig <- (t(residue) %*% residue)/2
		lower<-max(beta^2/u^(2/alpha))
		upper<-Inf
		sigma2 <-rtinvgamma(1,shape.sig,rate.sig,lower,upper)
		sigma2Samples[k] <- sigma2

		# Sample u
		T<- abs(beta/sqrt(sigma2))^alpha
		lambdaPrime <- lambda
		u<- rep(0, m)
		for (i in seq(m)) {
			u[i] <- rexp(1,lambdaPrime) + T[i]
		}
		uSamples[k, ] <- u

		# Update lambda
		shape.lamb =  r + m + (m/alpha)
		rate.lamb = delta + sum(abs(u))
		lambda <- rgamma(1, shape=shape.lamb, rate=rate.lamb)
		lambdaSamples[k] <- lambda
		}
	
	# Collect all quantities and prepare outputs
	beta.post=betaSamples[seq(n.burn+1, max.steps, n.thin),]
	sigma2.post=sigma2Samples[seq(n.burn+1, max.steps, n.thin)]
	lambda.post=lambdaSamples[seq(n.burn+1, max.steps, n.thin)]

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
	  lambda<-mean(lambda.post)
	} else if (posterior.summary.lambda == 'median') {
	  lambda<-median(lambda.post)
	} else{
	  stop('posterior.summary.lambda must be either mean or median.')
	}
	
	# Credible intervals for beta
	lowerbeta<-colQuantiles(beta.post,prob=(1 - beta.ci.level)/2)
	upperbeta<-colQuantiles(beta.post,prob=(1 + beta.ci.level)/2)
	
	# Assign names
	names(beta)<-names(lowerbeta)<-names(upperbeta)<-colnames(beta.post)<-colnames(x)
	
	# Credible intervals for lambda
	lambdaci<-credInt(lambda.post, cdf = NULL,conf=lambda.ci.level, type="twosided")      

	# Track stopping time
	stop.time <- Sys.time()
	time<-round(difftime(stop.time, start.time, units="min"),3)
	cat(c("Job finished at:",date()),fill=TRUE)
	cat("Computational time:", time, "minutes \n")
	
	# Return output
	return(list(time = format(time),
	            beta = beta,
	            lowerbeta = lowerbeta,
	            upperbeta = upperbeta,
	            lambda = lambda,
	            lambdaci = lambdaci,
	            beta.post = beta.post,
	            sigma2.post = sigma2.post,
	            lambda.post = lambda.post))
	}


