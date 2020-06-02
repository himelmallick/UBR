###################################################################################### 
# Helper functions for MH sampler for alpha - adapted from the BayesBridge R package #
######################################################################################

#############################################################
# Helper functions for Truncated Inverse Gamma Distribution #
#############################################################

# Truncated inverse gamma distribution
rtinvgamma<-function(n,a,b,L,U){
  l.prob<-pgamma(1/U,shape=a,rate=b)
  u.prob<-pgamma(1/L,shape=a,rate=b)
  1/qgamma(runif(n,l.prob,u.prob),shape=a,rate=b)
}

#######################
# Sampling from alpha #
#######################

# Posterior likelihood for alpha
llh.alpha <- function(alpha, beta,lambda, group=group)
{
  p = length(beta);
  gsize=group$gsize;
  G<-length(gsize);
  dd<-c();
  for (i in 1:G){
    q<-sum(abs(beta[group$index[[i]]]))^alpha
    dd<-c(dd,q)
  }
  
  sum((gsize/alpha) * log(lambda)) - sum(lgamma((gsize/alpha)+1)) - sum(lambda*dd);
}

# Posterior of alpha with a uniform prior
draw.alpha <- function(alpha, beta, lambda, group, pr.a, pr.b,ep=0.1)
{
  
  a.old = alpha;
  a.new=rbeta(1,pr.a,pr.b)
  
  log.accept = llh.alpha(a.new, beta,lambda,group) - llh.alpha(a.old, beta,lambda,group) +
    dbeta(a.new, pr.a, pr.b, log=TRUE) - dbeta(a.old, pr.a, pr.b, log=TRUE);  
  if (runif(1) < exp(log.accept)) alpha = a.new
  
  return(alpha);
}

# Posterior of alpha with a beta prior
draw.alpha.2 <- function(alpha, beta, lambda,group, ep=0.1)
{
  
  a.old = alpha;
  lb = 0.05
  ub = 0.95
  
  ## We might want to change these bounds if we have some belief about the value
  ## of alpha.
  l.new = max(c(lb, a.old-ep))
  r.new = min(c(ub, a.old+ep))
  d.new = r.new - l.new
  a.new = runif(1, l.new, r.new);
  
  d.old = min(c(ub, a.new+ep)) -  max(c(lb, a.new-ep))
  
  
  log.accept = llh.alpha(a.new,beta,lambda,group=group) - llh.alpha(a.old, beta,lambda,group=group);
  
  if (runif(1) < exp(log.accept)) alpha = a.new
  
  return(alpha);
}


##################################################
# Utility functions for creating group and index #
##################################################

# Create group information 
create.group<-function(group){
  size<-c()
  G<-length(unique(group))
  for (i in 1:G)
  {
    size[i]<-length(group[group==unique(group)[i]])
  }
  
  ss<-cumsum(size)
  
  x <- list(NA,NA)
  
  names(x)<-c("gindex","gsize")
  x$gsize<-size
  
  
  x$gindex <- list(NA)
  x$gindex[[1]]<-seq(1:x$gsize[1])
  
  for (i in 2:G){
    x$gindex[[i]]<-seq((ss[i-1]+1),ss[i],1)
  }
  
  return(x)
}

# Create index from grouping information
create.index<-function(group){
  W<-list()
  G<-length(group$gsize)
  for (g in 1:G){
    W0 <- replicate(group$gsize[[g]], group$gindex[[g]], FALSE)
    W<-append(W,W0);
  }
  return(W);
}


############################################################
# Truncated Multivariate Normal Sampler On L1-Norm Simplex #
############################################################

rtmvn.simp.1<- function(n,mu,sigma,bound,init){
  
  p <- ncol(sigma) # number of parameters, i.e. length of beta vector
  mu<-as.vector(mu)  
  sigma<-as.matrix(sigma)
  keep<-matrix(0,n,p)
  
  
  for (i in 1:n){
    for (j in 1:p){  # p is dimension of the multivariate normal
      sj<-sigma[-j,j] # remove jth row and extract rth column
      Sj <- sigma[-j,-j]  # p-1 by p-1 matrix by removing the jth column and jth row of Sigma   
      Sj<-solve(Sj)
      mean<-as.matrix(mu[j] + t(sj)%*%Sj%*%(init[-j]-mu[-j]))
      sd<-sqrt(sigma[j,j] - t(sj)%*%Sj%*%sj)
      upper <- bound - sum(abs(init[-j]))
      if (upper<0) upper<-0;
      lower<- - upper
      keep[i,j] <- tmvmixnorm::rtuvn(1,mean=mean, sd=sd,lower=lower,upper=upper)	 
      init[j]<-keep[i,j]
    }
  }
  return(keep)
}

# Sampling from conditional multivariable normal 
condNormal <- function(x.given, mu, sigma, given.ind, req.ind){
  # Returns conditional mean and variance of x[req.ind] 
  # Given x[given.ind] = x.given
  # where X is multivariate Normal with
  # mean = mu and covariance = sigma
  # 
  B <- sigma[req.ind, req.ind]
  C <- sigma[req.ind, given.ind, drop=FALSE]
  D <- sigma[given.ind, given.ind]
  CDinv <- C %*% solve(D)
  cMu <- c(mu[req.ind] + CDinv %*% (x.given - mu[given.ind]))
  cVar <- B - CDinv %*% t(C)
  list(condMean=cMu, condVar=cVar)
}

