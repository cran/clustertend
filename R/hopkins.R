# hopkins.R

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Package `clustertend` is deprecated.  Use package `hopkins` instead.")
}

## ---------------------------------------------------------------------------

#' @title Calculate the Hopkins' statistic
#'
#' Calculate the Hopkins' statistic of given data.
#'
#' Note:"Package \code{clustertend} is deprecated.  Use package \code{hopkins} instead.
#' 
#' Sample data must be preprocessed into dataframe or matrix form before given as the value of parameter "data".
#'
#' @param X Data (matrix or data.frame) to check clusterability.
#' 
#' @param n The number of rows to sample from X. The default is 1/10th the number of rows of X.
#' 
#' @return The value returned is actually 1-Hopkins statistic.
#' 
#' @author Luo YiLan, Zeng RuTong.
#' 
#' @examples
#' set.seed(1)
#' hopkins(iris[,-5], n=15)
#' 
#' @references 
#' Lawson, R.G. and Jurs, P.C.(1990).
#' New index for clustering tendency and its application to chemical problems.
#' Journal of Chemical Information and Computer Sciences. 30(1):36-41.
#' 
#' @importFrom stats dist runif
#' @export
#' 
hopkins <- function(X,n=ceiling(nrow(X)/10)) {
  .Deprecated(msg="Package `clustertend` is deprecated.  Use package `hopkins` instead.")

  if(is.data.frame(X))
    X <- as.matrix(X)
  if (!(is.matrix(X)))
    stop("X must be data.frame or matrix")
  
  if(n>=nrow(X))
    stop("n must be no larger than num of samples")
  c <- apply(X,2,min) # minimum value per colume
  d <- apply(X,2,max)
  p <- matrix(0,ncol=ncol(X),nrow=n) # n vectors of space
  for(i in 1:ncol(X)) {
    p[,i] <- runif(n,min=c[i],max=d[i])
  }

  #k <- round(runif(n,1,nrow(X)))
  k <- sample(1:nrow(X), n)
  
  q <- as.matrix(X[k,])
  distp <- rep(0,nrow(X))
  #distq=rep(0,nrow(X)-1)
  distq <- 0;
  minp <- rep(0,n)
  minq <- rep(0,n)
  for(i in 1:n) {
    distp[1] <- dist(rbind(p[i,],X[1,]))
    minqi <- dist(rbind(q[i,],X[1,]))
    for(j in 2:nrow(X)) {
      distp[j] <- dist(rbind(p[i,],X[j,]))
      error <- q[i,]-X[j,]
      if(sum(abs(error))!=0) {
        #distq[j] <- dist(rbind(q[i,],X[j,]))
        distq <- dist(rbind(q[i,],X[j,]))
        if(distq<minqi)
          minqi <- distq;
      }
    }
    minp[i] <- min(distp)
   # minq[i] <- apply(distq,1,min)
   minq[i] <- minqi;
  }
  list(H=(sum(minq)/(sum(minp)+sum(minq))))
  
}

