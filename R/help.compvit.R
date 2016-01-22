help.compvit <- function(x,y,trans,ab)
  # use viterbi algorithm to deconvolute the status
{
  K <- ncol(trans) #two states in our model
  Nrep = ncol(x)
  N <- nrow(x)
  dnx <- .help.factorial(x)
  dny <- .help.factorial(y)
  dn  <- .help.factorial(x+y)
  dens <- exp( lgamma(matrix(1,N,1) %*% colSums(ab))*Nrep -
    t(rowsum(lgamma(t(matrix(c(matrix(c(x+y)) %*% matrix(1,1,K) + 
                                 matrix(1,N*Nrep,1) %*% (colSums(ab))),nrow=N))),rep(1:K,each=Nrep)) ) +
    t(rowsum(lgamma(t(matrix(c(matrix(c(x)) %*% matrix(1,1,K) + 
                                 matrix(1,N*Nrep,1) %*% (ab[1,])),nrow=N))),rep(1:K,each=Nrep)) ) +
    t(rowsum(lgamma(t(matrix(c(matrix(c(y)) %*% matrix(1,1,K) + 
                                 matrix(1,N*Nrep,1) %*% (ab[2,])),nrow=N))),rep(1:K,each=Nrep)) ) -
    lgamma(matrix(1,N,1) %*% ab[1,])*Nrep -
    lgamma(matrix(1,N,1) %*% ab[2,])*Nrep -
    dnx %*% matrix(1,Nrep,K) - dny %*% matrix(1,Nrep,K) + dn %*% matrix(1,Nrep,K) )  
  logdens <- log(dens + .Machine$double.xmin) # to avoid inf
  
  #dynamic programming recursion
  trans <- log(trans + .Machine$double.xmin)
  plogl <- matrix(0,N,K)
  bcktr <- matrix(0,N-1,K)
  
  plogl[1,] <- logdens[1,] - log(K)
  for (t in 2:N){
    tmp <- matrix(plogl[t-1,,drop=F]) %*% matrix(1,1,K) + trans
    plogl[t,] <- apply(tmp,2,max)
    bcktr[t-1,] <- max.col(t(tmp))
    plogl[t,] <- plogl[t,] + logdens[t,]
    
  }
  
  class <- matrix(0,nrow=N)
  class[N] <- which.max(plogl[N,])
  
  for (t in (N-1):1){
    class[t] <- bcktr[t,class[t+1]]
  }
  list(class=class, logl=max(plogl[N,]))
}

