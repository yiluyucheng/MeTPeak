help.obj_plus <- function(t, post, wght, x, y, n, ab, F, h, m, d, D, N, Nrep){
# try to avoid repeated computation and can not afford the time consumption
  # m = nrow(F)
  # Nrep = dim(x)[2]
  # N = nrow(post)
  K = ncol(post)
  # wght = colSums(post)
  d = F %*% (matrix(c(ab))) - h
  D = diag(c(1/d))
  H_o = matrix(0,1,2*K-1) # initialize the vector for later assignment off the diagonal

  # compute the shared data plus parameters in advance
  n_a_b = t(matrix(c(matrix(c(n)) %*% matrix(1,1,K) + matrix(1,N*Nrep,1) %*% (colSums(ab))),nrow=N))
  x_a = t(matrix(c(matrix(c(x)) %*% matrix(1,1,K) + matrix(1,N*Nrep,1) %*% ab[1,]),nrow=N))
  y_b = t(matrix(c(matrix(c(y)) %*% matrix(1,1,K) + matrix(1,N*Nrep,1) %*% ab[2,]),nrow=N))

  # caculate the jacobian for multiple states (very complicated to write in a vector form avoiding for loop)
  J = matrix(1,2,1) %*% (digamma(colSums(ab)) * wght) * Nrep -
    matrix(1,2,1) %*% diag(rowsum(digamma(n_a_b),rep(1:K,each=Nrep)) %*% post) +
    rbind(diag(rowsum(digamma(x_a),rep(1:K,each=Nrep)) %*% post),
          diag(rowsum(digamma(y_b),rep(1:K,each=Nrep)) %*% post)) - 
    rbind( digamma(ab[1,]) * wght * Nrep, 
           digamma(ab[2,]) * wght * Nrep )

  # first compute the diagonal of the Hessian	
  H_d = matrix(1,2,1) %*% (trigamma(colSums(ab)) * wght) * Nrep -
    matrix(1,2,1) %*% diag(rowsum(trigamma(n_a_b),rep(1:K,each=Nrep)) %*% post) +
    rbind(diag(rowsum(trigamma(x_a),rep(1:K,each=Nrep)) %*% post),
          diag(rowsum(trigamma(y_b),rep(1:K,each=Nrep)) %*% post)) -
    rbind( trigamma(ab[1,]) * wght * Nrep,
           trigamma(ab[2,]) * wght * Nrep )
  # second compute the off diagonal elements of the Hessian
  H_o[seq(1,2*K-1,by=2)] =  trigamma(colSums(ab)) * wght * Nrep - 
    diag(rowsum(trigamma(n_a_b),rep(1:K,each=Nrep)) %*% post)
  H = diag(c(H_d))
  diag(H[-1,]) = c(H_o)
  diag(H[,-1]) = c(H_o)
  # likelihood 
  px = lgamma(matrix(1,N,1) %*% colSums(ab))*Nrep -
    t(rowsum(lgamma(n_a_b),rep(1:K,each=Nrep)) ) +
    t(rowsum(lgamma(x_a),rep(1:K,each=Nrep)) ) +
    t(rowsum(lgamma(y_b),rep(1:K,each=Nrep)) ) -
    lgamma(matrix(1,N,1) %*% ab[1,])*Nrep -
    lgamma(matrix(1,N,1) %*% ab[2,])*Nrep
  # return function values
  g = -t * matrix(c(J)) - t(F) %*% D %*% matrix(1,m,1)
  H = -t * H
  f = -t * sum(matrix(1,1,N) %*% (post*px)) - log(-t(d)) %*% matrix(1,m,1)
  return(list(f=f,g=g,H=H))

}
  
