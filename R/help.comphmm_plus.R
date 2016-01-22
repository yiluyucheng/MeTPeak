help.comphmm_plus <- function(x,y,trans,ab,F,h,Nit=100,plotflag=FALSE){

	m    = nrow(F)
	d    = F %*% (matrix(ab)) - h
	D    = diag(c(1/d))
	Nrep = ncol(x)
	K    = ncol(ab)
	N    = nrow(x)
	dnx  = .help.factorial(x);
	dny  = .help.factorial(y);
  n = x + y
	dn   = .help.factorial(n)
	logl = numeric(Nit) 
	#initialization
	dens <- matrix(0,N,K)
	forwrd <- matrix(0,N,K)
	bckwrd <- matrix(0,N,K)
	scale <- rbind(1,matrix(0,N-1))
	mu = 1e6
  PRECISION = 1e6
  t_INITIAL = 1
	NTTOL = 1e-6
	ALPHA = 0.1
	BETA = 0.7
  
  
  for (nit in seq(1,Nit)){
  	# initial the barrier step
  	t_step = t_INITIAL
    # 1: E-step to compute density
    dens <- exp( lgamma(matrix(1,N,1) %*% colSums(ab))*Nrep -
    t(rowsum(lgamma(t(matrix(c(matrix(c(n)) %*% matrix(1,1,K) + 
                                 matrix(1,N*Nrep,1) %*% (colSums(ab))),nrow=N))),rep(1:K,each=Nrep)) ) +
    t(rowsum(lgamma(t(matrix(c(matrix(c(x)) %*% matrix(1,1,K) + 
                                 matrix(1,N*Nrep,1) %*% (ab[1,])),nrow=N))),rep(1:K,each=Nrep)) ) +
    t(rowsum(lgamma(t(matrix(c(matrix(c(y)) %*% matrix(1,1,K) + 
                                 matrix(1,N*Nrep,1) %*% (ab[2,])),nrow=N))),rep(1:K,each=Nrep)) ) -
    lgamma(matrix(1,N,1) %*% ab[1,])*Nrep -
    lgamma(matrix(1,N,1) %*% ab[2,])*Nrep -
    dnx %*% matrix(1,Nrep,K) - dny %*% matrix(1,Nrep,K) + dn %*% matrix(1,Nrep,K) ) 
       #     dens <- apply(dens,2,function(x) {x[x==0] <- .Machine$double.xmin; return (x)})
    # 2: E-step forward recursion
    forwrd[1,] <- dens[1,]/K
    for (t in 2:N){
      forwrd[t,] <- (forwrd[t-1,] %*% trans) * dens[t,]
      #system scaling
      scale[t] <- sum(forwrd[t,])
      forwrd[t,] <- forwrd[t,] / scale[t]
    }
    #compute log-likelyhood
    logl[nit+1] <- log(sum(forwrd[N,])) + sum(log(scale))
            # message(sprintf('Iteration %d:\t%.4f', (nit), logl[nit+1]))
    
    #3: E-step backward recursion
    bckwrd[N,] <- matrix(1,1,K)
    for (t in (N-1):1){
      bckwrd[t,] <- (bckwrd[t+1,] * dens[t+1,]) %*% t(trans)
      bckwrd[t,] <- bckwrd[t,] / scale[t]
    }
    
    #4: M-step reestimate transition matrix
    trans <- trans * (t(forwrd[1:(N-1),]) %*% (dens[2:N,] * bckwrd[2:N,] ))
    trans <- trans / (apply(trans,1,sum)) %*% matrix(1,1,K) 
    #5: M-step reestimate parameters
    postprob <- forwrd * bckwrd
    postprob <- postprob / (apply(postprob,1,sum)) %*% matrix(1,1,K)
    wght = colSums(postprob)

    while (t_step <= PRECISION){
    	res_obj = help.obj_plus(t_step,postprob,wght,x,y,n,ab,F,h,m,d,D,N,Nrep)
    	
    	# if res_obj is singular, then break
    	errFlag = FALSE
    	errorHandler <- tryCatch({        # result from Gaussian elimination (H is singular)
    	  v = -solve(res_obj$H,res_obj$g)
    	}, simpleError = function(e) {  # only catch simpleErrors
    	  errFlag = TRUE
    	})
    	if (errFlag) {break}
    	
    	lambda = t(res_obj$g) %*% v
    	# perform backtracking forward line search
    	s = 1
    	while(min(h-F%*%(matrix(ab)+s*v))<0){ s = BETA*s }

    	res_obj_new = help.obj_plus(t_step,postprob,wght,x,y,n,ab+matrix(s*v,2,K),F,h,m,d,D,N,Nrep)
    	while (res_obj_new$f > res_obj$f + ALPHA*s*lambda){
    		s = BETA * s
    	 	res_obj_new = help.obj_plus(t_step,postprob,wght,x,y,n,ab+matrix(s*v,2,K),F,h,m,d,D,N,Nrep)
    	}
    	if (abs(lambda/2) < NTTOL){break}

    	ab = ab + matrix(s*v,2,K)
    	t_step = mu * t_step
    }
    
    # likelihood does not increase 
    if (abs((logl[nit+1]-logl[nit])/logl[nit+1])<1e-1){break}
    
 }   
  
  list(ab=ab,trans=trans,postprob=postprob,logl=logl[2:nit])

}