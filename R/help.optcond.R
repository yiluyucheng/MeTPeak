.help.optcond <- function(K=2,statsep=NULL,maxup=100){
# maxup: upbound of the parameters
# statsep : hard bound for seperating different states
  F = rbind(diag(-1,2*K,2*K),diag(1,2*K,2*K))
  h = rbind(matrix(0,2*K,1),matrix(maxup,2*K,1))
  if (is.null(statsep)){
    return( list(F=F,h=h) )
  } else{
    statsep = sort(statsep,decreasing=TRUE) # first component is the biggest methylation degree
    Ft = c() 
    for(k in 2:length(statsep)){
      Ft = rbind(Ft,
                 c(matrix(0,1,2*(k-1)),1-statsep[k-1],-statsep[k-1],matrix(0,2*(length(statsep)+1-k),1)),
                 c(matrix(0,1,2*(k-1)),statsep[k]-1,statsep[k],matrix(0,2*(length(statsep)+1-k),1)) )
    }    
    F = rbind(F,
              c(statsep[1]-1,statsep[1],matrix(0,2*length(statsep),1)),
              Ft,
              c(matrix(0,1,2*(K-1)),1-statsep[K-1],-1*statsep[K-1])) 
    h = rbind(h,matrix(0,2*length(statsep)))
    return(list(F=F,h=h) )
  }
}