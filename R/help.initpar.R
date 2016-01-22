.help.initpar <- function(K=2,meanval=NULL){
	if (is.null(meanval)){
		meanval = seq(1/(K+1),1-1/(K+1),by=1/(K+1))
	}else{
		K = length(meanval)
	}
	meanval = sort(meanval,decreasing=TRUE)# make sure the first component is the biggest
	a = sort(runif(K,1,8),decreasing=TRUE) # randomly generate the variable a
	b = a/meanval - a
	return (ab=data.matrix(rbind(a,b)))
}