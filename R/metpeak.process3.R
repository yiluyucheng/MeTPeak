.metpeak.process3 <- function(IP,INPUT,batch_id,minimal_counts_in_fdr=10){
#   using betabinomial.hmm to call peaks per gene
  print("Peak Calling Stage!")
  Ng = unique(batch_id)
  nip <- ncol(IP)
  nin <- ncol(INPUT)
  INPUT_mean <- rowMeans(INPUT)
  IP_mean <- rowMeans(IP)
  if (nip > nin) {
    avg_input <- round(matrix(rep(INPUT_mean,nip-nin),ncol=nip-nin))
    INPUT <- cbind(INPUT,avg_input) 
  }
  else if (nip < nin){
    avg_ip <- matrix(rep(IP_mean,nin-nip),ncol=nin-nip)
    IP <- cbind(IP,avg_ip) 
  }

# initialize the HMM parameters
K = 3;
# opt = .help.optcond(K=K,c(0.7,0.4))
opt = .help.optcond(K=K)
F = opt$F
h = opt$h
# ab = .help.initpar(K=K) # automatically generate the inital parameters
ab = rbind(c(12,8,1),c(2,6,5))
ip_mean = mean(IP_mean) #averaged reads per nucleid
  
  # initialize the variables
  pvalues <- rep(1,nrow(IP))
  for (ii in Ng){
    print(ii)
    flag <- batch_id==ii
    ip=as.matrix(IP[flag,])
    input=as.matrix(INPUT[flag,])    
    
    if (sum(ip)<ip_mean*nip){
      res = list(ab,matrix(1/K,K,K),matrix(0,nrow(ip),K) )
      res[[3]][,K] = matrix(1,nrow(ip))
      cl = list(class = matrix(3,nrow(ip)) )
      
    }else{
#     
      # using C++ to do the computation
      res = cmpHmm(ip,input+1,matrix(1/K,K,K),ab,F,h) #res[[2]] = trans; to avoid 0, input + 1;
      cl = help.compvit(ip,input+1,res[[2]],res[[1]])
      
    }
    
    # find the most significant region of a gene and treat it as a potential peak region
    maxID = which.max(res[[1]][1,]/colSums(res[[1]])) # The most significant class res[[1]] = ab
#     plot the result for future view usage
#     dot <- 1:nrow(ip); 
#     matplot(cbind(ip,input),col = 1:(ncol(ip)*2));
#     points(dot[cl$class==maxID],rep(max(ip)/2,sum(cl$class==maxID))) # test part    
#     anno = .get.gene.anno(ii,ANNOTATION,ANNOTATION_BATCH_ID)
#     title(anno$gene)

    # determine the peak region within a gene  
    peak_loci <- which(cl$class==maxID)
    peak_loci_end <- length(peak_loci) # no peak condition included
    if (peak_loci_end > 0){
      peak_id <- c(0,which( (peak_loci[-1] - peak_loci[-peak_loci_end] ) > 1 ),peak_loci_end)
      for (jj in 1:(length(peak_id)-1)){
        jjpeak <- peak_loci[peak_id[jj]+1]:peak_loci[peak_id[jj+1]]
        res[[3]][jjpeak,maxID] <- median(res[[3]][jjpeak,maxID])  #res[[3]] = res$postprob
      }
    }  
    pvalues[flag] = 1 - res[[3]][,maxID]  
  }

  log.fdr=log(p.adjust(pvalues,method='fdr'))

  # with significant number of reads only
  ID=which( (IP_mean+INPUT_mean) > minimal_counts_in_fdr) #should be vector not matrix
  log.fdr_sig=log(p.adjust(pvalues[ID], method = "fdr"))
  log.fdr[ID]=log.fdr_sig
  # fold enrichment
  log.fc=log(IP_mean/(INPUT_mean+1))
  
  # output result
  PW=list(log.p=log(pvalues),log.fdr=log.fdr,log.fc=log.fc)
  return(PW)
  
}