#Main procedure
ssd = function(p1,p2,k,ratio=1,alpha=0.05,beta=0.2,cc=0.02,d=0.2,r=0.3,m,scheme='M2',Niter=500) {
  #Original method
  orig_sol = function(p1,p2,lambda0){
             p=0.5*(p1+p2)
             temp_orig=sum(((p1-p2)^2)/p)
             ceiling(2*lambda0/temp_orig) 
  }
  #Minimum difference method
  min_sol = function(p1,p2,lambda0,cc){
            p=0.5*(p1+p2)
            temp_min=abs(p1-p2)
            temp_min[temp_min<cc]=cc
            temp_min=sum(temp_min^2/p)
            ceiling(2*lambda0/temp_min) 
  }
  #Correction method
  corr_D = function(p1,p2,tj){
           k=length(p1)
           p=0.5*(p1+p2)
           zj=rnorm(2*k,mean=0,sd=1)
           vB=2*(p1-p2)*sqrt(tj)*zj[1:k]/p
           vC=(p1-p2)^2*sqrt(tj)*zj[(k+1):(2*k)]/(2*p)
           sum(vB)-sum(vC) 
  }   
  corr_sol = function(p1,p2,n0,lambda0,Niter){
             k=length(p1)
             p=0.5*(p1+p2)
             a=sum(((p1-p2)^2)/p)
             A=2*lambda0*(1/sqrt(n0))/(a^2)
             tj=p1*(1-p1)+p2*(1-p2)
             ntemp=replicate(Niter,corr_D(p1,p2,tj))
             ceiling(orig_sol(p1,p2,lambda0)+A*mean(ntemp)) 
  }            
  #Parametric bootstrap method
  gen_boot = function(p1,p2,n0,lambda0){
             sample1=rmultinom(1,n0,p1)
             sample2=rmultinom(1,n0,p2)
             estp1=sample1/n0
             estp2=sample2/n0
             estp=0.5*(estp1+estp2)
             while (min(estp)==0){
                    sample1=rmultinom(1,n0,p1)
                    sample2=rmultinom(1,n0,p2)
                    estp1=sample1/n0
                    estp2=sample2/n0
                    estp=0.5*(estp1+estp2) 
             }
             orig_sol(estp1,estp2,lambda0) 
  }
  boot_sol = function(p1,p2,n0,Niter,lambda0){
             boottemp=replicate(Niter,gen_boot(p1,p2,n0,lambda0))
             c(ceiling(mean(boottemp)),ceiling(median(boottemp)), ceiling(quantile(boottemp,0.75)), ceiling(quantile(boottemp,0.80))) 
  }
  #Unbalanced Parametric bootstrap method
  ugen_boot = function(p1,p2,n1,n2,lambda0){
              sample1=rmultinom(1,n1,p1)
              sample2=rmultinom(1,n2,p2)
              estp1=sample1/n1
              estp2=sample2/n2
              estp=0.5*(estp1+estp2)
              while (min(estp)==0){
                    sample1=rmultinom(1,n1,p1)
                    sample2=rmultinom(1,n2,p2)
                    estp1=sample1/n1
                    estp2=sample2/n2
                    estp=0.5*(estp1+estp2) 
              }
              unbalance_sol(estp1,estp2,n1,n2,lambda0) 
  }
  uboot_sol = function(p1,p2,n1,n2,Niter,lambda0){
              boottemp=replicate(Niter,ugen_boot(p1,p2,n1,n2,lambda0))
              c(ceiling(mean(boottemp)),ceiling(median(boottemp)), ceiling(quantile(boottemp,0.75)), ceiling(quantile(boottemp,0.80))) 
  }
  m1_sol = function(r,d,k,lambda){
           ceiling( 2*lambda/(r*d*k) ) 
  }
  unbalance_sol = function(p1,p2,n1,n2,lambda0) {
                  r=n1/n2
                  p=0.5*(p1+p2)
                  temp_orig=sum(((p1-p2)^2)/p)
                  nn1=ceiling((1+r)*lambda0/temp_orig)
                  nn2=ceiling(nn1/r) 
                  c(nn1,nn2)
  }
  #Main part
  if (!missing(p1)) {
         if ((min(p1)<0)|(max(p1)>1)|(sum(p1)!=1)){
            stop("The specified p1 is not valid")
         }
      }
      if (!missing(p2)) {
         if ((min(p2)<0)|(max(p2)>1)|(sum(p2)!=1)){
            stop("The specified p2 is not valid")
         }
      }
      if (!missing(k)) {
         if ((k<0)|(floor(k)!=k)){
            stop("The specified k is not valid")
         }
      }
      if (!missing(ratio)) {
         if ((ratio<0)|(ratio>1)){
            stop("The specified ratio is not valid")
         }
      }
      if (!missing(alpha)) {
         if ((alpha<0)|(alpha>1)){
            stop("The specified alpha is not valid")
         }
      }
      if (!missing(beta)) {
         if ((beta<0)|(beta>1)){
            stop("The specified beta is not valid")
         }
      }
      if (!missing(d)) {
         if ((d<0)|(d>1)){
            stop("The specified d is not valid")
         }
      }
      if (!missing(r)) {
         if (r<0){
            stop("The specified r is not valid")
         }
      }
      if (!missing(scheme)) {
         if (!(scheme %in% c('M1','M2'))) {
            stop("The specified scheme is not valid")
         }
      }
      if (!missing(Niter)) {
         if ((Niter<0)|(floor(Niter)!=Niter)){
            stop("The specified Niter is not valid")
         }
      }
      if (!missing(cc)) {
         if ((cc<0)|(cc>1)){
            stop("The specified cc is not valid")
         }
      }
      if (!missing(m)) {
         if ((m<0)|(floor(m)!=m)){
            stop("The specified m is not valid")
         }
      }
      if (!missing(m)&!missing(ratio)) {
         n2=ceiling(m/ratio)
      }
      kk=0
      if (!missing(p1)&!missing(p2)) {
         if (length(p1)!=length(p2)) {
            stop("The specified p1 and p2 are not valid")
         } else {
            kk=length(p1)
         }
      }
      if (!missing(k)) {
         kk=k
      }
      if (kk!=0) {
         lambda=1
         ncchi=qchisq(beta,df=kk-1,ncp=lambda)
         cchi=qchisq(1-alpha,df=kk-1)
         while (ncchi<cchi) {
               lambda=lambda+0.01
               ncchi=qchisq(beta,df=kk-1,ncp=lambda)
               cchi = qchisq(1-alpha, df=kk-1) 
         }
      } else {
         stop("Either p1 and p2 or k needs to be specified")
      }
      if (scheme=='M2') {
         if (!missing(p1)&!missing(p2)) {
            if (ratio==1) {
               cat("The calculated sample sizes under scheme M2 are following:\n")
               cat(paste("Original Method:",orig_sol(p1,p2,lambda),"\n",sep=' '))
               cat(paste("Minimum Difference Method:",min_sol(p1,p2,lambda,cc),paste("(cc=",cc,")",sep=''),"\n",sep=' '))
               if (!missing(m)) {      
	                cat(paste("Correction Method:",corr_sol(p1,p2,m,lambda,Niter),"\n",sep=' '))
                  boot_result=boot_sol(p1,p2,m,Niter,lambda)
	                cat(paste("Bootstrap Mean Method:",boot_result[1],"\n",sep=' '))
	                cat(paste("Bootstrap Median Method:",boot_result[2],"\n",sep=' '))
	                cat(paste("Bootstrap 75% Percentile Method:",boot_result[3],"\n",sep=' '))
	                cat(paste("Bootstrap 80% Percentile Method:",boot_result[4],"\n",sep=' '))
               }
            }
            if (ratio!=1) {
               if (!missing(m)) {
                  cat("The calculated sample sizes under scheme M2 are following:\n")
                  nnn=unbalance_sol(p1,p2,m,n2,lambda)
                  cat(paste("Original Method: n1=",nnn[1],";","n2=",nnn[2],"\n",sep=' ') )
                  uboot_result=uboot_sol(p1,p2,m,n2,Niter,lambda)
	                cat(paste("Bootstrap Mean Method: n1=",uboot_result[1],";","n2=",ceiling(uboot_result[1]/ratio),"\n",sep=' '))
	                cat(paste("Bootstrap Median Method: n1=",uboot_result[2],";","n2=",ceiling(uboot_result[2]/ratio),"\n",sep=' '))
	                cat(paste("Bootstrap 75% Quantile Method: n1=",uboot_result[3],";","n2=",ceiling(uboot_result[3]/ratio),"\n",sep=' '))
	                cat(paste("Bootstrap 80% Quantile Method: n1=",uboot_result[4],";","n2=",ceiling(uboot_result[4]/ratio),"\n",sep=' '))
               } else {
                  stop("m needs to be specified")
               }
            }
         } else {
           stop("p1 and p2 need to be specified")
         }
      }
      if (scheme=='M1') {
         if (!missing(k)) {
            cat("The calculated sample sizes under scheme M1 are following:\n")
            cat(paste("Scheme M1:",m1_sol(r,d,k,lambda),"\n",sep=' '))
         } else {
            stop("The parameter setting is not valid")
         }
      }
}



