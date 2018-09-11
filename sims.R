rm(list=ls())
require(MASS);require("clue")
source("EPalgosim.R") 
source("functions.R")

n1 = 200; n2= 200 # community size
lo = 0.2; hi = 0.8;  # see model setup in paper
n1.lo = n1*0.50; n1.hi = n1*0.50
H = (3:8)/2
nsim = 100

b.true = c(rep(1,n1),rep(2,n2))
N1 <- 1:n1; N2 <- n1 + 1:n2
lambda = matrix(NA, nrow = length(b.true), ncol = max(b.true))
E <- rep(NA, length(H)) # 2 x expected no. of edges in the network
PA.err = matrix(NA, nrow=nsim, ncol=length(H))   # community detection error
DC.err = matrix(NA, nrow=nsim, ncol=length(H))   
E.PA = matrix(NA, nrow=nsim, ncol=length(H))    # popularity estimation error
E.DC = matrix(NA, nrow=nsim, ncol=length(H))    
aic.PA = matrix(NA, nrow=nsim, ncol=length(H))  # AIC
aic.DC = matrix(NA, nrow=nsim, ncol=length(H))

ptm <- proc.time()
for (ihom in 1:length(H)){
  h = H[ihom]		# homophily parameter
  foo1 <- sqrt(h)/sqrt(1+h); foo2 <- 1/sqrt(1+h)
  print(proc.time()-ptm); print(h)
  
  # ##### Construct lambda #####
  lambda[1:n1.hi,1] <- hi; lambda[(n1.hi+1:n1.lo),1] <- lo
  lambda[1:n1.hi,2] <- lo; lambda[(n1.hi+1:n1.lo),2] <- hi
  lambda[(n1+1:n1.hi),1] <- lo; lambda[(n1+n1.hi+1:n1.lo),1] <- hi
  lambda[(n1+1:n1.hi),2] <- hi; lambda[(n1+n1.hi+1:n1.lo),2] <- lo
  lambda[N1,1] <- foo1*lambda[N1,1]; lambda[N1,2] <- foo2*lambda[N1,2]
  lambda[N2,2] <- foo1*lambda[N2,2]; lambda[N2,1] <- foo2*lambda[N2,1]
  P.true <- P_sim(lambda,b.true)
  mu.true <- mu.hat(P.true,b.true)
  E[ihom] <- sum(P.true) # This is TWICE the expected network degree
  
  for (isim in 1:nsim){
    set.seed(1729+isim)
    A = generate(lambda,b=b.true)  # generate network
    b.can = EPalgo(A,eps=0) # EP algorithm (no perturbation)
    Q.PA.can = rep(NA, ncol(b.can))	# array to store Q values
    Q.DC.can = rep(NA, ncol(b.can))	# array to store Q values
    for (i in 1:ncol(b.can)){
      #check if any cluster is empty
      foo = rep(NA, max(b.true))
      for (clus in 1:max(b.true)) {foo[clus]=sum(b.can[,i]==clus)}
      if (min(foo)==0) {stop('Empty groups are not allowed')} 
      Q.PA.can[i] = Q.PA(A, b=b.can[,i])   # fit PABM
      Q.DC.can[i] = Q.DC(A, b=b.can[,i])   # fit DCBM
    } # end of i for loop
    foo1 = order(-Q.PA.can)[1] 
    b.PA = b.can[,foo1]   # community assignment that maximises Q.PA
    foo2 = order(-Q.DC.can)[1]
    b.DC = b.can[,foo2]   # community assignment that maximises Q.DC
    
    ##### Analysis: 1. Community detection #####
    x = as.cl_hard_partition(b.true)
    y1 = as.cl_hard_partition(b.PA)
    y2 = as.cl_hard_partition(b.DC)
    PA.err[isim,ihom] = 1-cl_agreement(x,y1,method="manhattan")
    DC.err[isim,ihom] = 1-cl_agreement(x,y2,method="manhattan")
    
    ##### Analysis: 2. Estimation accuracy #####
    if (table(b.true,b.PA)[1,1]<table(b.true,b.PA)[1,2]) {   # check if label rotation is needed
      b.PA <- ifelse(b.PA==1,2,1)}; table(b.true,b.PA) # label rotation, if needed
    if (table(b.true,b.DC)[1,1]<table(b.true,b.DC)[1,2]) {   # check if label rotation is needed
      b.DC <- ifelse(b.DC==1,2,1)}; table(b.true,b.DC) # label rotation, if needed
    P.PA.hat = P_PA(A,b.PA)     # prob matrix fitted by PABM by comm detection
    P.PA.true = P_PA(A,b.true)  # prob matrix fitted by PABM from true assignment
    P.DC.hat = P_DC(A,b.DC)     # prob matrix fitted by DCBM by comm detection
    P.DC.true = P_DC(A,b.true)  # prob matrix fitted by DCBM from true assignment
    mu.PA.hat <- mu.hat(P.PA.hat,b.PA)    # popularity matrix fitted by PABM by comm detection
    mu.PA.true <- mu.hat(P.PA.true,b.true) # popularity matrix fitted by PABM by comm detection
    mu.DC.hat <- mu.hat(P.DC.hat,b.DC)    # popularity matrix fitted by DCBM by comm detection
    mu.DC.true <- mu.hat(P.DC.true,b.true) # popularity matrix fitted by DCBM by comm detection
    E.PA[isim,ihom] <- sum((mu.PA.hat - mu.true)^2)/E[ihom]
    E.DC[isim,ihom] <- sum((mu.DC.hat - mu.true)^2)/E[ihom]
    
    ##### Analysis: 3. AIC/BIC (model fit adjusted for model size) #####
    n=n1+n2
    aic.PA[isim,ihom] <- (4*n-2)-2*log.l(A,P.PA.hat)
    aic.DC[isim,ihom] <- 2*(n+1)-2*log.l(A,P.DC.hat)
    print(isim)}} # isim, ihom

proc.time()-ptm
