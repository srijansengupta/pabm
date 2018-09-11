AdjList<-function(A){ # function that takes an adjacency matrix and gives out a nested list of neighbors
  N<-nrow(A)          # no. of nodes
  B<-vector("list",N) # create an empty list of length n
  for (i in 1:N){
    foo <- which(A[i,]==1)
    B[[i]]<-c(B[[i]],foo)  # the ith entry is the list of neighbors of node i.
  }
return(B)}

##### function for estimated popularity #####
mu.hat<-function(Phat,bhat){  
  N<-length(bhat)
  K<-max(bhat)
  mu.hat<-matrix(NA,nrow=N,ncol=K)
  for (k in 1:K){
    nodes = which(bhat == k)
    for (i in 1:N){  
    mu.hat[i,k] = sum(Phat[i,nodes])
      }}
  return(mu.hat)}

##### function for log-likelihood #####
log.l<-function(A,Phat){  #formula: \sum_{i<j} [A_{ij}log(P_{ij})-P_{ij}] - 1/2 \sum_i P_{ii}
  A[lower.tri(A,diag=TRUE)]=0
  foo=A*log(Phat)
  foo1=sum(foo,na.rm=TRUE)
  foo2=sum(diag(Phat)) # this is 0.5p_{ii}, see functions P_PA, P_DC 
  Phat[lower.tri(Phat,diag=TRUE)]=0
  foo3=sum(Phat)
  return(foo1-foo3-foo2)}
######################

##### Functions for PABM #####

#####################################################################
##### fn to calculate popularity M and block interaction O ##########
#####################################################################
f.PA<-function(A,b){	
K<-max(b)       # no. of communities
N<-nrow(A)      # no. of nodes
M<-matrix(NA,nrow=N,ncol=K)  # popularity matrix
O<-matrix(NA,nrow=K,ncol=K)  # community interaction matrix
for (i in 1:N){		# calculate M
  for (r in 1:K){
  nodes = which(b == r)
  M[i,r] = sum(A[i,nodes])
  }}
for (r in 1:K){		# calculate O
  for (s in r:K){
  nodes1 = which(b == r)
  nodes2 = which(b == s)
  O[r,s] = sum(A[nodes1,nodes2])
  O[s,r] = O[r,s]
  }}
list(M=M, O=O)}

#####################################
########## PABM Likelihood ##########
#####################################
Q.PA <- function(A, b){
foo<-f.PA(A,b)
O=foo$O; M = foo$M
s1 = sum(M*log(M),na.rm=TRUE) # na.rm = TRUE ignores M=0 cases as log(0) = NA
s2 = sum(O*log(O),na.rm=TRUE) # na.rm = TRUE ignores O=0 cases as log(0) = NA
return(2*s1-s2)}

#####################################
########## PABM AIC ##########
#####################################
AIC.PA <- function(A, b){
  n<- nrow(A)
  E<-sum(A)/2
  return(4*n-2-Q.PA(A,b)+2*E)}

BIC.PA <- function(A, b){
  n<- nrow(A)
  E<-sum(A)/2
  return((4*n-2)*log(n*(n+1)/2)-Q.PA(A,b)+2*E)}


#############################################################
############### MLE initiated with random clusters ########## 
#############################################################
PAmle <- function(A, k, n.iter,b){
  N = nrow(A)
    q.com = rep(NA, n.iter+1)  		# vessel for values of complete data likelihood
    b1 <- b
    q.com[1] = Q.PA(A=A,b=b1)
    
    for (iter in 1:n.iter){
      b0 = b1
      b.new = LabelSwitchPA(A,b0)
      b1 = b.new
      q.com[iter+1] = Q.PA(A,b1)
      Q1 = q.com[iter+1]
        if (q.com[iter+1] <= q.com[iter]) {
          b.max = b1;Q.max = Q1 ; break}
      }
list(q.last=Q.max,iter=iter,b.hat=b.max,lik.iter = q.com)}

############################################################
############## obtain MLE of probability matrix ########## 
############################################################
P_PA <- function(A, b){
  N<-ncol(A)
  foo<-f.PA(A,b)
  M<-foo$M; O<-foo$O
  lambda <- matrix(NA, nrow=N, ncol = max(b))
  P = matrix(NA,nrow=N,ncol=N)
  for (i in 1:N){
    s <- b[i]
    for (r in 1:max(b)){
      lambda[i,r] <- M[i,r]/sqrt(O[s,r])} }
  for (i in 1:N){  	# populate diagonal entries
    r = b[i]
    P[i,i] = (1/2)*(lambda[i,r]^2)}	# divide by 2 --- see notes section 3.3
    for (i in 1:(N-1)) {
    for (j in (i+1):N) {
      p = lambda[i,b[j]]*lambda[j,b[i]]
      P[i,j] = p
      P[j,i] = P[i,j]
    }}
  return(P)}














##### Functions for DCBM #####

##############################################################
####### Calculate interaction O and community stubs d ########
##############################################################
f.DC<-function(A,b){	
K<-max(b)
O<-matrix(NA,nrow=K,ncol=K)  # interaction matrix
for (i in 1:K){		# calculate O
for (j in i:K){
nodes1 = which(b == i)
nodes2 = which(b == j)
O[i,j] = sum(A[nodes1,nodes2])
O[j,i] = O[i,j]
}}
list(O=O, d=rowSums(O))}

##############################################
########## DCBM Likelihood ###################
##############################################
Q.DC <- function(A, b){
K <- max(b)
q <- matrix(0, nrow = K, ncol = K)
foo<-f.DC(A,b)
O<-foo$O; d<-foo$d 
for (i in 1:K){
for (j in 1:K){
if (O[i,j]>0){
q[i,j] = O[i,j]*log(O[i,j]/(d[i]*d[j]))}
}}
return(sum(q))}	# formula for Q

#####################################
########## DCBM AIC ##########
#####################################
AIC.DC <- function(A, b){
  n<- nrow(A)
  deg<-colSums(A); foo3<-sum(deg*log(deg),na.rm=TRUE)
  E<-sum(A)/2
  return(2*n+2-2*Q.DC(A,b)-4*foo3+4*E)}

BIC.DC <- function(A, b){
  n<- nrow(A)
  deg<-colSums(A); foo3<-sum(deg*log(deg),na.rm=TRUE)
  E<-sum(A)/2
  return((2*n+2)*log(n*(n+1)/2)-2*Q.DC(A,b)-4*foo3+4*E)}



############################################################
############## MLE initiated with given clusters ########## 
############################################################
DCmle <- function(A, k, n.iter,b){
N = nrow(A)
q.com = rep(NA, n.iter+1)			# vessel for values of complete data likelihood
b1 <- b
q.com[1] = Q.DC(A=A,b=b1)

for (iter in 1:n.iter){
b0 = b1
b.new = LabelSwitchDC(A,b0)

b1 = b.new
q.com[iter+1] = Q.DC(A,b1)
Q1 = q.com[iter+1]
  if (q.com[iter+1] <= q.com[iter]) {
    b.max = b1
    Q.max = Q1 
    break}
  }
list(q.last=Q.max,iter=iter,b.hat=b.max,lik.iter = q.com)}

############################################################
############## obtain MLE of probability matrix ########## 
############################################################
P_DC<-function(A,b){    
  N<-ncol(A)
  deg = rowSums(A)		# node degrees
  P = matrix(0,nrow=N,ncol=N)
  foo<-f.DC(A,b)
  O<-foo$O; d<-foo$d
  for (i in 1:N){		# populate diagonal entries
    r = b[i]
    P[i,i] = (1/2)*((deg[i])^2)*O[r,r]/((d[r])^2)}	# divide by 2 --- see notes
  for (i in 1:(N-1)){		# calculate P
    for (j in (i+1):N){
      r = b[i]; s = b[j]
      P[i,j] = deg[i]*deg[j]*O[r,s]/(d[r]*d[s])
      P[j,i] = P[i,j]
    }}
  return(P)}

######################

#### functions useful in simulations, to keep the simualtion code succinct
generate <- function(lambda, b){
  if (nrow(lambda) != length(b)) stop('dimension mismatch between lambda and b')
  N = nrow(lambda)
  k = ncol(lambda)
  b = as.factor(b)
  A = matrix(0, nrow = N, ncol = N)
  for (i in 1:(N-1)) {
    for (j in (i+1):N) {
      p = lambda[i,b[j]]*lambda[j,b[i]]
      A[i,j] = rbinom(n = 1, size = 1, prob = min(max(p,0),1))
      A[j,i] = A[i,j]}}
  return(A)}
#####################################################################
P_sim <- function(lambda, b){
  if (nrow(lambda) != length(b)) stop('dimension mismatch between lambda and b')
  N = nrow(lambda)
  k = ncol(lambda)
  b = as.factor(b)
  P = matrix(0, nrow = N, ncol = N)
  for (i in 1:(N-1)) {
    for (j in (i+1):N) {
      p = lambda[i,b[j]]*lambda[j,b[i]]
      P[i,j] = p
      P[j,i] = P[i,j]
    }}
  return(P)}
#####################################################################
P_sim_DC <- function(m,rho,h,b){
  N = length(b); x <- 2/(m+1)
  k = max(b)
  b = as.factor(b)
  P = matrix(0, nrow = N, ncol = N)
  omega <- matrix(1,nrow=k,ncol=k);diag(omega)<-h;omega<-rho*omega
  theta <- rbinom(n = N, size = 1, prob = 0.5)
  theta <- x + (m-1)*x*theta
  for (i in 1:(N-1)) {
    for (j in (i+1):N) {
      P[i,j] = theta[i]*theta[j]*omega[b[i],b[j]]
      P[j,i] = P[i,j]
    }}
  return(P)}
#####################################################################
generate_DC <- function(P){
  N <- nrow(P)
  A = matrix(0, nrow = N, ncol = N)
  for (i in 1:(N-1)) {
    for (j in (i+1):N) {
      p = min(P[i,j],1)
      A[i,j] = rbinom(n = 1, size = 1, prob = min(max(p,0),1))
      A[j,i] = A[i,j]}}
  return(A)}
