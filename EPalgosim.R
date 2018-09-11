EPalgo<-function(A,eps=0){

##### perturbed adj matrix ##### (Le pg 15, Amini 2013)
tau = eps*(mean(colSums(A))/nrow(A))
A = A + tau

foo<-eigen(A, symmetric = TRUE)
val = abs(foo$values)			# pick the top 2
id = order(-val)				# eigenvalues of A and put
id_vec = id[1:2]				# their eigenvectors into a 2*N matrix
# columns of foo$vectors = eigenvectors of A
# we want a 2xn matrix whose rows are the leading eigenvectors
X = t(foo$vectors[,id_vec])		
y = X[,1:2]

comms = 1:2
u <- list(comms)
v = expand.grid(rep(u, 2))
v = as.matrix(v)
# initialize with the parallelogram
epts = y%*%t(v)	# extreme pts are the columns of this matrix
b.can = t(v)	# candidate configurations.
row.names(b.can)=NULL

ptm<-proc.time()
for (i in 3:ncol(X)){
b.can1 = rbind(b.can,rep(1,ncol(b.can)))
b.can2 = rbind(b.can,rep(2,ncol(b.can)))
b.can = cbind(b.can1,b.can2)
foo = X[,1:i]%*%b.can
hull = chull(t(foo))
epts = foo[,hull]
b.can = b.can[,hull]}	# next i = next row of X
proc.time()-ptm

##### remove invalid candidates
k = max(b.can)
foo = b.can
foo1 = NA
for (i in 1:ncol(b.can)){
foo2 = rep(NA,k)
for (clus in 1:k){foo2[clus]=sum(b.can[,i]==clus)}
if (min(foo2)==0){foo1 = c(foo1,i)}
}
if (length(foo1)>1) {foo1 = foo1[-1]
b.can = b.can[,-foo1]}

###### remove eqv candidates
foo1 = NA
for (i in 2:ncol(b.can)){
for (j in 1:i){
foo4 = abs(b.can[,i] - b.can[,j])
if (mean(foo4) == 1){ # this means b.can[,i] and b.can[,j] are exactly 1 apart
foo1 = c(foo1,i)
break}
}}
if (length(foo1)>1){
  foo1 = foo1[-1]
  b.can = b.can[,-foo1]
}
return(b.can)}	# end of function
