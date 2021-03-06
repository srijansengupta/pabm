rm(list = ls())
library("igraph")
library("clue")

##### Data input and pre-processing #####
foo=read.table("BritishMPcomm.txt",header=T)
mat=read.table("politicsuk-retweets.mtx",header=F)
mat[1,]  # no of nodes, and no of edges (all 5 communities)
mat = mat[-1,]
A = matrix(0,nrow(foo),nrow(foo))
for (k in 1:nrow(mat)){ #adjacency matrix with 5 parties
  edgei = mat[k,1]
  edgej = mat[k,2]
  i = which(foo[,1]==edgei)
  j = which(foo[,1]==edgej)
  A[i,j]<-1;A[j,i]<-1
}
b.true <- foo[,2]
levels(b.true) <- 1:5
b.true <- as.numeric(b.true)
table(b.true) # membership sizes of the 5 parties
nodes1 = which(b.true==1); nodes2 = which(b.true==2)
nodes = c(nodes1, nodes2)
A <- A[nodes,nodes] # Select only the two big communities
G<-graph.adjacency(A, mode="undirected")
is.connected(G, mode="strong") # check if the 2-comm network is connected: nope
foo1 <- clusters(G, mode="strong") 
nodes2 <- which(foo1$membership==1) # largest connected component (basically 31 nodes are disconnected)
A <- A[nodes2,nodes2] # form network with 329 nodes
G<-graph.adjacency(A, mode="undirected") 
is.connected(G, mode="strong") # now the network is connected
b.true <- b.true[nodes2]
id.BritMP <- foo[nodes2,]
write.table(id.BritMP,"id-BritMP.csv",sep=",")  # csv file with id and affiliation for each MP
nrow(A); sum(A)/2; max(b.true) # nodes, edges, number of communities
rm(list=setdiff(ls(), c("A","b.true"))) # delete unncessary objects
 
source("EPalgosim.R") 
source("functions.R")
##### get extreme points (with and without regularization) #####
Eps <- (0:19)/20
b.can <- NULL
for (ieps in 1:length(Eps)){
  b.can<-c(b.can,list(EPalgo(A,eps=Eps[ieps]))) # eps=0 for unperturbed version
  print(ieps)
}

##### estimate communities using PABM and DCBM modularity #####
##### 1. UNPERTURBED #####
b.can.unpert <- b.can[[1]]     # matrix of candidate assignments for eps = 0
n.can.unpert <- dim(b.can.unpert)[2]  # total number of candidate assignments
Q.PA.unpert = rep(NA, n.can.unpert)  # array to store modularity values
Q.DC.unpert = rep(NA, n.can.unpert)  # array to store modularity values

##### modularity calculation #####
ptm<-proc.time()
for (ican in 1:n.can.unpert){
  Q.PA.unpert[ican] = Q.PA(A, b=b.can.unpert[,ican])
  Q.DC.unpert[ican] = Q.DC(A, b=b.can.unpert[,ican])
  if(ican%%50==0) print(ican)}          # report every 50 passes
proc.time()-ptm
index.PA<-which(Q.PA.unpert==max(Q.PA.unpert))
index.DC<-which(Q.DC.unpert==max(Q.DC.unpert))
b.PA.unpert = b.can.unpert[,index.PA]
b.DC.unpert = b.can.unpert[,index.DC]

##### Community detection error (unperturbed) #####
x = as.cl_hard_partition(b.true) # true community assignments
y.PA = as.cl_hard_partition(b.PA.unpert)
y.DC = as.cl_hard_partition(b.DC.unpert)
1 - round(cl_agreement(x,y.PA,method="manhattan"),4) #PABM comm detection error
1 - round(cl_agreement(x,y.DC,method="manhattan"),4) #DCBM comm detection error

##### number of misclustered nodes #####
b.PA <- b.PA.unpert;b.DC <- b.DC.unpert
if (table(b.true,b.PA)[1,1]<table(b.true,b.PA)[1,2]) {   # check if label rotation is needed
  b.PA <- ifelse(b.PA==1,2,1)}; table(b.true,b.PA) # label rotation, if needed
if (table(b.true,b.DC)[1,1]<table(b.true,b.DC)[1,2]) {   # check if label rotation is needed
  b.DC <- ifelse(b.DC==1,2,1)}; table(b.true,b.DC) # label rotation, if needed

##### Model fitting error (unperturbed) #####
P.PA.hat = P_PA(A,b.PA)     # prob matrix fitted by PABM by comm detection
P.PA.true = P_PA(A,b.true)  # prob matrix fitted by PABM from true assignment
mu.PA.hat <- mu.hat(P.PA.hat,b.PA)    # popularity matrix fitted by PABM by comm detection
mu.PA.true <- mu.hat(P.PA.true,b.true) # popularity matrix fitted by PABM by comm detection
P.DC.hat = P_DC(A,b.DC)     # prob matrix fitted by DCBM by comm detection
P.DC.true = P_DC(A,b.true)  # prob matrix fitted by DCBM from true assignment
mu.DC.hat <- mu.hat(P.DC.hat,b.DC)    # popularity matrix fitted by DCBM by comm detection
mu.DC.true <- mu.hat(P.DC.true,b.true) # popularity matrix fitted by DCBM by comm detection

M.true <- f.PA(A, b.true)$M
E <- sum(A)/2       # no. of edges
F1.PA <- sum((mu.PA.hat - M.true)^2)/(2*E)
F2.PA <- sum((mu.PA.true - M.true)^2)/(2*E)
F1.DC <- sum((mu.DC.hat - M.true)^2)/(2*E)
F2.DC <- sum((mu.DC.true - M.true)^2)/(2*E)

F1.PA;F1.DC # model fitting error
F2.PA;F2.DC # model fitting error


##### 2. PERTURBED #####
# combine candidates for various epsilons
b.can.pert <- b.can[[2]]     # matrix of combined candidate assignments
for (ieps in 3:length(Eps)){
  b.can.pert <- cbind(b.can.pert, b.can[[ieps]])}
# identify duplicates and remove them
dup <- duplicated(b.can.pert,MARGIN=2)*1
can <- which(dup==0)
b.can.pert <- b.can.pert[,can]
n.can.pert <- dim(b.can.pert)[2]  # total number of candidate assignments
Q.PA.pert = rep(NA, n.can.pert)  # array to store PABM modularity values
Q.DC.pert = rep(NA, n.can.pert)  # array to store DCBM modularity values

##### modularity calculation #####
ptm<-proc.time()
for (ican in 1:n.can.pert){
  Q.PA.pert[ican] = Q.PA(A, b=b.can.pert[,ican])
  Q.DC.pert[ican] = Q.DC(A, b=b.can.pert[,ican])
  if(ican%%50==0) print(ican)}          # report every 50 passes
proc.time()-ptm
index.PA<-which(Q.PA.pert==max(Q.PA.pert))[1] # there is a tie
index.DC<-which(Q.DC.pert==max(Q.DC.pert))[1] # there is a tie
b.PA.pert = b.can.pert[,index.PA]
b.DC.pert = b.can.pert[,index.DC]

##### Community detection error (unperturbed) #####
x = as.cl_hard_partition(b.true) # true community assignments
y.PA = as.cl_hard_partition(b.PA.pert)
y.DC = as.cl_hard_partition(b.DC.pert)
1 - round(cl_agreement(x,y.PA,method="manhattan"),4) #PABM comm detection error
1 - round(cl_agreement(x,y.DC,method="manhattan"),4) #DCBM comm detection error

##### number of misclustered nodes #####
b.PA <- b.PA.pert;b.DC <- b.DC.pert
if (table(b.true,b.PA)[1,1]<table(b.true,b.PA)[1,2]) {   # check if label rotation is needed
  b.PA <- ifelse(b.PA==1,2,1)}; table(b.true,b.PA) # label rotation, if needed
if (table(b.true,b.DC)[1,1]<table(b.true,b.DC)[1,2]) {   # check if label rotation is needed
  b.DC <- ifelse(b.DC==1,2,1)}; table(b.true,b.DC) # label rotation, if needed

##### Model fitting error (unperturbed) #####
P.PA.hat = P_PA(A,b.PA)     # prob matrix fitted by PABM by comm detection
P.PA.true = P_PA(A,b.true)  # prob matrix fitted by PABM from true assignment
mu.PA.hat <- mu.hat(P.PA.hat,b.PA)    # popularity matrix fitted by PABM by comm detection
mu.PA.true <- mu.hat(P.PA.true,b.true) # popularity matrix fitted by PABM by comm detection
P.DC.hat = P_DC(A,b.DC)     # prob matrix fitted by DCBM by comm detection
P.DC.true = P_DC(A,b.true)  # prob matrix fitted by DCBM from true assignment
mu.DC.hat <- mu.hat(P.DC.hat,b.DC)    # popularity matrix fitted by DCBM by comm detection
mu.DC.true <- mu.hat(P.DC.true,b.true) # popularity matrix fitted by DCBM by comm detection

M.true <- f.PA(A, b.true)$M
E <- sum(A)/2       # no. of edges
F1.PA <- sum((mu.PA.hat - M.true)^2)/(2*E)
F2.PA <- sum((mu.PA.true - M.true)^2)/(2*E)
F1.DC <- sum((mu.DC.hat - M.true)^2)/(2*E)
F2.DC <- sum((mu.DC.true - M.true)^2)/(2*E)

F1.PA;F1.DC # model fitting error
F2.PA;F2.DC # model fitting error