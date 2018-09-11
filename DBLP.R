rm(list = ls())
library("igraph")
library("clue")

##### Data input and pre-processing #####
a = read.delim("author_label.txt", header = FALSE, col.names = c("a.id","a.lab","a.name"))
a.list = unique(a$a.id)
a <- a[,1:2]       # 0: Database; 1: Data Mining; 2: AI; 3: Information Retrieval
pa = read.delim("paper_author.txt", header = FALSE, col.names = c("p.id","a.id")) # paper author
pa = subset(pa, a.id %in% a.list)
pc = read.delim("paper_conf.txt", header = F, col.names = c("p.id","c.id"))  # paper at what conference?
conf = rep(NA, dim(pa)[1])
for (i in 1:dim(pa)[1]){
  id = pa[i,]$p.id
  conf[i] = subset(pc, p.id == id)$c.id}
pac = cbind(pa, conf) # paper - author -conference
n.a = length(unique(pa$a.id)); n.p <- length(unique(pa$p.id)); n.c = length(unique(pac$conf)) 
papers <- NULL # list whose ith element is a list of papers written by ith author
confs <- NULL # list whose ith element is a list of conferences where ith author presented
comm <- rep(NA, n.a)
for (i in 1:n.a){
  foo <- a[i,1]
  comm[i] <- subset(a, a.id == foo)$a.lab
  papers <- c(papers,list(subset(pa, a.id == foo)$p.id))
  confs <- c(confs,list(subset(pac, a.id == foo)$conf))}
A <- matrix(0, nrow = n.a, ncol = n.a)
for (i in 1:(n.a-1)){
  if (length(confs[[i]])>0) { # if ith author has presented at at least one conference
    confsi <- confs[[i]] # conferences at which ith author presented
    for (j in (i+1):n.a){
      if (length(confs[[j]])>0) { # if jth author has presented at at least one conference
        confsj <- confs[[j]] # conferences at which jth author presented
        if (length(intersect(confsi,confsj))>0) { # common conference
          A[i,j] = 1; A[j,i] = 1
        }}}}}
areas = c(0,3) #database and IR
foo <- which(comm %in% areas)
A <- A[foo,foo]
b.true <- comm[foo]
b.true <- ifelse(b.true==areas[1],1,2)
G<-graph.adjacency(A, mode="undirected")
is.connected(G, mode="strong") # check if the 2-comm network is connected: yes!
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


