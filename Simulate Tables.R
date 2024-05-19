library(mvtnorm)
library(MASS)
library(matrixStats)
library("FactoMineR")
library("cluster")
library("dbscan")
library("flashClust")
library(sn)
rm(list = ls()) 
options(digits=5)

########### Functions to calculate the high-dimensional dissimilarity delta #########################
delta <- function(both) { 
  myN=nrow(both)
  Newzbar=matrix(0,nrow=myN,ncol=myN)
  con=1/(myN-2)
  for (ix in 1:(myN-1))
  { xx=both[ix,];        mux=mean(xx);   vx=sqrt(var(xx))
  for (iy in (ix+1):myN)
  {   yy=both[iy,];      muy=mean(yy);   vy=sqrt(var(yy))
  part=0
 # cat('mux=',mux,'muy=',muy,'vx=',vx,'vy=',vy,fill=TRUE)
  for (it in 1:myN)  if ((it!=ix) & (it!=iy)) {   
    zz=both[it,]
    muz=mean(zz)
    vz=sqrt(var(zz))
    dxz=sqrt(((mux-muz)^2+(vx-vz)^2)) 
    dyz=sqrt(((muy-muz)^2+(vy-vz)^2))
    part=part+abs(dxz-dyz)
  }
  Newzbar[ix,iy]<- con*part
  Newzbar[iy,ix]<-Newzbar[ix,iy]
  }
  }
  return(Newzbar)  
}
########### Functions to calculate the high-dimensional dissimilarity Rho #########################
Rho0 <- function(ipdmat) {
  myN=nrow(ipdmat)
  R0=matrix(0,nrow=myN,ncol=myN)
  con=1/((myN-2)*sqrt(px))
  for (ix in 1:(myN-1))
  { 
    for (iy in (ix+1):myN)
    {   
      part1=ipdmat[ix,]
      p1=part1[-ix]
      iy1=iy-1
      Q1=p1[-iy1]
      part2=ipdmat[iy,]
      p2=part2[-ix]
      Q2=p2[-iy1]
      tt=mean(abs(Q1-Q2))
      R0[ix,iy]=con*mean(abs(Q1-Q2))
      R0[iy,ix]= R0[ix,iy]    
    }
  }
  return(R0)  
}

############# HDCPC clustering ###########################
HDCPC_voting <- function(mydata,k,index= c("delta","rho"))
{
  if(index== "delta")
  {dis_matrix =delta(mydata)}
  if(index=="rho")
  {ipdmat <- as.matrix(dist(mydata,method = "euclidean"))
  dis_matrix=Rho0(ipdmat)}
  
  num_obs = nrow(mydata)
  save = dis_matrix
  save_cluster = list()
  agreement = matrix(0,nrow = num_obs,ncol = num_obs)
  
  for(i in rep(1:num_obs))
  { 
    count_index = rep(1:num_obs)
    cluster_index = rep(k,num_obs)
    focus = dis_matrix[i, ]
    sorted_focus = sort(focus)
    sorted_indices <- order(focus)
    indices = count_index[sorted_indices]
    sorted_withindices = as.matrix(rbind(focus,count_index))[,sorted_indices]
    indent = c(sorted_focus[-1],0)
    differencesR = (indent-sorted_focus)[-length((indent-sorted_focus))]
    cutoff = sort(order(differencesR,decreasing = T)[1:(k-1)])
    loc=1
    ijk=1
    for(jj in cutoff)
    {
      cluster_index[indices[loc:jj]]=ijk
      ijk=ijk+1
      loc = jj+1
    }
    save_cluster[[i]]=cluster_index
  }
  
  for (row_i in rep(1:(num_obs-1)))
  {
    for (col_j in rep((row_i+1):num_obs))
    {
      agreement[col_j,row_i]=agreement[row_i,col_j]=randindex(save_cluster[[row_i]],save_cluster[[col_j]])
    }
  }
  return(save_cluster[[which.min(rowSums(agreement))]])
}


########### Minimal Spanning Tree Clustering #################################
msTreePrim <- function(nodes, arcs, start.node = 1) {
  
  # Duplicate the arcs
  arcs <- rbind(arcs, matrix(c(arcs[, 2], arcs[, 1], arcs[, 3]), ncol = 3))
  
  # Initialize with empty tree and start.node
  tree.arcs <- matrix(ncol = 3)[-1, ]
  tree.nodes <- nodes[nodes == start.node]
  
  stages <- 0  # initialize counter
  stages.arcs <- c()  # vector to store stage number in wich each arc was added
  
  # Iterate until every node was checked
  while (length(tree.nodes) < length(nodes)) {
    
    # Arcs leaving tree nodes and reaching unadded nodes
    k <- which(arcs[, 1] %in% tree.nodes &
                 arcs[, 2] %in% nodes[-which(nodes %in% tree.nodes)])
    validArcs <- matrix(arcs[k, ], ncol = 3)
    
    # Valid arcs with minimum cost
    l <- which(validArcs[, 3] == min(validArcs[, 3]))
    min.arc <- matrix(validArcs[l, ], ncol = 3)
    
    # Save one arc in minimum tree
    tree.arcs <- rbind(tree.arcs, min.arc[1, ])
    
    # Update tree nodes
    tree.nodes <- c(tree.nodes, min.arc[1, 2])
    
    stages <- stages + 1  # counter
    # Save in which stage an arc was added to the tree and update
    stages.arcs <- c(stages.arcs,
                     rep(stages, nrow(tree.arcs) - length(stages.arcs)))
    
  }
  colnames(tree.arcs) <- c("ept1", "ept2", "weight")
  output <- list("tree.nodes" = tree.nodes, "tree.arcs" = tree.arcs,
                 "stages" = stages, "stages.arcs" = stages.arcs)
  return(output)
}
arcmatrix <- function(data,diss){
  N=nrow(data)
  arc = c()
  L0 = diss(data)
  for (i in rep(1:(N-1)))
  {
    for (j in rep((i+1):N))
    {
      arc = rbind(arc,c(i,j,L0[i,j]))
    }
  }
  return(arc)
}
arcmatrix_rho <- function(data,diss=Rho0){
  N=nrow(data)
  arc = c()
  ipdmat <- as.matrix(dist(data,method = "euclidean"))
  L0 = diss(ipdmat)
  for (i in rep(1:(N-1)))
  {
    for (j in rep((i+1):N))
    {
      arc = rbind(arc,c(i,j,L0[i,j]))
    }
  }
  return(arc)
}
arcmatrix_delta <- function(data,diss=delta){
  N=nrow(data)
  arc = c()
  L0 = diss(data)
  for (i in rep(1:(N-1)))
  {
    for (j in rep((i+1):N))
    {
      arc = rbind(arc,c(i,j,L0[i,j]))
    }
  }
  return(arc)
}
arcmatrix_dist <- function(data){
  N=nrow(data)
  arc = c()
  L0 = as.matrix(dist(data,diag = 1,upper = 1))
  for (i in rep(1:(N-1)))
  {
    for (j in rep((i+1):N))
    {
      arc = rbind(arc,c(i,j,L0[i,j]))
    }
  }
  return(arc)
}
mstree_cluster_rho <- function(data,dis=Rho0,k){
  N=nrow(data)
  min_tree = msTreePrim(rep(1:N),arcmatrix_rho(data,dis))
  cutoff = order(min_tree$tree.arcs[,3],decreasing=TRUE)[1:(k-1)]
  cutoff=sort(cutoff, decreasing = TRUE)
  min_tree$tree.nodes
  cluster_tag = rep(k,N)
  for (ks in rep(1:(k-1)))
  {
    kks=k-ks
    cluster_tag[min_tree$tree.nodes[1:cutoff[ks]]]=kks
  }
  return(cluster_tag)
}
mstree_cluster_delta <- function(data,dis=delta,k){
  N=nrow(data)
  min_tree = msTreePrim(rep(1:N),arcmatrix_delta(data,dis=delta))
  cutoff = order(min_tree$tree.arcs[,3],decreasing=TRUE)[1:(k-1)]
  cutoff=sort(cutoff, decreasing = TRUE)
  min_tree$tree.nodes
  cluster_tag = rep(k,N)
  for (ks in rep(1:(k-1)))
  {
    kks=k-ks
    cluster_tag[min_tree$tree.nodes[1:cutoff[ks]]]=kks
  }
  return(cluster_tag)
}
mstree_dist_cluster <- function(data,k){
  N=nrow(data)
  aa=arcmatrix_dist(data)
  min_tree = msTreePrim(rep(1:N),aa)
  cutoff = order(min_tree$tree.arcs[,3],decreasing=TRUE)[1:(k-1)]
  cutoff=sort(cutoff, decreasing = TRUE)
  min_tree$tree.nodes
  cluster_tag = rep(k,N)
  for (ks in rep(1:(k-1)))
  {
    kks=k-ks
    cluster_tag[min_tree$tree.nodes[1:cutoff[ks]]]=kks
  }
  return(cluster_tag)
}

############## Rand Index ###################################
randindex <- function(a,b)
{
  sumsum = 0
  n=dim(as.data.frame(a))[1]
  coef = n*(n-1)/2
  for (i in rep(1:(n-1)))
  {
    for (j in rep((i+1):n))
    {
      first = as.numeric(a[i]==a[j])
      second = as.numeric(b[i]==b[j])
      sum = as.numeric((first+second)==1)
      sumsum=sumsum+sum
    }
  }
  return(sumsum/coef)
}

############## k-Means Clustering ###############################
Euclideanmeanscluster<- function(mydata,k,max=10){
  #max is the maximum number of groups returned
  max_conve = 100
  iteration_time = 0
  #num_obs: number of observations
  num_obs = nrow(mydata)
  dis_matrix = as.matrix(dist(mydata,diag=T,upper=T))
  #initialize tags
  center  = c()
  first_center = sample(rep(1:num_obs),1)
  second_center = which.max(dis_matrix[first_center,])
  center= c(first_center,second_center)
  for( jjj in rep(1:(k-2)))
  {
    dis_matrix[center,]
    smallest_values <- apply(dis_matrix[center,], 2, min)
    center = c(center,which.max(smallest_values))
  }
  cluster_tag = rep(1:num_obs)
  for (jjjj in rep(1:num_obs)[-center])
  {
    assignjjjj = center[which.min(dis_matrix[jjjj,center])]
    cluster_tag[jjjj] = assignjjjj
  }
  for (mmax in rep(1:max_conve))
  {
    #num_groups: number of cluster groups
    num_groups = nrow(as.data.frame(unique(cluster_tag)))
    iteration_time = iteration_time + 1
    converge_save  = cluster_tag
    update_tag = rep(0,num_obs)
    for (i in rep(1:num_obs))
    {
      dis_eachgroup = rep(0,num_groups)   
      count = 0
      for (group in unique(cluster_tag))
      {
        count= count+1
        dis_eachgroup[count] = sum(dis_matrix[which(cluster_tag==group),i])
        save = as.data.frame(table(cluster_tag))
        colnames(save) = c("x","freq")
        save[save$x==cluster_tag[i],]$freq = save[save$x==cluster_tag[i],]$freq-1
        dis_eachgroup[count] = (sum(dis_matrix[which(cluster_tag==group),i]))/(save[save$x==group,]$freq)
      }
      update_tag[i] =  unique(cluster_tag)[which.min(dis_eachgroup)]
      
      
    }
    
    for (jk in rep(1:nrow(as.matrix((update_tag)))))
    { update_tag[jk] = update_tag[update_tag[jk]]}
    if(identical(update_tag,converge_save)||(nrow(as.data.frame(unique(update_tag)))<=k))
    {
      break;
    }
    cluster_tag = update_tag
    converge_save = update_tag
    
  }
  num_groups = nrow(as.data.frame(unique(cluster_tag)))
  while(num_groups > max)
  {
    num_groups = nrow(as.data.frame(unique(cluster_tag)))
    save_forpairs = matrix(99, nrow = num_groups,ncol = num_groups)
    for (ci in rep(1:(num_groups-1)))
    {
      for (bi in rep((ci+1):num_groups))
      {
        sum_avgdis = 0
        group1 = (unique(cluster_tag))[ci]
        group2 = (unique(cluster_tag))[bi]
        whichgroup1 = which(cluster_tag==group1)
        whichgroup2 = which(cluster_tag==group2)
        howmanygroup1 = dim(as.data.frame(which(cluster_tag==group1)))[1]
        howmanygroup2 = dim(as.data.frame(which(cluster_tag==group2)))[1]
        coeff = 1/(howmanygroup1*howmanygroup2)
        for (di in whichgroup1)
        {
          for (ei in whichgroup2)
          {
            sum_avgdis = sum_avgdis + dis_matrix[di,ei]
          }
        }
        save_forpairs[ci,bi] = sum_avgdis*coeff
      }
    }
    combine = which(save_forpairs == min(save_forpairs),arr.ind = T)
    whichcombine_group1 = which(cluster_tag==((unique(cluster_tag))[combine[1]]))
    cluster_tag[whichcombine_group1] = ((unique(cluster_tag))[combine[2]])
    num_groups = nrow(as.data.frame(unique(cluster_tag)))
  }
  return(cluster_tag)
}

deltameanscluster<- function(mydata,k,max=100){
  #max is the maximum number of groups returned
  iteration_time = 0
  max_conve = 100
  num_obs = nrow(mydata)
  dis_matrix =delta(mydata)
  center  = c()
  first_center = sample(rep(1:num_obs),1)
  second_center = which.max(dis_matrix[first_center,])
  center= c(first_center,second_center)
  for( jjj in rep(1:(k-2)))
  {
    dis_matrix[center,]
    smallest_values <- apply(dis_matrix[center,], 2, min)
    center = c(center,which.max(smallest_values))
  }
  cluster_tag = rep(1:num_obs)
  for (jjjj in rep(1:num_obs)[-center])
  {
    assignjjjj = center[which.min(dis_matrix[jjjj,center])]
    cluster_tag[jjjj] = assignjjjj
  }
  
  for (mmax in rep(1:max_conve))
  {
    #num_groups: number of cluster groups
    num_groups = nrow(as.data.frame(unique(cluster_tag)))
    iteration_time = iteration_time + 1
    converge_save  = cluster_tag
    update_tag = rep(0,num_obs)
    for (i in rep(1:num_obs))
    {
      dis_eachgroup = rep(0,num_groups)   
      count = 0
      for (group in unique(cluster_tag))
      {
        count= count+1
        dis_eachgroup[count] = sum(dis_matrix[which(cluster_tag==group),i])
        save = as.data.frame(table(cluster_tag))
        colnames(save) = c("x","freq")
        save[save$x==cluster_tag[i],]$freq = save[save$x==cluster_tag[i],]$freq-1
        dis_eachgroup[count] = (sum(dis_matrix[which(cluster_tag==group),i]))/(save[save$x==group,]$freq)
      }
      update_tag[i] =  unique(cluster_tag)[which.min(dis_eachgroup)]
      
      
    }
    
    for (jk in rep(1:nrow(as.matrix((update_tag)))))
    { update_tag[jk] = update_tag[update_tag[jk]]}
    if(identical(update_tag,converge_save)||(nrow(as.data.frame(unique(update_tag)))<=k))
    {
      break;
    }
    cluster_tag = update_tag
    converge_save = update_tag
    
  }
  num_groups = nrow(as.data.frame(unique(cluster_tag)))
  while(num_groups > max)
  {
    num_groups = nrow(as.data.frame(unique(cluster_tag)))
    # num_group_pairs = num_groups*(num_groups-1)/2
    save_forpairs = matrix(99, nrow = num_groups,ncol = num_groups)
    for (ci in rep(1:(num_groups-1)))
    {
      for (bi in rep((ci+1):num_groups))
      {
        sum_avgdis = 0
        group1 = (unique(cluster_tag))[ci]
        group2 = (unique(cluster_tag))[bi]
        whichgroup1 = which(cluster_tag==group1)
        whichgroup2 = which(cluster_tag==group2)
        howmanygroup1 = dim(as.data.frame(which(cluster_tag==group1)))[1]
        howmanygroup2 = dim(as.data.frame(which(cluster_tag==group2)))[1]
        coeff = 1/(howmanygroup1*howmanygroup2)
        for (di in whichgroup1)
        {
          for (ei in whichgroup2)
          {
            sum_avgdis = sum_avgdis + dis_matrix[di,ei]
          }
        }
        save_forpairs[ci,bi] = sum_avgdis*coeff
      }
    }
    combine = which(save_forpairs == min(save_forpairs),arr.ind = T)
    whichcombine_group1 = which(cluster_tag==((unique(cluster_tag))[combine[1]]))
    cluster_tag[whichcombine_group1] = ((unique(cluster_tag))[combine[2]])
    num_groups = nrow(as.data.frame(unique(cluster_tag)))
 
  }
  return(cluster_tag)
}
rhomeanscluster <- function(mydata,k,max=100){
  #max is the maximum number of groups returned
  iteration_time = 0
  #num_obs: number of observations
  num_obs = nrow(mydata)
  ipdmat <- as.matrix(dist(mydata,method = "euclidean"))
  dis_matrix = Rho0(ipdmat)
  max_conve = 100
  #initialize tags
  center  = c()
  first_center = sample(rep(1:num_obs),1)
  second_center = which.max(dis_matrix[first_center,])
  center= c(first_center,second_center)
  for( jjj in rep(1:(k-2)))
  {
    dis_matrix[center,]
    smallest_values <- apply(dis_matrix[center,], 2, min)
    center = c(center,which.max(smallest_values))
  }
  cluster_tag = rep(1:num_obs)
  for (jjjj in rep(1:num_obs)[-center])
  {
    assignjjjj = center[which.min(dis_matrix[jjjj,center])]
    cluster_tag[jjjj] = assignjjjj
  }
  for (mmax in rep(1:max_conve))
  {
    #num_groups: number of cluster groups
    num_groups = nrow(as.data.frame(unique(cluster_tag)))
    iteration_time = iteration_time + 1
    converge_save  = cluster_tag
    update_tag = rep(0,num_obs)
    for (i in rep(1:num_obs))
    {
      dis_eachgroup = rep(0,num_groups)   
      count = 0
      for (group in unique(cluster_tag))
      {
        count= count+1
        dis_eachgroup[count] = sum(dis_matrix[which(cluster_tag==group),i])
        save = as.data.frame(table(cluster_tag))
        colnames(save) = c("x","freq")
        save[save$x==cluster_tag[i],]$freq = save[save$x==cluster_tag[i],]$freq-1
        dis_eachgroup[count] = (sum(dis_matrix[which(cluster_tag==group),i]))/(save[save$x==group,]$freq)
      }
      update_tag[i] =  unique(cluster_tag)[which.min(dis_eachgroup)]
      
      
    }
    
    for (jk in rep(1:nrow(as.matrix((update_tag)))))
    { update_tag[jk] = update_tag[update_tag[jk]]}
    if(identical(update_tag,converge_save)||(nrow(as.data.frame(unique(update_tag)))<=k))
    {
      break;
    }
    cluster_tag = update_tag
    converge_save = update_tag
    
  }
  num_groups = nrow(as.data.frame(unique(cluster_tag)))
  while(num_groups > max)
  {
    num_groups = nrow(as.data.frame(unique(cluster_tag)))
    # num_group_pairs = num_groups*(num_groups-1)/2
    save_forpairs = matrix(99, nrow = num_groups,ncol = num_groups)
    for (ci in rep(1:(num_groups-1)))
    {
      for (bi in rep((ci+1):num_groups))
      {
        sum_avgdis = 0
        group1 = (unique(cluster_tag))[ci]
        group2 = (unique(cluster_tag))[bi]
        whichgroup1 = which(cluster_tag==group1)
        whichgroup2 = which(cluster_tag==group2)
        howmanygroup1 = dim(as.data.frame(which(cluster_tag==group1)))[1]
        howmanygroup2 = dim(as.data.frame(which(cluster_tag==group2)))[1]
        coeff = 1/(howmanygroup1*howmanygroup2)
        for (di in whichgroup1)
        {
          for (ei in whichgroup2)
          {
            sum_avgdis = sum_avgdis + dis_matrix[di,ei]
          }
        }
        save_forpairs[ci,bi] = sum_avgdis*coeff
      }
    }
    combine = which(save_forpairs == min(save_forpairs),arr.ind = T)
    whichcombine_group1 = which(cluster_tag==((unique(cluster_tag))[combine[1]]))
    cluster_tag[whichcombine_group1] = ((unique(cluster_tag))[combine[2]])
    num_groups = nrow(as.data.frame(unique(cluster_tag)))
 
  }
  return(cluster_tag)
}
Euclideanmeanscluster_update<- function(mydata,k){
  save=as.data.frame(kmeans(mydata,centers=k)[1])
  return(t(save))
}

###### Estimate of Number of Clusters #############################
KL = function(mydata,dis="Euclidean",Kmax=10,base_fun){
  ipdmatt <- as.matrix(dist(mydata,method = "euclidean"))
  dis.list=c("Euclidean","rho","delta")
  loc = which(!is.na(match(dis.list,dis)))
  if(loc==2)
  {ipdmatt = Rho0(ipdmatt)}
  if(loc==3)
  {ipdmatt=delta(mydata)}
  num_obs = nrow(mydata)
  d=ncol(mydata)
  W1=sum(ipdmatt)/(2*num_obs)
  save_KL = rep(0,Kmax) 
  W_save = c(W1)
  for(i in rep(2:Kmax))
  {
    if(identical(base_fun,HDCPC_voting) )
    {cluster = base_fun(mydata,k=i,index=dis)}
    else{    cluster = base_fun(mydata,k=i)}
    
    W_next_save = 0
    for( j in unique(cluster))
    {
      Cj = sum(cluster==j)
      ipdmat_j <- ipdmatt[cluster==j,cluster==j]
      W_next_save=W_next_save+sum(ipdmat_j)/(2*Cj)
    }
    W_save=c(W_save,W_next_save)
  }
  Diff_k=rep(0,Kmax)
  W_save_update =rep(0,Kmax)
  for (ii in rep(1:Kmax))
  {
    W_save_update[ii] = ii^(2/d)*W_save[ii]
  }
  for(jj in rep(2:Kmax))
  {
    Diff_k[jj] = W_save_update[jj-1]-W_save_update[jj]
  }
  KL = rep(0,(Kmax-1))
  for( zz in rep(2:(Kmax-1)))
  {
    KL[zz] = abs(Diff_k[zz]/Diff_k[zz+1])
  }
  KL[KL==Inf]=0
  return(which.max(KL))
}



##################################  Main  ################################
set.seed(41253)
rhox=0
n1=n2=n3=n4=n5=8
N=CAPN=n1+n2+n3+n4+n5
sim=500

labval=rep(0,N)
pxlist= c(4,8,16,32,64,128,256,512,1024,2048,4096,8192)
performance=matrix(0,nrow = length(pxlist),ncol=10)
rownames(performance) = pxlist
colnames(performance)=c("New-R","New D","KRho","KDelta","kmns","hdbs","Eucl","msT_E","msT_rho","msT_delta")

for(Idis in 5:5) {
    for (pxx in rep(1:length(pxlist))) {
       px=pxlist[pxx]
       mu1=rep(1,px)
       mu2=rep(1.5,px)
       mu3=rep(2,px)
       mu4=rep(2.5,px)
       mu5=rep(3,px)
       sig1=matrix(rhox,nrow=px, ncol=px);  diag(sig1)=1
       sig2=matrix(rhox,nrow=px, ncol=px);  diag(sig2)=1
       sig3=matrix(rhox,nrow=px, ncol=px);  diag(sig3)=1
       sig4=matrix(rhox,nrow=px, ncol=px);  diag(sig4)=1
       sig5=matrix(rhox,nrow=px, ncol=px);  diag(sig5)=1
       for (IV in 1:sim) {
         # normal
           cat('Idis',Idis,'=px=',px,'sim=',IV,fill=TRUE)
         if (Idis==1) {
         samplex1= rmvnorm(n1, mean = mu1,  sig1) 
         samplex2= rmvnorm(n2, mean = mu2, sig2) 
         samplex3= rmvnorm(n3, mean = mu3, sig3) 
         samplex4= rmvnorm(n4, mean = mu4, sig4) 
         samplex5= rmvnorm(n5, mean = mu5, sig5) 
         comsam=rbind(samplex1, samplex2, samplex3, samplex4 ,samplex5)
         }
         if (Idis==2) {
           # abs normal
           samplex1= abs(rmvnorm(n1, mean = mu1,  sig1)) 
           samplex2= abs(rmvnorm(n2, mean = mu2, sig2))
           samplex3= abs(rmvnorm(n3, mean = mu3, sig3)) 
           samplex4= abs(rmvnorm(n4, mean = mu4, sig4)) 
           samplex5= abs(rmvnorm(n5, mean = mu5, sig5)) 
           comsam=rbind(samplex1, samplex2, samplex3, samplex4 ,samplex5)
         }
       else if(Idis==3){
         #Bernoulli  
         samplex1=matrix(rbinom(n1*px, 1, .1),nrow=n1, ncol=px)
         samplex2=matrix(rbinom(n2*px, 1, .25),nrow=n2, ncol=px)
         samplex3=matrix(rbinom(n3*px, 1, .50),nrow=n2, ncol=px)
         samplex4=matrix(rbinom(n4*px, 1, .75),nrow=n2, ncol=px)
         samplex5=matrix(rbinom(n5*px, 1, .90),nrow=n2, ncol=px)
         comsam=rbind(samplex1,samplex2,samplex3,samplex4,samplex5)
        }
        else if(Idis==4) {
          # Skew normal   
          samplex1=(rmsn(n1, xi=rep(1, px), Omega=sig1, alpha=rep(1, px)))
          samplex2=(rmsn(n2, xi=rep(1.5, px), Omega=sig2, alpha=rep(1, px)))
          samplex3=(rmsn(n3, xi=rep(2, px), Omega=sig3, alpha=rep(1, px)))
          samplex4=(rmsn(n4, xi=rep(2.5, px), Omega=sig4, alpha=rep(1, px)))
          samplex5=(rmsn(n5, xi=rep(3, px), Omega=sig5, alpha=rep(1, px)))
          comsam=rbind(samplex1, samplex2, samplex3, samplex4 ,samplex5)
         }
         else if(Idis==5){
           # t distribution
           DFF=7
           samplex1=rmvt(n1, sigma=sig1, df=DFF)+1 
           samplex2=rmvt(n2, sigma=sig2, df=DFF)+2
           samplex3=rmvt(n3, sigma=sig3, df=DFF)+3
           samplex4=rmvt(n4, sigma=sig4, df=DFF)+4
           samplex5=rmvt(n5, sigma=sig5, df=DFF)+5
           comsam=rbind(samplex1,samplex2,samplex3,samplex4,samplex5 )
         }
         else if(Idis==6){
           # Chisq-3 distribution
           samplex1=matrix(rchisq(n1*px,3),nrow=n1, ncol=px)
           samplex2=matrix(rchisq(n1*px,3),nrow=n2, ncol=px)+.5
           samplex3=matrix(rchisq(n1*px,3),nrow=n2, ncol=px)+1
           samplex4=matrix(rchisq(n1*px,3),nrow=n2, ncol=px)+1.5
           samplex5=matrix(rchisq(n1*px,3),nrow=n2, ncol=px)+2
           comsam=rbind(samplex1,samplex2,samplex3,samplex4,samplex5)
         }
         # normal with scale
         else if (Idis==7) {
           sig1=matrix(rhox,nrow=px, ncol=px);  diag(sig1)=1
           sig2=matrix(rhox,nrow=px, ncol=px);  diag(sig2)=1.5
           sig3=matrix(rhox,nrow=px, ncol=px);  diag(sig3)=2
           sig4=matrix(rhox,nrow=px, ncol=px);  diag(sig4)=2.5
           sig5=matrix(rhox,nrow=px, ncol=px);  diag(sig5)=3
           samplex1= rmvnorm(n1, mean = mu1,  sig1) 
           samplex2= rmvnorm(n2, mean = mu2, sig2) 
           samplex3= rmvnorm(n3, mean = mu3, sig3) 
           samplex4= rmvnorm(n4, mean = mu4, sig4) 
           samplex5= rmvnorm(n5, mean = mu5, sig5) 
           comsam=rbind(samplex1, samplex2, samplex3, samplex4 ,samplex5)
         }
         obslab=c(rep(1,n1),rep(2,n2),rep(3,n3),rep(4,n4),rep(5,n5))
  # shuffle the sample 
    mixed_index= sample(1:N)
    comsam=comsam[mixed_index,]
    obslab=obslab[mixed_index]
    
    ##################################################New-R ###################

    k_keep =  KL(comsam,dis="rho",Kmax=10,HDCPC_voting)
    plab=HDCPC_voting(comsam,k=k_keep,index="rho")
    rdx=randindex (obslab,plab)
    performance[pxx,1]=performance[pxx,1]+rdx
    
    ##################################################New-D ###################
    k_keep = KL(comsam,dis="delta",Kmax=10,HDCPC_voting)
    plab=HDCPC_voting(comsam,k=k_keep,index="delta")
    rdx=randindex (obslab,plab)
    performance[pxx,2]=performance[pxx,2]+rdx
    
    #########################################  rhomeanscluster ###################
    ipdmat <- as.matrix(dist(comsam,method = "euclidean"))
    k_keep = KL(comsam,dis="rho",Kmax=10,rhomeanscluster)
    plab=rhomeanscluster(comsam,k=k_keep,max=10)
    rdx=randindex (obslab,plab)
    performance[pxx,3]=performance[pxx,3]+rdx
 
    ########################## deltameanscluster #################
    k_keep = KL(comsam,dis="delta",Kmax=10,deltameanscluster)
    plab= deltameanscluster(comsam,k=k_keep,max=10)
    rdx=randindex (obslab,plab)
    performance[pxx,4]=performance[pxx,4]+rdx
 
    ############# R kmeans  #############
    k_keep = KL(comsam,dis="Euclidean",Kmax=10,Euclideanmeanscluster_update)
    RKM=kmeans(comsam, centers=k_keep, iter.max = 20, nstart = 2,algorithm = "Hartigan-Wong", trace=FALSE)
    plab=RKM$cluster
    rdx=randindex (obslab,plab)
    performance[pxx,5]=performance[pxx,5]+rdx
 
    ################ hdbx package  ###########
    hdbx <- hdbscan(comsam, minPts = 5)
    plab=hdbx$cluster
    rdx=randindex (obslab,plab)
    performance[pxx,6]=performance[pxx,6]+rdx
 
    
    #####################   Euclideanmeanscluster   #############
    k_keep = KL(comsam,dis="Euclidean",Kmax=10,Euclideanmeanscluster_update)
    plab=Euclideanmeanscluster(comsam,k=k_keep,max=10)
    rdx=randindex (obslab,plab)
    performance[pxx,7]=performance[pxx,7]+rdx
 
    
    ##############  mstree_dist_cluster       ###########
    k_keep = KL(comsam,dis="Euclidean",Kmax=10,mstree_dist_cluster)
    plab=mstree_dist_cluster(comsam,k=k_keep)
    rdx=randindex(obslab,plab)
    performance[pxx,8] = performance[pxx,8]+rdx
 
    
    ##############   mstree_cluster_rho   ###########
    k_keep = KL(comsam,dis="rho",Kmax=10,mstree_cluster_rho)
    plab=mstree_cluster_rho(comsam,dis=Rho0,k=k_keep)
    rdx=randindex(obslab,plab)
    performance[pxx,9] = performance[pxx,9]+rdx
    
    
    ######################  mstree_cluster_delta     ######
    k_keep = KL(comsam,dis="delta",Kmax=10,mstree_cluster_delta)
    plab=mstree_cluster_delta(comsam,dis=delta,k=k_keep)
    rdx=randindex(obslab,plab)
    performance[pxx,10] =  performance[pxx,10]+rdx
 
  } #for sim
       show(performance)
  } # for p
 
  performance = performance/sim
  show(round(performance,digits=3))
} #Idis



 