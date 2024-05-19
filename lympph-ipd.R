library(spls)
data(lymphoma)
lym=as.data.frame(lymphoma)
#str(lym)
xx=lym$x
yy=lym$y
show(yy)
n1=sum(yy==0)
n2=sum(yy==1)
n3=sum(yy==2)
N=dim(lym)[1]
p=dim(lym)[2]
cat('p=',p,'N=',N,'n1=',n1,'n2=',n2,'n3=',n3,fill=TRUE)
 
 
m1=n1*(n1-1)/2
m2=n2*(n2-1)/2
m3=n3*(n3-1)/2


xg1=lym[yy==0,]
xg2=lym[yy==1,]
xg3=lym[yy==2,]


xmat1 <- as.matrix(dist(xg1,method = "euclidean"))

ipdx1=rep(0,m1)
count=1
for(i in 1:(n1-1)){for(j in (i+1):n1){ipdx1[count]=xmat1[i,j]; count <- count+1 }}



xmat2 <- as.matrix(dist(xg2,method = "euclidean"))

ipdx2=rep(0,m2)
count=1
for(i in 1:(n2-1)){for(j in (i+1):n2){ipdx2[count]=xmat2[i,j]; count <- count+1 }}


xmat3 <- as.matrix(dist(xg3,method = "euclidean"))

ipdx3=rep(0,m3)
count=1
for(i in 1:(n3-1)){for(j in (i+1):n3){ipdx3[count]=xmat3[i,j]; count <- count+1 }}


hist(ipdx1,ylab='Freq.',xlab='IPD',main="DLBCL IPD Distribution")
hist(ipdx2,ylab='Freq.',xlab='IPD',main="FL IPD Distribution")
hist(ipdx3,ylab='Freq.',xlab='IPD',main="CLL IPD Distribution")



ipd=c(ipdx1,ipdx2,ipdx3)
ipdx1bar=sum(ipdx1)/m1
ipdx2bar=sum(ipdx2)/m2
ipdx3bar=sum(ipdx3)/m3



ipd=c(ipdx1,ipdx2,ipdx3bar)
cat('ipdx1bar=',ipdx1bar,'ipdx2bar=',ipdx2bar,'ipdx3bar=',ipdx3bar,fill=TRUE)
minipd=min(ipd)
maxipd=max(ipd)
deltaint=100
unit=(maxipd-minipd)/deltaint
delta=seq(minipd,maxipd,unit)
cdfx1=cdfx2=cdfx3=cdfx4=cdfx5=rep(0,deltaint)
for (i in 1:deltaint){
  cdfx1[i]=sum(ipdx1<=delta[i])/m1
  cdfx2[i]=sum(ipdx2<=delta[i])/m2
  cdfx3[i]=sum(ipdx3<=delta[i])/m3
  
}
plot(cdfx1,type="l",lty=2,  col="red", lwd=1,ylab="ECDF", xlab = "IPD Distance")
title(main = "Lymphoma Data")
lines(cdfx2,lty=3,lwd=2,col="blue")
lines(cdfx3,lty=4,lwd=2,col="green")
legend("topleft",legend=c("CLL","FL","DLBCL"), col=c("green","blue","red"), lty=c(2,3,4),lwd=1, ncol=1)



