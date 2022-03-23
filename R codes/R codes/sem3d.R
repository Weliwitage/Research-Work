
library(MASS)
library(spatstat)
library(Matrix)
setwd("C:/Users/jimaj/Desktop/Project_2020")
v<-box3(xrange = c(0, 5), yrange = c(0, 5), zrange = c(0, 10), unitname = NULL)

p3d<-read.csv("dp3d.csv")
c13d<-read.csv("d3c1.csv")

x<-c(p3d$V1,c13d$V1)
y<-c(p3d$V2,c13d$V2)
z<-c(p3d$V3,c13d$V3)

data_3d<-pp3(x,y,z,box3=v)
plot(data_3d,type="p",pch=20,main="Simulated 3D point process ")
### edge corrected lisa function

tc.pdl.pp3<-function(pp){
  
  
  box<-pp$domain
  e<- (sqrt(5)/(10))/(sqrt(intensity(pp)))
  dt<-seq(e+0.004,0.25*min(box$xrange[2],box$yrange[2],box$zrange[2]),length=25) 
  df<-data.frame(cbind(pp$data$x,pp$data$y,pp$data$z))
  ed<- dist(df, method = "euclidean", diag = FALSE, upper = FALSE, p = 2)
  ed<-as.matrix(ed)
  
  a<-volume(pp$domain)
  n<-length(pp$data$x)
  
  
  par(new=FALSE)
  
  edge.Trans.pp3 <- function(X, trim = spatstat.options("maxedgewt")) {
    stopifnot(verifyclass(X, "pp3"))
    stopifnot(verifyclass(X$domain, "box3"))
    dx <- outer(X$data$x, X$data$x, "-")
    dy <- outer(X$data$y, X$data$y, "-")
    dz <- outer(X$data$z, X$data$z, "-")
    width <- diff(X$domain$xrange)
    depth <- diff(X$domain$yrange)
    height <- diff(X$domain$zrange)
    weight <- volume(X$domain) / ((width - abs(dx)) * (depth - abs(dy)) * (height - abs(dz)))
    if (length(weight) > 0) {
      weight <- pmin.int(weight, trim)
    }
    matrix(weight, ncol = npoints(X))
  }
  ##########################################
  
  
  ec<-edge.Trans.pp3(pp)
  
  
  
  Y<-matrix(data = NA,nrow = length(dt),ncol = n)
  
  for (k in 1:n) {
    fi<-matrix(data = ed[-k,k],nrow = length(dt),ncol = n-1, byrow=TRUE)
    fi<-fi-dt
    
    f<-matrix(data = 0,nrow = length(dt),ncol = n-1, byrow=TRUE)
    good<-which(-e <= fi  & fi<=e )
    f[good]<-(3/(4*e))*(1- (fi[good]/e)^2)
    
    
    Y[,k]<- ((n-1)/(4*pi*dt^2*a))*rowSums(f*matrix(data= c(ec[-k,k]),ncol = n-1,nrow = length(dt),byrow = TRUE))
    
    
    
    plot(dt,Y[,k],type="l",pch=20,ylim = c(0,30),xlim=c(0,dt[length(dt)]+0.004),ylab = "pdf",xlab = "r")
    par(new=TRUE)
  }
  Y
}


###############

#t number of iterations


pp<-data_3d
t<-20

#edge corrected lisa functions of data 
Y<-tc.pdl.pp3(pp)

#expected value under homo.poisson 
hp<-intensity(pp)^2+(intensity(pp)/volume(pp$domain))

d<-Y-hp
d
#d1 L1 distances between lisa and expected
d1<-colSums(abs(d))
#d1

########l2
#d1<-sqrt(colSums((d^2)))

#l initial classification, 
l<-c()

for (i in 1:length(d1)) {
  
  
  if (d1[i]> mean(d1))
  {
    l[i]<-1
  }
  
  else 
    l[i]<-0
}

l



X<-d1
W<-pp$domain
a<-volume(W)
n<-length(X)


## simulation of clutter and feature points
#selected window of feature
wf<-box3(xrange = c(3, 4), yrange = c(2,3), zrange = c(5,6), unitname = NULL)

#spp<-rpoispp3(300/a, domain = v, nsim=1, drop=TRUE)
#dspp<-as.data.frame(cbind(spp$data$x,spp$data$y,spp$data$z))
#write.csv(dspp,file = "dspp3.csv")
spp<-read.csv("dspp3.csv")
#dscp<-as.data.frame(cbind(scp$x,scp$y))
#write.csv(dscp,file = "dscp.csv")
#scp<-rpoispp3(50,domain = wf, nsim=1, drop=TRUE)
scp<-read.csv("dscp3_40_113.csv")
spp
scp

sx<-c(spp$V1,scp$V1)
sy<-c(spp$V2,scp$V2)
sz<-c(spp$V3,scp$V3)
sd<-pp3(sx,sy,sz,box3=v)
par(new=FALSE)
plot(sd,type="p",pch=20,main="Simulated 3D mixture process for s- step in SEM ")

sy<-tc.pdl.pp3(sd)

dp<-sy-hp
ddp<-colSums(abs(dp))
spn<-length(spp$X)
scn<-length(scp$X)
denp<-density(ddp[1:spn])
denc<-density(ddp[spn+1:scn])

rX<- range(X)
rp<-range(denp$x)
rc<-range(denc$x)

rX
rp
rc
subset(X,rp[1]<X & X<rp[2])

subset(X,rc[1]<X & X<rc[2])



#l1 and l2 are intensities of clutter and featires, n1 and n2 are number of points of clutter and features
p<-matrix(data=NA,ncol=1,nrow = t)
l1<-matrix(data=NA,ncol=1,nrow = t)
l2<-matrix(data=NA,ncol=1,nrow = t)
n1<-matrix(data=NA,ncol=1,nrow = t)
n2<-matrix(data=NA,ncol=1,nrow = t)
#initial values 
p[1]<-0.7
l1[1]<-length(subset(l,l==0))/a
l2[1]<-length(subset(l,l==1))/a
n1[1]<-length(subset(l,l==0))
n2[1]<- length(subset(l,l==1)) 
t1<-t-1

cdenp<-cbind(denp$x,denp$y)
cdenc<-cbind(denc$x,denc$y)

for (i in 1:t1) {
  
  
  f1<-matrix(data=NA,ncol=1,nrow = n)
  f2<-matrix(data=NA,ncol=1,nrow = n)
  
  for (j in 1:n) {
    
    if (X[j]>rp[1]&X[j]<rp[2])
    {
      
      newd<-(c(denp$x,X[j]))
      s<-sort(newd)
      w<-which(s==X[j])
      V1<-s[w-1]
      V2<-s[w+1]
      
      v1<-subset(cdenp,cdenp[,1]==V1)
      v2<-subset(cdenp,cdenp[,1]==V2)
      
      
      y1j<-(v1[,2]*(v2[,1]-X[j])+v2[,2]*(X[j]-v1[,1]))/(v2[,1]-v1[,1])
      
    }
    else
    {
      y1j<-0
    }
    
    
    if (X[j]>rc[1]&X[j]<rc[2])
    {
      
      newd<-(c(denc$x,X[j]))
      s<-sort(newd)
      w<-which(s==X[j])
      V1<-s[w-1]
      V2<-s[w+1]
      
      v1<-subset(cdenc,cdenc[,1]==V1)
      v2<-subset(cdenc,cdenc[,1]==V2)
      
      
      y2j<-(v1[,2]*(v2[,1]-X[j])+v2[,2]*(X[j]-v1[,1]))/(v2[,1]-v1[,1])
      
    }
    else
    {
      y2j<-0
    }
    
    
    
    
    
    f1[j]<- (p[i]*y1j)/((p[i]*y1j)+((1-p[i])*y2j))
    f2[j]<-1-f1[j]
    
  }
  
  f1
  f2
  
  
  p[i+1]<-sum(f1)/n
  z<-matrix(data=NA,ncol=1,nrow = n)
  for (j in 1:n) {
    if(  f1[j] > f2[j] )
    {
      z[j]<-0
    }
    
    else
      
      z[j]<-1
    
  }
  z 
  
  n1[i+1]<-length(subset(z,z==0))
  n2[i+1]<-length(subset(z,z==1))
  
  l1[i+1]<-n1[i+1]/a
  l2[i+1]<-n2[i+1]/a
  
  
  
  
  
}

p
n1
n2
l1
l2
z

par(mfrow=c(5,1))

plot(p,type = "b",col="red",pch=19,ylab = "mixing probaility")
plot(n1,type = "b",col="red",pch=19,ylab = "n_c")
plot(n2,type = "b",col="red",pch=19,ylab="n_f")
plot(l1,type = "b",col="red",pch=19,ylab="l_c")
plot(l2,type = "b",col="red",pch=19,ylab="l_f")



D<-cbind(z,data_3d$data$x,data_3d$data$y,data_3d$data$z)
Dp<-subset(D,z==0)
Dc<-subset(D,z==1)
par(mfrow=c(1,1))
plot(pp3(Dc[,2],Dc[,3],Dc[,4],v),col="red",type="p",pch=20,main = "3D SEM for L1")



D<-cbind(z,data_3d$data$x,data_3d$data$y,data_3d$data$z,Z=c(rep(0,92),rep(1,106)))

subset(D,D[,1]==1 & D[,5]==1)
subset(D,D[,1]==1 & D[,5]==0)
subset(D,D[,1]==0 & D[,5]==0)
subset(D,D[,1]==0 & D[,5]==1)
par(new=TRUE)
plot(pp3(Dc[,2],Dc[,3],Dc[,4],v),col="red")
###################
######################################################

par(mfrow=c(1,3))

plot(pp3(scp$V1,scp$V2,scp$V3,box3=v),type="p",pch=20,main="Simulated 3D feature points for s-step in SEM")
plot(pp3(spp$V1,spp$V2,spp$V3,box3=v),type="p",pch=20,main="Simulated 3D clutter points for s-step in SEM")
