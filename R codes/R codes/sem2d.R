
library(MASS)
library(spatstat)
library(Matrix)
setwd("C:/Users/jimaj/Desktop/Project_2020")
area<-owin(c(0,2),c(2,4),unitname = NULL)
dp_2<-read.csv("dp_2.csv")
dcp_2<-read.csv("dcp_2.csv")
x<-c(dp_2$V1,dcp_2$V1)
y<-c(dp_2$V2,dcp_2$V2)
data_2<-ppp(x,y,area)
par(new=FALSE)
plot(data_2,type="p",pch=20,main="Simulated 2D point pattern ")
###

### edge corrected lisa function

tc.pdf.ppp<-function(pp){
  
  W<-pp$window
  df<-data.frame(cbind(pp$x,pp$y))
  ed<-as.matrix(pairdist.ppp(pp))
  a<-area.default(pp$window)
  n<-pp$n
  e<- (sqrt(5)/(10))/(sqrt(intensity(pp)))
  dt<-seq(e+0.004,0.25*min(pp$window$xrange[2],pp$window$yrange[2]),length=25)
  
  par(new=FALSE)
  ec<-edge.Trans(pp)
  
  Y<-matrix(data = NA,nrow = length(dt),ncol = n)
  Y
  r<-c()
  for (k in 1:n) {
    
    
    fi<-matrix(data = ed[-k,k],nrow = length(dt),ncol = n-1, byrow=TRUE)
    fi<-fi-dt
    f<-matrix(data = 0,nrow = length(dt),ncol = n-1, byrow=TRUE)
    good<-which(-e <= fi  & fi<=e )
    f[good]<-(3/(4*e))*(1- (fi[good]/e)^2)
    
    #r<- ((n-1)/(2*pi*dt*a))*rowSums(((2*pi*f)%*% diag(1/c(ec[-k,k])))%*%diag(c(ed[-k,k])))
    r<- ((n-1)/(2*pi*dt*a))*rowSums(f%*% diag(c(ec[-k,k])))
    Y[,k]<-r
    plot(dt,r,type="l",pch=20,ylim = c(0,300000),xlim=c(0,dt[length(dt)]+0.004),ylab = "pdf",xlab = "r")
    par(new=TRUE)
  }
  Y
  
}
###############

#t number of iterations


pp<-data_2
t<-20

#edge corrected lisa functions of data 
Y<-tc.pdf.ppp(pp)

#expected value under homo.poisson 
hp<-intensity(pp)^2+(intensity(pp)/area.default(pp$window))

d<-Y-hp
d
#d1 L1 distances between lisa and expected
d1<-colSums(abs(d))
#d1

###### d2

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

N<-cbind(data_2$x,data_2$y,d1,l)
N

X<-d1
W<-pp$window
a<-area.default(pp$window)
n<-length(X)


## simulation of clutter and feature points
#selected window of feature
#wf<-owin(c(0.2,0.3),c(2.8,2.9),unitname = NULL)
#spp<-rpoispp(300/area(W),win = W,drop = TRUE)
#dspp<-as.data.frame(cbind(spp$x,spp$y))
#write.csv(dspp,file = "dspp.csv")
spp<-read.csv("dspp.csv")
#dscp<-as.data.frame(cbind(scp$x,scp$y))
#write.csv(dscp,file = "dscp.csv")
#scp<-rpoispp(35/area(wf),win = wf,drop = TRUE)
scp<-read.csv("dscp_35.csv")
spp
scp

sx<-c(spp$V1,scp$V1)
sy<-c(spp$V2,scp$V2)
sd<-ppp(sx,sy,area)
par(new=FALSE)
plot(sd,type="p",pch=20,main="Simulated feature-clutter for s-step ")

sy<-tc.pdf.ppp(sd)

dp<-sy-hp
ddp<-sqrt(colSums((dp^2)))
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
p[1]<-0.6
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


D<-cbind(z,data_2$x,data_2$y)
Dp<-subset(D,z==0)
  Dc<-subset(D,z==1)
  par(mfrow=c(1,1))
plot(ppp(Dp[,2],Dp[,3],area),type="p",pch=20,main="SEM with L2 and 0.6 ")
par(new=TRUE)
plot(ppp(Dc[,2],Dc[,3],area),cols="red",type="p",pch=20,main = "2D SEM")


D<-cbind(z,data_2$x,data_2$y,Z=c(rep(0,105),rep(1,98)))

subset(D,D[,1]==1 & D[,4]==1)
subset(D,D[,1]==1 & D[,4]==0)
subset(D,D[,1]==0 & D[,4]==0)
subset(D,D[,1]==0 & D[,4]==1)
###################
######################################################









