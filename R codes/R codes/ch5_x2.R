
####    Author: W.Jithmi Jannidi
####    M.Sc Thesis : Local Indicators of Spatial Association of 2d & 3d Point Patterns

####    R code for the example discussed in the chapter 5 , see page 66
####    simulate point process with more than one cluster


#####
library(MASS)
library(spatstat)
library(Matrix)
setwd("C:/Users/jimaj/Desktop/Project_2020")
area<-owin(c(0,1),c(0,1),unitname = NULL)
#pp<-rpoispp(300,win = area,drop = TRUE)
 #plot(pp,type="p",pch=20)
# cp<-rpoispp(10000,win =owin(c(0.5,0.6),c(0.5,0.6)),drop = TRUE)
# par(new=TRUE)
#  plot(cp,type="p",pch=1,cols="red")
# 
#dp_2<-as.data.frame(cbind(pp$x,pp$y))
#  dp_2
 #write.csv(dp_2,file = "dp_25.csv")
#  dcp_2<-as.data.frame(cbind(cp$x,cp$y))
#  dcp_2
#  write.csv(dcp_2,file = "dcp_25.csv")
# 
#  cp1<-rpoispp(20000,win =owin(c(0.2,0.25),c(0.3,0.4)),drop = TRUE)
#  dcp_1<-as.data.frame(cbind(cp1$x,cp1$y))
#  dcp_1
# write.csv(dcp_1,file = "dcp_15.csv")
 dp_2<-read.csv("dp_25.csv")
dcp_2<-read.csv("dcp_25.csv")
dcp_1<-read.csv("dcp_15.csv")



# x<-c(dcp_1$V1,dcp_2$V1)
# y<-c(dcp_1$V2,dcp_2$V2)
# data_2<-ppp(x,y,area)
# 
# par(mfrow=c(1,3))
# # par(new=FALSE)
# plot(data_2,type="p",pch=19,main="Origina feature of X_2 ")
# 
# plot(ppp(dp_2$V1,dp_2$V2,area),type="p",pch=19,main="Origina clutter of X_2 ")
x<-c(dp_2$V1, dcp_1$V1,dcp_2$V1)
y<-c(dp_2$V2, dcp_1$V2,dcp_2$V2)
data_2<-ppp(x,y,area)
plot(data_2,type="p",pch=".",main=" point process X_2 ") ## simulated 2d pp with two features


# 
# plot(ppp(dp_2$V1,dp_2$V2,area),type="p",pch=19,main="Simulated 2D point pattern ")
# par(new=TRUE)
# plot(ppp(dcp_1$V1,dcp_1$V2,area),type="p",pch=19,main=" ",cols = "red")
# par(new=TRUE)
# plot(ppp(dcp_2$V1,dcp_2$V2,area),type="p",pch=19,main=" ",cols = "blue")




#######################2d lisa

# tc.pdf.ppp<-function(pp){
#   plot(pp)
#   W<-pp$window
#   
#   
#   df<-data.frame(cbind(pp$x,pp$y))
#   #ed<- dist(df, method = "euclidean", diag = FALSE, upper = FALSE, p = 2)
#   ed<-as.matrix(pairdist.ppp(pp))
#   a<-area.default(pp$window)
#   n<-pp$n
#   e<- (sqrt(5)/(10))/(sqrt(intensity(pp)))
#   dt<-seq(e+0.004,0.25*min(pp$window$xrange[2],pp$window$yrange[2]),length=25)
#   r<-c()
#   par(new=FALSE)
#   ec<-edge.Trans(pp)
#   ### matrix of local pdf of each points <- Y
#   Y<-matrix(data = NA,nrow = length(dt),ncol = n)
#   Y
#   
#   
#   for (k in 1:n) {
#     for (j in 1:length(dt)) {
#       fi<-c(subset(ed[,k],ed[,k]!=0)-dt[j])
#       
#       f<-c()
#       nn<-n-1
#       for (i in 1:nn)
#         
#       {
#         if (-e <= fi[i]  & fi[i]<=e )
#         {
#           f[i]<-(3/(4*e))*(1- (fi[i]/e)^2)
#         }
#         
#         
#         
#         else
#         {
#           f[i]<- 0
#         }
#         
#       }
#       r[j]<- ((n-1)/(2*pi*dt[j]*a))*sum(c(subset(ec[,k],ec[,k]!=1))*f)
#       Y[j,k]<-r[j]
#       
#     }
#     
#     plot(dt,r,type="l",pch=20,ylim = c(0,300000),xlim=c(0,dt[length(dt)]+0.004),ylab = "edge corrected pdf",xlab = "r")
#     par(new=TRUE)
#     
#   }
#   
#   return(list(Y))
# }
# 
# 
# 
# 
 # Y2<-tc.pdf.ppp(data_2)[[1]]
 # Y2<-as.data.frame(Y2)
 # # 
 # Y2<-write.csv(Y2,file = "Y25.csv")
Y2<-read.csv("Y25.csv")
Y2<-as.matrix(Y2[,1:data_2$n+1])

h2<-hclust( dist(t(Y2)))
h2
par(mfrow=c(1,1))
plot(h2,xlab = "")

memb<-cutree(h2,k=2)
memb


t<-cbind(data_2$x,data_2$y,c(rep(1,325),rep(2,204)),c(as.numeric(memb)))
t1<-subset(t, t[,4]==1)
t2<-subset(t,t[,4]==2)
#t3<-subset(t,t[,4]==3)
par(mfrow=c(1,2))
plot(ppp(t1[,1],t1[,2],area),type="p",pch=".",main="estimated clutter ")


plot(ppp(t2[,1],t2[,2],area),type="p",pch=".",main="estimated feature")




subset(t,t[,3]==1 & t[,4]==1)
subset(t,t[,3]==1 & t[,4]==2)
subset(t,t[,3]==2 & t[,4]==2)
subset(t,t[,3]==2 & t[,4]==1)




################################   nnclean function EM method

x<-c(dp_2$V1,dcp_1$V1,dcp_2$V1)
y<-c(dp_2$V2,dcp_1$V2,dcp_2$V2)
data_2<-ppp(x,y,area)
plot(data_2,type="p",pch=".",main="Simulated 2D point pattern ")


n2<-nnclean(data_2,k=5)
d2<-cbind(n2$x,n2$y,n2$marks$class)
noise2<-subset(d2,d2[,3]==1)
fea<-subset(d2,d2[,3]==2)
fea
t<-cbind(c(rep(1,325),rep(2,204)),d2[,3])

par(mfrow=c(1,2))
plot(ppp(noise2[,1],noise2[,2],area),type="p",pch=".",main = "estimated clutter")

plot(ppp(fea[,1],fea[,2],area),type="p",pch=".",main = "estimated feature")

subset(t,t[,1]==1 & t[,2]==1)
subset(t,t[,1]==1 & t[,2]==2)
subset(t,t[,1]==2 & t[,2]==2)
subset(t,t[,1]==2 & t[,2]==1)

fea
noise2


########################## SEM


pp<-data_2
t<-20

#edge corrected lisa functions of data 
Y<-Y2

#expected value under homo.poisson 
hp<-intensity(pp)^2+(intensity(pp)/area.default(pp$window))

d<-Y-hp
d
#d1 L1 distances between lisa and expected
#d1<-colSums(abs(d))


###### d2

d1<-sqrt(colSums((d^2)))

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
wf<-owin(c(0.2,0.3),c(0.8,0.9),unitname = NULL)
#spp<-rpoispp(300/area(W),win = W,drop = TRUE)
#dspp<-as.data.frame(cbind(spp$x,spp$y))
#write.csv(dspp,file = "x2semc.csv")
spp<-read.csv("x2semc.csv")
#dscp<-as.data.frame(cbind(scp$x,scp$y))
#write.csv(dscp,file = "x2semf.csv")
#scp<-rpoispp(225/area(wf),win = wf,drop = TRUE)
scp<-read.csv("x2semf.csv")
spp
scp


plot(ppp(spp$V1,spp$V2,area),type="p",pch=20,main="Simulated clutter for s-step ")
plot(ppp(scp$V1,scp$V2,area),type="p",pch=20,main="Simulated feature for s-step ")

sx<-c(spp$V1,scp$V1)
sy<-c(spp$V2,scp$V2)
sd<-ppp(sx,sy,area)
par(new=FALSE)
plot(sd,type="p",pch=20,main="Simulated feature-clutter for s-step ")

sy<-tc.pdf.ppp(sd)

dp<-sy-hp
ddp<-sqrt(colSums((dp^2)))
#ddp<-colSums(abs(dp))
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
p[1]<-0.3
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
par(mfrow=c(1,2))
plot(ppp(Dp[,2],Dp[,3],area),type="p",pch=20,main="Estimated clutter  ")

plot(ppp(Dc[,2],Dc[,3],area),type="p",pch=20,main = "Estimated feature")


D<-cbind(z,data_2$x,data_2$y,Z=c(rep(0,325),rep(1,204)))

subset(D,D[,1]==1 & D[,4]==1)
subset(D,D[,1]==1 & D[,4]==0)
subset(D,D[,1]==0 & D[,4]==0)
subset(D,D[,1]==0 & D[,4]==1)






















