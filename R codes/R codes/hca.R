
####    Author: W.Jithmi Jannidi
####    M.Sc Thesis : Local Indicators of Spatial Association of 2d & 3d Point Patterns

####    R code for Hierarchical cluster analysis (hca) section 4.3.1.a
####    examples of simulated point processes ex.4.2.1.b 2d pp and ex.4.2.2.a 3d pp
####    areas covered in this code : edge corrected local product density function, hca

library(MASS)
library(spatstat)

## the function tc.pfd.ppp is created to calculate translation edge corrected product density function of a 2d pp
##  input is any 2d point pattern =pp
## output Y= matrix of each points's local edge corrected product density function at each distance of dt
tc.pdf.ppp<-function(pp){
  
  plot(pp)    #plot the point pattern (pp)
  W<-pp$window   #window of pp
  
  
  df<-data.frame(cbind(pp$x,pp$y))   #store x,y coordinates in a data frame
  #ed<- dist(df, method = "euclidean", diag = FALSE, upper = FALSE, p = 2)
  ed<-as.matrix(pairdist.ppp(pp))  # store pairwise distances between points in a matrix
  a<-area.default(pp$window)   # area of the window
  n<-pp$n    # total number of points
  e<- (sqrt(5)/(10))/(sqrt(intensity(pp)))   # bandwidth of the kernel function , see 4.2 equation
  dt<-seq(e+0.004,0.25*min(pp$window$xrange[2],pp$window$yrange[2]),length=25)  # create a sequence of 25 distances lying in the window 
  r<-c()  # edge corrected local pdf 
  par(new=FALSE)
  ec<-edge.Trans(pp) # Translation  edge correction weights 
  ### Y is a matrix of dt=number of distance values by n=number of point to store local edge corrected pdf of each points
  Y<-matrix(data = NA,nrow = length(dt),ncol = n)
  Y
  
  # kth point's jth distance value
  for (k in 1:n) {
    for (j in 1:length(dt)) {
      # fi input value of function fe in eq. 4.7 see page 26
      fi<-c(subset(ed[,k],ed[,k]!=0)-dt[j]) # distance between pairwise distance and distance from created sequence see. eq.4.7  != for not comparing with point itself 
      f<-c() #Epanechnikov kernel function see eq 4.2
      nn<-n-1
      for (i in 1:nn)
        
      {
        if (-e <= fi[i]  & fi[i]<=e )
        {
          f[i]<-(3/(4*e))*(1- (fi[i]/e)^2)
        }
        
        
        
        else
        {
          f[i]<- 0
        }
        
      }
      r[j]<- ((n-1)/(2*pi*dt[j]*a))*sum(c(subset(ec[,k],ec[,k]!=1))*f) # local product density value of kth point for jth distance
      Y[j,k]<-r[j]
      
    }
    
    plot(dt,r,type="l",pch=20,ylim = c(0,300000),xlim=c(0,dt[length(dt)]+0.004),ylab = "edge corrected pdf",xlab = "r")
    par(new=TRUE)
    
  }
  
  return(list(Y))
}


area<-owin(c(0,2),c(2,4),unitname = NULL)
dp_2<-read.csv("dp_2.csv")
dcp_2<-read.csv("dcp_2.csv")
x<-c(dp_2$V1,dcp_2$V1)
y<-c(dp_2$V2,dcp_2$V2)
data_2<-ppp(x,y,area)
plot(data_2,type="p",pch=20,main="Simulated 2D point pattern ")



tc.pdf.ppp(data_2)

Y<-tc.pdf.ppp(data_2)[[1]]
Y

h2<-hclust( dist(t(Y)))    #´hierarchical cluster analysis on edge corrected local product density functions
h2
plot(h2)

memb<-cutree(h2,k=2)  #cut the hierarchical clustering tree into k number of groups
memb


t<-cbind(data_2$x,data_2$y,c(rep(1,105),rep(2,98)),memb) # x coordinates,y coordinates, labels of points 1 for parent 2 for feature, labels from hca
t1<-subset(t, t[,4]==1) #points identified as parents from hca
t2<-subset(t,t[,4]==2)  # points identified as features from hca
plot(ppp(t1[,1],t1[,2],area),type="p",pch=20,main="2D HCA ")

par(new=TRUE)
plot(ppp(t2[,1],t2[,2],area),type="p",pch=20,main="2D HCA",cols = "red")

subset(t,t[,3]==1 & t[,4]==1) #correctly identified parent points
subset(t,t[,3]==1 & t[,4]==2) #parent points identified as features
subset(t,t[,3]==2 & t[,4]==2)  #correctly identified feature points
subset(t,t[,3]==2 & t[,4]==1) #feature points identified as parents


######################### let's do the same procedure for a 3d point pattern

v<-box3(xrange = c(0, 5), yrange = c(0, 5), zrange = c(0, 10), unitname = NULL)

p3d<-read.csv("dp3d.csv")
c13d<-read.csv("d3c1.csv")


x<-c(p3d$V1,c13d$V1)
y<-c(p3d$V2,c13d$V2)
z<-c(p3d$V3,c13d$V3)

data_3d<-pp3(x,y,z,box3=v)
plot(data_3d,type="p",pch=20,main="Simulated 3D point process ")


## the function tc_pfd.pp3 is created to calculate translation edge corrected product density function of a 3d pp
##  input is any 3d point pattern =pp
## output Y= matrix of each points's local edge corrected product density function at each distance of dt

tc_pdf.pp3<-function(pp){
  box<-pp$domain
  e<- (sqrt(5)/(10))/(sqrt(intensity(pp)))
  dt<-seq(e+0.004,0.25*min(box$xrange[2],box$yrange[2],box$zrange[2]),length=25) 
  
  # edge.Trans.pp3 - function for translation edge correction for 3d point pattern 
  # maxedgewt trimmed edge correction weights,not exceed to this value
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
      weight <- pmin.int(weight, trim) # return single min
    }
    matrix(weight, ncol = npoints(X))
  }
  #########################################
  
  ec<-edge.Trans.pp3(pp)
  ec
  
  df<-data.frame(cbind(pp$data$x,pp$data$y,pp$data$z))
  #ed <- 3d euclidean distances
  ed<- dist(df, method = "euclidean", diag = FALSE, upper = FALSE, p = 2)
  ed<-as.matrix(ed)
  # a<- volume of the box
  a<-volume(pp$domain)
  n<-length(pp$data$x)
  r<-c()
  par(new=FALSE)
  
  ### matrix of local pdf of each points - Y
  Y<-matrix(data = NA,nrow = length(dt),ncol = n)
  Y
  
  
  for (k in 1:n) {
    for (j in 1:length(dt)) {
      fi<-c(subset(ed[,k],ed[,k]!=0)-dt[j])
      
      f<-c()
      nn<-n-1
      for (i in 1:nn)
        
      {
        if (-e <= fi[i]  & fi[i]<=e )
        {
          f[i]<-(3/(4*e))*(1- (fi[i]/e)^2)
        }
        
        
        
        else
        {
          f[i]<- 0
        }
        
      }
      r[j]<- ((n-1)/(4*pi*dt[j]^2*a))*sum(f*c(subset(ec[,k],ec[,k]!=ec[k,k])))
      
      Y[j,k]<-r[j]
      
    }
    
    plot(dt,r,type="l",pch=20,ylim = c(0,30),xlim=c(0,dt[length(dt)]+0.004),ylab = "edge corrected 3D pdf",xlab = "r")
    par(new=TRUE)
  }
  
  return(list(Y))
}

tc_pdf.pp3(data_3d)



Y<-tc_pdf.pp3(data_3d)[[1]]
Y

h3<-hclust( dist(t(Y)))
h3
plot(h3)

memb<-cutree(h3,k=2)
memb


t<-cbind(data_3d$data$x,data_3d$data$y,data_3d$data$z,c(rep(1,92),rep(2,106)),memb)
t
t1<-subset(t, t[,5]==1)
t2<-subset(t, t[,5]==2)
plot(pp3(t1[,1],t1[,2],t1[,3],box3=v),type="p",pch=20,main="3D HCA",col="red")

par(new=TRUE)
plot(pp3(t2[,1],t2[,2],t2[,3],box3=v),type="p",pch=20,main="3D HCA",col="red")

subset(t,t[,4]==2 & t[,5]==2)
subset(t,t[,4]==2 & t[,5]==1)
subset(t,t[,4]==1 & t[,5]==1)
subset(t,t[,4]==1 & t[,5]==2)



