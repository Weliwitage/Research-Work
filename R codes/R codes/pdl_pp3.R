####    Author: W.Jithmi Jannidi
####    M.Sc Thesis : Local Indicators of Spatial Association of 2d & 3d Point Patterns

####    R code for generating local product density functions of 3d pp 
####      with and without edge correction





library(MASS)
library(spatstat)
library(Matrix)

v<-box3(xrange = c(0, 5), yrange = c(0, 5), zrange = c(0, 10), unitname = NULL)

setwd("C:/Users/jimaj/Desktop/Project_2020")
p3d<-read.csv("dp3d.csv")
c13d<-read.csv("d3c1.csv")


x<-c(p3d$V1,c13d$V1)
y<-c(p3d$V2,c13d$V2)
z<-c(p3d$V3,c13d$V3)

data_3d<-pp3(x,y,z,box3=v)
plot(data_3d,type="p",pch=20,main="Simulated 3D point process ")

###

###pdl.pp3 product density lisa functions of 3d pp

pdl.pp3<-function(pp){
  
  box<-pp$domain
  e<- (sqrt(5)/(10))/(sqrt(intensity(pp)))
  dt<-seq(e+0.004,0.25*min(box$xrange[2],box$yrange[2],box$zrange[2]),length=25) 
  df<-data.frame(cbind(pp$data$x,pp$data$y,pp$data$z))
  ed<- dist(df, method = "euclidean", diag = FALSE, upper = FALSE, p = 2)
  ed<-as.matrix(ed)
  
  a<-volume(pp$domain)
  n<-length(pp$data$x)
  
  
  
  Y<-matrix(data = NA,nrow = length(dt),ncol = n)
  
  for (k in 1:n) {
    fi<-matrix(data = ed[-k,k],nrow = length(dt),ncol = n-1, byrow=TRUE)
    fi<-fi-dt
    
    f<-matrix(data = 0,nrow = length(dt),ncol = n-1, byrow=TRUE)
    good<-which(-e <= fi  & fi<=e )
    f[good]<-(3/(4*e))*(1- (fi[good]/e)^2)
    
    
    Y[,k]<- ((n-1)/(4*pi*dt^2*a))*rowSums(f)
    
    
    
    plot(dt,Y[,k],type="l",pch=20,ylim = c(0,30),xlim=c(0,dt[length(dt)]+0.004),ylab = "pdf",xlab = "r")
    par(new=TRUE)
  }
  Y
}

pdl.pp3(data_3d)


#### tc.pdl.pp3 product density lisa functions of a 3d pp with transform edge correction

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

tc.pdl.pp3(data_3d)



