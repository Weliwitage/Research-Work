####    Author: W.Jithmi Jannidi
####    M.Sc Thesis : Local Indicators of Spatial Association of 2d & 3d Point Patterns

####    R code for generating local product density functions of 2d pp 
####      with and without edge correction




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

###pdl.ppp product density lisa functions of 2d pp


pdl.ppp<-function(pp){
  
  W<-pp$window
  #  df<-data.frame(cbind(pp$x,pp$y))
  ed<-as.matrix(pairdist.ppp(pp))
  a<-area.default(pp$window)
  n<-pp$n
  e<- (sqrt(5)/(10))/(sqrt(intensity(pp)))
  dt<-seq(e+0.004,0.25*min(pp$window$xrange[2],pp$window$yrange[2]),length=25)
  
  
  
  Y<-matrix(data = NA,nrow = length(dt),ncol = n)

  for (k in 1:n) {
    fi<-matrix(data = ed[-k,k],nrow = length(dt),ncol = n-1, byrow=TRUE)
    fi<-fi-dt
   
    f<-matrix(data = 0,nrow = length(dt),ncol = n-1, byrow=TRUE)
    good<-which(-e <= fi  & fi<=e )
    f[good]<-(3/(4*e))*(1- (fi[good]/e)^2)
    
    
    Y[,k]<- ((n-1)/(2*pi*dt*a))*rowSums(f)
    
    
    
    plot(dt,Y[,k],type="l",pch=20,ylim = c(0,300000),xlim=c(0,dt[length(dt)]+0.004),ylab = "pdf",xlab = "r")
    par(new=TRUE)
  }
  Y
  }
  
pdl.ppp(data_2)

###tc.pdl.ppp product density lisa functions with transform edge corrections of 2d pp


tc.pdl.ppp<-function(pp){
  
  W<-pp$window
  #  df<-data.frame(cbind(pp$x,pp$y))
  ed<-as.matrix(pairdist.ppp(pp))
  a<-area.default(pp$window)
  n<-pp$n
  e<- (sqrt(5)/(10))/(sqrt(intensity(pp)))
  dt<-seq(e+0.004,0.25*min(pp$window$xrange[2],pp$window$yrange[2]),length=25)
  
  par(new=FALSE)
  ec<-edge.Trans(pp)
  
  
  Y<-matrix(data = NA,nrow = length(dt),ncol = n)
  
  for (k in 1:n) {
    fi<-matrix(data = ed[-k,k],nrow = length(dt),ncol = n-1, byrow=TRUE)
    fi<-fi-dt
    
    f<-matrix(data = 0,nrow = length(dt),ncol = n-1, byrow=TRUE)
    good<-which(-e <= fi  & fi<=e )
    f[good]<-(3/(4*e))*(1- (fi[good]/e)^2)
    
    
    Y[,k]<- ((n-1)/(2*pi*dt*a))*rowSums(f*matrix(data= c(ec[-k,k]),ncol = n-1,nrow = length(dt),byrow = TRUE))
    
    
    
    plot(dt,Y[,k],type="l",pch=20,ylim = c(0,300000),xlim=c(0,dt[length(dt)]+0.004),ylab = "pdf",xlab = "r")
    par(new=TRUE)
  }
  Y
}

tc.pdl.ppp(data_2)



