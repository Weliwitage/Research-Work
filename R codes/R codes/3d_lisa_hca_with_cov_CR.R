
####    Author: W.Jithmi Jannidi
####    M.Sc Thesis : Local Indicators of Spatial Association of 2d & 3d Point Patterns

####    R code for Hierarchical cluster analysis (hca) with covariance section 4.3.1.b 3d case
####    simulated point processes ex.4.2.2.a 3d pp 
####    areas covered in this code : edge corrected local product density function, hca


library(MASS)
library(spatstat)
library(Matrix)
#setwd("C:/Users/jimaj/Desktop/Project_2020")
v<-box3(xrange = c(0, 5), yrange = c(0, 5), zrange = c(0, 10), unitname = NULL)

p3d<-read.csv("dp3d.csv")
c13d<-read.csv("d3c1.csv")

x<-c(p3d$V1,c13d$V1)
y<-c(p3d$V2,c13d$V2)
z<-c(p3d$V3,c13d$V3)

data_3d<-pp3(x,y,z,box3=v)
plot(data_3d,type="p",pch=20,main="Simulated 3D point process ")

###


hcpp3<-function(pp){
  v<-box3(xrange = c(0, 5), yrange = c(0, 5), zrange = c(0, 10), unitname = NULL)
  
    
  
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
  
  
  
  #}
  ### let the estimator of lambha<- l1, lambha^2<-l2 and lambha^3<-l3
  
  l1<-(n-2)/a
  l2<-((n-2)*(n-3))/(a^2)
  l3<-((n-2)*(n-3)*(n-4))/(a^3)
  
  
  
  g<-function(x,y,r){
    s<-ed[x,y]-dt[r]
    if (-e <= s  & s<=e ){
      f<-(3/(4*e))*(1- (s/e)^2)
    }
    else{
      f<-0
    }
    
    G<-ec[x,y]*f
    
    return(G)
  }
  
  
  
  ##### while loop gives uniform point pattern which has  at least one data point 
  M<-rpoispp3(0.4,domain=box)
  while (length(M$data$x)==0) {
    M<-rpoispp3(0.4,domain=box)
  }
  
  x<-c(pp$data$x,M$data$x)
  y<-c(pp$data$y,M$data$y)
  z<-c(pp$data$z,M$data$z)
  pm<-pp3(x,y,z,box3=v)
  edpm<-edge.Trans.pp3(pm)
  
  Mec<-edpm[1:length(pp$data$x),length(pp$data$x)+1:length(M$data$x)]
  
  cd<-crossdist.pp3(pp,M) 
  
  ##g function for MC
  
  #CR: use distances precomputed with crossdist 
  Mg<-function(i,m,r){
    
    # dim<-dist(cbind(c(pp[i]$x,M[m]$x),c(pp[i]$y,M[m]$y)))
    s<-cd[i,m]-dt[r]
    if (-e <= s  & s<=e ){
      f<-(3/(4*e))*(1- (s/e)^2)
    }
    else{
      f<-0
    }
    
    G<-Mec[i,m]*f
    
    return(G)
  }
  
  
  Mn<-length(M$data$x)
  # MII =eq 3.2
  
  MII<-function(i,l,t1){
    h<-c()
    for (j in 1:Mn) {
      
      h[j]<-(Mg(i,j,t1)-Mg(l,j,t1))^2
      
    }
    
    im<-(a/Mn)* sum(h)
    return(im)
  }
  
  
  # MIII =eq 3.3
  
  MIII<-function(i,l,t1,t2){
    h<-c()
    for (j in 1:Mn) {
      
      h[j]<- (Mg(i,j,t1)-Mg(l,j,t1))*(Mg(i,j,t2)-Mg(l,j,t2))
      
    }
    
    im<-(a/Mn)* sum(h)
    return(im)
  }
  
  #  factor1<-(l1/a)/(2*pi)^2
  factor2<-(l3+5*l2/a +2*l1/a^2)/(4*pi^2)^2
  
  d<-function(r,c){
    
    d3<-Y[,r]-Y[,c] 
    d1<-t(d3)
    Si<-matrix(data=NA,ncol = length(dt),nrow = length(dt))
    # Si is a symmetric matrix, so only one side off diagonal elements are calculated by for loop
    
    for (i  in 1: length(dt)) {
      for (j in i: length(dt)){
        if(i==j){
          Si[i,j]<-factor2/dt[i]^4*MII(r,c,i)  # CR: Multiplication factor is now precomputed, first summand should disappear as g with translation edge correction is symmetric in the two points (please check) 
        }
        else{
          Si[i,j]<- factor2/(dt[i]^2*dt[j]^2)* MIII(r,c,i,j) #CR: Same as above
          Si[j,i]<-Si[i,j]
        }
      }
    }
    
    Si
    
    d2<-ginv(Si)
    dd<-d1 %*% d2 %*% d3
    return(dd)
  }
  
  
  
  
  D<-matrix(data=0,nrow = n,ncol = n)
  
  ### diagonal elements are zeros and symmetric matrix, only one side off diagonal elements are calculated by for loop
  
  start.time <- Sys.time()
  for (i in 1:(n-1))
  {
    for (j in (i+1):n) {
      
      #  if(i==j)
      # {
      #    D[i,j]<-0
      # }
      #  else{
      D[i,j]<-d(i,j)
      D[j,i]<-D[i,j]
      # }
    }
  }
  end.time <- Sys.time()
  end.time-start.time
  D
  return(D)
  ################hclust 
  
  
  # par(new=FALSE)
  #  hc<-hclust(dist(DM))
  # return(list(DM,plot(hc),plot(pp)))
  
}


#pp<-rpoispp3(lambda=2)
#plot(pp)
start.time <- Sys.time()

matD<-hcpp3(data_3d)

end.time <- Sys.time()
ttt<-end.time - start.time

ttt<-write.csv(ttt,"tI3.csv")
md<-write.csv(matD,"matDI3.csv")
####################################################################
########################## hca




