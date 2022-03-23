
####    Author: W.Jithmi Jannidi
####    M.Sc Thesis : Local Indicators of Spatial Association of 2d & 3d Point Patterns

####    R code for Hierarchical cluster analysis (hca) with covariance section 4.3.1.b 2d case
####    simulated point processes ex.4.2.1.b 2d pp 
####    areas covered in this code : edge corrected local product density function, hca


# install required R packages
library(MASS)
library(spatstat)
library(Matrix)
#setwd("C:/Users/jimaj/Desktop/Project_2020")


area<-owin(c(0,2),c(2,4),unitname = NULL)  # creating a observation window for 2d point pattern
dp_2<-read.csv("dp_2.csv")                  # reading x,y coordinates of simulated point pattern
dcp_2<-read.csv("dcp_2.csv")
x<-c(dp_2$V1,dcp_2$V1)
y<-c(dp_2$V2,dcp_2$V2)
data_2<-ppp(x,y,area)
par(new=FALSE)
plot(data_2,type="p",pch=20,main="Simulated 2D point pattern ")



hcppp<-function(pp){
  
               #### step 1 : generate edge correction weights fro a given point process
  
  W<-pp$window  #info of window
  #  df<-data.frame(cbind(pp$x,pp$y))
  ed<-as.matrix(pairdist.ppp(pp))  # store pairwise distances between points in a matrix
  a<-area.default(pp$window)  # area of the window
  n<-pp$n   #total number of point of the point process
  e<- (sqrt(5)/(10))/(sqrt(intensity(pp)))  #
  dt<-seq(e+0.004,0.25*min(pp$window$xrange[2],pp$window$yrange[2]),length=25)
  
  par(new=FALSE)
  ec<-edge.Trans(pp)
  
            #### step 2: generate edge corrected local product density functions
  
  ## Y - matrix of local edge corrected product density functions, read hca.R file for more details
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
  
  #}
  ### let the estimator of lambha<- l1, lambha^2<-l2 and lambha^3<-l3
  ##  we estimate lamba,lamba^2 , lamba^3 using unbiased estimators see page 43 
  
  l1<-(n-2)/a
 l2<-((n-2)*(n-3))/(a^2)
 l3<-((n-2)*(n-3)*(n-4))/(a^3)

  # function g in covariance equation 4.25
 
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
  
  
  
  ##### M -simulated point pattern  for MC intergration  , while loop gives uniform point pattern which has  at least one data point 
  
  M<-rpoispp(100,win =W)
  while (M$n==0) {
    M<-rpoispp(100,win =W)
  }
  
  plot(M)
  Mec<-edge.Trans(pp,M) # edge correction 
  
  cd<-crossdist(pp,M) #pairwise distances
  
  ##g function for MC in eq.4.28, eq.4.29 Mg
  
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
  
  
  
  # MII =eq 4.28
  
  MII<-function(i,l,t1){
    h<-c()
    for (j in 1:M$n) {
      
      h[j]<-(Mg(i,j,t1)-Mg(l,j,t1))^2
      
    }
    
    im<-(a/M$n)* sum(h)
    return(im)
  }
  
  
  # MIII =eq 4.29
  
  MIII<-function(i,l,t1,t2){
    h<-c()
    for (j in 1:M$n) {
      
      h[j]<- (Mg(i,j,t1)-Mg(l,j,t1))*(Mg(i,j,t2)-Mg(l,j,t2))
      
    }
    
    im<-(a/M$n)* sum(h)
    return(im)
  }
  
  #  factor1<-(l1/a)/(2*pi)^2
  factor2<-(l3+5*l2/a +2*l1/a^2)/(2*pi)^2
  
  # dd is the matrix dd in eq. 4.21, di,j is the  distance between all local product densities of ith and jth points 
  d<-function(r,c){
    
    d3<-Y[,r]-Y[,c] #see eq 4.21
    d1<-t(d3)
    Si<-matrix(data=NA,ncol = length(dt),nrow = length(dt)) #covariance matrix of eq. 4.21
    # Si is a symmetric matrix, so only one side off diagonal elements are calculated by for loop
    
    for (i  in 1: length(dt)) {
      for (j in i: length(dt)){
        if(i==j){
          Si[i,j]<-factor2/dt[i]^2*MII(r,c,i)  # CR: Multiplication factor is now precomputed, first sum should disappear as g with translation edge correction is symmetric in the two points (please check) 
        }
        else{
          Si[i,j]<- factor2/(dt[i]*dt[j])* MIII(r,c,i,j) #CR: Same as above
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
  
  
  
}


##### end

start.time <- Sys.time()

#matD<-hcppp(data_2)

end.time <- Sys.time()
ttt<-end.time - start.time
ttt<-write.csv(ttt,"t2new.csv")
md<-write.csv(matD,"matD2new.csv")
####################################################################


############# check the r code with small point process

area<-owin(c(0,2),c(2,4),unitname = NULL) # crearing observation window for 2d point pattern
pp<-rpoispp(5,win = area,drop = TRUE) #generate random poison pp with intensity 25 as parent process
plot(pp)

hcppp(pp)
