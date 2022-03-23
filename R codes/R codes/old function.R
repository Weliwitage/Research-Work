

#old function
#dOld<-function(r,c){

#  d1<-t(Y[,r]-Y[,c])
#  d3<-Y[,r]-Y[,c]

#  Si<-matrix(data=NA,ncol = length(dt),nrow = length(dt))
# Si is a symmetric matrix, so only one side off diagonal elements are calculated by for loop
#  for (i  in 1: length(dt)) {

#    for (j in i: length(dt)){
#      if(i==j){
#        Si[i,j]<- ( (l1/a)/(2*pi*dt[i])^2)*(g(r,c,i)-g(c,r,i))^2 + ((l3+5*l2/a +2*l1/a^2)/(2*pi*dt[i]^2))*MII(r,c,i)
#      }
#     else{
#        Si[i,j]<- ((l1/a)/((2*pi)^2*dt[i]*dt[j]))*((g(r,c,i)-g(c,r,i))*(g(r,c,j)-g(c,r,j))) +((l3+5*l2/a +2*l1/a^2)/((2*pi)^2 *dt[i]*dt[j]))* MIII(r,c,i,j)
#       Si[j,i]<-Si[i,j]
#      }
#    }
#  }

#  Si

#    d2<-ginv(Si)
#   dd<-d1 %*% d2 %*% d3
#    return(dd)
# }

