
####    Author: W.Jithmi Jannidi
####    M.Sc Thesis : Local Indicators of Spatial Association of 2d & 3d Point Patterns

####    R code for nnclean function /EM method for clustering of pp
####    examples of simulated point processes ex.4.2.1.b 2d pp and ex.4.2.2.a 3d pp

library(spatstat)
#2d pp

area<-owin(c(0,2),c(2,4),unitname = NULL)
dp_2<-read.csv("dp_2.csv")
dcp_2<-read.csv("dcp_2.csv")
x<-c(dp_2$V1,dcp_2$V1)
y<-c(dp_2$V2,dcp_2$V2)
data_2<-ppp(x,y,area)
plot(data_2,type="p",pch=20,main="Simulated 2D point pattern ")


n2<-nnclean(data_2,k=5)
d2<-cbind(n2$x,n2$y,n2$marks$class)
noise2<-subset(d2,d2[,3]==1)
fea<-subset(d2,d2[,3]==2)
fea
t<-cbind(c(rep(1,105),rep(2,98)),d2[,3])
subset(t,t[,1]==2 & t[,2]==2)
#plot(ppp(noise2[,1],noise2[,2],area),type="p",pch=20)
plot(ppp(fea[,1],fea[,2],area),type="p",pch=20,main="k=5")

par(new=FALSE)
par(mfrow=c(2,2))
par(mfrow=c(1,1))



### 3d  0.4,10

v<-box3(xrange = c(0, 5), yrange = c(0, 5), zrange = c(0, 10), unitname = NULL)

p3d<-read.csv("dp3d.csv")
c13d<-read.csv("d3c1.csv")

x<-c(p3d$V1,c13d$V1)
y<-c(p3d$V2,c13d$V2)
z<-c(p3d$V3,c13d$V3)

data_3d<-pp3(x,y,z,box3=v)
plot(data_3d,type="p",pch=20,main="Simulated 3D point process ")

n3<-nnclean(data_3d,k=5)
d3<-cbind(n3$data$x,n3$data$y,n3$data$z, n3$data$label)
subset(d3,d3[,4]==2)
tt<-cbind(c(rep(1,92),rep(2,106)),d3[,4])
subset(tt,tt[,1]==2 & tt[,2]==2)


noise3<-subset(d3,d3[,4]==1)
fea3<-subset(d3,d3[,4]==2)
#plot(ppp(noise3[,1],noise3[,2],area),type="p",pch=20)
plot(pp3(fea3[,1],fea3[,2],fea3[,3],box3=v),type="p",pch=20,main="k=5")

