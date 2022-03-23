##### This R code is related to simulation of 2d and 3d point patterns used in 
##### M.Sc Thesis study by W.Jithmi Jannidi (Author) 

    ###### By referring this R code, one can learn how to simulate a clustered point process and plot it, in both
    ######    2d and 3d cases using spatstat R package 

library(spatstat) # spatstat- R package for spatial point pattern(spp)

area<-owin(c(0,2),c(2,4),unitname = NULL) # crearing observation window for 2d point pattern
# pp<-rpoispp(25,win = area,drop = TRUE) #generate random poison pp with intensity 25 as parent process
# plot(pp)
# cp<-rpoispp(5000,win =owin(c(1.5,1.6),c(2.5,2.7)),drop = TRUE) #generate random pp with intensity 5000 as feature process which lies inside parent pp window
# plot(cp)

# dp_2<-as.data.frame(cbind(pp$x,pp$y)) #x,y coordinates as a data frame
# dp_2
# write.csv(dp_2,file = "dp_2.csv") #store x,y coordinates in csv format
# dcp_2<-as.data.frame(cbind(cp$x,cp$y))
# dcp_2
# write.csv(dcp_2,file = "dcp_2.csv")


dp_2<-read.csv("dp_2.csv")  #reading coordinates of points from csv files 
dcp_2<-read.csv("dcp_2.csv")



x<-c(dp_2$V1,dcp_2$V1) #x coordinates of both parent and feature patterns
y<-c(dp_2$V2,dcp_2$V2) #y coordinates of both parent and feature patterns
data_2<-ppp(x,y,area) #creating a clustered point pattern by combining parent and feature patterns
par(new=FALSE)
plot(data_2,type="p",pch=20,main="Simulated 2D point pattern ") #plot of clustered point pattern
##########################

v<-box3(xrange = c(0, 5), yrange = c(0, 5), zrange = c(0, 10), unitname = NULL)  # creating volume box for 3d point pattern

# disribution of parent with lambha 0.4
# 
# p3d<-rpoispp3(0.4, domain = v, nsim=1, drop=TRUE)
# p3d
# plot(p3d,type="p",pch=20)
# dp3d<-as.data.frame(cbind(p3d$data$x,p3d$data$y,p3d$data$z))  #x,y,z coordinates of 3d parent point process
# write.csv(dp3d,file = "dp3d.csv")
 

# here 2 cluster centers  chosen by hand 
#c1 (2,3,8)
#c2 (3,2,3)


#select cluster volume as unit 1 and  intensity is 10
#c1<-rpoispp3(10, domain = box3(xrange = c(2-1, 2+1), yrange = c(3-1,3+1), zrange = c(8-1,8+1), unitname = NULL), nsim=1, drop=TRUE)
#c2<-rpoispp3(10, domain = box3(xrange = c(3-1, 3+1), yrange = c(2-1,2+1), zrange = c(3-1,3+1), unitname = NULL), nsim=1, drop=TRUE)

# plot(c1)
# d3c1<-as.data.frame(cbind(c1$data$x,c1$data$y,c1$data$z)) 
# write.csv(d3c1,file = "d3c1.csv")


#d3c2<-as.data.frame(cbind(c2$data$x,c2$data$y,c2$data$z))  #x,y,z coordinates of 3d feature point process
#write.csv(d3c2,file = "d3c2.csv")

p3d<-read.csv("dp3d.csv")
c13d<-read.csv("d3c1.csv")
#d3c2<-read.csv("d3c2.csv")

x<-c(p3d$V1,c13d$V1)
y<-c(p3d$V2,c13d$V2)
z<-c(p3d$V3,c13d$V3)

data_3d<-pp3(x,y,z,box3=v)
plot(data_3d,type="p",pch=20,main="Simulated 3D point process ")   #simulated 3d point process




