####################################################################
##########              A Simulation Study of            ###########
##########    Principal Manifolds Estimation Algorithm   ###########
##########               Kun (Michael) Meng              ###########
########## Department of Biostatistics, Brown University ###########
####################################################################

# install.packages("plot3D")
library("plot3D")

################### Section 1, d=1 and D=2 #########################
####################################################################

### Case I.

# Simulated data
I=1000                              # Sample size
t=rnorm(I,mean = 0,sd=1)
sd.noise=0.15                       # Standard deviation of the noise 
e1=rnorm(I,mean = 0,sd=sd.noise)    # noise
e2=rnorm(I,mean = 0,sd=sd.noise)
X=matrix(NA,ncol = 2,nrow = I)
manifold=function(tau){ return(c(tau,sin(tau+pi/2))) }
for(i in 1:I){ X[i,]=manifold(t[i]) }
data.points=X+cbind(e1,e2) 

# Apply PME()
ptm <- proc.time()
result=PME(x.obs=data.points, d=1)
proc.time() - ptm

# Plot the fitted manifold
f=result$embedding.map
t.test=seq(from=-100,to=100,by=0.05)
t.length=length(t.test)
x.test=matrix(NA,ncol=2,nrow = t.length)
for(i in 1:t.length){ x.test[i,]=f(t.test[i]) }
plot(data.points[,1], data.points[,2],
     pch=20, col="grey",
     main="Principal Manifold Estimation",
     xlab=" ", ylab=" ")
lines(x.test[,1],x.test[,2],col="red",type = "l",lwd=3)


### Case II

manifold=function(t){ return(c(cos(t),sin(t)))}
I=1000
t=runif(I,min = 0,max = 1.5*pi)
X=manifold(t)
sd.noise=0.1
e1=rnorm(I,mean = 0, sd=sd.noise)
e2=rnorm(I,mean = 0, sd=sd.noise)
data.points=X+cbind(e1,e2)  

ptm <- proc.time()
result=PME(x.obs=data.points, d=1)
proc.time() - ptm

f=result$embedding.map
t.test=seq(from=-100,to=100,by=0.05)
t.length=length(t.test)
x.test=matrix(NA,ncol=2,nrow = t.length)
for(i in 1:t.length){ x.test[i,]=f(t.test[i]) }
plot(data.points[,1], data.points[,2],
     pch=20, col="grey",
     main="Principal Manifold Estimation",
     xlab=" ", ylab=" ")
lines(x.test[,1],x.test[,2],col="red",type = "l",lwd=3)


### Case III

I=1000                             
t=runif(I, min = -3*pi, max = 3*pi)
sd.noise=0.1                       
e1=rnorm(I,mean = 0,sd=sd.noise)    
e2=rnorm(I,mean = 0,sd=sd.noise)
X=matrix(NA,ncol = 2,nrow = I)
manifold=function(tau){ return(c(tau,sin(tau))) }
for(i in 1:I){ X[i,]=manifold(t[i]) }
data.points=X+cbind(e1,e2) 

ptm <- proc.time()
result=PME(x.obs=data.points, d=1)
proc.time() - ptm

f=result$embedding.map
t.test=seq(from=-100,to=100,by=0.05)
t.length=length(t.test)
x.test=matrix(NA,ncol=2,nrow = t.length)
for(i in 1:t.length){ x.test[i,]=f(t.test[i]) }
plot(data.points[,1], data.points[,2],
     pch=20, col="grey",
     main="Principal Manifold Estimation",
     xlab=" ", ylab=" ")
lines(x.test[,1],x.test[,2],col="red",type = "l",lwd=3)



################### Section 2, d=1 and D=3 #########################
####################################################################

### Case I

I=1000
manifold=function(t){ return(c(t,t^2,t^3)) }
t=seq(from=-1,to=1,length.out = I)
X=matrix(0,nrow = length(t),ncol = 3)
for(i in 1:length(t)){ X[i,]=manifold(t[i]) }
noise=0.1
e1=rnorm(I,mean=0,sd=noise)
e2=rnorm(I,mean=0,sd=noise)
e3=rnorm(I,mean=0,sd=noise)
data.points=X+cbind(e1,e2,e3)

ptm <- proc.time()
result=PME(x.obs=data.points, d=1)
proc.time() - ptm

f=result$embedding.map
t.test=seq(from=-10,to=10,by=0.005)
t.length=length(t.test)
x.test=matrix(NA,ncol=3,nrow = t.length)
for(i in 1:t.length){ x.test[i,]=f(t.test[i]) }
index=(x.test[,1]>=min(data.points[,1]))&x.test[,1]<=max(data.points[,1])&(x.test[,2]>=min(data.points[,2]))&x.test[,2]<=max(data.points[,2])&(x.test[,3]>=min(data.points[,3]))&x.test[,3]<=max(data.points[,3])
scatter3D(data.points[,1], data.points[,2], data.points[,3], 
          pch = 20, box=TRUE, cex = 0.2, colkey = FALSE, 
          border="black", shade=0.8, 
          bty = "g", ticktype = "detailed",
          main="Principal Manifold Estimation")
scatter3D(x.test[index,1], x.test[index,2],x.test[index,3], 
          pch = 20, box=TRUE, cex = 0.2, colkey = FALSE, col = "red", 
          border="black", shade=0.8, main=" ",add = TRUE)


### Case II

I=1000
manifold=function(t){ return(c(t,cos(t),sin(t))) }
t=seq(from=0,to=3*pi,length.out = I)
X=matrix(0,nrow = length(t),ncol = 3)
for(i in 1:length(t)){ X[i,]=manifold(t[i]) }
noise=0.05
e1=rnorm(I,mean=0,sd=noise)
e2=rnorm(I,mean=0,sd=noise)
e3=rnorm(I,mean=0,sd=noise)
data.points=X+cbind(e1,e2,e3)

ptm <- proc.time()
result=PME(x.obs=data.points, d=1)
proc.time() - ptm

f=result$embedding.map
t.test=seq(from=-20,to=20,by=0.01)
t.length=length(t.test)
x.test=matrix(NA,ncol=3,nrow = t.length)
for(i in 1:t.length){ x.test[i,]=f(t.test[i]) }
index=(x.test[,1]>=min(data.points[,1]))&x.test[,1]<=max(data.points[,1])&(x.test[,2]>=min(data.points[,2]))&x.test[,2]<=max(data.points[,2])&(x.test[,3]>=min(data.points[,3]))&x.test[,3]<=max(data.points[,3])
scatter3D(data.points[,1], data.points[,2], data.points[,3], 
          pch = 20, box=TRUE, cex = 0.2, colkey = FALSE, 
          border="black", shade=0.8, 
          bty = "g", ticktype = "detailed",
          main="Principal Manifold Estimation")
scatter3D(x.test[index,1], x.test[index,2],x.test[index,3], 
          pch = 20, box=TRUE, cex = 0.2, colkey = FALSE, col = "red", 
          border="black", shade=0.8, main=" ",add = TRUE)



################### Section 3, d=2 and D=3 #########################
####################################################################

### Case I

manifold=function(t){ return(c(t[1],t[2],norm.euclidean(t)^2)) }
I=1000
noise=0.05
e1=rnorm(I,mean=0,sd=noise)
e2=rnorm(I,mean=0,sd=noise)
e3=rnorm(I,mean=0,sd=noise)
input.size=1
t.1=runif(I,min = -input.size,max=input.size)
t.2=runif(I,min = -input.size,max=input.size)
t=cbind(t.1,t.2)
data.points=t(apply(t,1,manifold))+cbind(e1,e2,e3)

ptm <- proc.time()
result=PME(x.obs=data.points, d=2)
proc.time() - ptm

f=result$embedding.map
t.plot.1=seq(from=-2, to=2, length.out = 200)
t.plot.2=seq(from=-2, to=2, length.out = 200)
t.length=length(t.plot.1)
surf.plot=matrix(0,ncol=3,nrow=1)
for(i in 1:t.length){
  for(j in 1:t.length){
    surf.plot=rbind(surf.plot,f(c(t.plot.1[i],t.plot.2[j])))
  }
}
x.test=surf.plot[-1,]
index=(x.test[,1]>=min(data.points[,1]))&x.test[,1]<=max(data.points[,1])&(x.test[,2]>=min(data.points[,2]))&x.test[,2]<=max(data.points[,2])&(x.test[,3]>=min(data.points[,3]))&x.test[,3]<=max(data.points[,3])
scatter3D(data.points[,1], data.points[,2], data.points[,3], 
          pch = 20, box=FALSE, cex = 0.5, colkey = FALSE, 
          border="black", shade=0.8, 
          ticktype = "detailed",
          main="Principal Manifold Estimation")
scatter3D(x.test[index,1], x.test[index,2],x.test[index,3], 
          pch = 20, box=FALSE, cex = 0.1, colkey = FALSE, col = "grey", 
          border="black", shade=0.8, main=" ",add = TRUE)



