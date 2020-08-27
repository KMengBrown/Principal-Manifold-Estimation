####################################################################
##########          PME, HS, and ISOMAP (d=1, D=3)       ###########
##########               Kun (Michael) Meng              ###########
########## Department of Biostatistics, Brown University ###########
####################################################################

# Packages
library("princurve")


# Simulations
####################################################################

MSD.PME=vector()
MSD.HS=vector()
MSD.ISOMAP=vector()

N=10
for(kkk in 1:N){
  
  print(kkk)
  
  # Data
  I=1000
  manifold=function(t){ return(c(t,cos(t),sin(t))) }
  t=seq(from=0,to=3*pi,length.out = I)
  X=matrix(0,nrow = length(t),ncol = 3)
  for(i in 1:length(t)){ X[i,]=manifold(t[i]) }
  noise=0.05
  e1=rnorm(I,mean=0,sd=noise)
  e2=rnorm(I,mean=0,sd=noise)
  e3=rnorm(I,mean=0,sd=noise)
  X=X+cbind(e1,e2,e3)
  
  # PME
  print("Principal Manifold Estimation")
  ptm <- proc.time()
  result=PME(x.obs=X, d=1)
  proc.time() - ptm  
  MSD.PME[kkk]=min(result$MSD)
  
  # HS
  print("HS")
  fit_smooth=principal_curve(X, thresh = 0, trace = FALSE, maxit = 500, smoother = "smooth_spline")
  fit_lowess=principal_curve(X, thresh = 0, trace = FALSE, maxit = 500, smoother = "lowess")
  fit_periodic=principal_curve(X, thresh = 0, trace = FALSE, maxit = 500, smoother = "periodic_lowess")
  MSD.HS[kkk]=min(c(fit_smooth$dist/I, fit_lowess$dist/I, fit_periodic$dist/I))
  
  # ISOMAP
  print("ISOMAP")
  Y=X
  dissimilarity.matrix=as.matrix(dist(Y))                             # The (i,j)th element of this matrix is the Euclidean 
  isomap.initial=isomap(dissimilarity.matrix,ndim = 1,k=10)           # distrance between mu[i,] and mu[j,].
  tnew=isomap.initial$points      
  D=3
  d=1
  lambda=4-d  
  smoothing.para=0
  I=dim(X)[1]
  theta.hat=rep(1,I)
  W=diag(theta.hat)                                                   # The matrix W
  w=smoothing.para
  I=length(theta.hat)
  T=cbind(rep(1,I),tnew)                                              # The matrix T
  E=matrix(NA,ncol=I,nrow=I)                                          
  for(j in 1:I){
    E.prepare=function(t){ eta.kernel(t-tnew[j,],lambda) }
    E[,j]=apply(tnew,1,E.prepare)                                     # The matrix E
  }
  M1=cbind(2*E%*%W%*%E+2*w*E,2*E%*%W%*%T,T)
  M2=cbind(2*t(T)%*%W%*%E,2*t(T)%*%W%*%T,matrix(0,ncol=d+1,nrow=d+1))
  M3=cbind(t(T),matrix(0,ncol=d+1,nrow=d+1),matrix(0,ncol=d+1,nrow=d+1))
  M=rbind(M1,M2,M3)                                                   # The coefficient matrix of the linear equations
  b=rbind(2*E%*%W%*%Y,2*t(T)%*%W%*%Y,matrix(0,nrow=d+1,ncol=D))       # The nonhomogeneous term of the linear equations
  sol=ginv(M)%*%b                                                     # Solve the linear equations
  eta.func=function(t){
    eta.func.prepare=function(tau){ return(eta.kernel(t-tau,lambda)) }
    return(matrix(apply(tnew,1,eta.func.prepare),ncol=1))
  }
  fnew=function(t){                                                   # fnew() gives the first step from f0().
    return(as.vector(t(sol[1:I,])%*%eta.func(t)+t(sol[(I+1):(I+d+1),])%*%matrix(c(1,t),ncol=1))) 
  }
  S=0
  proj.ind=vector()
  for(i in 1:I){
    #print(as.character(i))
    proj.ind[i]=projection(X[i,], fnew, tnew[i,])
    S=S+dist.euclidean(X[i,], fnew(proj.ind[i]))^2
  }
  MSD.ISOMAP[kkk]=S/I
  
}

MSD.comparison=cbind(MSD.PME, MSD.HS, MSD.ISOMAP)
colnames(MSD.comparison)=c("PME", "HS", "ISOMAP")
MSD.comparison




