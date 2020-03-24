#########################################################################
##########    Principal Manifold Estimation (PME) Algorithm   ###########
##########                  Kun (Michael) Meng                ###########
##########    Department of Biostatistics, Brown University   ###########
#########################################################################


########### Introduction ###########################################
####################################################################

# We propose an R-software function for d-dimensional estimating principle manifolds, 
# where d = 1, 2, 3, and manifolds embedded into D-dimensional Euclidean space with D 
# strictly larger than d.

# The proposed function is "PME()." To construct this function, we need some functions 
# for "high dimensional mixture density estimation (HDMDE)."

# This manuscript - whatever its name is - is organized as follows.
#   (i)   In Section 1, we define some basic functions.
#   (ii)  In Section 2, we construct a function for HDMDE. In this section, we first construct 
#         the "weight.seq" function. Then we apply this function to construct "hdmde" function.
#   (iii) Based on this "hdmde" function, "PME" estimation function is finally given in Section 3.
# The theoretical foundation of this R-software function is the following paper.
#
# K. Meng and A. Eloyan, Principal Manifolds: A framework Using Sobolev Spaces and Model
#                        Complexity Selection Using Mixture Densities

########### Section 0, Packages ####################################
####################################################################

#install.packages("MASS")
library("MASS")
#install.packages("Matrix")
library("Matrix")
#install.packages("vegan")
library("vegan")
#install.packages("plot3D")
library("plot3D")



########### Section 1, Some Basic Functions ########################
####################################################################

## Subsection 1.1, Functions for Euclidean metrics

# Norm function in an Euclidean space of any dimension
norm.euclidean=function(x){ return(norm(matrix(x,ncol=1),type = "F")) }
# Distance function in an Euclidean space of any dimension
dist.euclidean=function(x,y){ return(norm.euclidean(x-y)) }
## Subsection 1.2, Kernels for minimization in a semi-normed space of Sobolev type


## Subsection 1.2 Kernels

# Smoothing kernel for density estimation 
# (We applied Gaussian kernel.)
ker=function(x,mu,sigma){
  yseq <- sapply((x-mu)/sigma,dnorm)
  return((sigma^{-length(x)})*prod(yseq))
}

# Reproducing Kernels associated with Sobolev space D^{-2}L^2(R^d)
eta.kernel=function(t,lambda){
  if(lambda%%2==0){
    if(norm.euclidean(t)==0){
      y=0
    }else{
      y=(norm.euclidean(t)^lambda)*log(norm.euclidean(t))
    }
  }else{
    y=norm.euclidean(t)^lambda
  }
  return(y)
}


## Subsection 1.3, Projection Index function
projection=function(x,f,initial.guess){
  DD=function(t){ return(dist.euclidean(x,f(t))) }
  est=nlm(DD,p=initial.guess)
  return(est$estimate)
}




##### Section 2, High Dimensional Mixture Density Estimation #######
####################################################################

## Subsection 2.1 
# When \mu's and \sigma are given, the following function estimates \hat{\theta}'s.

weight.seq=function(x.obs, mu, sigma, epsilon=0.001, max.iter=1000){ 
  
  # "x.obs" is the data set of interest. 
  #         There are n observations, and each observation is a D-dimensional point.
  #         x.obs is a n-by-D matrix.
  # "mu" is a vector of the knots in a mixture density estimation.
  # "sigma" is the bandwidth of this density estimation.
  # "epsilon" is a predetermined tolerance of the Euclidean distance between thetas in two consecutive steps.
  # "max.iter" is a predetermined upper bound of the number of steps in this iteration.
  
  n.D=dim(x.obs)                  
  n=n.D[1]
  D=n.D[2]
  N=dim(mu)[1]                     
  
  A=matrix(NA,ncol = N,nrow = n)
  for(j in 1:N){
    A.prepare=function(x){ return(ker(x,mu[j,],sigma))  }
    A[,j]=apply(x.obs,1,A.prepare)   # A[i,j] is \psi_sigma (x_i-mu_j).
  }
  
  
  theta.old=rep(1/N,N)               # The initial guess for weights, say \theta_j's.
  abs.diff=10*epsilon                # Absolute value of the difference between "theta.new" and "theta.old".
  count=0                            # Counting the number of steps in this iteration.
  lambda.hat.old=c(n,rep(-1,D))      # The initial guess of the Lagrangian multipliers
  
  while((abs.diff>epsilon)&(count<=max.iter)){   # The iteration for computing desired theta's
    
    W=t(t(A)*theta.old)              # \theta_j^{(k)} \times \psi_\sigma(x_i-mu_j)
    W=W/apply(W,1,sum)               # W[i,j] is the posterior probability of Z=j|X=x_i, say w_{i,j}(\theta.old).
    w=apply(W,2,sum)                 # w[j] = \sum_{i=1}^n w_{i,j}
    
    flambda=function(lambda){        # This function is for computing Lagrangian multipliers.
      
      denom.temp=apply(t(t(cbind(rep(1,dim(mu)[1]),mu))*lambda),1,sum)        # The denominator sequence: 
      # \lambda_1+\lambda_2^T \mu_j, j=1,2,...,N.
      num.temp=mu*w
      
      f1<-sum(w/denom.temp)          # \sum_{j=1}^N \frac{ w_ij }{ \lambda_1+\lambda_2^T \mu_j }
      f2=apply(num.temp*(as.vector(1/denom.temp)),2,sum)
      f=dist.euclidean(f1,1)+dist.euclidean(f2,apply(x.obs,2,mean))
      return(f) 
      
    }
    lambda.hat=nlm(flambda,lambda.hat.old,iterlim=1000)$estimate              # The lagrangian multipliers.
    # We set the Lagrangian multipliers in the previous step 
    # as the initial guess in this step.
    theta.new=w/(apply(t(t(cbind(rep(1,dim(mu)[1]),mu))*lambda.hat),1,sum))   # The new theta's computed from the old theta's
    abs.diff=dist.euclidean(theta.new,theta.old)                              # The Euclidean distance between the old and new theta vectors.
    if(is.na(abs.diff)){ abs.diff=0; theta.new=theta.old }                    # It helps us avoid "NA trouble".
    
    theta.old=pmax(theta.new,0)      # Set the new theta as the old theta for the next iteration step.   
    theta.old=pmin(theta.old,1)      # pmax() and pmin() guarantee that theta_j's are in [0,1]. 
    count=count+1
    lambda.hat.old=lambda.hat
  }
  
  theta.hat=pmax(theta.new,0)        # The iteration has stopped.  
  theta.hat=pmin(theta.hat,1)        # Again, we guarantee that theta_j's are in [0,1].
  return(theta.hat)                  # It returns the estimation of weights \theta_j's.
}



## Subsection 2.2
# Based on the function weight.seq() in Subsection 2.1, we give the following
# high dimensional mixture density estiamtion function.

hdmde=function(x.obs, N0, alpha, max.comp){
  
  # "x.obs" is the data set of interest. 
  #         There are n observations, and each observation is a D-dimensional point.
  #         x.obs is a n-by-D matrix.
  # "N0" is a predetermined lower bound for N - the number of density components
  # "alpha" is the predetermined confidence level.
  # "max.comp" is the upper bound of the number of components in the desired mixture density.
  
  zalpha=qnorm(1-alpha/2)  
  n.D=dim(x.obs)                      # "x.obs" is a n by D matrix. Each row of it denotes a data point in R^D.
  n=n.D[1]
  D=n.D[2]
  N=N0
  
  km=kmeans(x.obs, N, nstart = 100)   # Use k-means clustering to get N clusters of x.obs.
  mu=km$centers                       # Set the centers of the N clusters as the means of the density components.
  
  sigma.vec=rep(NA,N)                 # The following block estimates \sigma_N.
  for(j in 1:N){                      
    index.temp=which(km$cluster==j)   
    xi.j=x.obs[index.temp,]           
    sig.prepare=function(x){ return((dist.euclidean(x,mu[j,]))^2) }
    s=apply(xi.j,1,sig.prepare)
    sigma.vec[j]=mean(s)
  }
  sig=sqrt(mean(sigma.vec)/dim(x.obs)[2])
  
  theta.hat=weight.seq(x.obs,mu,sig)  # It gives an estimation of theta_j's with weight.seq().
  
  # The following block gives an approximation to the underlying density function of interest.
  # This estimation is of the form of weights times scaled kernels.
  f.test=function(x){
    fun.prepare=function(t){ return(ker(x,t,sig)) }
    comp.vec=apply(mu,1,fun.prepare)
    return(sum(theta.hat*comp.vec))
  }
  
  p.old=apply(x.obs,1,f.test)         # The second/negative term p_N in Delta.hat.
  
  test.rejection=1
  while((test.rejection==1)&(N<=min(n,max.comp))){
    
    N <- N+1
    
    ##################################################
    # The following is a repetition of the codes above.
    km=kmeans(x.obs,N, nstart = 100)              
    mu=km$centers                   
    
    sigma.vec=rep(NA,N)              
    for(j in 1:N){                      
      index.temp=which(km$cluster==j)   
      xi.j=matrix(x.obs[index.temp,],nrow = length(index.temp))        
      sig.prepare=function(x){ return(dist.euclidean(x,mu[j,])^2) }
      s=apply(xi.j,1,sig.prepare)
      sigma.vec[j]=mean(s)
    }
    sig=sqrt(mean(sigma.vec)/dim(x.obs)[2])
    
    theta.hat=weight.seq(x.obs,mu,sig) 
    
    f.test=function(x){
      fun.prepare=function(t){ return(ker(x,t,sig)) }
      comp.vec=apply(mu,1,fun.prepare)
      return(sum(theta.hat*comp.vec))
    }
    # The part above is a repetition.
    ##################################################
    
    p.new=apply(x.obs,1,f.test)        # The first/positive term p_{N+1} in Delta.hat.
    delta.hat=p.new-p.old              # Delta.hat
    sigma.hat.sq=mean((delta.hat-mean(delta.hat))^2)
    Z.I.N=sqrt(dim(x.obs)[1])*mean(delta.hat)/sqrt(sigma.hat.sq)
    
    if ((Z.I.N<=zalpha) & (Z.I.N>=-zalpha) & (!is.na(Z.I.N))) { test.rejection=0 }
    p.old=p.new
  }
  
  f=f.test
  
  resp=list(estimating.pdf=f,theta.hat=theta.hat,mu=mu,N=N,k.means.result=km,sigma=sig,Z.I.N=Z.I.N)
  
  # Output:
  # "estimating.pdf" is the derived mixture density approximating an underlying desnity.
  # "theta.hat" is a vector of weights for knots of this mixture density.
  # "mu" is a vector of knots of this mixture density.
  # "N" is the number of knots of this mixture density.
  # "sigma" is the variance shared by the components of this mixture density. 
  # "Z.I.N" is the statistic determining the size of N.
  # "k.means.result" gives the result list of "kmeans(obs.x,N)".
  
  return(resp)
}




############ Section 3, Principal Manifold Estimation ######################
############################################################################

PME=function(x.obs, d, N0=0, tuning.para.seq=exp((-15:5)), alpha=0.05, max.comp=100, epsilon=0.05, max.iter=100, print.MSDs=TRUE){
  
  # "x.obs" is the data set of interest. 
  #         There are n observations, and each observation is a D-dimensional point.
  #         x.obs is a n-by-D matrix.
  # "d" is the intrinsic dimension of the underlying manifold
  # "N0" is a predetermined lower bound for N - the number of density components, default value is 30*d
  # "tuning.para.seq" is a vector of tuning parameter candidates, its default value is exp((-15:5)).
  #                   If you would like to fit a manifold for a specific lambda, set tuning.prar.seq=c(lambda).
  # "alpha" is the pre-determined confidence level, which determines the number of the components in a mixture density.
  # "max.comp" is the upper bound of the number of components in the desired mixture density.
  # "epsilon" is the tolerance for distance between manifolds f.old and f.new.
  # "max.iter" is the upper bound of the number of steps in this iteration.
  
  # Remark: The larger N0 is, the less time consuming the function is.
  
  
  dimension.size=dim(x.obs)
  D=dimension.size[2]                                                   # "D" is the dimension of the input space.
  n=dimension.size[1]                                                   # "n" is the number of observed D-dimensional data points.
  lambda=4-d                                                            # "lambda" determines the form of reproducing kernels
  
  if(N0==0){ N0=20*D }
  
  est=hdmde(x.obs, N0, alpha, max.comp)                                 # "hdmde" gives \hat{Q}_N.
  theta.hat=est$theta.hat  
  centers=est$mu
  sigma=est$sigma
  W=diag(theta.hat)                                                     # The matrix W
  X=est$mu
  I=length(theta.hat)
  
  dissimilarity.matrix=as.matrix(dist(X))                               # The (i,j)th element of this matrix is the Euclidean 
  isomap.initial=isomap(dissimilarity.matrix,ndim = d, k=10)            # distrance between mu[i,] and mu[j,].
  t.initial=isomap.initial$points                                       # Give the initial projection indices by ISOMAP.
  
  MSE.seq=vector()
  SOL=list()
  TNEW=list()
  
  for (tuning.ind in 1:length(tuning.para.seq)) {
    
    print(paste("The tuning parameter is lambda[", as.character(tuning.ind), "] = ", as.character(tuning.para.seq[tuning.ind]), "."))
    
    w=tuning.para.seq[tuning.ind]
    tnew=t.initial
    T=cbind(rep(1,I),tnew)                                             # The matrix T
    
    E=matrix(NA,ncol=I,nrow=I)                                          
    for(j in 1:I){
      E.prepare=function(t){ eta.kernel(t-tnew[j,],lambda) }
      E[,j]=apply(tnew,1,E.prepare)                                     # The matrix E
    }
    
    # This block gives the first step of iteration.
    ###############################################
    M1=cbind(2*E%*%W%*%E+2*w*E,2*E%*%W%*%T,T)
    M2=cbind(2*t(T)%*%W%*%E,2*t(T)%*%W%*%T,matrix(0,ncol=d+1,nrow=d+1))
    M3=cbind(t(T),matrix(0,ncol=d+1,nrow=d+1),matrix(0,ncol=d+1,nrow=d+1))
    M=rbind(M1,M2,M3)                                                   # The coefficient matrix of the linear equations
    
    b=rbind(2*E%*%W%*%X,2*t(T)%*%W%*%X,matrix(0,nrow=d+1,ncol=D))       # The nonhomogeneous term of the linear equations
    sol=ginv(M)%*%b                                                     # Solve the linear equations
    
    eta.func=function(t){
      eta.func.prepare=function(tau){ return(eta.kernel(t-tau,lambda)) }
      return(matrix(apply(tnew,1,eta.func.prepare),ncol=1))
    }
    
    fnew=function(t){                                                    
      return(as.vector(t(sol[1:I,])%*%eta.func(t)+t(sol[(I+1):(I+d+1),])%*%matrix(c(1,t),ncol=1))) 
    }
    
    # ISOMAP gives the initial projection indices. Then the initial projection indices give the initial manifold f0.
    # The new projection indices are derived by projecting mu_j's onto f0. 
    
    f0=fnew
    
    X.initial.guess=cbind(X,tnew)                                         # The "tnew" here is derived from ISOMAP.
    projection.index.f0=function(x.init){ projection(x.init[1:D],f0,x.init[(D+1):(D+d)]) }  
    # The first D columns of x.init corresponds to X and the last d columns corresponds to tnew.
    # projection() is applied to X[j,] with initial guess tnew[j,], which is the projection index for X[j,] onto the old manifold f0.
    tnew=matrix(t(apply(X.initial.guess,1,projection.index.f0)),nrow=I)   # This method can help us avoid the chaos from improper initial guess.
    
    # Sum of the squared distances between x_i and its projection onto manifold f.
    # The first D columns of "x.prin" corresponds to points in the input space
    # and the last d columns of "x.prin" corresponds to the projection indices of these points onto f.
    SSD.prepare=function(x.prin,f){ return(dist.euclidean(x.prin[1:D],f(x.prin[(D+1):(D+d)]))^2) }
    
    X.projection.index=cbind(X,tnew)                                      # "tnew" here is the projection index onto fnew, rather than f0. 
    SSD.prepare.again=function(x.init){ return(SSD.prepare(x.init,fnew)) }
    SSD.new=sum(as.vector(apply(X.projection.index,1,SSD.prepare.again))) # "SSD" stands for "sum of squared distances."
    # This block gives the first step of iteration.
    ###############################################
    
    # The iteration for PME is given by the following loop.
    count=1
    SSD.ratio=10*epsilon                                                  # A quantity measuring the distance between f0 and fnew.
    
    while((SSD.ratio > epsilon)&(SSD.ratio<=5)&(count<=(max.iter-1))){
      
      SSD.old=SSD.new
      f0=fnew                                                             # Set the new manifold in the previous step as the old manifold in the next step.
      
      # Repetition 
      #############
      T=cbind(rep(1,I),tnew)                                             
      
      E=matrix(NA,ncol=I,nrow=I)                                          
      for(j in 1:I){
        E.prepare=function(t){ eta.kernel(t-tnew[j,],lambda) }
        E[,j]=apply(tnew,1,E.prepare)                                     
      }
      
      # This block gives the first step of iteration.
      ###############################################
      M1=cbind(2*E%*%W%*%E+2*w*E,2*E%*%W%*%T,T)
      M2=cbind(2*t(T)%*%W%*%E,2*t(T)%*%W%*%T,matrix(0,ncol=d+1,nrow=d+1))
      M3=cbind(t(T),matrix(0,ncol=d+1,nrow=d+1),matrix(0,ncol=d+1,nrow=d+1))
      M=rbind(M1,M2,M3)                      
      
      b=rbind(2*E%*%W%*%X,2*t(T)%*%W%*%X,matrix(0,nrow=d+1,ncol=D))      
      sol=ginv(M)%*%b                                               
      
      eta.func=function(t){
        eta.func.prepare=function(tau){ return(eta.kernel(t-tau,lambda)) }
        return(matrix(apply(tnew,1,eta.func.prepare),ncol=1))
      }
      
      fnew=function(t){                                                   
        return(as.vector(t(sol[1:I,])%*%eta.func(t)+t(sol[(I+1):(I+d+1),])%*%matrix(c(1,t),ncol=1))) 
      }
      
      
      t.old <- tnew
      # Repetition 
      #############
      
      X.initial.guess=cbind(X,tnew)                                       # The "tnew" here is the projection index to f0.
      projection.index.f0=function(x.init){ projection(x.init[1:D],f0,x.init[(D+1):(D+d)]) }                                                                
      tnew=matrix(t(apply(X.initial.guess,1,projection.index.f0)),nrow=I) # The "tnew" here is the projection index to fnew, rather than f0.
      
      X.projection.index=cbind(X,tnew)
      SSD.prepare.again=function(x.init){ return(SSD.prepare(x.init,fnew)) }
      SSD.new=sum(as.vector(apply(X.projection.index,1,SSD.prepare.again)))
      
      SSD.ratio=abs(SSD.new-SSD.old)/SSD.old
      count <- count+1
      
      print(paste("SSD.ratio is ",as.character(round(SSD.ratio, 4))," and this is the ",as.character(count),"th step of iteration."))
    }
    
    
    # For a fixed tuning parameter value, the corresponding MSD is computed by the following chunk.
    km=est$k.means.result
    data.initial=matrix(0,nrow = 1, ncol = D+d)
    for(i in 1:I){
      index.temp=which(km$cluster==i)
      length.temp=length(index.temp)
      X.i=x.obs[index.temp,]
      t.temp=matrix(rep(tnew[i,1],length.temp))
      for(j in 1:d){ t.temp=cbind(t.temp,rep(tnew[i,j],length.temp)) }
      t.temp=matrix(t.temp[,-1],nrow = length.temp)
      data.initial=rbind(data.initial,cbind(X.i,t.temp))
    }
    data.initial=data.initial[-1,]
    proj.para.prepare=function(data.init){ return(projection(data.init[1:D],fnew,data.init[(D+1):(D+d)])) }
    proj.para=matrix(t(apply(data.initial,1,proj.para.prepare)),ncol = d)
    proj.points=t(apply(proj.para,1,fnew))
    diff.data.fit=apply(data.initial[,1:D]-proj.points,1,norm.euclidean)
    MSE=mean(diff.data.fit^2)
    
    MSE.seq[tuning.ind]=MSE
    print(paste("When lambda = ", as.character(w), ", MSD = ", as.character(MSE), "."))
    SOL[[tuning.ind]]=sol
    TNEW[[tuning.ind]]=tnew
    
    # To reduce the computational burden, if the MSD in the k-th step of the for-loop is
    # smaller than that in the next 4 steps of this for-loop (k+1, k+2, k+3, k+4), 
    # we stop this for-loop. 
    if(tuning.ind>=4){
      if((MSE.seq[tuning.ind]>MSE.seq[tuning.ind-1])&(MSE.seq[tuning.ind-1]>MSE.seq[tuning.ind-2])&(MSE.seq[tuning.ind-2]>MSE.seq[tuning.ind-3])){break}
    }
    
  }
  
  # The following chunk gives the f_\lambda with the optimal \lambda.
  optimal.ind=min(which(MSE.seq==min(MSE.seq)))
  sol.opt=SOL[[optimal.ind]]
  tnew.opt=TNEW[[optimal.ind]]
  eta.func=function(t){
    eta.func.prepare=function(tau){ return(eta.kernel(t-tau,lambda)) }
    return(matrix(apply(tnew.opt,1,eta.func.prepare),ncol=1))
  }
  f.optimal=function(t){                                                   
    return(as.vector(t(sol.opt[1:I,])%*%eta.func(t)+t(sol.opt[(I+1):(I+d+1),])%*%matrix(c(1,t),ncol=1))) 
  }
  
  
  if(print.MSDs==TRUE){
    plot(log(tuning.para.seq[1:length(MSE.seq)]), MSE.seq, 
         xlab = "Log Lambda", ylab = "MSD", 
         type = "l")
    lines(log(tuning.para.seq[1:length(MSE.seq)]), MSE.seq, 
          type = "p", pch=20, col="orange", cex=2)
    abline(v=log(tuning.para.seq[optimal.ind]), lwd=1.5, col="darkgreen", lty=2)
    
    print(paste("The optimal tuning parameter is ", 
                as.character(tuning.para.seq[optimal.ind]), 
                ", and the MSD of the optimal fit is ",
                as.character(MSE.seq[optimal.ind]), "."))
  }
  
  resp=list(embedding.map=f.optimal, 
            MSD=MSE.seq,  
            knots=centers,
            weights.of.knots=theta.hat,
            coe.kernel=sol.opt[1:I,],
            coe.poly=sol.opt[(I+1):(I+d+1),],
            SOL=SOL,
            TNEW=TNEW,
            T.parameter=sol.opt,
            Coef=tnew.opt)
  return(resp)
  
  # Output:
  # embedding.map: The fitted embedding map R^d -> R^D. 
  # MSD: A vector of mean squared distances. 
  #      Each component of this vector corresponds to a tuning parameter candidate.
  # knots: Knots in the discrete measure \hat{Q}_N
  # weights.of.knots: A vector of weights for the knots of \hat{Q}_N. 
  #                   The k-th component is the weight for the k-th knot.
  # Lists T.parameter and Coef are quantities determining the analytic formula of f.optimal.
  # coe.kernel and coe.poly are quantities for the Interior Identification function.
  
}

##################################################################
## Section 3 completes the principal manifold estiamtion function.
##################################################################
