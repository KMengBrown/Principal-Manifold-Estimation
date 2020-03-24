# Principal-Manifold-Estimation

Description: This repository is devoted to the principal manifold estimation (PME) algorithm. The PME algorithm comes from the paper "Principal Manifolds: A Framework and Model Complexity Selection" by Kun Meng and Ani Eloyan (hereafter ME). This paper is available at https://arxiv.org/pdf/1711.06746.pdf. 

Depends R (>= 3.0)

In addition to this "README" file, there four files in this repository.

1. "Principal_Manifold_Estimation.R": This file contains the functions "hdmde" and "PME." The "hdmde" function implements the high-dimensional mixture density estimation (HDMDE) presented in Algorithm 1 of ME, and the "PME" function implements the principal manifold estimation presented in Algorithm 2 of ME. 

The syntax of PME:
```r
result=PME(x.obs=data_points, d=intrinsic_dimension, N0=0, tuning.para.seq=exp((-15:5)), alpha=0.05, max.comp=100, epsilon=0.05, max.iter=100, print.MSDs=TRUE)
# "data_points" is an n-by-D matrix, each row of which denotes a D-dimensional data point.
# "intrinsic_dimension": is the dimension of the potentially interested underlying manifold.
# More details are provided in the "Principal_Manifold_Estimation.R" file.
```

2. "Examples_for_PME.R": This file presents six simulated examples with dimension pairs (d=1, D=2), (d=1,D=3), and (d=2,D=3) to demonstrate the performance of PME.

```r
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
```

3. "Graphical_Output_of_Examples.pdf": This file shows the graphical output of "Examples_for_PME.R."

4. "Comparison_between_PME_and_existing_method.R": This file reproduces the simulation study in Figure 6 (b) of ME, comparing PME, HS, and ISOMAP. We measure the performances of these three methods using mean squared distance (MSD). In this file, for simplicity, we do only ten simulations. The corresponding output is the ten MSDs of each method, say 30 numbers in a 10-by-3 matrix. HS stands for the HS principal curve algorithm (T. Hastie and W. Stuetzle. Principal curves. Journal of the American Statistical Association, 84(406): 502-516, 1989). In this file, we apply the "principal_curve" function to perform HS. ISOMAP is based on (J. B. Tenenbaum, V. De Silva, and J. C. Langford. A global geometric framework for nonlinear dimensionality reduction. science, 290(5500):2319-2323, 2000.) and is implemented using the "isomap" function. 

Kun (Michael) Meng,

Ph.D. Student,
Department of Biostatistics, 
Brown University
