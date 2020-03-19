# Principal-Manifold-Estimation

This repository is devoted to the principal manifold estimation (PME) algorithm. The PME algorithm comes from the paper "Principal Manifolds: A Framework Using Sobolev Spaces and Model Complexity Selection Using Mixture Densities" by Kun Meng and Ani Eloyan. This paper is available at https://arxiv.org/pdf/1711.06746.pdf and will be updated soon. 

In addition to this "README" file, there four files in this repository.

1. "Principal_Manifold_Estimation.R": This file contains functions "hdmde" and "PME." The "hdmde" function conducts high-dimensional mixture density estimation (HDMDE), and the "PME" function conducts principal manifold estimation. "hdmde" is contained in "PME" and works for "PME."

2. "Examples_for_PME.R": This file applies six examples with dimension pairs (d=1, D=2), (d=1,D=3), and (d=2,D=3) to demonstrate the performance of PME.

3. "Graphical_Output_of_Example.pdf": After running the entire "Examples_for_PME.R" file, you will get this graphical output.

4. "Comparison_between_PME_and_existing_method.R": This file reproduces a simulation study comparing PME, HS, and ISOMAP. The performances of these three methods are measured using mean squared distance (MSD). In this file, we do ten simulations. The corresponding output is the ten MSDs of each method. HS stands for the HS principal curve algorithm (T. Hastie and W. Stuetzle. Principal curves. Journal of the American Statistical Association, 84(406): 502?516, 1989). In this file, we apply the "principal_curve" function to perform HS. ISOMAP is based on (J. B. Tenenbaum, V. De Silva, and J. C. Langford. A global geometric framework for nonlinear dimensionality reduction. science, 290(5500):2319?2323, 2000.) and is implemented using the "isomap" function. 

Kun (Michael) Meng,\\ 
Ph.D. Student,\\
Department of Biostatistics, \\
Brown University