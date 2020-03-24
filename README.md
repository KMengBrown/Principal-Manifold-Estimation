# Principal-Manifold-Estimation

This repository is devoted to the principal manifold estimation (PME) algorithm. The PME algorithm comes from the paper "Principal Manifolds: A Framework and Model Complexity Selection" by Kun Meng and Ani Eloyan (hereafter ME). This paper is available at https://arxiv.org/pdf/1711.06746.pdf. 

In addition to this "README" file, there four files in this repository.

1. "Principal_Manifold_Estimation.R": This file contains the functions "hdmde" and "PME." The "hdmde" function implements the high-dimensional mixture density estimation (HDMDE) presented in Algorithm 1 of ME, and the "PME" function implements the principal manifold estimation presented in Algorithm 2 of ME. 

2. "Examples_for_PME.R": This file presents six simulated examples with dimension pairs (d=1, D=2), (d=1,D=3), and (d=2,D=3) to demonstrate the performance of PME.

3. "Graphical_Output_of_Examples.pdf": This file shows the graphical output of "Examples_for_PME.R."

4. "Comparison_between_PME_and_existing_method.R": This file reproduces the simulation study in Figure 6 (b) of ME, comparing PME, HS, and ISOMAP. We measure the performances of these three methods using mean squared distance (MSD). In this file, for simplicity, we do only ten simulations. The corresponding output is the ten MSDs of each method, say 30 numbers in a 10-by-3 matrix. HS stands for the HS principal curve algorithm (T. Hastie and W. Stuetzle. Principal curves. Journal of the American Statistical Association, 84(406): 502-516, 1989). In this file, we apply the "principal_curve" function to perform HS. ISOMAP is based on (J. B. Tenenbaum, V. De Silva, and J. C. Langford. A global geometric framework for nonlinear dimensionality reduction. science, 290(5500):2319-2323, 2000.) and is implemented using the "isomap" function. 

Kun (Michael) Meng,

Ph.D. Student,
Department of Biostatistics, 
Brown University