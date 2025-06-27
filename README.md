# libc_MG_poisson
Implementation of geometric and algebraic multigrid for Poisson equation in 1D, 2D and 3D

* GMG: Geometric multigrid
* AMG: Algebraic multigrid

Author: Pengliang Yang, Harbin Institute of Technology, China

E-mail: ypl.2100@gmail.com

Contents:
========

* sxamg-master: a AMG toolbox sxamg developed by refactorization of FASP by Dr. Liu.

* demo_mumps_poisson_2d: direct solve Poisson equation with MUMPS.

* demo_gmg_amg_poisson: implemented my own GMG while calling sxamg for AMG

* demo_mg_poisson_1/2/3d: my own implementation of both GMG and AMG. All of them are working nicely.

* demo_cg_poisson_3d: 3D poisson equation solved using conjugate gradient (CG) method 
