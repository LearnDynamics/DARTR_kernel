# DARTR_kernel
Data Adaptive RKHS Tikhonov Regularization for learning kernels in operators


### Problem Statement:
We fit to the function data $(u_i,f_i)$ an operator $R_\phi[u_i] = f_i$, where  
    $$ R_\phi[u](x) = \int \phi(|y|)g[u](x,x+y) dy 
                   = \sum_r \phi(r) [ g[u](x,x+r)+ g[u](x,x-r) ] dr.$$

Examples include 

- Linear integral operator with $g[u](x,y) =  u(x+y) $  
- Nonlinear operator with $g[u](x,y) = div(u(x+y))u(x) $
- nonlocal operator with $g[u](x,y) =  u(x+y)-u(x)$   


### About the code

- start from [run_demo.m], change the settings there
- The current code is only for vector estimator: viewing the kernel as a vector on the grid points. Code for estimator on basis functions available at request.  