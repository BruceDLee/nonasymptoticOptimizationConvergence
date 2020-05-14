# nonasymptoticOptimizationConvergence

A working draft of code to implement the procedure for computing finite step performance guarantees of first order
optimization algorithms on strongly convex functions with Lipschitz continuous gradients presented in [1]


The code requires YALMIP and MOSEK

The key comparisons are as follows:

-testNonasymptoticWorstCaseTraj: computes the bound and a trajectory that achieves the bound for a sample algorithm.

-comparePEP: compares the solution generated by the approach in [1] with that from [2].

-compareAsymptotic: compares the asymptotic performance bound from [1] with that from [4].

-compareFiniteAsymptotic: generates Fig 1. in [1]. Demonstrates that the asymptotic convergence rate bound matches the decay
rate of the finite step performance bound.

**References**

[1] Lee, B. and Seiler, P. "Finite step performance of first-order methods using interpolation conditions without function evaluations," arxiv, 2020.   
[2] A. B. Taylor, J. M. Hendrickx, and F. Glineur, “Smooth strongly convex interpolation and exact worst-case performance 
of first-order methods,” Mathematical Programming, vol. 161, pp. 307–345, (1–2) 2017.   
[3] A. B. Taylor, J. M. Hendrickx, and F. Glineur, “Performance estimation toolbox (PESTO): Automated worst-case analysis 
of first-order optimization methods,” in 2017 IEEE 56th Annual Conference on Decision and Control (CDC), 2017, pp. 1278–1283.  
[4] L. Lessard, B. Recht, and A. Packard, “Analysis and design of optimization algorithms via integral quadratic constraints,”
SIAM Journal on Optimization, vol. 26, no. 1, pp. 57–95, 2016.
