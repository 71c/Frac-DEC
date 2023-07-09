# Code for paper "Structure-Preserving Discretization of Fractional Vector Calculus using Discrete Exterior Calculus"

This repository contains the code used in the paper.

The files used in the paper are:
- `matlab/test_fractional_gradient_discretization_v3.m`, `matlab/test_fractional_curl_discretization_v3.m`, and `matlab/test_fractional_divergence_discretization_v3.m`: used to generate the plots in Section 4.1
- `matlab/test_dd_0.m`: used to generate the plots in Section 4.2
- `mathematica/numerical_experiment_calculations.nb`: does the analytical calculations for the numerical experiments in Section 4.1. This uses `mathematica/ToMatlab.m`, which is a slight modification to the code by Harri Ojanen available [here](https://library.wolfram.com/infocenter/MathSource/577/), which converts Mathematica code to Matlab code. The code is modified to add support for hypergeometric function, gamma function, incomplete gamma function, and sinint and cosint functions.
