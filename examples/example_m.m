%*******************************************************************************
%
%*******************************************************************************
addpath("../interfaces");

options         = 0;

solvers         = 0;

max_coef        = 34;
n_quad          = 8;
n_stokes        = 1;
n_derivs        = 3;
n_layers        = 4;
n_theta_0s      = 1;
n_kernels       = 0;
n_kernel_quad   = 16;
kernels         = zeros(0, 1, 'int32')
n_out_levels    = 2;
n_out_thetas    = 3;
n_out_phis      = 1;

xrtm = xrtm_int_mex(options, solvers, max_coef, n_quad, n_stokes, n_derivs, n_layers, n_theta_0s, n_kernels, n_kernel_quad, kernels, n_out_levels, n_out_thetas);
