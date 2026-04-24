#! /usr/bin/env python


#*******************************************************************************
#
#*******************************************************************************
import numpy as np
import sys

sys.path.insert(0, '../interfaces')

# Module that defines the XRTM class which should not be modfied by the user.
import xrtm


#*******************************************************************************
# Define inputs.
#*******************************************************************************
options       = ['calc_derivs', 'delta_m', 'n_t_tms', 'output_at_levels', 'source_solar']

solvers       = ['eig_add']

max_coef      = 34
n_quad        = 8
n_stokes      = 1
n_derivs      = 3
n_layers      = 4
n_theta_0s    = 1
n_kernel_quad = 16
kernels       = ['lambertian']
n_out_levels  = 2
n_out_thetas  = 3
n_out_phis    = 1

F_0           = 1.
theta_0       = 35.

out_levels    = [0, 4]

out_thetas    = [0., 30., 60.]
phi           = [45.]
out_phis      = [phi, phi, phi]

ltau          = [.02, .05, 1., .1]

omega         = [1., 1., .9, 1.]

n_ray_coef    = 3
ray_coef      = [1.000000e+00,
                 0.000000e+00,
                 4.798741e-01]

n_aer_coef    = 34
aer_coef      = [1.000000e+00,
                 1.865569e+00,
                 1.789985e+00,
                 1.220838e+00,
                 7.472409e-01,
                 4.017337e-01,
                 2.173326e-01,
                 1.054020e-01,
                 5.737447e-02,
                 2.570752e-02,
                 1.527185e-02,
                 6.202491e-03,
                 4.278587e-03,
                 1.529611e-03,
                 1.276447e-03,
                 3.964385e-04,
                 4.036524e-04,
                 1.112891e-04,
                 1.338887e-04,
                 3.468251e-05,
                 4.611093e-05,
                 1.204792e-05,
                 1.637357e-05,
                 4.577401e-06,
                 5.975423e-06,
                 1.849954e-06,
                 2.241820e-06,
                 7.774087e-07,
                 8.673507e-07,
                 3.351400e-07,
                 3.476180e-07,
                 1.472730e-07,
                 1.448956e-07,
                 6.591328e-08]

n_coef        = [n_ray_coef, n_ray_coef, n_aer_coef, n_ray_coef]

coef          = np.zeros((n_layers, 1, max_coef), dtype = np.double)
coef[0,0,0:n_ray_coef] = np.reshape(ray_coef, (1, 1, n_ray_coef))
coef[1,0,0:n_ray_coef] = np.reshape(ray_coef, (1, 1, n_ray_coef))
coef[2,0,0:n_aer_coef] = np.reshape(aer_coef, (1, 1, n_aer_coef))
coef[3,0,0:n_ray_coef] = np.reshape(ray_coef, (1, 1, n_ray_coef))

albedo        = .2


#*******************************************************************************
# Create an XRTM instance.
#*******************************************************************************
try:
    model = xrtm.xrtm(options, solvers, max_coef, n_quad, n_stokes, n_derivs,
            n_layers, n_theta_0s, n_kernel_quad, kernels, n_out_levels, n_out_thetas)
except xrtm.error as e:
    print(str(e) + '\nERROR: xrtm.init()')
    exit()


#*******************************************************************************
# Set inputs.
#
# Inputs must be set before the first model run.  For subsequent runs only the
# inputs that change need to be set.  For example calculating the radiance
# across the O2-A band spectrum, assuming constant scattering properites, would
# require only updating ltau and omega for each point.
#*******************************************************************************
try:
    model.set_fourier_tol(.0001)
except xrtm.error as e:
    print(str(e) + '\nERROR: xrtm.set_fourier_tol()')
    exit()

try:
    model.set_out_levels(out_levels)
except xrtm.error as e:
    print(str(e) + '\nERROR: xrtm.set_out_levels()')
    exit()

try:
    model.set_out_thetas(out_thetas)
except xrtm.error as e:
    print(str(e) + '\nERROR: xrtm.set_out_thetas()')
    exit()

try:
    model.set_F_iso_top(0.)
except xrtm.error as e:
    print(str(e) + '\nERROR: xrtm.set_F_iso_top()')
    exit()

try:
    model.set_F_iso_bot(0.)
except xrtm.error as e:
    print(str(e) + '\nERROR: xrtm.set_F_iso_bot()')
    exit()

try:
    model.set_F_0(F_0)
except xrtm.error as e:
    print(str(e) + '\nERROR: xrtm.set_F_0()')
    exit()

try:
    model.set_theta_0(theta_0)
except xrtm.error as e:
    print(str(e) + '\nERROR: xrtm.set_theta_0()')
    exit()

try:
    model.set_phi_0(0.)
except xrtm.error as e:
    print(str(e) + '\nERROR: xrtm.set_phi_0()')
    exit()


# Set optical property inputs
try:
    model.set_ltau_n(ltau)
except xrtm.error as e:
    print(str(e) + '\nERROR: xrtm.set_ltau_n()')
    exit()

try:
    model.set_omega_n(omega)
except xrtm.error as e:
    print(str(e) + '\nERROR: xrtm.set_omega_n()')
    exit()

try:
    model.set_coef_n(n_coef, coef)
except xrtm.error as e:
    print(str(e) + '\nERROR: xrtm.set_coef_n()')
    exit()


# Alternatively the inputs can be set one layer at a time.
for i in range(0, n_layers):
    try:
        model.set_ltau_1(i, ltau[i])
    except xrtm.error as e:
        print(str(e) + '\nERROR: xrtm.set_ltau_n()')
        exit()

    try:
        model.set_omega_1(i, omega[i])
    except xrtm.error as e:
        print(str(e) + '\nERROR: xrtm.set_omega_n()')
        exit()

    try:
        model.set_coef_1(i, n_coef[i], coef[i][0:1,0:n_coef[i]])
    except xrtm.error as e:
        print(str(e) + '\nERROR: xrtm.set_coef_1()')
        exit()

# Set surface albedo
try:
    model.set_kernel_ampfac(0, albedo)
except xrtm.error as e:
    print(str(e) + '\nERROR: xrtm.set_kernel_ampfac()')
    exit()


#*******************************************************************************
# Set linearized inputs.
#*******************************************************************************
try:
    model.set_ltau_l_11(2, 0, 1.)
except xrtm.error as e:
    print(str(e) + '\nERROR: xrtm.set_ltau_l_11()')
    exit()

try:
    model.set_omega_l_11(2, 1, 1.)
except xrtm.error as e:
    print(str(e) + '\nERROR: xrtm.set_omega_l_11()')
    exit()

try:
    model.set_kernel_ampfac_l_1(0, 2, 1.)
except xrtm.error as e:
    print(str(e) + '\nERROR: xrtm.set_kernel_ampfac_l_1()')
    exit()

try:
    model.update_varied_layers()
except xrtm.error as e:
    print(str(e) + '\nERROR: model.update_varied_layers()')
    exit()


#*******************************************************************************
# Run the model for radiances and associated derivatives.  If this is the
# initial run and all the required inputs have not been initialized then XRTM
# will print(a appropriate message and return < 0.
#*******************************************************************************
try:
    I_p, I_m, K_p, K_m = model.radiance('eig_add', n_out_phis, out_phis)
except xrtm.error as e:
    print(str(e) + '\nERROR: model.radiance()')
    exit()


#*******************************************************************************
# Output results.
#*******************************************************************************
for i in range(0, n_out_levels):
    print('level: %d' % i)
    print('     intensity:')

    for j in range(0, n_out_thetas):
        print('          theta = %9.2E, I_p = %13.6E, I_m = %13.6E' % \
              (out_thetas[j], I_p[i,j,0,0], I_m[i,j,0,0]))
    for j in range(0, n_derivs):
        print('     derivative: %d' % j)
        for k in range(0, n_out_thetas):
            print('          theta = %9.2E, K_p = %13.6E, K_m = %13.6E' % \
                  (out_thetas[k], K_p[i,j,k,0,0], K_m[i,j,k,0,0]))
print()


#*******************************************************************************
# Let's get the inputs and print them just for a test.
#*******************************************************************************
print(model.get_fourier_tol())
print(model.get_out_levels())
print(model.get_out_thetas())
print(model.get_F_iso_top())
print(model.get_F_iso_bot())
print(model.get_F_0())
print(model.get_theta_0())
print(model.get_phi_0())
print()

print(model.get_ltau(0))
print(model.get_ltau(2))
print(model.get_omega(0))
print(model.get_omega(2))
print(model.get_coef(0, 0, n_ray_coef - 1))
print(model.get_coef(2, 0, n_aer_coef - 1))
print()

print(model.get_kernel_ampfac(0))
print()

print(model.get_ltau_l(2, 0))
print(model.get_omega_l(2, 1))
print(model.get_kernel_ampfac_l(0, 2))
print()


#*******************************************************************************
# Delete xrtm instance.
#*******************************************************************************
del model
