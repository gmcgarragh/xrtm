#! /usr/bin/env julia


#*******************************************************************************
#
#*******************************************************************************
# Module that defines the XRTM interface.
push!(LOAD_PATH, "../interfaces")
using XRTM

using Printf


#*******************************************************************************
# Define inputs.
#*******************************************************************************
options       = ["calc_derivs", "delta_m", "n_t_tms", "output_at_levels", "source_solar"]

solvers       = ["eig_add"]

max_coef      = 34
n_quad        = 8
n_stokes      = 1
n_derivs      = 3
n_layers      = 4
n_theta_0s    = 1
n_kernel_quad = 16
kernels       = ["lambertian"]
n_out_levels  = 2
n_out_thetas  = 3
n_out_phis    = 1

F_0           = 1.
theta_0       = 35.

out_levels    = Array{Int32}([0, 4])

out_thetas    = [0., 30., 60.]
phi           = 45.
out_phis      = reshape([phi, phi, phi], 1, 3)

ltau          = [.02, .05, 1., .1]

omega         = [1., 1., .9, 1.]

n_ray_coef    = 3
ray_coef      = reshape([1.000000e+00,
                         0.000000e+00,
                         4.798741e-01], n_ray_coef, 1)

n_aer_coef    = 34
aer_coef      = reshape([1.000000e+00,
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
                         6.591328e-08], n_aer_coef, 1)

n_coef        = Array{Int32}([n_ray_coef, n_ray_coef, n_aer_coef, n_ray_coef])

coef          = zeros(max_coef, 1, n_layers)
coef[1:n_ray_coef,1,1] = ray_coef
coef[1:n_ray_coef,1,2] = ray_coef
coef[1:n_aer_coef,1,3] = aer_coef
coef[1:n_ray_coef,1,4] = ray_coef

albedo        = .2


#*******************************************************************************
# Create an XRTM instance.
#*******************************************************************************
xrtm = XRTM.create(options, solvers, max_coef, n_quad, n_stokes, n_derivs, n_layers,
                   n_theta_0s, n_kernel_quad, kernels, n_out_levels, n_out_thetas)


#*******************************************************************************
# Set inputs.
#
# Inputs must be set before the first model run.  For subsequent runs only the
# inputs that change need to be set.  For example calculating the radiance
# across the O2-A band spectrum, assuming constant scattering properites, would
# require only updating ltau and omega for each point.
#*******************************************************************************
XRTM.set_fourier_tol(xrtm, .0001)
XRTM.set_out_levels(xrtm, out_levels)
XRTM.set_out_thetas(xrtm, out_thetas)
XRTM.set_F_iso_top(xrtm, 0.)
XRTM.set_F_iso_bot(xrtm, 0.)
XRTM.set_F_0(xrtm, F_0)
XRTM.set_theta_0(xrtm, theta_0)
XRTM.set_phi_0(xrtm, 0.)

# Set optical property inputs
XRTM.set_ltau_n(xrtm, ltau)
XRTM.set_omega_n(xrtm, omega)
XRTM.set_coef_n(xrtm, n_coef, coef)

# Alternatively the inputs can be set one layer at a time.
for i in 1:n_layers
    XRTM.set_ltau_1(xrtm, i - 1, ltau[i])
    XRTM.set_omega_1(xrtm, i - 1, omega[i])
    XRTM.set_coef_1(xrtm, i - 1, n_coef[i], coef[1:n_coef[i],:,i])
end

# Set surface albedo
XRTM.set_kernel_ampfac(xrtm, 0, albedo)


#*******************************************************************************
# Set linearized inputs.
#*******************************************************************************
XRTM.set_ltau_l_11(xrtm, 2, 0, 1.)
XRTM.set_omega_l_11(xrtm, 2, 1, 1.)
XRTM.set_kernel_ampfac_l_1(xrtm, 0, 2, 1.)

XRTM.update_varied_layers(xrtm)


#*******************************************************************************
# Run the model for radiances and associated derivatives.  If this is the
# initial run and all the required inputs have not been initialized then XRTM
# will print(a appropriate message and return < 0.
#*******************************************************************************
I_p, I_m, K_p, K_m = XRTM.radiance(xrtm, "eig_add", n_out_phis, out_phis)


#*******************************************************************************
# Output results.
#*******************************************************************************
for i in 1:n_out_levels
    @printf("level: %d\n", i)
    @printf("     intensity:\n")

    for j in 1:n_out_thetas
        @printf("          theta = %9.2E, I_p = %13.6E, I_m = %13.6E\n",
                out_thetas[j], I_p[1,1,j,i], I_m[1,1,j,i])
    end

    for j in 1:n_derivs
        @printf("     derivative: %d\n", j)
        for k in 1:n_out_thetas
            @printf("          theta = %9.2E, K_p = %13.6E, K_m = %13.6E\n",
                    out_thetas[k], K_p[1,1,k,j,i], K_m[1,1,k,j,i])
        end
    end
end

println()


#*******************************************************************************
# Let's get the inputs and print them just for a test.
#*******************************************************************************
println(XRTM.get_fourier_tol(xrtm))
println(XRTM.get_out_levels(xrtm))
println(XRTM.get_out_thetas(xrtm))
println(XRTM.get_F_iso_top(xrtm))
println(XRTM.get_F_iso_bot(xrtm))
println(XRTM.get_F_0(xrtm))
println(XRTM.get_theta_0(xrtm))
println(XRTM.get_phi_0(xrtm))
println()

println(XRTM.get_ltau(xrtm, 0))
println(XRTM.get_ltau(xrtm, 2))
println(XRTM.get_omega(xrtm, 0))
println(XRTM.get_omega(xrtm, 2))
println(XRTM.get_coef(xrtm, 0, 0, n_ray_coef - 1))
println(XRTM.get_coef(xrtm, 2, 0, n_aer_coef - 1))
println()

println(XRTM.get_kernel_ampfac(xrtm, 0))
println()

println(XRTM.get_ltau_l(xrtm, 2, 0))
println(XRTM.get_omega_l(xrtm, 2, 1))
println(XRTM.get_kernel_ampfac_l(xrtm, 0, 2))
println()


#*******************************************************************************
# Destroy the xrtm.  (Currently the XRTM instance is not garbage collected.)
#*******************************************************************************
XRTM.destroy(xrtm)
