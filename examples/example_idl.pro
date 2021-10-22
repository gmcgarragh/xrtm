pro example_idl


;*******************************************************************************
;
;*******************************************************************************
HOME = getenv('HOME')

PATH = '../interfaces:' + '<IDL_DEFAULT>'

pref_set, 'IDL_PATH',     PATH, /commit
pref_set, 'IDL_DLM_PATH', PATH, /commit


;*******************************************************************************
; Define inputs.
;*******************************************************************************
options         = ['calc_derivs', 'delta_m', 'n_t_tms', 'output_at_levels', 'source_solar']

solvers         = ['eig_add']

max_coef        = 34
n_quad          = 8
n_stokes        = 1
n_derivs        = 3
n_layers        = 4
n_theta_0s      = 1
n_kernel_quad   = 16
kernels         = ['lambertian']
n_out_levels    = 2
n_out_thetas    = 3
n_out_phis      = 1

F_0             = 1.
theta_0         = 35.

out_levels      = [0l, 4l]

out_thetas      = [0.d, 30.d, 60.d]
phi             = 45.d
out_phis        = dblarr(n_out_phis, n_out_thetas)
out_phis(*,*)   = phi

ltau            = [.02d, .05d, 1.d, .1d]

omega           = [1.d, 1.d, .9d, 1.d]

n_ray_coef      = 3L

ray_coef        = [1.000000d+00, $
                   0.000000d+00, $
                   4.798741d-01]
n_aer_coef      = 34L

aer_coef        = [1.000000d+00, $
                   1.865569d+00, $
                   1.789985d+00, $
                   1.220838d+00, $
                   7.472409d-01, $
                   4.017337d-01, $
                   2.173326d-01, $
                   1.054020d-01, $
                   5.737447d-02, $
                   2.570752d-02, $
                   1.527185d-02, $
                   6.202491d-03, $
                   4.278587d-03, $
                   1.529611d-03, $
                   1.276447d-03, $
                   3.964385d-04, $
                   4.036524d-04, $
                   1.112891d-04, $
                   1.338887d-04, $
                   3.468251d-05, $
                   4.611093d-05, $
                   1.204792d-05, $
                   1.637357d-05, $
                   4.577401d-06, $
                   5.975423d-06, $
                   1.849954d-06, $
                   2.241820d-06, $
                   7.774087d-07, $
                   8.673507d-07, $
                   3.351400d-07, $
                   3.476180d-07, $
                   1.472730d-07, $
                   1.448956d-07, $
                   6.591328d-08]

n_coef          = [n_ray_coef, n_ray_coef, n_aer_coef, n_ray_coef]

coef            = dblarr(max_coef, 1, 4)
coef[*,*,*]     = 0.d
coef[0:n_ray_coef-1,0,0] = ray_coef
coef[0:n_ray_coef-1,0,1] = ray_coef
coef[0:n_aer_coef-1,0,2] = aer_coef
coef[0:n_ray_coef-1,0,3] = ray_coef

albedo          = .2


;*******************************************************************************
; Create an XRTM instance.
;*******************************************************************************
xrtm_create, xrtm, options, solvers, max_coef, n_quad, n_stokes, n_derivs, $
             n_layers, n_theta_0s, n_kernel_quad, kernels, n_out_levels, $
             n_out_thetas


;*******************************************************************************
; Set inputs.
;
; Inputs must be set before the first model run.  For subsequent runs only the
; inputs that change need to be set.  For example calculating the radiance
; across the O2-A band spectrum, assuming constant scattering properites, would
; require only updating ltau and omega for each point.
;*******************************************************************************
xrtm_set_fourier_tol, xrtm, .0001d
xrtm_set_out_levels, xrtm, out_levels
xrtm_set_out_thetas, xrtm, out_thetas
xrtm_set_F_iso_top, xrtm, 0.
xrtm_set_F_iso_bot, xrtm, 0.
xrtm_set_F_0, xrtm, F_0
xrtm_set_theta_0, xrtm, theta_0
xrtm_set_phi_0, xrtm, 0.

; Set optical property inputs
xrtm_set_ltau_n, xrtm, ltau
xrtm_set_omega_n, xrtm, omega
xrtm_set_coef_n, xrtm, n_coef, coef

; Alternatively the inputs can be set one layer at a time.
for i = 0, n_layers - 1 do begin
    xrtm_set_ltau_1, xrtm, i, ltau[i]
    xrtm_set_omega_1, xrtm, i, omega[i]
    xrtm_set_coef_1, xrtm, i, n_coef[i], reform(coef[0:n_coef[i]-1,*,i], n_coef[i], 1)
endfor

; Set surface albedo
xrtm_set_kernel_ampfac, xrtm, 0, albedo


;*******************************************************************************
; Set linearized inputs.
;*******************************************************************************
xrtm_set_ltau_l_11, xrtm, 2, 0, 1.
xrtm_set_omega_l_11, xrtm, 2, 1, 1.
xrtm_set_kernel_ampfac_l_1, xrtm, 0, 2, 1.

xrtm_update_varied_layers, xrtm


;*******************************************************************************
; Run the model for radiances and associated derivatives.  If this is the
; initial run and all the required inputs have not been initialized then XRTM
; will print a appropriate message and return < 0.
;*******************************************************************************
xrtm_radiance, xrtm, 'eig_add', n_out_phis, out_phis, I_p, I_m, K_p, K_m


;*******************************************************************************
; Output results.
;*******************************************************************************
for i = 0, n_out_levels - 1 do begin
    print, i, format = '("level: ", I1)'
    print,    format = '("     intensity:")'

    for j = 0, n_out_thetas - 1 do begin
        print, out_thetas[j], I_p[0,0,j,i], I_m[0,0,j,i], format = $
               '("          theta = ", E9.2, ", I_p = ", E13.6, ", I_m = ", E13.6)'
    endfor
    for j = 0, n_derivs - 1 do begin
        print, j, format = '("     derivative: ", I1)'
        for k = 0, n_out_thetas - 1 do begin
            print, out_thetas[k], K_p[0,0,k,j,i], K_m[0,0,k,j,i], format = $
                   '("          theta = ", E9.2, ", K_p = ", E13.6, ", K_m = ", E13.6)'
        endfor
    endfor
endfor
print


;*******************************************************************************
; Let's get the inputs and print them just for a test.
;*******************************************************************************
xrtm_get_fourier_tol, xrtm, x & print, x
xrtm_get_out_levels,  xrtm, x & print, x
xrtm_get_out_thetas,  xrtm, x & print, x
xrtm_get_F_iso_top,   xrtm, x & print, x
xrtm_get_F_iso_bot,   xrtm, x & print, x
xrtm_get_F_0,         xrtm, x & print, x
xrtm_get_theta_0,     xrtm, x & print, x
xrtm_get_phi_0,       xrtm, x & print, x
print

xrtm_get_ltau,  xrtm, 0, x & print, x
xrtm_get_ltau,  xrtm, 2, x & print, x
xrtm_get_omega, xrtm, 0, x & print, x
xrtm_get_omega, xrtm, 2, x & print, x
xrtm_get_coef,  xrtm, 0, 0, n_ray_coef - 1, x & print, x
xrtm_get_coef,  xrtm, 2, 0, n_aer_coef - 1, x & print, x
print

xrtm_get_kernel_ampfac, xrtm, 0, x & print, x
print

xrtm_get_ltau_l,           xrtm, 2, 0, x & print, x
xrtm_get_omega_l,          xrtm, 2, 1, x & print, x
xrtm_get_kernel_ampfac_l,  xrtm, 0, 2, x & print, x
print


;*******************************************************************************
; Destroy the model.
;*******************************************************************************
xrtm_destroy, xrtm


end
