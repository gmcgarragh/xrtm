!*******************************************************************************
!
!*******************************************************************************
program example_f90

     use xrtm_int_f90

     implicit none

     integer              :: i
     integer              :: j
     integer              :: k

     integer              :: error

     integer              :: options
     integer              :: solvers
     integer              :: max_coef
     integer              :: n_quad
     integer              :: n_stokes
     integer              :: n_derivs
     integer              :: n_layers
     integer              :: n_theta_0s
     integer              :: n_kernels
     integer              :: n_kernel_quad
     integer              :: kernels(1)
     integer              :: n_out_levels
     integer              :: n_out_thetas
     integer              :: n_out_phis

     integer, allocatable :: out_levels(:)

     integer              :: n_ray_coef
     integer              :: n_aer_coef

     integer, allocatable :: n_coef(:)

     real(8)              :: F_0
     real(8)              :: theta_0

     real(8), allocatable :: out_thetas(:)
     real(8), allocatable :: out_phis(:,:)

     real(8), allocatable :: ltau(:)
     real(8), allocatable :: omega(:)
     real(8)              :: ray_coef(3,1)
     real(8), allocatable :: aer_coef(:,:)
     real(8), allocatable :: coef(:,:,:)

     real(8)              :: albedo

     real(8), allocatable :: I_p(:,:,:,:)
     real(8), allocatable :: I_m(:,:,:,:)
     real(8), allocatable :: K_p(:,:,:,:,:)
     real(8), allocatable :: K_m(:,:,:,:,:)

     ! Structure the represents the XRTM instance the members of which should
     ! never be modfied by the user.
     type(xrtm_type)      :: xrtm


     !**************************************************************************
     ! Define inputs.
     !**************************************************************************
     options       = 0
     options       = ior(options, XRTM_OPTION_CALC_DERIVS)
     options       = ior(options, XRTM_OPTION_DELTA_M)
     options       = ior(options, XRTM_OPTION_N_T_TMS)
     options       = ior(options, XRTM_OPTION_OUTPUT_AT_LEVELS)
     options       = ior(options, XRTM_OPTION_SOURCE_SOLAR)

     solvers       = 0
     solvers       = ior(solvers, XRTM_SOLVER_EIG_ADD)

     max_coef      = 34
     n_quad        = 8
     n_stokes      = 1
     n_derivs      = 3
     n_layers      = 4
     n_theta_0s    = 1
     n_kernels     = 1
     n_kernel_quad = 16
     kernels(1)    = XRTM_KERNEL_LAMBERTIAN
     n_out_levels  = 2
     n_out_thetas  = 3
     n_out_phis    = 1

     allocate(out_levels(n_out_levels))
     allocate(out_thetas(n_out_thetas))
     allocate(out_phis(n_out_thetas, n_out_phis))
     allocate(ltau(n_layers))
     allocate(omega(n_layers))
     allocate(aer_coef(max_coef, 1))
     allocate(n_coef(n_layers))
     allocate(coef(max_coef, 1, n_layers))

     F_0           = 1.
     theta_0       = 35.

     out_levels    = (/ 0, 4 /)

     out_thetas    = (/ 0., 30., 60. /)
     out_phis      = reshape((/ 45., 45., 45./), shape(out_phis))

     ltau          = (/ .02, .05, 1., .1 /)

     omega         = (/ 1., 1., .9, 1. /)

     n_ray_coef    = 3

     ray_coef      = reshape((/ 1.000000e+00, &
                                0.000000e+00, &
                                4.798741e-01 /), shape(ray_coef))
     n_aer_coef    = 34

     aer_coef      = reshape((/ 1.000000e+00, &
                                1.865569e+00, &
                                1.789985e+00, &
                                1.220838e+00, &
                                7.472409e-01, &
                                4.017337e-01, &
                                2.173326e-01, &
                                1.054020e-01, &
                                5.737447e-02, &
                                2.570752e-02, &
                                1.527185e-02, &
                                6.202491e-03, &
                                4.278587e-03, &
                                1.529611e-03, &
                                1.276447e-03, &
                                3.964385e-04, &
                                4.036524e-04, &
                                1.112891e-04, &
                                1.338887e-04, &
                                3.468251e-05, &
                                4.611093e-05, &
                                1.204792e-05, &
                                1.637357e-05, &
                                4.577401e-06, &
                                5.975423e-06, &
                                1.849954e-06, &
                                2.241820e-06, &
                                7.774087e-07, &
                                8.673507e-07, &
                                3.351400e-07, &
                                3.476180e-07, &
                                1.472730e-07, &
                                1.448956e-07, &
                                6.591328e-08 /), shape(aer_coef))

     n_coef(1)     = n_ray_coef
     n_coef(2)     = n_ray_coef
     n_coef(3)     = n_aer_coef
     n_coef(4)     = n_ray_coef

     coef(:,:,1)   = ray_coef
     coef(:,:,2)   = ray_coef
     coef(:,:,3)   = aer_coef
     coef(:,:,4)   = ray_coef

     albedo        = .2


     !**************************************************************************
     ! Create an XRTM instance.
     !**************************************************************************
     call xrtm_create_f90(xrtm, options, solvers, max_coef, n_quad, n_stokes, &
                          n_derivs, n_layers, n_theta_0s, n_kernels, n_kernel_quad, &
                          kernels, n_out_levels, n_out_thetas, error)
     if (error /= 0) then
          write (0, *) 'error calling xrtm_create_f90()'
          stop
     endif


     !**************************************************************************
     ! Set inputs.
     !
     ! Inputs must be set before the first model run.  For subsequent runs only
     ! the inputs that change need to be set.  For example calculating the
     ! radiance across the O2-A band spectrum, assuming constant scattering
     ! properites, would require only updating ltau and omega for each point.
     ! *************************************************************************
     call xrtm_set_fourier_tol_f90(xrtm, .0001d0, error)
     if (error /= 0) then
          write (0, *) 'error calling xrtm_set_fourier_tol_f90()'
          stop
     endif

     call xrtm_set_out_levels_f90(xrtm, out_levels, error)
     if (error /= 0) then
          write (0, *) 'error calling xrtm_set_out_levels_f90()'
          stop
     endif

     call xrtm_set_out_thetas_f90(xrtm, out_thetas, error)
     if (error /= 0) then
          write (0, *) 'error calling xrtm_set_out_thetas,_f90()'
          stop
     endif

     call xrtm_set_F_0_f90(xrtm, F_0, error)
     if (error /= 0) then
          write (0, *) 'error calling xrtm_set_F_0_f90()'
          stop
     endif

     call xrtm_set_F_iso_top_f90(xrtm, 0.d0, error)
     if (error /= 0) then
          write (0, *) 'error calling xrtm_set_F_iso_top_f90()'
          stop
     endif

     call xrtm_set_F_iso_bot_f90(xrtm, 0.d0, error)
     if (error /= 0) then
          write (0, *) 'error calling xrtm_set_F_iso_bot_f90()'
          stop
     endif

     call xrtm_set_theta_0_f90(xrtm, theta_0, error)
     if (error /= 0) then
          write (0, *) 'error calling xrtm_set_theta_0_f90()'
          stop
     endif

     call xrtm_set_phi_0_f90(xrtm, 0.d0, error)
     if (error /= 0) then
          write (0, *) 'error calling xrtm_set_phi_0_f90()'
          stop
     endif


     ! Set optical property inputs
     call xrtm_set_ltau_n_f90(xrtm, ltau, error)
     if (error /= 0) then
          write (0, *) 'error calling xrtm_set_ltau_n_f90()'
          stop
     endif

     call xrtm_set_omega_n_f90(xrtm, omega, error)
     if (error /= 0) then
          write (0, *) 'error calling xrtm_set_omega_n_f90()'
          stop
     endif

     call xrtm_set_coef_n_f90(xrtm, n_coef, coef, error)
     if (error /= 0) then
          write (0, *) 'error calling xrtm_set_coef_n_f90()'
          stop
     endif

     ! Alternatively the inputs can be set one layer at a time.
     do i = 1, n_layers
          call xrtm_set_ltau_1_f90(xrtm, i - 1, ltau(i), error)
          if (error /= 0) then
               write (0, *) 'error calling xrtm_set_ltau_1_f90()'
               stop
          endif

          call xrtm_set_omega_1_f90(xrtm, i - 1, omega(i), error)
          if (error /= 0) then
               write (0, *) 'error calling xrtm_set_omega_1_f90()'
               stop
          endif

          call xrtm_set_coef_1_f90(xrtm, i - 1, n_coef(i), coef(:,:,i), error)
          if (error /= 0) then
               write (0, *) 'error calling xrtm_set_coef_1_f90()'
               stop
          endif
     enddo

     ! Set surface albedo
     call xrtm_set_kernel_ampfac_f90(xrtm, 0, albedo, error)
     if (error /= 0) then
          write (0, *) 'error calling xrtm_set_kernel_ampfac_f90()'
          stop
     endif


     !**************************************************************************
     ! Set linearized inputs.
     !**************************************************************************
     call xrtm_set_ltau_l_11_f90(xrtm, 2, 0, 1.d0, error)
     if (error /= 0) then
          write (0, *) 'error calling xrtm_set_ltau_l_11_f90()'
          stop
     endif

     call xrtm_set_omega_l_11_f90(xrtm, 2, 1, 1.d0, error)
     if (error /= 0) then
          write (0, *) 'error calling xrtm_set_omega_l_11_f90()'
          stop
     endif

     call xrtm_set_kernel_ampfac_l_1_f90(xrtm, 0, 2, 1.d0, error)
     if (error /= 0) then
          write (0, *) 'error calling xrtm_set_kernel_ampfac_l_1_f90()'
          stop
     endif

     call xrtm_update_varied_layers_f90(xrtm, error)
     if (error /= 0) then
          write (0, *) 'error calling xrtm_update_varied_layers_f90()'
          stop
     endif


     !**************************************************************************
     ! Allocate output arrays.
     !**************************************************************************
     allocate(I_p(n_stokes, n_out_phis, n_out_thetas,           n_out_levels))
     allocate(I_m(n_stokes, n_out_phis, n_out_thetas,           n_out_levels))
     allocate(K_p(n_stokes, n_out_phis, n_out_thetas, n_derivs, n_out_levels))
     allocate(K_m(n_stokes, n_out_phis, n_out_thetas, n_derivs, n_out_levels))


     !**************************************************************************
     ! Run the model for radiances and associated derivatives.  If this is the
     ! initial run and all the required inputs have not been initialized then
     ! XRTM will print a appropriate message and return < 0.
     !**************************************************************************
     call xrtm_radiance_f90(xrtm, XRTM_SOLVER_EIG_ADD, n_out_phis, out_phis, &
                            I_p, I_m, K_p, K_m, error)
     if (error /= 0) then
          write (0, *) 'error calling xrtm_radiance_f90()'
          stop
     endif


     !**************************************************************************
     ! Output results.
     !**************************************************************************
     do i = 1, n_out_levels
          print '("level: ", I1)', i
          print '("     intensity:")'
          do j = 1, n_out_thetas
               print '("          theta = ", ES9.2, ", I_p = ", ES13.6, ", I_m = ", ES13.6)', &
                     out_thetas(j), I_p(1,1,j,i), I_m(1,1,j,i)
          enddo
          do j = 1, n_derivs
               print '("     derivative: ", I1)', j
               do k = 1, n_out_thetas
                    print '("          theta = ", ES9.2, ", K_p = ", ES13.6, ", K_m = ", ES13.6)', &
                          out_thetas(k), K_p(1,1,k,j,i), K_m(1,1,k,j,i)
               enddo
          enddo
     enddo


     !**************************************************************************
     ! Deallocate memory and destroy the model.
     !**************************************************************************
     deallocate(out_levels)

     deallocate(n_coef)

     deallocate(out_thetas)
     deallocate(out_phis)

     deallocate(ltau)
     deallocate(omega)
     deallocate(aer_coef)
     deallocate(coef)

     deallocate(I_p)
     deallocate(I_m)
     deallocate(K_p)
     deallocate(K_m)

     call xrtm_destroy_f90(xrtm, error)


end program example_f90
