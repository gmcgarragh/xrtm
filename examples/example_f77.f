c***********************************************************************
c
c***********************************************************************
      program example_f77

      implicit none

      include '../interfaces/xrtm_int_f77.inc'

      integer i
      integer j
      integer k

      integer error

      integer options
      integer solvers
      integer max_coef
      parameter (max_coef = 34)
      integer n_quad
      parameter (n_quad = 8)
      integer n_stokes
      parameter (n_stokes = 1)
      integer n_derivs
      parameter (n_derivs = 3)
      integer n_layers
      parameter (n_layers = 4)
      integer n_kernels
      parameter (n_kernels = 1)
      integer n_kernel_quad
      parameter (n_kernel_quad = 16)
      integer n_out_levels
      parameter (n_out_levels = 2)
      integer n_out_thetas
      parameter (n_out_thetas = 3)
      integer n_out_phis
      parameter (n_out_phis = 1)

      integer out_levels(n_out_levels)
      data out_levels /0, 4/

      real*8  F_0
      data F_0 /1./
      real*8  theta_0
      data theta_0 /35./

      real*8  out_thetas(n_out_thetas)
      data out_thetas /0., 30., 60./
      real*8  out_phis(n_out_thetas, n_out_phis)
      data out_phis /45., 45., 45./

      real*8  ltau(n_layers)
      data ltau /.02, .05, 1., .1/
      real*8  omega(n_layers)
      data omega /1., 1., .9, 1./

      integer n_ray_coef
      parameter (n_ray_coef = 3)
      real*8  ray_coef(n_ray_coef,1)
      data ray_coef /1.000000e+00,
     &               0.000000e+00,
     &               4.798741e-01/

      integer n_aer_coef
      parameter (n_aer_coef = max_coef)
      real*8  aer_coef(n_aer_coef,1)
      data aer_coef /1.000000e+00,
     &               1.865569e+00,
     &               1.789985e+00,
     &               1.220838e+00,
     &               7.472409e-01,
     &               4.017337e-01,
     &               2.173326e-01,
     &               1.054020e-01,
     &               5.737447e-02,
     &               2.570752e-02,
     &               1.527185e-02,
     &               6.202491e-03,
     &               4.278587e-03,
     &               1.529611e-03,
     &               1.276447e-03,
     &               3.964385e-04,
     &               4.036524e-04,
     &               1.112891e-04,
     &               1.338887e-04,
     &               3.468251e-05,
     &               4.611093e-05,
     &               1.204792e-05,
     &               1.637357e-05,
     &               4.577401e-06,
     &               5.975423e-06,
     &               1.849954e-06,
     &               2.241820e-06,
     &               7.774087e-07,
     &               8.673507e-07,
     &               3.351400e-07,
     &               3.476180e-07,
     &               1.472730e-07,
     &               1.448956e-07,
     &               6.591328e-08/

      integer n_coef(n_layers)
      data n_coef /n_ray_coef, n_ray_coef, n_aer_coef, n_ray_coef/
      real*8  coef(max_coef,1,n_layers)

      integer kernels(n_kernels)
      data kernels /XRTM_KERNEL_LAMBERTIAN/

      real*8  albedo
      data albedo /.2/

      real*8  I_p(n_stokes, n_out_phis,
     &            n_out_thetas,           n_out_levels)
      real*8  I_m(n_stokes, n_out_phis,
     &            n_out_thetas,           n_out_levels)
      real*8  K_p(n_stokes, n_out_phis,
     &            n_out_thetas, n_derivs, n_out_levels)
      real*8  K_m(n_stokes, n_out_phis,
     &            n_out_thetas, n_derivs, n_out_levels)

      byte xrtm(N_BYTES_XRTM_TYPE)


c     ******************************************************************
c     *
c     ******************************************************************
      options = 0
      options = ior(options, XRTM_OPTION_CALC_DERIVS)
      options = ior(options, XRTM_OPTION_DELTA_M)
      options = ior(options, XRTM_OPTION_N_T_TMS)

      solvers = 0
      solvers = ior(solvers, XRTM_SOLVER_EIG_ADD)

      do i = 1, n_ray_coef
           coef(i,1,1) = ray_coef(i,1)
      enddo

      do i = 1, n_ray_coef
           coef(i,1,2) = ray_coef(i,1)
      enddo

      do i = 1, n_aer_coef
           coef(i,1,3) = aer_coef(i,1)
      enddo

      do i = 1, n_ray_coef
           coef(i,1,4) = ray_coef(i,1)
      enddo


c     ******************************************************************
c     *
c     ******************************************************************
      call xrtm_create_f77(xrtm, options, solvers, max_coef, n_quad,
     &                     n_stokes, n_derivs, n_layers, n_kernels,
     &                     n_kernel_quad, kernels, n_out_levels,
     &                     n_out_thetas, error)
      if (error /= 0) then
           write (0, *) 'error calling xrtm_create_f77()'
           stop
      endif


c     ******************************************************************
c     *
c     ******************************************************************
      call xrtm_set_fourier_tol_f77(xrtm, .0001d0, error)
      if (error .ne. 0) then
           write (0, *) 'error calling xrtm_set_fourier_tol_f77()'
           stop
      endif

      call xrtm_set_F_0_f77(xrtm, F_0, error)
      if (error .ne. 0) then
           write (0, *) 'error calling xrtm_set_F_0_f77()'
           stop
      endif

      call xrtm_set_theta_0_f77(xrtm, theta_0, error)
      if (error .ne. 0) then
           write (0, *) 'error calling xrtm_set_theta_0_f77()'
           stop
      endif

      call xrtm_set_phi_0_f77(xrtm, 0.d0, error)
      if (error .ne. 0) then
           write (0, *) 'error calling xrtm_set_phi_0_f77()'
           stop
      endif

      call xrtm_set_out_levels_f77(xrtm, out_levels, error)
      if (error .ne. 0) then
           write (0, *) 'error calling xrtm_set_out_levels_f77()'
           stop
      endif

      call xrtm_set_out_thetas_f77(xrtm, out_thetas, error)
      if (error .ne. 0) then
           write (0, *) 'error calling xrtm_set_out_thetas,_f77()'
           stop
      endif


      call xrtm_set_ltau_n_f77(xrtm, ltau, error)
      if (error .ne. 0) then
           write (0, *) 'error calling xrtm_set_ltau_n_f77()'
           stop
      endif

      call xrtm_set_omega_n_f77(xrtm, omega, error)
      if (error .ne. 0) then
           write (0, *) 'error calling xrtm_set_omega_n_f77()'
           stop
      endif

      call xrtm_set_coef_n_f77(xrtm, n_coef, coef, error)
      if (error .ne. 0) then
           write (0, *) 'error calling xrtm_set_coef_n_f77()'
           stop
      endif

      do i = 1, n_layers
           call xrtm_set_ltau_1_f77(xrtm, i - 1, ltau(i), error)
           if (error .ne. 0) then
                write (0, *) 'error calling xrtm_set_ltau_1_f77()'
                stop
           endif

           call xrtm_set_omega_1_f77(xrtm, i - 1, omega(i), error)
           if (error .ne. 0) then
                write (0, *) 'error calling xrtm_set_omega_1_f77()'
                stop
           endif

           call xrtm_set_coef_1_f77(xrtm, i - 1, n_coef(i),
     &                              coef(1,1,i), error)
           if (error .ne. 0) then
                write (0, *) 'error calling xrtm_set_coef_1_f77()'
                stop
           endif
      enddo

      call xrtm_set_kernel_ampfac_f77(xrtm, 0, albedo, error)
      if (error .ne. 0) then
           write (0, *) 'error calling xrtm_set_kernel_ampfac_f77()'
           stop
      endif


c     ******************************************************************
c     *
c     ******************************************************************
      call xrtm_set_ltau_l_11_f77(xrtm, 2, 0, 1.d0, error)
      if (error .ne. 0) then
           write (0, *) 'error calling xrtm_set_ltau_l_11_f77()'
           stop
      endif

      call xrtm_set_omega_l_11_f77(xrtm, 2, 1, 1.d0, error)
      if (error .ne. 0) then
           write (0, *) 'error calling xrtm_set_omega_l_11_f77()'
           stop
      endif

      call xrtm_set_kernel_ampfac_l_1_f77(xrtm, 0, 2, 1.d0, error)
      if (error .ne. 0) then
           write (0, *) 'error calling xrtm_set_kernel_ampfac_l_1_f77()'
           stop
      endif

      call xrtm_update_varied_layers_f77(xrtm, error)
      if (error .ne. 0) then
           write (0, *) 'error calling xrtm_update_varied_layers_f77()'
           stop
      endif


c     ******************************************************************
c     *
c     ******************************************************************
      call xrtm_radiance_f77(xrtm, XRTM_SOLVER_EIG_ADD, n_out_phis,
     &                       out_phis, I_p, I_m, K_p, K_m, error)
      if (error .ne. 0) then
           write (0, *) 'error calling xrtm_radiance_f77()'
           stop
      endif


c     ******************************************************************
c     *
c     ******************************************************************
      do i = 1, n_out_levels
           print '("level: ", I1)', i
           print '("     intensity:")'
           do j = 1, n_out_thetas
                print '("          theta = ", ES9.2,
     &                ", I_p = ", ES13.6, ", I_m = ", ES13.6)',
     &                out_thetas(j), I_p(1,1,j,i), I_m(1,1,j,i)
           enddo
           do j = 1, n_derivs
                print '("     derivative: ", I1)', j
                do k = 1, n_out_thetas
                     print '("          theta = ", ES9.2,
     &                    ", K_p = ", ES13.6, ", K_m = ", ES13.6)',
     &                    out_thetas(k), K_p(1,1,k,j,i), K_m(1,1,k,j,i)
                enddo
           enddo
      enddo


c     ******************************************************************
c     *
c     ******************************************************************
      call xrtm_destroy_f77(xrtm, error)


      end
