!*******************************************************************************
!
!    Copyright (C) 2007-2020 Greg McGarragh <greg.mcgarragh@colostate.edu>
!
!    This source code is licensed under the GNU General Public License (GPL),
!    Version 3.  See the file COPYING for more details.
!
!*******************************************************************************

subroutine call_soi_f(n_four, n_elem, n_coef, n_quad, n_layers, qx, F_0, &
                      theta_0, phi_0, umus, n_umus, phis, n_phis, planet_r, &
                      levels_z, n_kernels, n_kernel_quad, kernels, ampfac, &
                      params, n_coef_layer, chi, omega, ltau, I_p, I_m, d_tau, &
                      delta_m, n_t_tms, psa, epsilon, info, n_umus2)

      use soi_rt_model_az
      use soi_surface_module


      implicit none


      integer, parameter          :: max_params = 3

      integer, intent(in)         :: n_four
      integer, intent(in)         :: n_elem
      integer, intent(in)         :: n_coef
      integer, intent(in)         :: n_quad
      integer, intent(in)         :: n_layers

      real(8), intent(in)         :: qx(n_quad)

      real(8), intent(in)         :: F_0

      real(8), intent(in)         :: theta_0

      real(8), intent(in)         :: phi_0

      real(8), intent(in)         :: umus(n_umus)
      integer, intent(in)         :: n_umus

      real(8), intent(in)         :: phis(n_phis)
      integer, intent(in)         :: n_phis

      real(8), intent(in)         :: planet_r
      real(8), intent(in), target :: levels_z(n_layers+1)

      integer, intent(in)         :: n_kernels
      integer, intent(in)         :: n_kernel_quad
      integer, intent(in)         :: kernels(n_kernels)
      real(8), intent(in)         :: ampfac(n_kernels)
      real(8), intent(in)         :: params(max_params, n_kernels)

      integer, intent(in)         :: n_coef_layer(n_layers)
      real(8), intent(in)         :: chi(n_coef, n_elem, n_layers)

      real(8), intent(in)         :: omega(n_layers)

      real(8), intent(in)         :: ltau(n_layers)

      real(8), intent(out)        :: I_p(n_phis, n_umus2, n_layers+1)
      real(8), intent(out)        :: I_m(n_phis, n_umus2, n_layers+1)

      real(8), intent(in)         :: d_tau

      integer, intent(in)         :: delta_m
      integer, intent(in)         :: n_t_tms
      integer, intent(in)         :: psa

      real(8), intent(in)         :: epsilon

      integer, intent(inout)      :: info

      integer, intent(in)         :: n_umus2


      !*************************************************************************
      !
      !*************************************************************************
      integer, parameter     :: max_kernels = 9

      logical                :: delta_m2
      logical                :: azimuth_avg
      logical                :: psa2

      integer                :: i
      integer                :: ii
      integer                :: j
      integer                :: jj

      integer                :: i_kernel

      integer                :: n_t_tms2

      integer                :: kernel_index(max_kernels)    = (/1,7,4,5,2,3,6,8,9/)

      integer                :: kernel_n_params(max_kernels) = (/0,0,0,2,2,3,0,3,2/)

      real(8)                :: a
      real(8)                :: b

      real(8)                :: pi

      real(8), pointer       :: zheight(:)

      real(8)                :: omega2(n_layers)
      real(8)                :: ltau2(n_layers)

      real(8)                :: chi2(n_coef, n_layers)

      type (surface_type)    :: surface


      !**************************************************************************
      !
      !**************************************************************************
      if (n_umus == 0) then
           write (0, *) 'ERROR: soi does not do output at quadrature points'
           info = 1
           return
      endif

      if (n_umus > 1) then
           write (0, *) 'ERROR: soi can only handle one user zenith angle'
           info = 1
           return
      endif

      if (n_phis > 1) then
           write (0, *) 'ERROR: soi can only handle one user azimuth angle'
           info = 1
           return
      endif


      !*************************************************************************
      !
      !*************************************************************************
      a = qx(1)
      a = planet_r

      info = 0

      pi = 4.d0 * atan(1.0d0)


      !*************************************************************************
      !
      !*************************************************************************
      if (delta_m == 0) then
           delta_m2 = .false.
      else
           delta_m2 = .true.
      endif

      if (n_four == 0) then
           azimuth_avg = .true.
      else
           azimuth_avg = .false.
      endif

      if (n_t_tms == 0) then
           n_t_tms2 = 0
      else
           n_t_tms2 = 1
      endif

      if (psa == 0) then
           psa2 = .false.
           nullify(zheight)
      else
           psa2 = .true.
           zheight => levels_z
      endif

      do i = 1, n_layers
           ii = n_layers - i + 1
           omega2(ii) = omega(i)
           ltau2 (ii) = ltau (i)
      enddo

      do j = 1, n_layers
           jj = n_layers - j + 1
           do i = 1, n_coef_layer(j)
                chi2(i, jj) = chi(i, 1, j) / (2. * (i-1) + 1.)
           enddo

           do i = i, n_coef
                chi2(i, jj) = 0.
           enddo
      enddo

      surface%use_reflected_direct  = .true.
      surface%use_reflected_diffuse = .true.
      surface%use_surface_emission  = .false.
      surface%n_brdf_kernels        = n_kernels

      do i = 1, n_kernels
           i_kernel                      = kernel_index(kernels(i) + 1)
           surface%brdf_kernel(i)        = i_kernel - 1
           surface%kernel_amp_par(i)     = ampfac(i)
           surface%n_kernel_dist_par(i)  = kernel_n_params(i_kernel)
           do j = 1, kernel_n_params(i)
                surface%kernel_dist_par(j, i) = params(j, i)
           enddo
      enddo

      surface%n_brdf_quadratures = n_kernel_quad * 2


      !*************************************************************************
      !
      !*************************************************************************
      a = acos(umus(1)) * 180.d0 / pi

      call soi_solar_rt(ltau2, omega2, chi2, a, phis(1) - phi_0, &
                        theta_0, F_0, surface, b, zheight=zheight, &
                        delta_scaling=delta_m2, verbose=.false., &
                        nstreams=n_quad * 2, azimuth_avg=azimuth_avg, &
                        converge=epsilon, max_tau=d_tau, &
                        truncation=.false., adding=.false., &
                        tms=n_t_tms2, method=1, check_pf=.false., &
                        ps_correct=psa2, fourier_tol=epsilon, &
                        use_saved_surf=.false., use_saved_chap=.false., &
                        calc_ss=.false.)

      call soi_deallocate()

      I_p(1, 1, 1) = b
      I_m(1, 1, 1) = 0.
      I_p(1, 1, 2) = 0.
      I_m(1, 1, 2) = 0.


      return

end subroutine call_soi_f
