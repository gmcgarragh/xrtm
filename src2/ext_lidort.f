c*******************************************************************************
c
c    Copyright (C) 2007-2020 Greg McGarragh <greg.mcgarragh@colostate.edu>
c
c    This source code is licensed under the GNU General Public License (GPL),
c    Version 3.  See the file COPYING for more details.
c
c*******************************************************************************

      subroutine call_lidort_f
     &     (n_four, n_elem, n_coef, n_quad, n_derivs, n_layers, qx, F_0,
     &      theta_0, phi_0, ulevels, utaus, n_ulevels, umus, n_umus,
     &      phis, n_phis, planet_r,levels_z, n_kernels, n_kernel_quad,
     &      kernels, ampfac, ampfac_l, params, params_l, n_coef_layer,
     &      chi, chi_l,omega, omega_l, ltau, ltau_l, I_p,I_m, K_p, K_m,
     &      mean_p,mean_m, mean_p_l, mean_m_l, flux_p, flux_m, flux_p_l,
     &      flux_m_l, delta_m, n_t_tms, psa, quad_output, radiance,
     &      mean, flux, utau_output, derivs, epsilon, info, n_quad2)

      implicit none


      integer max_params
      parameter (max_params = 4)

      integer n_four
      integer n_elem
      integer n_coef
      integer n_quad
      integer n_derivs
      integer n_layers

      real*8 qx(n_quad)

      real*8 F_0

      real*8 theta_0

      real*8 phi_0

      integer ulevels(n_ulevels)
      real*8 utaus(n_ulevels)
      integer n_ulevels

      real*8 umus(n_umus)
      integer n_umus

      real*8 phis(n_phis)
      integer n_phis

      real*8 planet_r
      real*8 levels_z(n_layers+1)

      integer n_kernels
      integer n_kernel_quad
      integer kernels(n_kernels)
      real*8 ampfac(n_kernels)
      real*8 ampfac_l(n_derivs, n_kernels)
      real*8 params(max_params, n_kernels)
      real*8 params_l(max_params, n_derivs, n_kernels)

      integer n_coef_layer(n_layers)
      real*8 chi  (n_coef, n_elem, n_layers)
      real*8 chi_l(n_coef, n_elem, n_derivs, n_layers)

      real*8 omega(n_layers)
      real*8 omega_l(n_derivs, n_layers)

      real*8 ltau(n_layers)
      real*8 ltau_l(n_derivs, n_layers)

      real*8 I_p(n_phis, n_quad2, n_ulevels)
      real*8 I_m(n_phis, n_quad2, n_ulevels)

      real*8 K_p(n_phis, n_quad2, n_derivs, n_ulevels)
      real*8 K_m(n_phis, n_quad2, n_derivs, n_ulevels)

      real*8 mean_p(n_ulevels)
      real*8 mean_m(n_ulevels)

      real*8 mean_p_l(n_derivs, n_ulevels)
      real*8 mean_m_l(n_derivs, n_ulevels)

      real*8 flux_p(n_ulevels)
      real*8 flux_m(n_ulevels)

      real*8 flux_p_l(n_derivs, n_ulevels)
      real*8 flux_m_l(n_derivs, n_ulevels)

      integer delta_m
      integer n_t_tms
      integer psa
      integer quad_output
      integer radiance
      integer mean
      integer flux
      integer utau_output

      integer*1 derivs(n_derivs + 8, n_layers)

      real*8 epsilon

      integer info

      integer n_quad2


c     ******************************************************************
c     *
c     ******************************************************************
      logical flag

      integer max_kernels_xrtm
      parameter (max_kernels_xrtm = 10)
      integer max_kernels_lidort
      parameter (max_kernels_lidort = 9)

      character*10 kernel_names(max_kernels_lidort + 2)
      data kernel_names /
     &     'Lambertian',
     &     'Ross-thin ',
     &     'Ross-thick',
     &     'Li-sparse ',
     &     'Li-dense  ',
     &     'Hapke     ',
     &     'Roujean   ',
     &     'Rahman    ',
     &     'Cox-Munk  ',
     &     'GISS-soil ',
     &     'GISS-snow '/

      integer i
      integer i1
      integer i2
      integer j
      integer jj
      integer k
      integer k2
      integer l
      integer l2

      integer i_kernel

      integer i_layer
      integer i_deriv
      integer i_brdf

      integer status_inputcheck
      integer status_calculation

      integer quad_index(n_quad2)

      integer kernel_index(max_kernels_xrtm)
      data kernel_index /1,0,7,4,5,2,3,6,8,9/

      integer kernel_n_params(max_kernels_lidort)
      data kernel_n_params /0,0,0,2,2,3,0,3,2/

      byte    deriv_type(n_derivs)
      integer layers_index(n_derivs)
      integer derivs_index(n_derivs)
      integer brdf_index(n_derivs)

      real*8 a

      real*8 pi

      real*8 angles(n_quad2)

      real*8 brdf_factor(n_derivs)


c     ******************************************************************
c     *
c     ******************************************************************
      include 'LIDORT.PARS'
      include 'LIDORT_INPUTS.VARS'
      include 'LIDORT_RESULTS.VARS'
      include 'LIDORT_L_INPUTS.VARS'
      include 'LIDORT_L_RESULTS.VARS'

      include 'LIDORT_L_SOLUTION.VARS'

      L_U_WPOS = 0.
      L_U_WNEG = 0.


c     ******************************************************************
c     *
c     ******************************************************************
c      if (delta_m .ne. 0) then
c           do i = 1, n_layers
c                if (chi(n_coef+1, 1, i) .eq. 0.) then
c                     write (0, *) 'ERROR: lidort fails if using'     //
c     &                    ' delta-m scaling and the truncation'      //
c     &                    ' factor (last + 1 phase function momment)'//
c     &                    ' is = zero'
c                     info = 1
c                     return
c                endif
c           enddo
c      endif


c     ******************************************************************
c     *
c     ******************************************************************
      a = qx(1)

      info = 0

      pi = 4.d0 * atan(1.0d0)


c     ******************************************************************
c     *
c     ******************************************************************
      do_fullrad_mode         = .true.

      if (n_t_tms .eq. 0) then
           do_sscorrection    = .false.
      else
           do_sscorrection    = .true.
      endif
      do_dbcorrection         = .false.

      save_layer_msst         = .false.

      do_double_convtest      = .true.

      do_direct_beam          = .false.
      do_classical_solution   = .true.

      if (psa .eq. 0) then
           do_plane_parallel  = .true.
      else
           do_plane_parallel  = .false.
      endif
      do_refractive_geometry  = .false.
      do_chapman_function     = .true.

      do_rayleigh_only        = .false.
      do_isotropic_only       = .false.
      if (n_four .eq. 1) then
           do_no_azimuth      = .true.
      else
           do_no_azimuth      = .false.
      endif
      do_all_fourier          = .false.

      if (delta_m .eq. 0) then
           do_deltam_scaling  = .false.
      else
           do_deltam_scaling  = .true.
      endif

      do_solution_saving      = .false.
      do_bvp_telescoping      = .false.

      do_upwelling            = .true.
      do_dnwelling            = .true.

      if (quad_output .ne. 0) then
           do_quad_output     = .true.
           do_user_streams    = .false.
      else
           do_quad_output     = .false.
           do_user_streams    = .true.
      endif
      if (utau_output .eq. 0) then
           do_user_taus       = .false.
           do_lbound_taus     = .true.
      else
           do_user_taus       = .true.
           do_lbound_taus     = .false.
      endif

      do_additional_mvout = .false.
      do_mvout_only       = .false.

      if (radiance .ne. 0) then
           if (mean .eq. 0 .and. flux .eq. 0) then
                do_additional_mvout = .false.
           else
                do_additional_mvout = .true.
           endif
      else
           if (mean .eq. 0 .and. flux .eq. 0) then
                do_mvout_only       = .false.
           else
                do_mvout_only       = .true.
           endif
      endif

      do_debug_write          = .false.

      do_write_input          = .false.
      do_write_scenario       = .false.
      do_write_fourier        = .false.
      do_write_results        = .false.

      input_write_filename    = 'lidort_input.out'
      scenario_write_filename = 'lidort_scenario.out'
      fourier_write_filename  = 'lidort_fourier.out'
      results_write_filename  = 'lidort_results.out'

      lidort_error_filename   = 'lidort_error.out'
      lidort_error_init       = .true.

      do_fdtest               = .false.
      fdepsilon               = 0.

      nstreams                = n_quad

      nlayers                 = n_layers

      nmoments_input          = n_coef - 1

      flux_factor             = F_0

      lidort_accuracy         = epsilon

      zenith_tolerance        = 1.e-8

      nbeams                  = 1

      beam_szas(1)            = theta_0

c     double precision sza_local_input(0:maxlayers, maxbeams)

      earth_radius = planet_r

      rfindex_parameter = 0.

      n_user_relazms = n_phis
      do i = 1, n_phis
           user_relazms(i) = phis(i) - phi_0
      enddo

      a = 180.d0 / pi
      if (quad_output .ne. 0) then
           n_user_streams = n_quad
           do i = 1, n_quad
                quad_index(i) = n_quad - i + 1
                user_angles_input(i) = acos(qx(quad_index(i))) * a
           enddo
      else
           if (n_umus .eq. 0) then
                n_user_streams = n_quad
                do i = 1, n_quad
                     quad_index(i) = n_quad - i + 1
                     user_angles_input(i) = acos(qx(quad_index(i))) * a
                enddo
           else
                n_user_streams = n_umus
                do i = 1, n_umus
                     angles(i) = acos(umus(i)) * a
                enddo

                if (n_umus .eq. 1) then
                     quad_index(1) = 1
                else
                     call lidort_indexx(n_umus, angles, quad_index)
                endif

                do i = 1, n_umus
                     user_angles_input(i) = angles(quad_index(i))
                enddo
           endif
      endif

      n_out_usertaus = n_ulevels
      if (utau_output .eq. 0) then
           do i = 1, n_ulevels
                lbound_taus_input(i) = ulevels(i)
           enddo
      else
           do i = 1, n_ulevels
                user_taus_input(i)   = utaus(i)
           enddo
      endif

      if (psa .ne. 0) then
           do i = 0, n_layers
                height_grid(i) = levels_z(i+1)
           enddo
      endif
c     double precision pressure_grid   (0:maxlayers)
c     double precision temperature_grid(0:maxlayers)

c     integer finegrid(maxlayers)

      do i = 1, n_layers
           omega_total_input(i) = omega(i)
           deltau_vert_input(i) = ltau(i)
      enddo

      do j = 1, n_layers
           do i = 0, n_coef_layer(j) - 1
                phasmoms_total_input(i, j) = chi(i+1, 1, j)
           enddo

           do i = i, n_coef - 1
                phasmoms_total_input(i, j) = 0.
           enddo
      enddo

      if (n_kernels .eq. 1 .and. kernels(1) .eq. 0) then
           do_lambertian_surface = .true.
           lambertian_albedo     = ampfac(1)
      else
           do_lambertian_surface = .false.

           n_brdf_kernels        = n_kernels

           do i = 1, n_kernels
                i_kernel = kernel_index(kernels(i) + 1)
                if (i_kernel .eq. 0) then
                     print *, 'ERROR: unsupported kernel for lidort'
                     info = 1
                     return
                endif
                brdf_names(i) = kernel_names(i_kernel)

                if (i_kernel .ne. 1) then
                     lambertian_kernel_flag(i) = .false.
                else
                     lambertian_kernel_flag(i) = .true.
                endif

                which_brdf(i)              = i_kernel

                n_brdf_parameters(i)       = kernel_n_params(i_kernel)

                do j = 1, kernel_n_params(i_kernel)
                     brdf_parameters(i, j) = params(j, i)
                enddo

                brdf_factors(i)            = ampfac(i)
           enddo

           nstreams_brdf = n_kernel_quad * 2
      endif

      do_surface_emission       = .false.
      surfbb                    = 0.
      fp_surfbb                 = 0.


c     ******************************************************************
c     *
c     ******************************************************************
      if (n_derivs .eq. 0) then
           do_simulation_only        = .false.
           do_linearization          = .false.

           do_atmos_linearization    = .false.

           do_surface_linearization  = .false.
           do_surfbb_linearization   = .false.
      else
           do_simulation_only        = .false.
           do_linearization          = .true.

           do_atmos_linearization    = .true.

           do_surface_linearization  = .false.
           do_surfbb_linearization   = .false.

           n_totalatmos_wfs          = n_derivs

           do i = 1, n_layers
                layer_vary_flag  (i) = .false.
                layer_vary_number(i) = 0
           enddo

           do i = 1, n_kernels
                do_kernel_factor_wfs(i) = .false.
           enddo

           i_brdf = 0

           do j = 1, n_derivs
                do i = 1, n_layers
                     if (derivs(j, i) .ne. 0) then
                          do_atmos_linearization = .true.

                          layer_vary_flag  (i) = .true.

                          layer_vary_number(i) = layer_vary_number(i) + 1

                          jj = layer_vary_number(i)

                          if (omega(i) .ne. 0.) then
                               l_omega_total_input(jj, i) =
     &                              omega_l(j, i) / omega(i)
                          endif

                          if (ltau(i)  .ne. 0.) then
                               l_deltau_vert_input(jj, i) =
     &                               ltau_l(j, i) / ltau(i)
                          endif

                          do k = 0, n_coef_layer(i) - 1
                               if (chi(k+1, 1, i) .ne. 0.) then
                                    l_phasmoms_total_input(jj, k, i) =
     &                                  chi_l(k+1, 1, j, i) / chi(k+1, 1, i)
                               else
                                    l_phasmoms_total_input(jj, k, i) = 0.
                               endif
                          enddo

                          do k = k, n_coef - 1
                               l_phasmoms_total_input(jj, k, i) = 0.
                          enddo

                          deriv_type(j)   = 1
                          layers_index(j) = i
                          derivs_index(j) = layer_vary_number(i)

                          goto 666
                     endif
                enddo

                if (derivs(j, i) .ne. 0) then

                     do_surface_linearization = .true.

                     do k = 1, n_kernels
                          if (ampfac_l(j, k) .ne. 0.) then
                               if (do_kernel_factor_wfs(k) .eqv. .true.) then
                                    print *, 'ERROR: in lidort each kernel amplitude factor can only be linearized once'
                                    info = 1
                                    return
                               endif

                               do_kernel_factor_wfs(k) = .true.
                               brdf_factor(i_brdf + 1) = ampfac(k)
                               l = 0
                               goto 420
                          endif

                          do l = 1, kernel_n_params(k)
                               if (params_l(l, j, k) .ne. 0.) then
                                    if (do_kernel_params_wfs(k,l) .eqv. .true.) then
                                         print *, 'ERROR: in lidort each kernel paramter can only be linearized once'
                                         info = 1
                                         return
                                    endif
                                    do_kernel_params_wfs(k,l) = .true.
                                    brdf_factor(i_brdf + 1) = params(l, k)
                                    goto 420
                               endif
                          enddo
                     enddo

420                  do k2 = k, n_kernels
                          if (ampfac_l(j, k2) .ne. 0. .and. k2 .ne. k) then
                               print *, 'ERROR: lidort only supports one linearized brdf value per deriv'
                               info = 1
                               return
                          endif

                          do l2 = l + 1, kernel_n_params(k2)
                               if (params_l(l2, j, k2) .ne. 0. .and. l2 .ne. l) then
                                    print *, 'ERROR: lidort only supports one linearized brdf value per deriv'
                                    info = 1
                                    return
                               endif

                               l = 0
                          enddo
                     enddo

c                    logical do_kparams_derivs(max_brdf_kernels)

c                    integer n_totalbrdf_wfs
c                    integer n_kernel_factor_wfs
c                    integer n_kernel_params_wfs

                     i_brdf = i_brdf + 1

                     deriv_type(j) = 2
                     brdf_index(j) = i_brdf
                else
                     deriv_type(j) = 0
                endif


666             do i = i + 1, n_layers + 1
                     if (derivs(j, i) .ne. 0) then
                          print *, 'ERROR: lidort only supports one linearized layer per deriv, i_layer,i_deriv =', i, j
                          info = 1
                          return
                     endif
                enddo
           enddo

           atmoswf_names(1)   = 'atmos wf'
           surfacewf_names(1) = 'surface wf'
      endif


c     ******************************************************************
c     *
c     ******************************************************************
      if (n_derivs .eq. 0) then
           call lidort_master(status_inputcheck, status_calculation)
           if (status_inputcheck .eq. LIDORT_SERIOUS) then
                print *, 'ERROR: lidort_master(), '//
     &                   'status_inputcheck = ', status_inputcheck
                info = 1
                return
           endif
           if (status_calculation .eq. LIDORT_SERIOUS) then
                print *, 'ERROR: lidort_master(), '//
     &                   'status_calculation = ', status_calculation
                info = 1
                return
           endif
      else
           call lidort_l_master(status_inputcheck, status_calculation)
           if (status_inputcheck .eq. LIDORT_SERIOUS) then
                print *, 'ERROR: lidort_l_master(), '//
     &                   'status_inputcheck = ', status_inputcheck
                info = 1
                return
           endif
           if (status_calculation .eq. LIDORT_SERIOUS) then
                print *, 'ERROR: lidort_l_master(), '//
     &                   'status_calculation = ', status_calculation
                info = 1
                return
           endif
      endif

c     call lidort_status(status_inputcheck, status_calculation)


c     ******************************************************************
c     *
c     ******************************************************************
      flag = .false.
      if ((n_user_streams .eq. 1 .and. user_angles_input(1) .eq. zero)
     &     .or. do_mvout_only) then
           flag = .true.
      endif

      if (radiance .ne. 0) then
           do i = 1, n_user_streams
                i1 = quad_index(i)
                do k = 1, n_phis
                     if (flag) then
                          i2 = (i - 1) * n_phis + 1
                     else
                          i2 = (i - 1) * n_phis + k
                     endif
                     do l = 1, n_out_usertaus
                          I_p(k, i1, l) = intensity(l, i2, 1)
                          I_m(k, i1, l) = intensity(l, i2, 2)
                     enddo
                enddo
           enddo

           if (n_derivs .gt. 0) then
                do k = 1, n_totalatmos_wfs
                     if (deriv_type(k) .eq. 0) then
                          do i = 1, n_user_streams
                               i1 = quad_index(i)
                               do j = 1, n_phis
                                    do l = 1, n_out_usertaus
                                         K_p(j, i1, k, l) = 0.
                                         K_m(j, i1, k, l) = 0.
                                    enddo
                               enddo
                          enddo
                     else if (deriv_type(k) .eq. 1) then
                          i_layer = layers_index(k)
                          i_deriv = derivs_index(k)
                          do i = 1, n_user_streams
                               i1 = quad_index(i)
                               do j = 1, n_phis
                                    if (flag) then
                                         i2 = (i - 1) * n_phis + 1
                                    else
                                         i2 = (i - 1) * n_phis + j
                                    endif
                                    do l = 1, n_out_usertaus
                                         K_p(j, i1, k, l) =
     &                                   atmoswf(i_deriv, i_layer, l, i2, 1)
                                         K_m(j, i1, k, l) =
     &                                   atmoswf(i_deriv, i_layer, l, i2, 2)
                                    enddo
                               enddo
                          enddo
                     else if (deriv_type(k) .eq. 2) then
                          i_brdf = brdf_index(k)
                          do i = 1, n_user_streams
                               i1 = quad_index(i)
                               do j = 1, n_phis
                                    if (flag) then
                                         i2 = (i - 1) * n_phis + 1
                                    else
                                         i2 = (i - 1) * n_phis + j
                                    endif
                                    do l = 1, n_out_usertaus
                                         K_p(j, i1, k, l) =
     &                                   surfacewf(i_brdf, l, i2, 1) / brdf_factor(i_brdf)
                                         K_m(j, i1, k, l) =
     &                                   surfacewf(i_brdf, l, i2, 2) / brdf_factor(i_brdf)
                                    enddo
                               enddo
                          enddo
                     endif
                enddo
           endif
      endif

      if (mean .ne. 0) then
           call lidort_copy_mean_value
     &          (n_derivs, n_layers, n_out_usertaus, n_totalatmos_wfs,
     &           deriv_type, layers_index, derivs_index, brdf_index,
     &           mean_intensity, mint_atmoswf, mint_surfacewf,
     &           mean_p, mean_m, mean_p_l, mean_m_l)
      endif

      if (flux .ne. 0) then
           call lidort_copy_mean_value
     &          (n_derivs, n_layers, n_out_usertaus, n_totalatmos_wfs,
     &           deriv_type, layers_index, derivs_index, brdf_index,
     &           flux_integral, flux_atmoswf, flux_surfacewf,
     &           flux_p, flux_m, flux_p_l, flux_m_l)
      endif

c     print *, fourier_saved(1)

      return

      end



c***********************************************************************
c
c***********************************************************************
      subroutine lidort_copy_mean_value
     &     (n_derivs, n_layers, n_out_usertaus, n_totalatmos_wfs,
     &      deriv_type, layers_index, derivs_index, brdf_index,
     &      lidort_mean, lidort_atmos_l, lidort_surface_l,
     &      mean_p, mean_m, mean_p_l, mean_m_l)

      implicit none

      include 'LIDORT.PARS'

      integer n_derivs
      integer n_layers
      integer n_out_usertaus
      integer n_totalatmos_wfs

      byte    deriv_type(n_derivs)
      integer layers_index(n_derivs)
      integer derivs_index(n_derivs)
      integer brdf_index(n_derivs)

      real*8 lidort_mean(MAX_OUT_USERTAUS, MAXBEAMS, MAX_DIRECTIONS)
      real*8 lidort_atmos_l  (MAX_ATMOSWFS, MAXLAYERS, MAX_OUT_USERTAUS,
     &                        MAXBEAMS, MAX_DIRECTIONS)
      real*8 lidort_surface_l(MAX_SURFACEWFS, MAX_OUT_USERTAUS,
     &                        MAXBEAMS, MAX_DIRECTIONS)

      real*8 mean_p(n_layers+1)
      real*8 mean_m(n_layers+1)

      real*8 mean_p_l(n_derivs, n_layers+1)
      real*8 mean_m_l(n_derivs, n_layers+1)


      integer k
      integer l
      integer i_layer
      integer i_deriv
      integer i_brdf


      do l = 1, n_out_usertaus
           mean_p(l) = lidort_mean(l, 1, 1)
           mean_m(l) = lidort_mean(l, 1, 2)
      enddo

      do k = 1, n_totalatmos_wfs
           if (deriv_type(k) .eq. 0) then
                i_layer = layers_index(k)
                i_deriv = derivs_index(k)
                do l = 1, n_out_usertaus
                     mean_p_l(k, l) =
     &                    lidort_atmos_l(i_deriv, i_layer, l, 1, 1)
                     mean_m_l(k, l) =
     &                    lidort_atmos_l(i_deriv, i_layer, l, 1, 2)
                enddo
           else
                i_brdf = brdf_index(k)
                do l = 1, n_out_usertaus
                     mean_p_l(k, l) =
     &                    lidort_surface_l(i_brdf, l, 1, 1)
                     mean_m_l(k, l) =
     &                    lidort_surface_l(i_brdf, l, 1, 2)
                enddo
          endif
      enddo

      return

      end
