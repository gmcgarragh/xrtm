c*******************************************************************************
c
c    Copyright (C) 2007-2020 Greg McGarragh <greg.mcgarragh@colostate.edu>
c
c    This source code is licensed under the GNU General Public License (GPL),
c    Version 3.  See the file COPYING for more details.
c
c*******************************************************************************

      subroutine call_vlidort_f
     &     (n_four, n_elem, n_coef, n_quad, n_stokes, n_derivs,n_layers,
     &      qx, F_0, theta_0, phi_0, ulevels, utaus, n_ulevels, umus,
     &      n_umus, phis, n_phis, planet_r, levels_z,levels_b,surface_b,
     &      n_kernels, n_kernel_quad, kernels,ampfac, ampfac_l, params,
     &      params_l, n_coef_layer, chi, chi_l, omega, omega_l, ltau,
     &      ltau_l, I_p,I_m,K_p, K_m, mean_p, mean_m, mean_p_l,mean_m_l,
     &      flux_p, flux_m, flux_p_l, flux_m_l, delta_m, n_t_tms, psa,
     &      quad_output, radiance, mean, flux, thermal, utau_output,
     &      derivs, epsilon, info, n_quad2)

      implicit none


      integer max_params
      parameter (max_params = 3)

      integer n_four
      integer n_coef
      integer n_elem
      integer n_quad
      integer n_stokes
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

      real*8 levels_b(n_layers+1)
      real*8 surface_b

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

      real*8 I_p(n_stokes, n_phis, n_quad2, n_ulevels)
      real*8 I_m(n_stokes, n_phis, n_quad2, n_ulevels)

      real*8 K_p(n_stokes, n_phis, n_quad2, n_derivs, n_ulevels)
      real*8 K_m(n_stokes, n_phis, n_quad2, n_derivs, n_ulevels)

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
      integer thermal
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

      character*10 kernel_names(max_kernels_lidort + 3)
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
     &     'GissCoxMnk',
     &     'Gisssoil  ',
     &     'Gisssnow  '/

      integer i
      integer i1
      integer i2
      integer j
      integer jj
      integer k
      integer k2
      integer l
      integer l2
      integer m

      integer i_kernel

      integer i_layer
      integer i_deriv
      integer i_brdf

      integer status_inputcheck
      integer status_calculation

      integer quad_index(n_quad+n_umus)

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

      real*8 angles(n_quad+n_umus)

      real*8 brdf_factor(n_derivs)


c     ******************************************************************
c     *
c     ******************************************************************
      include 'VLIDORT.PARS'
      include 'VLIDORT_INPUTS.VARS'
      include 'VLIDORT_RESULTS.VARS'
      include 'VLIDORT_L_INPUTS.VARS'
      include 'VLIDORT_L_RESULTS.VARS'

      include 'VLIDORT_L_SOLUTION.VARS'

      L_UPAR_DN_2 = 0.
      L_UPAR_UP_2 = 0.


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

      if (thermal .ne. 0) then
           write (0, *) 'ERROR: vlidort does not support thermal emission'
           info = 1
           return
      endif


c     ******************************************************************
c     *
c     ******************************************************************
      a = n_four
      a = F_0
      if (thermal .ne. 0) then
           a = levels_b(1)
      endif
      a = quad_output

      info = 0

      pi = 4.d0 * atan(1.0d0)


c     ******************************************************************
c     *
c     ******************************************************************
      do_fullrad_mode         = .true.

c     do_sscorrection         = .false.
      if (n_t_tms .eq. 0) then
           do_sscorr_nadir    = .false.
      else
           do_sscorr_nadir    = .true.
      endif
      do_sscorr_outgoing      = .false.
      do_dbcorrection         = .false.

c     save_layer_msst         = .false.

      do_double_convtest      = .true.

      if (psa .eq. 0) then
           do_plane_parallel  = .true.
      else
           do_plane_parallel  = .false.
      endif
      do_refractive_geometry  = .false.
      do_chapman_function     = .true.

      do_rayleigh_only        = .false.

      if (delta_m .eq. 0) then
           do_deltam_scaling  = .false.
      else
           do_deltam_scaling  = .true.
      endif

      do_solution_saving      = .false.
      do_bvp_telescoping      = .false.

      do_upwelling            = .true.
      do_dnwelling            = .true.

      do_quad_output          = .false.
      do_user_streams         = .true.
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
      do_info_write           = .false.

      do_write_input          = .false.
      do_write_scenario       = .false.
      do_write_fourier        = .false.
      do_write_results        = .false.

      input_write_filename    = 'vlidort_input.out'
      scenario_write_filename = 'vlidort_scenario.out'
      fourier_write_filename  = 'vlidort_fourier.out'
      results_write_filename  = 'vlidort_results.out'

      vlidort_error_filename  = 'vlidort_error.out'
      vlidort_error_init      = .true.

      do_fdtest               = .false.
      fdepsilon               = 0.

      nstokes                 = n_stokes

      nstreams                = n_quad

      nlayers                 = n_layers

      ngreek_moments_input    = n_coef - 1

      vlidort_accuracy        = epsilon

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
                call vlidort_indexx(n_umus, angles, quad_index)
           endif

           do i = 1, n_umus
                user_angles_input(i) = angles(quad_index(i))
           enddo
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
                do k = 1, 16
                     greekmat_total_input(i, j, k) = 0.
                enddo

                greekmat_total_input(i, j, 1) = chi(i+1, 1, j)

                if (n_elem .gt. 1) then
                     greekmat_total_input(i, j, 6)  =  chi(i+1, 2, j)
                     greekmat_total_input(i, j, 2)  = -chi(i+1, 5, j)
                     greekmat_total_input(i, j, 5)  = -chi(i+1, 5, j)
                     greekmat_total_input(i, j, 11) =  chi(i+1, 3, j)
                     greekmat_total_input(i, j, 16) =  chi(i+1, 4, j)
                     greekmat_total_input(i, j, 12) =  chi(i+1, 6, j)
                     greekmat_total_input(i, j, 15) = -chi(i+1, 6, j)
                endif
           enddo

           do i = i, n_coef - 1
                do k = 1, 16
                     greekmat_total_input(i, j, k) = 0.
                enddo
           enddo
      enddo

      if (n_kernels == 1 .and. kernels(1) == 0) then
           do_lambertian_surface = .true.
           lambertian_albedo     = ampfac(1)
      else
           do_lambertian_surface = .false.

           n_brdf_kernels        = n_kernels

           do i = 1, n_kernels
                i_kernel = kernel_index(kernels(i) + 1)

                brdf_names(i)              = kernel_names(i_kernel)

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

c     integer n_brdf_stokessq

      do_shadow_effect          = .false.

      do_coxmunk_dbms           = .false.

      if (thermal .eq. 0) then
           do_surface_emission  = .false.
      else
           do_surface_emission  = .true.
           surfbb               = surface_b
      endif

      fp_surfbb                 = 0.
      n_emiss_stokessq          = 0


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

           do_atmos_linearization    = .false.

           do_surface_linearization  = .false.
           do_surfbb_linearization   = .false.

           n_totalatmos_wfs          = n_derivs

           do i = 1, n_layers
                layer_vary_flag  (i) = .false.
                layer_vary_number(i) = 0
           enddo

           do i = 1, n_kernels
                do_kernel_factor_wfs(i) = .false.

                do j = 1, kernel_n_params(i)
                     do_kernel_params_wfs(i,j) = .false.
                enddo
           enddo

           i_brdf = 0;

           do j = 1, n_derivs
                do i = 1, n_layers
                     if (derivs(j, i) .ne. 0) then
                          do_atmos_linearization = .true.

                          layer_vary_flag  (i) = .true.

                          layer_vary_number(i) = layer_vary_number(i) + 1

                          jj = layer_vary_number(i);

                          if (omega(i) .ne. 0.) then
                               l_omega_total_input(jj, i) =
     &                              omega_l(j, i) / omega(i)
                          endif

                          if (ltau(i)  .ne. 0.) then
                               l_deltau_vert_input(jj, i) =
     &                               ltau_l(j, i) / ltau(i)
                          endif

                          do k = 0, n_coef_layer(i) - 1
                               do l = 1, 16
                                    l_greekmat_total_input(jj, k, i, l) = 0.
                               enddo

                               if (chi(k+1, 1, i) .ne. 0.) then
                                    l_greekmat_total_input(jj, k, i, 1)  =  chi_l(k+1, 1, j, i) / chi(k+1, 1, i)
                               endif
                               if (n_elem .gt. 1) then
                                    if (chi(k+1, 2, i) .ne. 0.) then
                                         l_greekmat_total_input(jj, k, i, 6)  = chi_l(k+1, 2, j, i) / chi(k+1, 2, i)
                                    endif
                                    if (chi(k+1, 5, i) .ne. 0.) then
                                         l_greekmat_total_input(jj, k, i, 2)  = chi_l(k+1, 5, j, i) / chi(k+1, 5, i)
                                    endif
                                    if (chi(k+1, 5, i) .ne. 0.) then
                                         l_greekmat_total_input(jj, k, i, 5)  = chi_l(k+1, 5, j, i) / chi(k+1, 5, i)
                                    endif
                                    if (chi(k+1, 3, i) .ne. 0.) then
                                         l_greekmat_total_input(jj, k, i, 11) = chi_l(k+1, 3, j, i) / chi(k+1, 3, i)
                                    endif
                                    if (chi(k+1, 4, i) .ne. 0.) then
                                         l_greekmat_total_input(jj, k, i, 16) = chi_l(k+1, 4, j, i) / chi(k+1, 4, i)
                                    endif
                                    if (chi(k+1, 6, i) .ne. 0.) then
                                         l_greekmat_total_input(jj, k, i, 12) = chi_l(k+1, 6, j, i) / chi(k+1, 6, i)
                                    endif
                                    if (chi(k+1, 6, i) .ne. 0.) then
                                         l_greekmat_total_input(jj, k, i, 15) = chi_l(k+1, 6, j, i) / chi(k+1, 6, i)
                                    endif
                               endif
                          enddo

                          do k = k, n_coef - 1
                               do l = 1, 16
                                    l_greekmat_total_input(jj, k, i, l) = 0.
                               enddo
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
                                    write (0, *) 'ERROR: in vlidort each kernel amplitude factor can only be linearized once'
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
                                         write (0, *) 'ERROR: in vlidort each kernel parameter can only be linearized once'
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
                               write (0, *) 'ERROR: vlidort only supports one linearized brdf value per deriv'
                               info = 1
                               return
                          endif

                          do l2 = l + 1, kernel_n_params(k2)
                               if (params_l(l2, j, k2) .ne. 0. .and. l2 .ne. l) then
                                    write (0, *) 'ERROR: vlidort only supports one linearized brdf value per deriv'
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
                          write (0, *) 'ERROR: vlidort only supports one linearized layer/surface per deriv'
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
           call vlidort_master_brdf(status_inputcheck,
     &                              status_calculation)
           if (status_inputcheck .ne. 0) then
                write (0, *) 'ERROR: vlidort_master_brdf(), '//
     &                   'status_inputcheck = ', status_inputcheck
                info = 1
                return
           endif
           if (status_calculation .ne. 0) then
                write (0, *) 'ERROR: vlidort_master_brdf(), '//
     &                   'status_calculation = ', status_calculation
                info = 1
                return
           endif
      else
           call vlidort_l_master_brdf(status_inputcheck,
     &                                status_calculation)
           if (status_inputcheck .ne. 0) then
                write (0, *) 'ERROR: vlidort_l_master_brdf(), '//
     &                   'status_inputcheck = ', status_inputcheck
                info = 1
                return
           endif
           if (status_calculation .ne. 0) then
                write (0, *) 'ERROR: vlidort_l_master_brdf(), '//
     &                   'status_calculation = ', status_calculation
                info = 1
                return
           endif
      endif

      call vlidort_status(status_inputcheck, status_calculation)


c     ******************************************************************
c     *
c     ******************************************************************
      flag = .false.
      if ((nstokes .eq. 1 .and. (n_user_streams .eq. 1 .and.
     &     user_angles_input(1) .eq. zero) .and. (.not. do_ssfull .or.
     &     .not. do_sscorr_nadir .or. .not. do_sscorr_outgoing))
     &     .or. do_mvout_only) then
           flag = .true.
      endif

      if (radiance .ne. 0) then
           do i = 1, n_user_streams
                i1 = quad_index(i)
                do j = 1, n_stokes
                     do k = 1, n_phis
                          if (flag) then
                               i2 = (i - 1) * n_phis + 1
                          else
                               i2 = (i - 1) * n_phis + k
                          endif
                          do l = 1, n_out_usertaus
                               I_p(j, k, i1, l) = stokes(l, i2, j, 1)
                               I_m(j, k, i1, l) = stokes(l, i2, j, 2)
                          enddo
                     enddo
                enddo
           enddo

           if (n_derivs .gt. 0) then
                do l = 1, n_totalatmos_wfs
                     if (deriv_type(l) .eq. 0) then
                          do i = 1, n_user_streams
                               i1 = quad_index(i)
                               do j = 1, n_stokes
                                    do k = 1, n_phis
                                         do m = 1, n_out_usertaus
                                              K_p(j, k, i1, l, m) = 0.
                                              K_m(j, k, i1, l, m) = 0.
                                         enddo
                                    enddo
                               enddo
                          enddo
                     else if (deriv_type(l) .eq. 1) then
                          i_layer = layers_index(l);
                          i_deriv = derivs_index(l);
                          do i = 1, n_user_streams
                               i1 = quad_index(i)
                               do j = 1, n_stokes
                                    do k = 1, n_phis
                                         if (flag) then
                                              i2 = (i - 1) * n_phis + 1
                                         else
                                              i2 = (i - 1) * n_phis + k
                                         endif
                                         do m = 1, n_out_usertaus
                                              K_p(j, k, i1, l, m) =
     &                                             atmoswf(i_deriv, i_layer, m, i2, j, 1)
                                              K_m(j, k, i1, l, m) =
     &                                             atmoswf(i_deriv, i_layer, m, i2, j, 2)
                                         enddo
                                    enddo
                               enddo
                          enddo
                     else if (deriv_type(l) .eq. 2) then
                          i_brdf = brdf_index(l)
                          do i = 1, n_user_streams
                               i1 = quad_index(i)
                               do j = 1, n_stokes
                                    do k = 1, n_phis
                                         if (flag) then
                                              i2 = (i - 1) * n_phis + 1
                                         else
                                              i2 = (i - 1) * n_phis + k
                                         endif
                                         do m = 1, n_out_usertaus
                                              K_p(j, k, i1, l, m) =
     &                                             surfacewf(i_brdf, m, i2, j, 1) / brdf_factor(i_brdf)
                                              K_m(j, k, i1, l, m) =
     &                                             surfacewf(i_brdf, m, i2, j, 2) / brdf_factor(i_brdf)
                                         enddo
                                    enddo
                               enddo
                          enddo
                     endif
                enddo
           endif
      endif

      if (mean .ne. 0) then
           call vlidort_copy_mean_value
     &          (n_derivs, n_out_usertaus, n_totalatmos_wfs,
     &           deriv_type, layers_index, derivs_index, brdf_index,
     &           mean_stokes, mint_atmoswf, mint_surfacewf,
     &           mean_p, mean_m, mean_p_l, mean_m_l)
      endif

      if (flux .ne. 0) then
           call vlidort_copy_mean_value
     &          (n_derivs, n_out_usertaus, n_totalatmos_wfs,
     &           deriv_type, layers_index, derivs_index, brdf_index,
     &           flux_stokes, flux_atmoswf, flux_surfacewf,
     &           flux_p, flux_m, flux_p_l, flux_m_l)
      endif

c     print *, fourier_saved(1);

      return

      end



c***********************************************************************
c
c***********************************************************************
      subroutine vlidort_copy_mean_value
     &     (n_derivs, n_out_usertaus, n_totalatmos_wfs, deriv_type,
     &      layers_index, derivs_index, brdf_index, lidort_mean,
     &      lidort_atmos_l, lidort_surface_l, mean_p, mean_m,
     &      mean_p_l, mean_m_l)

      implicit none

      include 'VLIDORT.PARS'

      integer n_derivs
      integer n_out_usertaus
      integer n_totalatmos_wfs

      byte    deriv_type(n_derivs)
      integer layers_index(n_derivs)
      integer derivs_index(n_derivs)
      integer brdf_index(n_derivs)

      real*8 lidort_mean(MAX_OUT_USERTAUS, MAXBEAMS,
     &                   MAXSTOKES, MAX_DIRECTIONS)
      real*8 lidort_atmos_l  (MAX_ATMOSWFS, MAXLAYERS, MAX_OUT_USERTAUS,
     &                        MAXBEAMS, MAXSTOKES, MAX_DIRECTIONS)
      real*8 lidort_surface_l(MAX_SURFACEWFS, MAX_OUT_USERTAUS,
     &                        MAXBEAMS, MAXSTOKES, MAX_DIRECTIONS)

      real*8 mean_p(n_out_usertaus)
      real*8 mean_m(n_out_usertaus)

      real*8 mean_p_l(n_derivs, n_out_usertaus)
      real*8 mean_m_l(n_derivs, n_out_usertaus)


      integer k
      integer l
      integer i_layer
      integer i_deriv
      integer i_brdf


      do l = 1, n_out_usertaus
           mean_p(l) = lidort_mean(l, 1, 1, 1)
           mean_m(l) = lidort_mean(l, 1, 1, 2)
      enddo

      do k = 1, n_totalatmos_wfs
           if (deriv_type(k) .eq. 0) then
                i_layer = layers_index(k);
                i_deriv = derivs_index(k);
                do l = 1, n_out_usertaus
                     mean_p_l(k, l) =
     &                    lidort_atmos_l(i_deriv, i_layer, l, 1, 1, 1)
                     mean_m_l(k, l) =
     &                    lidort_atmos_l(i_deriv, i_layer, l, 1, 1, 2)
                enddo
           else
                i_brdf = brdf_index(k)
                do l = 1, n_out_usertaus
                     mean_p_l(k, l) =
     &                    lidort_surface_l(i_brdf, l, 1, 1, 1)
                     mean_m_l(k, l) =
     &                    lidort_surface_l(i_brdf, l, 1, 1, 2)
                enddo
          endif
      enddo

      return

      end
