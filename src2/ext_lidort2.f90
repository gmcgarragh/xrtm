!*******************************************************************************
!
!    Copyright (C) 2007-2020 Greg McGarragh <greg.mcgarragh@colostate.edu>
!
!    This source code is licensed under the GNU General Public License (GPL),
!    Version 3.  See the file COPYING for more details.
!
!*******************************************************************************

subroutine call_lidort2_f(n_four, n_elem, n_coef, n_quad, n_derivs, n_layers, qx, F_0, &
                           theta_0, phi_0, ulevels, utaus, n_ulevels, umus, n_umus, &
                           phis, n_phis, planet_r,levels_z, n_kernels, n_kernel_quad, &
                           kernels, ampfac, ampfac_l, params, params_l, n_coef_layer, &
                           chi, chi_l,omega, omega_l, ltau, ltau_l, I_p,I_m, K_p, K_m, &
                           mean_p,mean_m, mean_p_l, mean_m_l, flux_p, flux_m, flux_p_l, &
                           flux_m_l, delta_m, n_t_tms, psa, quad_output, radiance, &
                           mean, flux, utau_output, derivs, epsilon, info, n_quad2)

     use brdf_sup_inputs_def_m
     use lidort_inputs_def_m
     use lidort_io_defs_m
     use lidort_masters_m


     implicit none


     integer, parameter   :: max_params = 4

     integer, intent(in)  :: n_four
     integer, intent(in)  :: n_elem
     integer, intent(in)  :: n_coef
     integer, intent(in)  :: n_quad
     integer, intent(in)  :: n_derivs
     integer, intent(in)  :: n_layers

     real(8), intent(in)  :: qx(n_quad)

     real(8), intent(in)  :: F_0

     real(8), intent(in)  :: theta_0

     real(8), intent(in)  :: phi_0

     integer, intent(in)  :: ulevels(n_ulevels)
     real(8), intent(in)  :: utaus(n_ulevels)
     integer, intent(in)  :: n_ulevels

     real(8), intent(in)  :: umus(n_umus)
     integer, intent(in)  :: n_umus

     real(8), intent(in)  :: phis(n_phis)
     integer, intent(in)  :: n_phis

     real(8), intent(in)  :: planet_r
     real(8), intent(in)  :: levels_z(n_layers+1)

     integer, intent(in)  :: n_kernels
     integer, intent(in)  :: n_kernel_quad
     integer, intent(in)  :: kernels(n_kernels)
     real(8), intent(in)  :: ampfac(n_kernels)
     real(8), intent(in)  :: ampfac_l(n_derivs, n_kernels)
     real(8), intent(in)  :: params(max_params, n_kernels)
     real(8), intent(in)  :: params_l(max_params, n_derivs, n_kernels)

     integer, intent(in)  :: n_coef_layer(n_layers)
     real(8), intent(in)  :: chi  (n_coef, n_elem, n_layers)
     real(8), intent(in)  :: chi_l(n_coef, n_elem, n_derivs, n_layers)

     real(8), intent(in)  :: omega(n_layers)
     real(8), intent(in)  :: omega_l(n_derivs, n_layers)

     real(8), intent(in)  :: ltau(n_layers)
     real(8), intent(in)  :: ltau_l(n_derivs, n_layers)

     real(8), intent(out) :: I_p(n_phis, n_quad2, n_ulevels)
     real(8), intent(out) :: I_m(n_phis, n_quad2, n_ulevels)

     real(8), intent(out) :: K_p(n_phis, n_quad2, n_derivs, n_ulevels)
     real(8), intent(out) :: K_m(n_phis, n_quad2, n_derivs, n_ulevels)

     real(8), intent(out) :: mean_p(n_ulevels)
     real(8), intent(out) :: mean_m(n_ulevels)

     real(8), intent(out) :: mean_p_l(n_derivs, n_ulevels)
     real(8), intent(out) :: mean_m_l(n_derivs, n_ulevels)

     real(8), intent(out) :: flux_p(n_ulevels)
     real(8), intent(out) :: flux_m(n_ulevels)

     real(8), intent(out) :: flux_p_l(n_derivs, n_ulevels)
     real(8), intent(out) :: flux_m_l(n_derivs, n_ulevels)

     integer, intent(in)  :: delta_m
     integer, intent(in)  :: n_t_tms
     integer, intent(in)  :: psa
     integer, intent(in)  :: quad_output
     integer, intent(in)  :: radiance
     integer, intent(in)  :: mean
     integer, intent(in)  :: flux
     integer, intent(in)  :: utau_output

     integer*1 derivs(n_derivs + 8, n_layers)

     real(8), intent(in)  :: epsilon

     integer, intent(out) :: info

     integer, intent(in)  :: n_quad2


     !**************************************************************************
     !
     !**************************************************************************
     logical flag

     integer max_kernels_xrtm
     parameter (max_kernels_xrtm = 10)
     integer max_kernels_lidort
     parameter (max_kernels_lidort = 9)

     character(10) :: kernel_names(max_kernels_lidort) = (/'Lambertian', &
                                                           'Ross-thin ', &
                                                           'Ross-thick', &
                                                           'Li-sparse ', &
                                                           'Li-dense  ', &
                                                           'Hapke     ', &
                                                           'Roujean   ', &
                                                           'Rahman    ', &
                                                           'Cox-Munk  '/)
     integer :: i
     integer :: i1
     integer :: i2
     integer :: j
     integer :: k
     integer :: l

     integer :: i_kernel

!    integer :: i_layer
!    integer :: i_deriv
!    integer :: i_brdf

!    integer :: status_inputcheck
!    integer :: status_calculation

     integer quad_index(n_quad2)

     integer :: kernel_index(max_kernels_xrtm) = (/1,0,7,4,5,2,3,6,8,9/)

     integer :: kernel_n_params(max_kernels_lidort) = (/0,0,0,2,2,3,0,3,2/)

     real(8) :: a

     real(8) :: pi

     real(8) :: angles(n_quad2)

     type(lidort_fixed_inputs)    :: fixed
     type(lidort_modified_inputs) :: modified
     type(brdf_sup_inputs)        :: brdf_sup
     type(lidort_sup_inout)       :: sup_inout
     type(lidort_outputs)         :: outputs


     !**************************************************************************
     !
     !**************************************************************************
     i = n_four
     i = n_elem
     i = n_coef
     i = n_quad
     i = n_derivs
     i = n_layers

     a = qx(n_quad)

     a = F_0

     a = theta_0

     a = phi_0

     i = ulevels(1)
     a = utaus(1)
     i = n_ulevels

     a = umus(1)
     i = n_umus

     a = phis(1)
     i = n_phis

     a = planet_r
     a = levels_z(1)

     i = n_kernels
     i = n_kernel_quad
     i = kernels(1)
     a = ampfac(1)
     a = ampfac_l(1, 1)
     a = params(1, 1)
     a = params_l(1, 1, 1)

     i = n_coef_layer(1)
     a = chi  (1, 1, 1)
     a = chi_l(1, 1, 1, 1)

     a = omega(1)
     a = omega_l(1, 1)

     a = ltau(1)
     a = ltau_l(1, 1)

     i = delta_m
     i = n_t_tms
     i = psa
     i = quad_output
     i = radiance
     i = mean
     i = flux
     i = utau_output

     i = derivs(1, 1)

     a = epsilon

     i = n_quad2


     !**************************************************************************
     !
     !**************************************************************************
!    if (delta_m .ne. 0) then
!         do i = 1, n_layers
!              if (chi(n_coef+1, 1, i) .eq. 0.) then
!                   write (0, *) 'ERROR: lidort fails if using'     // &
!                        ' delta-m scaling and the truncation'      // &
!                        ' factor (last + 1 phase function momment)'// &
!                        ' is = zero'
!                   info = 1
!                   return
!              endif
!         enddo
!    endif


     !**************************************************************************
     !
     !**************************************************************************
     info = 0

     pi = 4.d0 * atan(1.0d0)


     !**************************************************************************
     !
     !**************************************************************************
     fixed%bool%ts_do_fullrad_mode                = .true.

     fixed%bool%ts_do_thermal_emission            = .false.
     fixed%bool%ts_do_surface_emission            = .false.

     if (psa .eq. 0) then
          fixed%bool%ts_do_plane_parallel         = .true.
     else
          fixed%bool%ts_do_plane_parallel         = .false.
     endif

     if (n_kernels .eq. 1 .and. kernels(1) .eq. 0) then
          fixed%bool%ts_do_brdf_surface           = .false.
     else
          fixed%bool%ts_do_brdf_surface           = .true.
     endif

     fixed%bool%ts_do_upwelling                   = .true.
     fixed%bool%ts_do_dnwelling                   = .true.

     fixed%bool%ts_do_surface_leaving             = .false.
     fixed%bool%ts_do_sl_isotropic                = .false.

     fixed%bool%ts_do_water_leaving               = .false.
     fixed%bool%ts_do_fluorescence                = .false.

     fixed%bool%ts_do_tf_iteration                = .false.

     fixed%bool%ts_do_wladjusted_output           = .false.

     fixed%bool%ts_do_toa_illumination            = .false.
     fixed%bool%ts_do_boa_illumination            = .false.

     fixed%bool%ts_do_albtrn_media(1)             = .false.
     fixed%bool%ts_do_albtrn_media(2)             = .false.

     fixed%bool%ts_do_planetary_problem           = .false.


     !**************************************************************************
     !
     !**************************************************************************
     fixed%cont%ts_taylor_order                   = 2

     fixed%cont%ts_nstreams                       = n_quad

     fixed%cont%ts_nlayers                        = n_layers

!    fixed%cont%ts_nfinelayers

     fixed%cont%ts_n_thermal_coeffs               = 2

     fixed%cont%ts_lidort_accuracy                = epsilon

!    fixed%cont%ts_tf_maxiter
!    fixed%cont%ts_tf_criterion

     fixed%cont%ts_toa_illumination               = 0.
     fixed%cont%ts_boa_illumination               = 0.


     !**************************************************************************
     !
     !**************************************************************************
     fixed%sunrays%ts_flux_factor                 = F_0


     !**************************************************************************
     !
     !**************************************************************************
     fixed%userval%ts_n_user_levels               = n_ulevels


     !**************************************************************************
     !
     !**************************************************************************
     if (psa .ne. 0) then
          do i = 0, n_layers
               fixed%chapman%ts_height_grid(i) = levels_z(i+1)
          enddo
     endif
!    fixed%chapman%ts_pressure_grid   (0:maxlayers)
!    fixed%chapman%ts_temperature_grid(0:maxlayers)

!    fixed%chapman%ts_finegrid(maxlayers)

     fixed%chapman%ts_rfindex_parameter           = 0.


     !**************************************************************************
     !
     !**************************************************************************
     do i = 1, n_layers
          fixed% optical%ts_deltau_vert_input(i) = ltau(i)
     enddo

     do j = 1, n_layers
          do i = 0, n_coef_layer(j) - 1
               fixed%optical%ts_phasmoms_total_input(i, j) = chi(i+1, 1, j)
          enddo

          do i = i, n_coef - 1
               fixed%optical%ts_phasmoms_total_input(i, j) = 0.
          enddo
     enddo

!    fixed%optical%ts_phasfunc_input_up
!    fixed%optical%ts_phasfunc_input_dn

     if (n_kernels .eq. 1 .and. kernels(1) .eq. 0) then
          fixed%optical%ts_lambertian_albedo = ampfac(1)
     endif

!    fixed%optical%ts_surface_bb_input

     fixed%optical%ts_atmos_wavelength = 0.


     !**************************************************************************
     !
     !**************************************************************************
     fixed%write%ts_do_debug_write                = .false.

     fixed%write%ts_do_write_input                = .false.
     fixed%write%ts_do_write_scenario             = .false.
     fixed%write%ts_do_write_fourier              = .false.
     fixed%write%ts_do_write_results              = .false.

     fixed%write%ts_input_write_filename          = 'lidort_input.out'
     fixed%write%ts_scenario_write_filename       = 'lidort_scenario.out'
     fixed%write%ts_fourier_write_filename        = 'lidort_fourier.out'
     fixed%write%ts_results_write_filename        = 'lidort_results.out'


     !**************************************************************************
     !
     !**************************************************************************
     modified%mbool%ts_do_focorr                  = .false.

     modified%mbool%ts_do_focorr_external         = .false.

     modified%mbool%ts_do_focorr_nadir            = .false.
     modified%mbool%ts_do_focorr_outgoing         = .false.

     modified%mbool%ts_do_sscorr_truncation       = .false.

     modified%mbool%ts_do_sscorr_usephasfunc      = .false.

     modified%mbool%ts_do_double_convtest         = .true.

     modified%mbool%ts_do_solar_sources           = .true.

     modified%mbool%ts_do_refractive_geometry     = .false.
     modified%mbool%ts_do_chapman_function        = .true.

     modified%mbool%ts_do_rayleigh_only           = .false.
     modified%mbool%ts_do_isotropic_only          = .false.
     if (n_four .eq. 1) then
          modified%mbool%ts_do_no_azimuth         = .true.
     else
          modified%mbool%ts_do_no_azimuth         = .false.
     endif
     modified%mbool%ts_do_all_fourier             = .false.

     if (delta_m .eq. 0) then
          modified%mbool%ts_do_deltam_scaling     = .false.
     else
          modified%mbool%ts_do_deltam_scaling     = .true.
     endif

     modified%mbool%ts_do_solution_saving         = .false.
     modified%mbool%ts_do_bvp_telescoping         = .false.

     if (quad_output .ne. 0) then
          modified%mbool%ts_do_user_streams       = .false.
     else
          modified%mbool%ts_do_user_streams       = .true.
     endif

     modified%mbool%ts_do_additional_mvout        = .false.
     modified%mbool%ts_do_mvout_only              = .false.

     if (radiance .ne. 0) then
          if (mean .eq. 0 .and. flux .eq. 0) then
               modified%mbool%ts_do_additional_mvout = .false.
          else
               modified%mbool%ts_do_additional_mvout = .true.
          endif
     else
          if (mean .eq. 0 .and. flux .eq. 0) then
               modified%mbool%ts_do_mvout_only       = .false.
          else
               modified%mbool%ts_do_mvout_only       = .true.
          endif
     endif

     modified%mbool%ts_do_thermal_transonly       = .false.

     modified%mbool%ts_do_observation_geometry    = .false.

     modified%mbool%ts_do_external_wleave         = .false.


     !**************************************************************************
     !
     !**************************************************************************
     modified%mcont%ts_nmoments_input             = n_coef - 1


     !**************************************************************************
     !
     !**************************************************************************
     modified%msunrays%ts_nbeams                  = 1

     modified%msunrays%ts_beam_szas(1)            = theta_0


     !**************************************************************************
     !
     !**************************************************************************
     modified%muserval%ts_n_user_relazms = n_phis
     do i = 1, n_phis
          modified%muserval%ts_user_relazms(i) = phis(i) - phi_0
     enddo

     a = 180.d0 / pi
     if (quad_output .ne. 0) then
          modified%muserval%ts_n_user_streams = n_quad
          do i = 1, n_quad
               quad_index(i) = n_quad - i + 1
               modified%muserval%ts_user_angles_input(i) = acos(qx(quad_index(i))) * a
          enddo
     else
          if (n_umus .eq. 0) then
               modified%muserval%ts_n_user_streams = n_quad
               do i = 1, n_quad
                    quad_index(i) = n_quad - i + 1
                    modified%muserval%ts_user_angles_input(i) = acos(qx(quad_index(i))) * a
               enddo
          else
               modified%muserval%ts_n_user_streams = n_umus
               do i = 1, n_umus
                    angles(i) = acos(umus(i)) * a
               enddo

               if (n_umus .eq. 1) then
                    quad_index(1) = 1
               else
                    call lidort_indexx(n_umus, angles, quad_index)
               endif

               do i = 1, n_umus
                    modified%muserval%ts_user_angles_input(i) = angles(quad_index(i))
               enddo
          endif
     endif

     if (utau_output .eq. 0) then
          do i = 1, n_ulevels
               modified%muserval%ts_user_levels(i) = ulevels(i)
          enddo
     else
          print *, 'ERROR: lidort does not support utau output'
          info = 1
          return
     endif

!    modified%muserval%ts_geometry_specheight

!    modified%muserval%ts_n_user_obsgeoms

!    modified%muserval%ts_user_obsgeom_input


     !**************************************************************************
     !
     !**************************************************************************
     modified%mchapman%ts_earth_radius            = planet_r


     !**************************************************************************
     !
     !**************************************************************************
     do i = 1, n_layers
          modified%moptical%ts_omega_total_input(i) = omega(i)
     enddo


     !**************************************************************************
     !
     !**************************************************************************
     if (n_kernels .eq. 1 .and. kernels(1) .eq. 0) then

     else
          brdf_sup%bs_do_user_streams            = .true.

          brdf_sup%bs_do_surface_emission        = .false.

          brdf_sup%bs_do_solar_sources           = .true.
          brdf_sup%bs_do_user_obsgeoms           = .false.

          brdf_sup%bs_nstreams = n_kernel_quad * 2

          brdf_sup%bs_nbeams                     = 1

          brdf_sup%bs_beam_szas(1)               = theta_0

          brdf_sup%bs_n_user_relazms = n_phis
          do i = 1, n_phis
               brdf_sup%bs_user_relazms(i) = phis(i) - phi_0
          enddo

          brdf_sup%bs_n_user_streams = n_umus
          do i = 1, n_umus
               angles(i) = acos(umus(i)) * a
          enddo

          if (n_umus .eq. 1) then
               quad_index(1) = 1
          else
               call lidort_indexx(n_umus, angles, quad_index)
          endif

          do i = 1, n_umus
               brdf_sup%bs_user_angles_input(i) = angles(quad_index(i))
          enddo

!         brdf_sup%bs_n_user_obsgeoms
!         brdf_sup%bs_user_obsgeom_input

          brdf_sup%bs_n_brdf_kernels = n_kernels

          do i = 1, n_kernels
               i_kernel = kernel_index(kernels(i) + 1)
               if (i_kernel .eq. 0) then
                    print *, 'ERROR: unsupported kernel for lidort'
                    info = 1
                    return
               endif
               brdf_sup%bs_brdf_names(i) = kernel_names(i_kernel)

               brdf_sup%bs_which_brdf(i)              = i_kernel

               brdf_sup%bs_n_brdf_parameters(i)       = kernel_n_params(i_kernel)

               do j = 1, kernel_n_params(i_kernel)
                    brdf_sup%bs_brdf_parameters(i, j) = params(j, i)
               enddo

               if (i_kernel .ne. 1) then
                    brdf_sup%bs_lambertian_kernel_flag(i) = .false.
               else
                    brdf_sup%bs_lambertian_kernel_flag(i) = .true.
               endif

               brdf_sup%bs_brdf_factors(i)            = ampfac(i)
          enddo
     endif


     !**************************************************************************
     !
     !**************************************************************************
     if (n_derivs .eq. 0) then
     else
     endif


     !**************************************************************************
     !
     !**************************************************************************
     if (n_derivs .eq. 0) then
          call lidort_master(fixed, modified, sup_inout, outputs)
!         if (status_inputcheck .eq. LIDORT_SERIOUS) then
!              print *, 'ERROR: lidort_master(), status_inputcheck = ', &
!                       status_inputcheck
!              info = 1
!              return
!         endif
!         if (status_calculation .eq. LIDORT_SERIOUS) then
!              print *, 'ERROR: lidort_master(), status_calculation = ', &
!                       status_calculation
!              info = 1
!              return
!         endif
     else
          call lidort_l_master(fixed, modified, sup_inout, outputs)
!         if (status_inputcheck .eq. LIDORT_SERIOUS) then
!              print *, 'ERROR: lidort_l_master(), status_inputcheck = ', &
!                       status_inputcheck
!              info = 1
!              return
!         endif
!         if (status_calculation .eq. LIDORT_SERIOUS) then
!              print *, 'ERROR: lidort_l_master(), status_calculation = ', &
!                       status_calculation
!              info = 1
!              return
!         endif
     endif


     !**************************************************************************
     !
     !**************************************************************************
     K_p = 0.
     K_m = 0.

     if (mean .ne. 0) then
          mean_p   = 0.
          mean_m   = 0.

          mean_p_l = 0.
          mean_m_l = 0.
     endif

     if (mean .ne. 0) then
          flux_p   = 0.
          flux_m   = 0.

          flux_p_l = 0.
          flux_m_l = 0.
     endif


     !**************************************************************************
     !
     !**************************************************************************
     flag = .false.
!     if ((n_user_streams .eq. 1 .and. user_angles_input(1) .eq. zero) &
!          .or. do_mvout_only) then
!          flag = .true.
!     endif

     if (radiance .ne. 0) then
          do i = 1, modified%muserval%ts_n_user_streams
               i1 = quad_index(i)
               do k = 1, n_phis
                    if (flag) then
                         i2 = (i - 1) * n_phis + 1
                    else
                         i2 = (i - 1) * n_phis + k
                    endif
                    do l = 1, fixed%userval%ts_n_user_levels
                         I_p(k, i1, l) = outputs%main%ts_intensity(l, i2, 1)
                         I_m(k, i1, l) = outputs%main%ts_intensity(l, i2, 2)
                    enddo
               enddo
          enddo

!          if (n_derivs .gt. 0) then
!               do k = 1, n_totalatmos_wfs
!                    if (deriv_type(k) .eq. 0) then
!                         do i = 1, n_user_streams
!                              i1 = quad_index(i)
!                              do j = 1, n_phis
!                                   do l = 1, n_out_usertaus
!                                        K_p(j, i1, k, l) = 0.
!                                        K_m(j, i1, k, l) = 0.
!                                   enddo
!                              enddo
!                         enddo
!                    else if (deriv_type(k) .eq. 1) then
!                         i_layer = layers_index(k);
!                         i_deriv = derivs_index(k);
!                         do i = 1, n_user_streams
!                              i1 = quad_index(i)
!                              do j = 1, n_phis
!                                   if (flag) then
!                                        i2 = (i - 1) * n_phis + 1
!                                   else
!                                        i2 = (i - 1) * n_phis + j
!                                   endif
!                                   do l = 1, n_out_usertaus
!                                        K_p(j, i1, k, l) = &
!                                        atmoswf(i_deriv, i_layer, l, i2, 1)
!                                        K_m(j, i1, k, l) = &
!                                        atmoswf(i_deriv, i_layer, l, i2, 2)
!                                   enddo
!                              enddo
!                         enddo
!                    else if (deriv_type(k) .eq. 2) then
!                         i_brdf = brdf_index(k)
!                         do i = 1, n_user_streams
!                              i1 = quad_index(i)
!                              do j = 1, n_phis
!                                   if (flag) then
!                                        i2 = (i - 1) * n_phis + 1
!                                   else
!                                        i2 = (i - 1) * n_phis + j
!                                   endif
!                                   do l = 1, n_out_usertaus
!                                        K_p(j, i1, k, l) = &
!                                        surfacewf(i_brdf, l, i2, 1) / brdf_factor(i_brdf)
!                                        K_m(j, i1, k, l) = &
!                                        surfacewf(i_brdf, l, i2, 2) / brdf_factor(i_brdf)
!                                   enddo
!                              enddo
!                         enddo
!                    endif
!               enddo
!          endif
     endif


     return

end subroutine call_lidort2_f
