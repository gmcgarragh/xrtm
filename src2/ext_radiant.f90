!*******************************************************************************
!
!    Copyright (C) 2007-2020 Greg McGarragh <greg.mcgarragh@colostate.edu>
!
!    This source code is licensed under the GNU General Public License (GPL),
!    Version 3.  See the file COPYING for more details.
!
!*******************************************************************************

subroutine call_radiant_f(n_four, n_elem, n_coef, n_quad, n_derivs, n_layers, &
                          qx, lambda, F_0, theta_0, phi_0, utaus, n_ulevels, &
                          umus, n_umus, phis, n_phis, top_t, planet_r, levels_z, &
                          levels_t, surface_t, n_kernels, n_kernel_quad, kernels, &
                          ampfac, ampfac_l, params, params_l, n_coef_layer, chi, &
                          chi_l, omega, omega_l, ltau, ltau_l, I_p,I_m, K_p, K_m, &
                          flux_p, flux_m, quadrature, delta_m, n_t_tms, psa, radiance, &
                          flux, thermal, utau_output, derivs, epsilon, info, n_quad2)

     use binning_banks4
     use radiant2z3
     use radiant2z3_io_defs


     implicit none


     integer, parameter     :: max_params = 4

     integer, intent(in)    :: n_four
     integer, intent(in)    :: n_elem
     integer, intent(in)    :: n_coef
     integer, intent(in)    :: n_quad
     integer, intent(in)    :: n_derivs
     integer, intent(in)    :: n_layers

     real(8), intent(in)    :: qx(n_quad)

     real(8), intent(in)    :: lambda

     real(8), intent(in)    :: F_0

     real(8), intent(in)    :: theta_0

     real(8), intent(in)    :: phi_0

     real(8), intent(in)    :: utaus(n_ulevels)
     integer, intent(in)    :: n_ulevels

     real(8), intent(in)    :: umus(n_umus)
     integer, intent(in)    :: n_umus

     real(8), intent(in)    :: phis(n_phis)
     integer, intent(in)    :: n_phis

     real(8), intent(in)    :: top_t

     real(8), intent(in)    :: planet_r
     real(8), intent(in)    :: levels_z(n_layers+1)
     real(8), intent(in)    :: levels_t(n_layers+1)

     real(8), intent(in)    :: surface_t

     integer, intent(in)    :: n_kernels
     integer, intent(in)    :: n_kernel_quad
     integer, intent(in)    :: kernels(n_kernels)
     real(8), intent(in)    :: ampfac(n_kernels)
     real(8), intent(in)    :: ampfac_l(n_derivs, n_kernels)
     real(8), intent(in)    :: params(max_params, n_kernels)
     real(8), intent(in)    :: params_l(max_params, n_derivs, n_kernels)

     integer, intent(in)    :: n_coef_layer(n_layers)
     real(8), intent(in)    :: chi(n_coef, n_elem, n_layers)
     real(8), intent(in)    :: chi_l(n_coef, n_elem, n_derivs, n_layers)

     real(8), intent(in)    :: omega(n_layers)
     real(8), intent(in)    :: omega_l(n_derivs, n_layers)

     real(8), intent(in)    :: ltau(n_layers)
     real(8), intent(in)    :: ltau_l(n_derivs, n_layers)

     real(8), intent(out)   :: I_p(n_phis, n_quad2, n_layers+1)
     real(8), intent(out)   :: I_m(n_phis, n_quad2, n_layers+1)

     real(8), intent(out)   :: K_p(n_phis, n_quad2, n_derivs, n_layers+1)
     real(8), intent(out)   :: K_m(n_phis, n_quad2, n_derivs, n_layers+1)

     real(8), intent(out)   :: flux_p(n_layers+1)
     real(8), intent(out)   :: flux_m(n_layers+1)

     integer, intent(in)    :: quadrature

     integer, intent(in)    :: delta_m
     integer, intent(in)    :: n_t_tms
     integer, intent(in)    :: psa
     integer, intent(in)    :: radiance
     integer, intent(in)    :: flux
     integer, intent(in)    :: thermal
     integer, intent(in)    :: utau_output

     integer(1), intent(in) :: derivs(n_derivs + 8, n_layers)

     real(8), intent(in)    :: epsilon

     integer, intent(out)   :: info

     integer, intent(in)    :: n_quad2


     !**************************************************************************
     !
     !**************************************************************************
     integer, parameter     :: max_kernels_xrtm    = 10
     integer, parameter     :: max_kernels_radiant = 9

     integer                :: i
     integer                :: j
     integer                :: k

     integer                :: i_kernel

     integer                :: r

     integer                :: i_layer
     integer                :: i_deriv
     integer                :: i_brdf

     integer                :: kernel_index   (max_kernels_xrtm)    = (/1,0,7,4,5,2,3,6,8,9/)

     integer                :: kernel_n_params(max_kernels_radiant) = (/0,0,0,2,2,3,0,3,2/)

     byte                   :: deriv_type(n_derivs)
     integer                :: layers_index(n_derivs)
     integer                :: derivs_index(n_derivs)
     integer                :: brdf_index(n_derivs)

     real(8)                :: a

     real(8)                :: pi

     type (radiant_control) :: rt_con
     type (planetary_scene) :: scene
     type (jacobian_inputs) :: jac
     type (radiant_outputs) :: rt_out


     !**************************************************************************
     !
     !**************************************************************************
     if (n_umus > 1) then
          write (0, *) 'ERROR: radiant can only handle one user zenith angle'
          info = 1
          return
     endif

     if (n_phis > 1) then
          write (0, *) 'ERROR: radiant can only handle one user azimuth angle'
          info = 1
          return
     endif



     !**************************************************************************
     !
     !**************************************************************************
     a = qx(1)
     a = phi_0

     if (n_derivs > 0 .and. n_kernels > 0) then
          a = ampfac_l(1, 1)
          if (kernel_n_params(1) > 0) then
               a = params_l(1, 1, 1)
          endif
     endif

     info = 0

     pi = 4.d0 * atan(1.0d0)


     !**************************************************************************
     !
     !**************************************************************************
     r = radiant_datatype_ctrl2(1, rt_con, scene, jac, rt_out, &
                                n_quad * 2, n_coef - 1, n_layers, max(1, n_derivs), n_ulevels)
     if (r .ne. 0) then
          write (0, *) 'ERROR: radiant_datatype_ctrl2(1, ...), status = ', r
          info = 1
          return
     endif


     !**************************************************************************
     !
     !**************************************************************************
     rt_con%error_output_lvl        = 0
     rt_con%errfile_unit            = 6
     rt_con%dbgfile_unit            = 6
     rt_con%use_infile              = .false.
     rt_con%use_outfile             = .false.
     rt_con%user_defined_filenames  = .false.
     rt_con%errfile_name            = 'radiant_error.out'
     rt_con%dbgfile_name            = 'radiant_debug.out'
     rt_con%infile_name             = 'radiant_infile.out'
     rt_con%outfile_name            = 'radiant_outfile.out'

     rt_con%streams                 = n_quad * 2
     rt_con%quadrature              = quadrature

     if (F_0 .ne. 0. .and. thermal .eq. 0) then
          rt_con%sources                 = 1
     else if (F_0 .eq. 0. .and. thermal .ne. 0) then
          rt_con%sources                 = 3
     else
          rt_con%sources                 = 2
     endif

     if (n_umus .eq. 0) then
          rt_con%apply_user_zenith_angle = .false.
          rt_con%user_zenith_angle       = 0.
     else
          rt_con%apply_user_zenith_angle = .true.
          rt_con%user_zenith_angle       = acos(umus(1))*180.d0/pi
     endif

     if (utau_output .eq. 0) then
          rt_con%get_user_rad            = .false.
          rt_con%n_user_tautot           = n_ulevels
     else
          rt_con%get_user_rad            = .true.
          rt_con%n_user_tautot           = n_ulevels
          do i = 1, n_ulevels
               rt_con%user_tautot(i)     = utaus(i)
          enddo
     endif

     if (n_four .eq. 0) then
          rt_con%azimuthal_rad_only = .true.
     else
          rt_con%azimuthal_rad_only = .false.
     endif
     if (delta_m .eq. 0) then
          rt_con%delta_m            = .false.
     else
          rt_con%delta_m            = .true.
     endif
     if (n_t_tms .eq. 0) then
          rt_con%ss_cor             = .false.
     else
          rt_con%ss_cor             = .true.
     endif
     rt_con%get_rad_direct          = .false.
     rt_con%get_rad_diffuse         = .true.
     if (flux .eq. 0) then
          rt_con%get_fluxes         = .false.
     else
          rt_con%get_fluxes         = .true.
     endif

     rt_con%fourier_tol             = epsilon

     rt_con%bin%use_layer_bank      = .true.
     rt_con%bin%use_column_bank     = .false.
     rt_con%bin%reset_bank          = .false.

     rt_con%bin%upper_blk_bnd       = 0
     rt_con%bin%upper_blk_col       = 0
     rt_con%bin%lower_blk_bnd       = 0
     rt_con%bin%lower_blk_col       = 0


     !**************************************************************************
     !
     !**************************************************************************
     scene%fsun                        = F_0
     scene%sza                         = theta_0
     scene%itms(:)                     = 0.
     scene%ss_itms(:)                  = 0.

     if (thermal .eq. 0) then
          scene%ttop                   = 0.
          scene%temiss                 = 0.
          scene%tlev(:)                = 0.
          scene%tsurf                  = 0.
          scene%planck_type            = 1
          scene%wvn                    = 0.
          scene%wvnlo                  = 0.
          scene%wvnhi                  = 0.
     else
          scene%ttop                   = top_t
          scene%temiss                 = 1.
          do i = 1, n_layers + 1
               scene%tlev(i)           = levels_t(i)
          enddo
          scene%tsurf                  = surface_t
          scene%planck_type            = 1
          scene%wvn                    = 1. / (lambda * 1.e-6 * 1.e2);
          scene%wvnlo                  = 0.
          scene%wvnhi                  = 0.
     endif

     scene%numlay                      = n_layers
     do i = 1, n_layers
          scene%tau(i)                 = ltau(i)
          scene%omega(i)               = omega(i)

          do j = 0, n_coef_layer(i) - 1
               scene%pfmom(j, i)       = chi(j+1, 1, i)
          enddo

          do j = j, n_coef - 1
               scene%pfmom(j, i)       = 0.
          enddo
     enddo

     if (psa .eq. 0) then
          scene%sphere%use_pseudo_spherical = .false.
     else
          scene%sphere%use_pseudo_spherical = .true.
          scene%sphere%planet_radius        = planet_r
          do i = 1, n_layers + 1
               scene%sphere%zlev(i) = levels_z(i)
          enddo
     endif

     scene%sphere%use_refraction       = .false.
     scene%sphere%refrac_ind_par       = 0.
     scene%sphere%refrac_lay_grid(:)   = 0
     scene%sphere%plev(:)              = 0.

     scene%sphere%los_cor              = .false.

     scene%surf%use_reflected_direct   = .true.
     scene%surf%use_reflected_diffuse  = .true.

     if (thermal .eq. 0) then
          scene%surf%use_surface_emission   = .false.
     else
          scene%surf%use_surface_emission   = .true.
     endif

     scene%surf%n_brdf_kernels         = n_kernels

     do i = 1, n_kernels
          i_kernel = kernel_index(kernels(i) + 1)
          if (i_kernel .eq. 0) then
               print *, 'ERROR: unsupported kernel for radiant'
               info = 1
               return
          endif

          scene%surf%brdf_kernel(i)            = i_kernel - 1
          scene%surf%kernel_amp_par(i)         = ampfac(i)
          scene%surf%n_kernel_dist_par(i)      = kernel_n_params(i_kernel)

          do j = 1, kernel_n_params(i_kernel)
               scene%surf%kernel_dist_par(j,i) = params(j, i)
          enddo
     enddo
     scene%surf%n_brdf_quadratures             = n_kernel_quad * 2


     !**************************************************************************
     !
     !**************************************************************************
     if (n_derivs .eq. 0) then
          jac%numpar = 1
          jac%get_atmos_jacobian(:, :)      = .false.
          jac%l_tau(:, :)                   = 0.
          jac%l_omega(:, :)                 = 0.
          jac%l_pfmom(:, :, :)              = 0.
          jac%get_surf_amp_jacobian(:)      = .false.
          jac%get_surf_dist_jacobian(:, :)  = .false.
     else
          jac%numpar = n_derivs

          jac%get_surf_amp_jacobian(:)      = .false.
          jac%get_surf_dist_jacobian(:, :)  = .false.

          deriv_type(:) = 0

          i_brdf = 0;

          do j = 1, n_derivs
               do i = 1, n_layers
                    if (derivs(j, i) .ne. 0) then
                         jac%get_atmos_jacobian(j, i) = .true.

                         jac%l_tau(j, i)              = ltau_l(j, i)
                         jac%l_omega(j, i)            = omega_l(j, i)
                         do k = 0, n_coef_layer(i) - 1
                              jac%l_pfmom(k, j, i)    = chi_l(k+1, 1, j, i)
                         enddo

                         do k = k, n_coef - 1
                              jac%l_pfmom(k, j, i)    = chi_l(k+1, 1, j, i)
                         enddo

                         deriv_type(j)   = 1
                         layers_index(j) = i
                         derivs_index(j) = j
                    else
                         jac%get_atmos_jacobian(j, i) = .false.
                    endif
               enddo

               if (derivs(j, i) .ne. 0) then
                    do k = 1, n_kernels
                         if (ampfac_l(j, k) .ne. 0.) then
                              jac%get_surf_amp_jacobian(k) = .true.
                         endif
                    enddo

                    i_brdf = i_brdf + 1

                    deriv_type(j) = 2
                    brdf_index(j) = i_brdf
               endif
          enddo
     endif


     !**************************************************************************
     !
     !**************************************************************************
     rt_con%user_azimuth_angle = phis(1)

     call radiant(rt_con, scene, jac, rt_out)

     if (rt_out%radiant_status .ne. 0) then
          write (0, *) 'ERROR: radiant(), rt_out%radiant_status = ', &
                       rt_out%radiant_status
          info = 1
          return
     endif

     if (radiance .ne. 0) then
          if (utau_output .eq. 0) then
               if (.not. rt_con%apply_user_zenith_angle) then
                    do i = 1, n_quad
                         I_p(1, i, 1) = rt_out%radiance(3, n_quad - i + 1)
                         I_m(1, i, 1) = rt_out%radiance(0, n_quad - i + 1)
                         I_p(1, i, 2) = rt_out%radiance(2, n_quad - i + 1)
                         I_m(1, i, 2) = rt_out%radiance(1, n_quad - i + 1)
                    enddo
               else
                    I_p(1, 1, 1) = rt_out%radiance(3, n_quad + 1)
                    I_m(1, 1, 1) = rt_out%radiance(0, n_quad + 1)
                    I_p(1, 1, 2) = rt_out%radiance(2, n_quad + 1)
                    I_m(1, 1, 2) = rt_out%radiance(1, n_quad + 1)
               endif

               do i = 1, n_derivs
                    if (deriv_type(i) .eq. 1) then
                         i_layer = layers_index(i);
                         i_deriv = derivs_index(i);

                         if (.not. rt_con%apply_user_zenith_angle) then
                              do k = 1, n_quad
                                   K_p(1, k, i, 1) = rt_out%l_radiance(3, n_quad - k + 1, i_deriv, i_layer)
                                   K_m(1, k, i, 1) = 0.
                                   K_p(1, k, i, 2) = rt_out%l_radiance(2, n_quad - k + 1, i_deriv, i_layer)
                                   K_m(1, k, i, 2) = rt_out%l_radiance(1, n_quad - k + 1, i_deriv, i_layer)
                              enddo
                         else
                              K_p(1, 1, i, 1) = rt_out%l_radiance(3, n_quad + 1, i_deriv, i_layer)
                              K_m(1, 1, i, 1) = 0.
                              K_p(1, 1, i, 2) = rt_out%l_radiance(2, n_quad + 1, i_deriv, i_layer)
                              K_m(1, 1, i, 2) = rt_out%l_radiance(1, n_quad + 1, i_deriv, i_layer)
                         endif
                    else if (deriv_type(i) .eq. 2) then
                         i_brdf = brdf_index(i)

                         if (.not. rt_con%apply_user_zenith_angle) then
                              do k = 1, n_quad
                                   K_p(1, k, i, 1) = rt_out%l_radiance_surf(3, n_quad - k + 1, 4, i_brdf)
                                   K_m(1, k, i, 1) = 0.
                                   K_p(1, k, i, 2) = rt_out%l_radiance_surf(2, n_quad - k + 1, 4, i_brdf)
                                   K_m(1, k, i, 2) = rt_out%l_radiance_surf(1, n_quad - k + 1, 4, i_brdf)
                              enddo
                         else
                              K_p(1, 1, i, 1) = rt_out%l_radiance_surf(3, n_quad + 1, 4, i_brdf)
                              K_m(1, 1, i, 1) = 0.
                              K_p(1, 1, i, 2) = rt_out%l_radiance_surf(2, n_quad + 1, 4, i_brdf)
                              K_m(1, 1, i, 2) = rt_out%l_radiance_surf(1, n_quad + 1, 4, i_brdf)
                         endif
                    endif
               enddo
          else
               if (.not. rt_con%apply_user_zenith_angle) then
                    do i = 1, n_quad
                         do j = 1, n_ulevels
                              I_p(1, i, j) = rt_out%user_radiance(2, n_quad - i + 1, j)
                              I_m(1, i, j) = rt_out%user_radiance(1, n_quad - i + 1, j)
                         enddo
                    enddo
               else
                    do j = 1, n_ulevels
                         I_p(1, 1, j) = rt_out%user_radiance(2, n_quad + 1, j)
                         I_m(1, 1, j) = rt_out%user_radiance(1, n_quad + 1, j)
                    enddo
               endif

               do i = 1, n_derivs
                    if (deriv_type(i) .eq. 1) then
                         i_layer = layers_index(i);
                         i_deriv = derivs_index(i);

                         if (.not. rt_con%apply_user_zenith_angle) then
                              do k = 1, n_quad
                                   do j = 1, n_ulevels
                                        K_p(1, k, i, j) = rt_out%l_user_radiance(2, n_quad - k + 1, i_deriv, i_layer, j)
                                        K_m(1, k, i, j) = rt_out%l_user_radiance(1, n_quad - k + 1, i_deriv, i_layer, j)
                                   enddo
                              enddo
                         else
                              do j = 1, n_ulevels
                                   K_p(1, 1, i, j) = rt_out%l_user_radiance(2, n_quad + 1, i_deriv, i_layer, j)
                                   K_m(1, 1, i, j) = rt_out%l_user_radiance(1, n_quad + 1, i_deriv, i_layer, j)
                              enddo
                         endif
                    else if (deriv_type(i) .eq. 2) then
                         i_brdf = brdf_index(i)

                         if (.not. rt_con%apply_user_zenith_angle) then
                              do k = 1, n_quad
                                   do j = 1, n_ulevels
                                        K_p(1, k, i, j) = rt_out%l_user_radiance_surf(2, n_quad - k + 1, 4, i_brdf, j)
                                        K_m(1, k, i, j) = rt_out%l_user_radiance_surf(1, n_quad - k + 1, 4, i_brdf, j)
                                   enddo
                              enddo
                         else
                              do j = 1, n_ulevels
                                   K_p(1, 1, i, j) = rt_out%l_user_radiance_surf(2, n_quad + 1, 4, i_brdf, j)
                                   K_m(1, 1, i, j) = rt_out%l_user_radiance_surf(1, n_quad + 1, 4, i_brdf, j)
                              enddo
                         endif
                    endif
               enddo
          endif
     endif

     if (flux .ne. 0) then
          if (utau_output .eq. 0) then
               flux_p(1) = rt_out%flux(3)
               flux_m(1) = rt_out%flux(0)
               flux_p(2) = rt_out%flux(2)
               flux_m(2) = rt_out%flux(1)
          else
               do i = 1, n_ulevels
                    flux_p(i) = rt_out%user_flux(2, i)
                    flux_m(i) = rt_out%user_flux(1, i)
               enddo
          endif
     endif


     !**************************************************************************
     !
     !**************************************************************************
     r = radiant_datatype_ctrl2(2, rt_con, scene, jac, rt_out, &
                                n_quad * 2, n_coef, n_layers, 1, 2)
     if (r .ne. 0) then
          write (0, *) 'ERROR: radiant_datatype_ctrl2(2, ...), status = ', r
          info = 1
          return
     endif


     return

end subroutine call_radiant_f
