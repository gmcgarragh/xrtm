!*******************************************************************************
!
!    Copyright (C) 2007-2020 Greg McGarragh <greg.mcgarragh@colostate.edu>
!
!    This source code is licensed under the GNU General Public License (GPL),
!    Version 3.  See the file COPYING for more details.
!
!*******************************************************************************

subroutine call_2stream_f(n_four, n_elem, n_coef, n_quad, n_derivs, n_layers, &
                          qx, lambda, F_0, theta_0, phi_0, umus, n_umus, phis, &
                          n_phis, top_t, planet_r, levels_z, levels_t, surface_t, &
                          n_kernels, n_kernel_quad, kernels, ampfac, ampfac_l, &
                          params, params_l, g, g_l, omega, omega_l, ltau, ltau_l, &
                          I_p,I_m, K_p, K_m, flux_p, flux_m, quadrature, delta_m, &
                          n_t_tms, psa, radiance, flux, thermal, derivs, epsilon, &
                          info, n_quad2)

     use twostream_master_m


     implicit none


     integer, parameter     :: max_params = 3

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

     real(8), intent(in)    :: g(n_layers)
     real(8), intent(in)    :: g_l(n_derivs, n_layers)

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

     integer(1), intent(in) :: derivs(n_derivs + 8, n_layers)

     real(8), intent(in)    :: epsilon

     integer, intent(out)   :: info

     integer, intent(in)    :: n_quad2


     !**************************************************************************
     !
     !**************************************************************************
     integer       :: i

     real(8)       :: a

     real(8)       :: pi

     logical       :: do_upwelling
     logical       :: do_dnwelling
     logical       :: do_plane_parallel
     logical       :: do_solar_sources
     logical       :: do_thermal_emission
     logical       :: do_surface_emission
     logical       :: do_d2s_scaling
     logical       :: do_brdf_surface

     integer       :: nlayers
     integer       :: ntotal
     integer       :: nthreads
     integer       :: nbeams
     integer       :: n_user_streams
     integer       :: n_user_relazms
     integer       :: n_geometries

     real(8)       :: beam_szas   (1)
     real(8)       :: user_angles (max(n_quad, n_umus))
     real(8)       :: user_relazms(n_phis)

     real(8)       :: stream_value

     real(8)       :: lambertian_albedo(1)


     real(8)       :: brdf_f_0(0:1, 1)
     real(8)       :: brdf_f  (0:1 )
     real(8)       :: ubrdf_f (0:1, n_umus)

     real(8)       :: thermal_bb_input(0:n_layers)
     real(8)       :: surfbb
     real(8)       :: emissivity

     integer       :: thread

     real(8)       :: flux_factor

     real(8)       :: earth_radius

     real(8)       :: height_grid(0:n_layers)

     real(8)       :: deltau_input(n_layers, 1)
     real(8)       :: omega_input (n_layers, 1)
     real(8)       :: asymm_input (n_layers, 1)
     real(8)       :: d2s_scaling (n_layers, 1)

     real(8)       :: intensity_toa(1, 1)
     real(8)       :: intensity_boa(1, 1)

     integer       :: status_inputcheck
     integer       :: c_nmessages
     character*100 :: c_messages(100)
     character*100 :: c_actions(100)

     integer       :: status_execution
     character*100 :: e_message
     character*100 :: e_trace_1
     character*100 :: e_trace_2


     !**************************************************************************
     !
     !**************************************************************************
     a = n_four
     a = n_elem
     a = n_coef
     a = n_derivs
     a = lambda
     a = phi_0
     a = top_t
     a = levels_t(1)
     a = surface_t
     a = n_kernels
     a = n_kernel_quad
     a = kernels(1)
     a = ampfac_l(1, 1)
     a = params(1, 1)
     a = params_l(1, 1, 1)
     a = g_l(1, 1)
     a = omega_l(1, 1)
     a = ltau_l(1, 1)
     a = I_p(1, 1, 1)
     a = I_m(1, 1, 1)
     a = K_p(1, 1, 1, 1)
     a = K_m(1, 1, 1, 1)
     a = flux_p(1)
     a = flux_m(1)
     a = quadrature
     a = delta_m
     a = n_t_tms
     a = radiance
     a = flux
     a = thermal
     a = derivs(1, 1)
     a = epsilon


     info = 0

     pi = 4.d0 * atan(1.0d0)


     !**************************************************************************
     !
     !**************************************************************************
     do_upwelling         = .true.
     do_dnwelling         = .true.
     do_plane_parallel    = .true.
     do_solar_sources     = .true.
     do_thermal_emission  = .false.
     do_surface_emission  = .false.
     do_d2s_scaling       = .false.
     do_brdf_surface      = .false.

     nlayers              = n_layers
     ntotal               = 2 * n_layers
     nthreads             = 1
     nbeams               = 1

     if (n_umus .eq. 0) then
          n_user_streams = n_quad
     else
          n_user_streams = n_umus
     endif

     n_user_relazms       = 1
     n_geometries         = 1

     beam_szas(1)         = theta_0

     a = 180.d0 / pi
     if (n_umus .eq. 0) then
          do i = 1, n_quad
               user_angles(i) = acos(qx(i)) * a
          enddo
     else
          do i = 1, n_umus
               user_angles(i) = acos(umus(i)) * a
          enddo
     endif

     user_relazms(1)      = phis(1)
     stream_value         = qx(1)
     lambertian_albedo(1) = ampfac(1)
     brdf_f_0             = 0.
     brdf_f               = 0.
     ubrdf_f              = 0.
     thermal_bb_input     = 0.
     surfbb               = 0.
     emissivity           = 0.
     thread               = 1;
     flux_factor          = F_0
     earth_radius         = planet_r
     height_grid(0)       = 1.
     height_grid(1)       = 0.
     if (psa .ne. 0) then
          height_grid     = levels_z
     endif
     deltau_input(: , 1)  = ltau
     omega_input(: , 1)   = omega
     asymm_input(: , 1)   = g
     d2s_scaling          = 0.


     !**************************************************************************
     !
     !**************************************************************************
     call twostream_master(do_upwelling, do_dnwelling, do_plane_parallel, &
                           do_solar_sources, do_thermal_emission, do_surface_emission, &
                           do_d2s_scaling, do_brdf_surface, nlayers, ntotal, nthreads, &
                           nbeams, n_user_streams, n_user_relazms, n_geometries, &
                           beam_szas, user_angles, user_relazms, stream_value, &
                           lambertian_albedo, brdf_f_0, brdf_f, ubrdf_f, &
                           thermal_bb_input, surfbb, emissivity, &
                           thread, flux_factor, earth_radius, height_grid, &
                           deltau_input, omega_input, asymm_input, d2s_scaling, &
                           intensity_toa, intensity_boa, &
                           status_inputcheck, c_nmessages, c_messages, c_actions, &
                           status_execution,  e_message, e_trace_1, e_trace_2)
     if (status_inputcheck .ne. 0) then
          write (0, *) 'ERROR: twostream_master(), status_inputcheck = ', status_inputcheck
          do i = 1, c_nmessages
               print *, c_messages(i), c_actions(i)
          enddo
          info = 1
          return
     endif

     if (status_execution .ne. 0) then
          write (0, *) 'ERROR: twostream_master(), status_inputcheck = ', status_inputcheck
          info = 1
          return
     endif


     !**************************************************************************
     !
     !**************************************************************************
     print *, intensity_toa(1, 1), intensity_boa(1, 1)


     return

end subroutine call_2stream_f
