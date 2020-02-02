!*******************************************************************************
!
!    Copyright (C) 2007-2020 Greg McGarragh <greg.mcgarragh@colostate.edu>
!
!    This source code is licensed under the GNU General Public License (GPL),
!    Version 3.  See the file COPYING for more details.
!
!*******************************************************************************

!*******************************************************************************
! n_stokes = 2 or 3, lambertian only
!*******************************************************************************
subroutine call_lrad_f(n_four, n_elem, n_coef, n_quad, n_stokes, n_derivs, &
                       n_layers, qx, F_0, mu_0, phi_0, n_ulevels, umus, n_umus, &
                       phis, n_phis, levels_z, n_kernels, n_kernel_quad, kernels, &
                       ampfac, ampfac_l, params, params_l, n_coef_layer, chi0, &
                       chi0_l, chi, chi_l, omega0, omega0_l, omega, omega_l, &
                       ltau, ltau_l, I_p, I_m, K_p, K_m, n_t_tms, psa, derivs, &
                       epsilon, info, n_umus2)

     use lrad_driver_module


     implicit none


     integer, parameter     :: max_params = 3

     integer, intent(in)    :: n_four
     integer, intent(in)    :: n_elem
     integer, intent(in)    :: n_coef
     integer, intent(in)    :: n_quad
     integer, intent(in)    :: n_stokes
     integer, intent(in)    :: n_derivs
     integer, intent(in)    :: n_layers

     real(8), intent(in)    :: qx(n_quad)

     real(8), intent(in)    :: F_0

     real(8), intent(in)    :: mu_0

     real(8), intent(in)    :: phi_0

     integer, intent(in)    :: n_ulevels

     real(8), intent(in)    :: umus(n_umus)
     integer, intent(in)    :: n_umus

     real(8), intent(in)    :: phis(n_phis)
     integer, intent(in)    :: n_phis

     real(8), intent(in)    :: levels_z(0:n_layers)

     integer, intent(in)    :: n_kernels
     integer, intent(in)    :: n_kernel_quad
     integer, intent(in)    :: kernels(n_kernels)
     real(8), intent(in)    :: ampfac(n_kernels)
     real(8), intent(in)    :: ampfac_l(n_derivs, n_kernels)
     real(8), intent(in)    :: params(max_params, n_kernels)
     real(8), intent(in)    :: params_l(max_params, n_derivs, n_kernels)

     integer, intent(in)    :: n_coef_layer(n_layers)
     real(8), intent(in)    :: chi0  (0:n_coef-1, n_elem, n_layers)
     real(8), intent(in)    :: chi0_l(0:n_coef-1, n_elem, n_derivs, n_layers)
     real(8), intent(in)    :: chi   (0:n_coef-1, n_elem, n_layers)
     real(8), intent(in)    :: chi_l (0:n_coef-1, n_elem, n_derivs, n_layers)

     real(8), intent(in)    :: omega0(n_layers)
     real(8), intent(in)    :: omega0_l(n_derivs, n_layers)
     real(8), intent(in)    :: omega (n_layers)
     real(8), intent(in)    :: omega_l (n_derivs, n_layers)

     real(8), intent(in)    :: ltau(n_layers)
     real(8), intent(in)    :: ltau_l(n_derivs, n_layers)

     real(8), intent(out)   :: I_p(n_stokes, n_phis, n_umus2, n_ulevels)
     real(8), intent(out)   :: I_m(n_stokes, n_phis, n_umus2, n_ulevels)

     real(8), intent(out)   :: K_p(n_stokes, n_phis, n_umus2, n_derivs, n_ulevels)
     real(8), intent(out)   :: K_m(n_stokes, n_phis, n_umus2, n_derivs, n_ulevels)

     integer, intent(in)    :: n_t_tms
     integer, intent(in)    :: psa

     integer(1), intent(in) :: derivs(n_derivs + 8, n_layers)

     real(8), intent(in)    :: epsilon

     integer, intent(out)   :: info

     integer, intent(in)    :: n_umus2


     !**************************************************************************
     !
     !**************************************************************************
     logical          :: regular_ps
     logical          :: linearize
     logical          :: s_linearize

     integer          :: i
     integer          :: ii
     integer          :: j
     integer          :: k
     integer          :: kk
     integer          :: l

     integer          :: lmax

     integer          :: n_coef2

     integer          :: nphibrdf
     integer          :: surftype
     integer          :: nspars

     byte             :: deriv_type(n_derivs)
     integer          :: layers_index(n_derivs)

     double precision :: a
     double precision :: b
     double precision :: c

     double precision :: pi

     double precision :: theta0
     double precision :: theta

     double precision :: fscale(n_layers)

     double precision :: opd (n_layers)
     double precision :: ssa0(n_layers)
     double precision :: ssa (n_layers)
     double precision :: ssa2(n_layers)
     double precision :: coefs0(0:2*n_quad-1, n_layers, 6)
     double precision :: coefs (0:2*n_quad-1, n_layers, 6)
     double precision :: Zmat(n_layers, n_stokes)

     double precision :: L_opd (n_layers, n_derivs)
     double precision :: L_ssa0(n_layers, n_derivs)
     double precision :: L_ssa (n_layers, n_derivs)
     double precision :: L_ssa2(n_layers, n_derivs)
     double precision :: L_coefs0(0:2*n_quad-1, n_layers, 6, n_derivs)
     double precision :: L_coefs (0:2*n_quad-1, n_layers, 6, n_derivs)
     double precision :: L_Zmat(n_layers, n_stokes, n_derivs)
     double precision :: L_fscale(n_layers, n_derivs)

     double precision :: spars(3)

     double precision :: R1(n_stokes)
     double precision :: Ls_R1(n_stokes, n_derivs)
     double precision :: L_R1(n_stokes, n_layers, n_derivs)

     double precision :: R2(n_stokes)
     double precision :: Ls_R2(n_stokes, n_derivs)
     double precision :: L_R2(n_stokes, n_layers, n_derivs)

     double precision :: Icorr
     double precision :: Ls_Icorr(n_derivs)
     double precision :: L_Icorr(n_layers, n_derivs)


     !**************************************************************************
     !
     !**************************************************************************
     if (n_umus > 1) then
          write (0, *) 'ERROR: lrad can only handle one user zenith angle'
          info = 1
          return
     endif

     if (n_phis > 1) then
          write (0, *) 'ERROR: lrad can only handle one user azimuth angle'
          info = 1
          return
     endif


     !**************************************************************************
     !
     !**************************************************************************
     a = n_four
     a = qx(1)
     a = F_0
     a = phi_0
     if (n_kernels > 1) then
          a = kernels(1)
          a = ampfac_l(1,1)
          a = params(1,1)
          a = params_l(1,1,1)
     endif
     if (n_derivs > 1) then
          a = derivs(1,1)
     endif

     info = 0

     pi = 4.d0 * atan(1.0d0)


     !**************************************************************************
     !
     !**************************************************************************
     theta0 = acos(mu_0) * 180.d0/pi

     if (psa == 0) then
          regular_ps = .false.
     else
          regular_ps = .true.
     endif

     if (n_derivs == 0) then
          linearize   = .false.
          s_linearize = .false.
     else
          linearize   = .true.
          s_linearize = .false.
     endif

     lmax = 0

     coefs   = 0.
     L_coefs = 0.

     a = 2. * 2. * n_quad + 1.

     if (n_elem == 1) then
          coefs0   = 0.
          coefs    = 0.
     endif

     do i = 1, n_layers
          ii = n_layers - i + 1

          opd(ii)  = ltau(i)
          ssa0(ii) = omega0(i)
          ssa(ii)  = omega(i)

          b = 1.
          if (n_t_tms /= 0) then
               if (n_coef_layer(i) > n_quad * 2) then
                    b = 1. - chi0(n_quad * 2, 1, i) / a * omega0(i)
               endif
          endif

          ssa2(ii) = omega0(i) / b

          lmax = max(lmax, n_coef_layer(i) - 1)

          n_coef2 = min(n_quad * 2 - 1, n_coef_layer(i) - 1)

          do j = 0, n_coef2
               coefs0(j, ii, 1) =  chi0(j, 1, i)
               coefs (j, ii, 1) =  chi (j, 1, i)

               if (n_elem > 1) then
                    coefs0(j, ii, 2) =  chi0(j, 2, i)
                    coefs0(j, ii, 3) =  chi0(j, 3, i)
                    coefs0(j, ii, 4) =  chi0(j, 4, i)
                    coefs0(j, ii, 5) =  chi0(j, 5, i)
                    coefs0(j, ii, 6) = -chi0(j, 6, i)

                    coefs (j, ii, 2) =  chi (j, 2, i)
                    coefs (j, ii, 3) =  chi (j, 3, i)
                    coefs (j, ii, 4) =  chi (j, 4, i)
                    coefs (j, ii, 5) =  chi (j, 5, i)
                    coefs (j, ii, 6) = -chi (j, 6, i)
               endif
          enddo

!         if (n_t_tms == 0) then
               fscale(ii) = 0.
!         else
!              if (n_coef_layer(i) <= n_quad * 2) then
!                   fscale(ii) = 0.
!              else
!                   fscale(ii) = chi(n_quad * 2, 1, i) / a
!              endif
!         endif
     enddo

     L_opd    = 0.
     L_ssa0   = 0.
     L_ssa    = 0.
     L_ssa2   = 0.
     L_coefs0 = 0.
     L_coefs  = 0.

     if (n_derivs == 0) then
          linearize   = .false.
          s_linearize = .false.
     else
          linearize   = .true.
          s_linearize = .false.

          do j = 1, n_derivs
               do i = 1, n_layers
                    ii = n_layers - i + 1

                    if (derivs(j, i) /= 0) then
                         L_opd(ii, j)  = ltau_l(j, i)
                         L_ssa0(ii, j) = omega0_l(j, i)
                         L_ssa(ii, j)  = omega_l(j, i)

                         if (n_t_tms == 0) then
                              L_ssa2(ii, j) = omega0_l(i, j)
                         else
                              if (n_coef_layer(i) > n_quad * 2) then
                                   c = (-chi0_l(n_quad * 2, 1, j, i) * omega0(i) - chi0(n_quad * 2, 1, i) * omega0_l(j, i)) / a
                              endif

                              L_ssa2(ii, j) = (omega0_l(i, j) - ssa2(ii) * c) / b
                         endif

                         n_coef2 = min(n_quad * 2 - 1, n_coef_layer(i) - 1)

                         do k = 0, n_coef2
                              L_coefs0(k, ii, 1, j) =  chi0_l(k, 1, j, i)
                              L_coefs (k, ii, 1, j) =  chi_l (k, 1, j, i)

                              if (n_elem > 1) then
                                   L_coefs0(k, ii, 2, j) =  chi0_l(k, 2, j, i)
                                   L_coefs0(k, ii, 3, j) =  chi0_l(k, 3, j, i)
                                   L_coefs0(k, ii, 4, j) =  chi0_l(k, 4, j, i)
                                   L_coefs0(k, ii, 5, j) =  chi0_l(k, 5, j, i)
                                   L_coefs0(k, ii, 6, j) = -chi0_l(k, 6, j, i)

                                   L_coefs (k, ii, 2, j) =  chi_l (k, 2, j, i)
                                   L_coefs (k, ii, 3, j) =  chi_l (k, 3, j, i)
                                   L_coefs (k, ii, 4, j) =  chi_l (k, 4, j, i)
                                   L_coefs (k, ii, 5, j) =  chi_l (k, 5, j, i)
                                   L_coefs (k, ii, 6, j) = -chi_l (k, 6, j, i)
                              endif
                         enddo

!                        if (n_t_tms == 0) then
                                   L_fscale(ii, j) = 0.
!                        else
!                             if (n_coef_layer(i) <= n_quad * 2) then
!                                  L_fscale(ii, j) = 0.
!                             else
!                                  L_fscale(ii, j) = chi_l(n_quad * 2, 1, j, i) / a
!                             endif
!                        endif

                         deriv_type  (j) = 1
                         layers_index(j) = ii

                         goto 666
                    endif
               enddo

               if (derivs(j, i) /= 0) then
                    if (s_linearize) then
                         write (0, *) 'ERROR: in lrad the albedo can only be linearized once'
                         info = 1
                         return
                    endif

                    if (ampfac_l(j,1) /= 1.) then
                         write (0, *) 'ERROR: in lrad the albedo input derivative umust be equal to unity'
                         info = 1
                         return
                    endif

                    s_linearize = .true.

                    deriv_type  (j) = 2
               endif

666       do i = i + 1, n_layers + 1
                    if (derivs(j, i) /= 0) then
                         write (0, *) 'ERROR: lrad only supports one linearized layer per deriv, i_layer,i_deriv =', i, j
                         info = 1
                         return
                    endif
               enddo
          enddo
     endif


     nphibrdf = 2 * n_kernel_quad

     surftype = 1

     nspars   = 1

     spars(1) = ampfac(1)


     !**************************************************************************
     !
     !**************************************************************************
     I_p = 0.
     I_m = 0.

     if (n_derivs /= 0) then
          K_p = 0.
          K_m = 0.
     endif

     theta = acos(umus(1)) * 180.d0/pi

     call calc_gsf(theta, theta0, phis(1), lmax)

     do k = 1, n_layers
          kk = n_layers - k + 1
          call calculate_zmat(chi0(0:n_coef_layer(k)-1,:,k), Zmat(kk,:), n_stokes)

          do l = 1, n_derivs
               call calculate_zmat(chi0_l(0:n_coef_layer(k)-1,:,l,k), L_Zmat(kk,:,l), n_stokes)
          enddo
     enddo

     call comp_surf_first(n_stokes, surftype, s_linearize, theta, theta0, phis(1), nspars, spars, R1, Ls_R1)

     call Lrad_first_init(theta, theta0, phis(1), levels_z, regular_ps, .false.)

     call Lrad_first(opd, ssa2, zmat, levels_z, fscale, L_opd, L_ssa2, L_zmat, L_fscale, theta, theta0, regular_ps, .false., linearize, s_linearize, R1, Ls_R1, L_R1)

     call Lrad_second_init(theta0, levels_z, regular_ps, .false.)

     call comp_surf_second(n_stokes, n_quad, nphibrdf, nspars, surftype, s_linearize, theta, theta0, spars)

     call Lrad_second(opd, ssa, coefs, L_opd, L_ssa, L_coefs, phis(1), epsilon, n_quad, regular_ps, .false., linearize, s_linearize, R2, Ls_R2, L_R2, Icorr, Ls_Icorr, L_Icorr)

     i = 1

!    I_p(:, 1, i, 1) = R1 / pi
!    do k = 1, n_derivs
!         if (deriv_type  (k) == 1) then
!              K_p(:, 1, i, k, 1) = L_R1(:, layers_index(k), k) / pi
!         else
!              K_p(:, 1, i, k, 1) = Ls_R1(:, 1) / pi
!         endif
!    enddo

     I_p(:, 1, i, 1) = R2 / pi
     do k = 1, n_derivs
          if (deriv_type  (k) == 1) then
               K_p(:, 1, i, k, 1) = L_R2(:, layers_index(k), k) / pi
          else
               K_p(:, 1, i, k, 1) = Ls_R2(:, 1) / pi
          endif
     enddo

!    I_p(:, 1, i, 1) = (R1 + R2) / pi
!    do k = 1, n_derivs
!         if (deriv_type  (k) == 1) then
!              K_p(:, 1, i, k, 1) = (L_R1(:, layers_index(k), k) + L_R2(:, layers_index(k), k)) / pi
!         else
!              K_p(:, 1, i, k, 1) = (Ls_R1(:, 1) + Ls_R2(:, 1)) / pi
!         endif
!    enddo

!    I_p(:, 1, i, 1) = (R2(1) - Icorr) / pi
!    do k = 1, n_derivs
!         if (deriv_type  (k) == 1) then
!              K_p(:, 1, i, k, 1) = (L_R2(1, layers_index(k), k) - L_Icorr(layers_index(k), k)) / pi
!         else
!              K_p(:, 1, i, k, 1) = (Ls_R2(1, 1) - Ls_Icorr(layers_index(k))) / pi
!         endif
!    enddo


     return

end
