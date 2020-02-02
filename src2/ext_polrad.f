c*******************************************************************************
c
c    Copyright (C) 2007-2020 Greg McGarragh <greg.mcgarragh@colostate.edu>
c
c    This source code is licensed under the GNU General Public License (GPL),
c    Version 3.  See the file COPYING for more details.
c
c*******************************************************************************

      subroutine call_polrad_f
     &                    (n_four, n_quad, n_stokes, n_layers, qx, F_0, theta_0,
     &                     phi_0, n_ulevels, umus, n_umus, phis, n_phis, albedo,
     &                     n_coef_layer, chi, omega, ltau, I_p, I_m, delta_m,
     &                     epsilon, info, n_quad2)

      implicit none


      integer n_four
      integer n_quad
      integer n_stokes
      integer n_layers

      real*8 qx(n_quad)

      real*8 F_0

      real*8 theta_0

      real*8 phi_0

      integer n_ulevels

      real*8 umus(n_umus)
      integer n_umus

      real*8 phis(n_phis)
      integer n_phis

      real*8 albedo

      integer n_coef_layer(n_layers)
      real*8 chi(4, 4, 0:5100, 40)

      real*8 omega(n_layers)

      real*8 ltau(n_layers)

      real*8 I_p(n_stokes, n_phis, n_quad2, n_ulevels)
      real*8 I_m(n_stokes, n_phis, n_quad2, n_ulevels)

      integer delta_m

      real*8 epsilon

      integer info

      integer n_quad2


c     ******************************************************************
c     *
c     ******************************************************************
c     integer i

      integer NDgeom
      integer NDmu
      integer NDlay
      integer NDcoef
      parameter(NDgeom=25, NDmu=100, NDlay=40, NDcoef=5100)

      real*8 a

      real*8 pi

      real*8 T(4, 4, NDgeom)
      real*8 Tst(4, 4, NDgeom)
      real*8 R(4, 4, NDgeom)
      real*8 Rst(4, 4, NDgeom)
      real*8 T1(4, 4, NDgeom)
      real*8 T1st(4, 4, NDgeom)
      real*8 R1(4, 4, NDgeom)
      real*8 R1st(4, 4, NDgeom)
      real*8 RL(4, 4, NDgeom)
      real*8 DL(4, 4, NDgeom)
      real*8 R1L(4, 4, NDgeom)
      real*8 D1L(4, 4, NDgeom)
      real*8 Ri(4, 4, NDgeom)
      real*8 Di(4, 4, NDgeom)
      real*8 R1i(4, 4, NDgeom)
      real*8 D1i(4, 4, NDgeom)
c     real*8 RiL(4, 4, NDgeom)
c     real*8 DiL(4, 4, NDgeom)
c     real*8 R1iL(4, 4, NDgeom)
c     real*8 D1iL(4, 4, NDgeom)

      integer Nmu0
      integer Nmu
      integer Nzero
      integer Nfirst
      integer Nall
      integer Nabove
      integer Nbelow
      parameter(Nmu0=1, Nmu=2, Nzero=1,
     +          Nfirst=2, Nall=3, Nabove=1, Nbelow=2 )

      real*8 Rflux(4, 4, NDgeom, Nmu, Nbelow, Nall)
      real*8 Tflux(4, 4, NDgeom, Nmu, Nbelow, Nall)
      real*8 URU(  4, 4,              Nbelow, Nall)
      real*8 UTU(  4, 4,              Nbelow, Nall)
      real*8 RLflux(4, 4, NDgeom, Nmu, Nbelow, Nall)
      real*8 DLflux(4, 4, NDgeom, Nmu, Nbelow, Nall)
      real*8 URUL(  4, 4,              Nbelow, Nall)
      real*8 Riflux(4, 4, NDgeom, Nmu, Nbelow, Nall)
      real*8 Diflux(4, 4, NDgeom, Nmu, Nbelow, Nall)
      real*8 Tiflux(4, 4, NDgeom, Nmu, Nbelow, Nall)
      real*8 URUi(  4, 4,              Nbelow, Nall)
      real*8 UDUi(  4, 4,              Nbelow, Nall)
      real*8 UTUi(  4, 4,              Nbelow, Nall)
c     real*8 RiLflux(4, 4, NDgeom, Nmu, Nbelow, Nall)
c     real*8 DiLflux(4, 4, NDgeom, Nmu, Nbelow, Nall)
c     real*8 URUiL(  4, 4,              Nbelow, Nall)

      real*8 Svin(4)
      real*8 SvR(4, NDgeom)
      real*8 SvT(4, NDgeom)


c     ******************************************************************
c     *
c     ******************************************************************
      a = qx(1)
      a = F_0
      a = phi_0
      a = albedo

      a = n_umus

      a = n_phis

      a = I_p(1, 1, 1, 1)
      a = I_m(1, 1, 1, 1)

      a = delta_m


      info = 0

      pi = 4.d0 * atan(1.0d0)


c     ******************************************************************
c     *
c     ******************************************************************
      a = 180.d0 / pi

      call adding(omega, ltau, chi, n_coef_layer, n_layers,
     +            acos(umus(1))*a, theta_0, phis(1), 1,
     +            n_quad, n_four, 1, n_stokes, epsilon, 0, 0.,
     +            R, T, R1, T1, Rst, Tst, R1st, T1st,
     +            Ri, R1i, Di, D1i,
     +            Rflux, Tflux, URU, UTU,
     +            Riflux, Tiflux, Diflux, URUi, UTUi, UDUi)


c     ******************************************************************
c     *
c     ******************************************************************
      if (albedo .ne. 0.) then
          call addlam(albedo, R, R1, T, T1, Rflux, Tflux, URU, UTU,
     +           1, n_stokes, RL, R1L, DL, D1L, RLflux, DLflux, URUL)
      endif


c     ******************************************************************
c     *
c     ******************************************************************
      Svin(1) = F_0
      Svin(2) = 0.
      Svin(3) = 0.
      Svin(4) = 0.

      call Stokout(theta_0, n_stokes, 1, Svin, R, SvR)
      call Stokout(theta_0, n_stokes, 1, Svin, T, SvT)

      print *, SvR(1,1), SvT(1,1)
      stop

      return

      end
