      subroutine pentdag(a,b,c,d,e,f,u,n)
      implicit none
      save

c..solves for a vector u of length n in the pentadiagonal linear system
c.. a_i u_(i-2) + b_i u_(i-1) + c_i u_i + d_i u_(i+1) + e_i u_(i+2) = f_i
c..input are the a, b, c, d, e, and f and they are not modified

c..in its clearest incarnation, this algorithm uses three storage arrays
c..called p, q and r. here, the solution vector u is used for r, cutting
c..the extra storage down to two arrays.

c..declare the pass
      integer          n
      double precision a(n),b(n),c(n),d(n),e(n),f(n),u(n)

c..local variables
      integer          nmax,i
      parameter        (nmax=500)
      double precision p(nmax),q(nmax),bet,den

c..initialize elimination and backsubstitution arrays
      if (c(1) .eq. 0.0)  stop 'eliminate u2 trivially'
      bet  = 1.0d0/c(1)
      p(1) = -d(1) * bet
      q(1) = -e(1) * bet
      u(1) = f(1)  * bet

      bet = c(2) + b(2)*p(1)
      if (bet .eq. 0.0) stop 'singular 1 in pentdag'
      bet = -1.0d0/bet
      p(2) = (d(2) + b(2)*q(1)) * bet
      q(2) = e(2) * bet
      u(2) = (b(2)*u(1) - f(2)) * bet

c..reduce to upper triangular
      do i=3,n
       bet = b(i) + a(i) * p(i-2)
       den = c(i) + a(i)*q(i-2) + bet*p(i-1)
       if (den .eq. 0.0) stop 'singular 2 in pentdag'
       den = -1.0d0/den
       p(i) = (d(i) + bet*q(i-1)) * den
       q(i) = e(i) * den
       u(i) = (a(i)*u(i-2) + bet*u(i-1) - f(i)) * den
      enddo

c..backsubstitution
      u(n-1) = u(n-1) + p(n-1) * u(n)
      do i=n-2,1,-1
       u(i) = u(i) + p(i) * u(i+1) + q(i) * u(i+2)
      enddo
      return
      end



      subroutine pentdag1(a,b,c,d,e,p,q,n)
      implicit none

c..solves for a vector u of length n in the pentadiagonal linear system
c.. a_i u_(i-2) + b_i u_(i-1) + c_i u_i + d_i u_(i+1) + e_i u_(i+2) = f_i
c..input are the a, b, c, d, e, and f and they are not modified

c..in its clearest incarnation, this algorithm uses three storage arrays
c..called p, q and r. here, the solution vector u is used for r, cutting
c..the extra storage down to two arrays.

c..declare the pass
      integer          n
      double precision a(n),b(n),c(n),d(n),e(n),p(n),q(n)

c..local variables
      integer          i
      double precision bet,den

c..initialize elimination and backsubstitution arrays
      if (c(1) .eq. 0.0)  stop 'eliminate u2 trivially'
      bet  = 1.0d0/c(1)
      p(1) = -d(1) * bet
      q(1) = -e(1) * bet

      bet = c(2) + b(2)*p(1)
      if (bet .eq. 0.0) stop 'singular 1 in pentdag'
      bet = -1.0d0/bet
      p(2) = (d(2) + b(2)*q(1)) * bet
      q(2) = e(2) * bet

c..reduce to upper triangular
      do i=3,n
       bet = b(i) + a(i) * p(i-2)
       den = c(i) + a(i)*q(i-2) + bet*p(i-1)
       if (den .eq. 0.0) stop 'singular 2 in pentdag'
       den = -1.0d0/den
       p(i) = (d(i) + bet*q(i-1)) * den
       q(i) = e(i) * den
      enddo

      return
      end



      subroutine pentdag2(a,b,c,f,u,p,q,n)
      implicit none

c..solves for a vector u of length n in the pentadiagonal linear system
c.. a_i u_(i-2) + b_i u_(i-1) + c_i u_i + d_i u_(i+1) + e_i u_(i+2) = f_i
c..input are the a, b, c, d, e, and f and they are not modified

c..in its clearest incarnation, this algorithm uses three storage arrays
c..called p, q and r. here, the solution vector u is used for r, cutting
c..the extra storage down to two arrays.

c..declare the pass
      integer          n
      double precision a(n),b(n),c(n),f(n),u(n),p(n),q(n)

c..local variables
      integer          i
      double precision bet,den

c..initialize elimination and backsubstitution arrays
      bet  = 1.0d0/c(1)
      u(1) = f(1)  * bet

      bet = c(2) + b(2)*p(1)
      bet = -1.0d0/bet
      u(2) = (b(2)*u(1) - f(2)) * bet

c..reduce to upper triangular
      do i=3,n
       bet = b(i) + a(i) * p(i-2)
       den = c(i) + a(i)*q(i-2) + bet*p(i-1)
       den = -1.0d0/den
       u(i) = (a(i)*u(i-2) + bet*u(i-1) - f(i)) * den
      enddo

c..backsubstitution
      u(n-1) = u(n-1) + p(n-1) * u(n)
      do i=n-2,1,-1
       u(i) = u(i) + p(i) * u(i+1) + q(i) * u(i+2)
      enddo
      return
      end
