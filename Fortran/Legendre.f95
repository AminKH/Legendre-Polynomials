! Program by Amin Khiabani
! aminkhiabani@outlook.com

      RECURSIVE FUNCTION factorial(n) RESULT(res)
          integer :: n
          INTEGER(kind=8) res
          IF (n .EQ. 0) THEN
              res = 1
          ELSE
              res = n * factorial(n - 1)
          END IF
      END

      recursive function Legendref(n,x) result(Pn)
         integer :: n
         real(kind=8):: Pn , x

         if (n == 0) then
           Pn = 1.0
         else if (n == 1) then
           Pn = x
         else
           Pn = (-(n-1)*Legendref(n-2,x) + (2*n-1)*x*Legendref(n-1,x))/n
         end if
      end function Legendref


      subroutine LegendrefSeries(n,x,LFS)
            implicit none
            integer :: n, I
            real(kind=8) :: x
            real(kind=8), dimension(0:n) :: LFS

            LFS(0) = 1.0
            LFS(1) = x
            if(n==0) then
                  LFS(0) = LFS(0)
            elseif(n==1) then
                  LFS(1) = LFS(1)
            elseif(n>=2) then
                  do I = 2,n
                        LFS(I)= (-(I-1)*LFS(I-2) + (2*I-1)*x*LFS(I-1))/I
                  end do
            end if

      end subroutine

      function Legendre(n,x)
            implicit none
            integer :: n, I
            real(kind=8) :: x , Legendre
            real(kind=8),dimension(0:1) :: LFS
            LFS(0) = 1.0
            LFS(1) = x
            if(n==0) then
                  Legendre = LFS(0)
            elseif(n==1) then
                  Legendre = LFS(1)
            elseif(n>=2) then
                 do I = 2,n
                        Legendre = (-(I-1)*LFS(0) + (2*I-1)*x*LFS(1))/I
                        LFS(0) = LFS(1)
                        LFS(1) = Legendre
                  end do
            end if

      end function

      function ALegendreFnm(n,m,x)
            implicit none
            integer :: n,m,k,l
            real(kind=8) :: x,s
            real(kind=8) :: ALegendreFnm, Legendre
            real(kind=8),dimension(0:n+1) :: ALFS,BLFS,Pnm

            s = dsqrt(1.D0 - x*x)

            select case (n)
                  case (0)
                      ALegendreFnm = 1.D0
                  case (1)
                     select case (m)
                        case (0)
                              ALegendreFnm = x
                        case (1)
                              ALegendreFnm = s
                     end select

                  case default

                        ALFS(0) = (3.D0*x*x-1.D0)/2.D0
                        ALFS(1) = 3.D0*s*x
                        ALFS(2) = 3.D0*(1.D0-x*x)

                        if(n == 2) then
                              select case (m)
                              case (0)
                                    ALegendreFnm = ALFS(0)
                              case (1)
                                    ALegendreFnm = ALFS(1)
                              case (2)
                                    ALegendreFnm = ALFS(2)
                              end select
                        else

                             if(m == 0) then
                                    ALegendreFnm = Legendre(n,x)
                                    return
                              else
                                    BLFS(0) = x
                                    BLFS(1) = s
                                    l = 2
                                    do while(l< n)
                                          Pnm = 0.
                                          Pnm(0) = Legendre(l+1,x)
                                          k = 1
                                          do while(k<=m)
                                                if(k<=l-1) then
                                                     Pnm(k) = ((2*l+1)*x*ALFS(k)-(l+k)*BLFS(k))/(l-k+1)
                                                elseif(k>l-1 .and. k<=l+2)then
                                                     Pnm(k) = -((l-k+2)*x*Pnm(k-1)-(l+k)*ALFS(k-1))/s
                                                   !  Pnm(k) = 2.D0*(k+1)*x*Pnm(k-1)/s-(l-k)*(l+k+1)*Pnm(k-2)
                                                end if
                                                k = k + 1
                                          end do
                                          BLFS = ALFS
                                          ALFS = Pnm
                                          l = l + 1
                                    end do
                                    ALegendreFnm = Pnm(m)
                              end if
                       end if
            end select

      end function


subroutine FulNormALegenSerie(n,x,Pnm)
      implicit none
      integer :: n , l , k
      real(kind=8) :: x , s , Anm , Bnm
      real(kind=8) :: Legendre, ld
      real(kind=8),dimension(0:n+1) :: ALFS,BLFS,Pnm

      s = dsqrt(1.0 - x*x)
      Pnm = 0

      BLFS(1) = s*dsqrt(3.D0)

      ALFS(1) = 3.D0*s*x*dsqrt(5.D0/3.D0)
      ALFS(2) = 3.D0*(1.D0-x*x)*dsqrt(5.D0/12.D0)

      if(n == 1)  then
            Pnm = BLFS
      elseif(n==2) then
            Pnm = ALFS
      else
            l = 3
            do while(l<= n)
                  k = 1
                  do while(k <= l)
                        if(k < l-1) then
                              Anm = dsqrt(real((2*l-1)*(2*l+1),kind=8) &
                              /real((l-k)*(l+k),kind=8))

                              Bnm = -dsqrt(real(((l-1)*(l-1)-k*k),kind=8) &
                              /real(4*(l-1)*(l-1)-1,kind=8))

                              Pnm(k) = Anm*(x*ALFS(K) + Bnm*BLFS(K))

                        elseif(k >= l-1) then
                              if(k == l-1) then
                                    Pnm(k) = dsqrt(2.D0*real(k,kind=8)+3.D0)*x*ALFS(l-1)
                              elseif(k == l) then
                                    ld = real(k,kind=8)
                                    Pnm(k) = s*dsqrt(1.D0+1.D0/(2.D0*ld))*ALFS(l-1)
                              endif
                        end if
                        k = k + 1
                  end do
                  BLFS = ALFS
                  ALFS = Pnm
                  l = l + 1
            end do
      end if
      Pnm(0) = dsqrt(2D0*real(n,kind=8)+1.D0)*Legendre(n,x)

end subroutine


function normalGravity(latitute)

      implicit none

      real(kind=8) :: latitute , normalGravity
      real(kind=8) :: x, e4 , e6
      real(kind=8) :: Rlat, Term

      real(kind=8), dimension(4) ::  a2n

      real(kind=8), parameter ::  PI = 3.141592653589793238462643D0
      real(kind=8), parameter ::  DEGRAD = PI/180.D0
      real(kind=8), parameter ::  e2 = 0.00669438002290D0
      real(kind=8), parameter ::  k = 0.001931851353D0

      integer :: l

      Rlat = latitute*DEGRAD
      x = dsin(Rlat)

      Term= 0.D0

      a2n(1) = 0.5D0*e2 + k
      e4 = e2*e2
      a2n(2) = 3.D0*e4/8.D0 + 0.5D0*e2*k
      e6 = e4*e2
      a2n(3) = 5.D0*e6/16.D0 + 3.D0*e4*k/8.D0
      a2n(4) = 35.D0*e4*e4/128.D0 + 5.D0*e6*k/16.D0

      do l = 1 , 4
            Term = Term + a2n(l)*x**(2*l)
      end do

      normalGravity = 9.7803267715D0*(1.D0+Term)

end function

function CFgravity(latitute,Height)
  implicit none

      real(kind=8) :: latitute , Height , Rlat,CFgravity
      real(kind=8) :: E , beta, z ,dp2, rp,D, R
      real(kind=8) :: cosbeta, q0 , qp, W, rp2, bp
      real(kind=8) :: sinbeta , sinbeta2, E2, bp2
      real(kind=8) :: Omega2
      real(kind=8), parameter ::  PI = 3.141592653589793238462643D0
      real(kind=8), parameter ::  DEGRAD = PI/180.D0
      real(kind=8), parameter ::  a = 6378137.D0
      real(kind=8), parameter ::  b = 6356752.314D0
      real(kind=8), parameter ::  Omega = 0.7292115D-4
      real(kind=8), parameter ::  GM = 0.39860050000D+15

      Rlat = latitute*DEGRAD
      beta = datan(b*dtan(Rlat)/a)
      z = b*dsin(beta) + Height*dsin(Rlat)
      rp = a*dcos(beta) + Height*dcos(Rlat)
      dp2 = rp*rp - z*z
      rp2 = rp*rp + z*z
      E = dsqrt(a*a-b*b)
      E2 = E*E
      D = dp2/E2
      R = rp2/E2
      cosbeta = dsqrt(0.5D0+0.5D0*R-dsqrt(0.25D0+0.25D0*R*R-0.5D0*D))
      bp = dsqrt(rp2-E2*cosbeta*cosbeta)
      bp2 = bp*bp
      q0 = 0.5D0*((1.D0+3.D0*b*b/E2)*datan(E/b)-3.D0*b/E)
      qp = 3.D0*(1.D0+bp2/E2)*(1.D0-bp*datan(E/bp)/E)-1D0
      sinbeta = dsin(dacos(cosbeta))
      sinbeta2 = sinbeta*sinbeta
      W = dsqrt((bp*bp+E2*sinbeta2)/(bp*bp+E2))
      Omega2 = Omega*Omega

      CFGravity = (GM/(bp2+E2)+Omega2*a*a*E*qp*(0.5D0*sinbeta2-1.D0/6.D0)/ &
      ((bp2+E2)*q0)- Omega2*bp*cosbeta*cosbeta)/W

end function

function adaptedGravity(Latitute,Height)

      implicit none
      real(kind=8) :: Latitute, Height
      real(kind=8) :: x, x2, adaptedGravity
      real(kind=8), parameter :: ga = 9.7803267715D0
      real(kind=8), parameter ::  PI = 3.141592653589793238462643D0
      real(kind=8), parameter ::  DEGRAD = PI/180.D0

      x = dsin(Latitute*DEGRAD)
      x2 = x*x

      adaptedGravity = ga*(1+x2*(0.0052790414D0 + x2*(0.0000232718D0+ &
      x2*(0.0000001262D0+0.0000000007D0*x2))))-Height*(0.03087798D-4 - &
      x2*( 0.0000439D-4 + 0.00000020D-4*x2))- Height*Height*(-0.00007265D-8 &
       + 0.00000021D-8*x2)

end function



