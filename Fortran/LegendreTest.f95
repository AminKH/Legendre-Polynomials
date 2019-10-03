! Program by
! Amin Khiabani
! aminkhiabani@outlook.com

subroutine ellipsoidal(a,method,latitute,ellipsoid)

      ! returns array of r, ellipsoidal radius 
      ! and ratio a to r , Rn
      ! a Earth eaquadorial radius
      ! method of calculationg r
      ! Geodetic or Geocentric latitude
      

      integer :: method
      real(kind=8) :: a, latitute
      real(kind=8) ,dimension(2) :: ellipsoid
      real(kind=8) :: x , x2, w2

      real(kind=8), parameter ::  PI = 3.141592653589793238462643D0
      real(kind=8), parameter ::  DEGRAD = PI/180.D0
      real(kind=8), parameter ::  e2 = 0.00669438002290D0
      real(kind=8), parameter ::  ep2 = 0.00673949677548D0

      x = dsin(latitute*DEGRAD)
      x2 = x*x

      if(method == 0)then
            ellipsoid(1) = a/Dsqrt(1+ep2*x2)
      else if(method == 1) then
            w2 = 1.D0-e2*x2
            ellipsoid(1) = a*dsqrt(1.D0 + e2*(e2-2.D0)*x2)/dsqrt(w2)
      end if
      ellipsoid(2) = a/ellipsoid(1)

end subroutine

subroutine EXgravity(method,C_Coef)
      ! returns C coefficient of liner normal gravity equation
      implicit none

      integer :: method
      real(kind=8) ::  Rphi, x , x2, s , s2 , latitude
      real(kind=8) ::  Rn , r, Legendre, normalGravity
      real(kind=8) :: Omega2
      real(kind=8) ,dimension(2) :: ellipsoid

      integer,PARAMETER :: N = 5
      integer,PARAMETER :: LDA = 2*N
      real(kind=8) ,dimension(0:N-1) ::  Phi, Gama
      real(kind=8) ,dimension(N) ::C_Coef

      REAL(kind=8) :: Ac(LDA,LDA), B(LDA), WORK(LDA), RCOND
      INTEGER   IWORK(LDA), I, J, ITASK, IND, k

      real(kind=8), parameter ::  PI = 3.141592653589793238462643D0
      real(kind=8), parameter ::  DEGRAD = PI/180.D0
      real(kind=8), parameter ::  a = 6378137.D0
      real(kind=8), parameter ::  ep2 = 0.0067394967754D0
      real(kind=8), parameter ::  Omega = 0.7292115D-4
      real(kind=8), parameter ::  GM = 0.39860050000D+15

      data Phi / 0.D0, 30.D0, 45.D0, 60.D0 ,90.D0 /

      Omega2 = Omega*Omega

      Gama(0) = 9.7803267715D0
      Gama(1) = normalGravity(30.D0)
      Gama(2) = 9.806199203D0
      Gama(3) = normalGravity(60.D0)
      Gama(4) = 9.8321863685D0

      C_Coef = 0.D0

      do i = 0, N-1
            call ellipsoidal(a,method,Phi(i),ellipsoid)
            Rphi = Phi(i)*DEGRAD
            x = dsin(Rphi)
            s = dcos(Rphi)
            x2 = x*x
            s2 = s*s
            r = ellipsoid(1)
            Rn = ellipsoid(2)
            B(i+1) = Gama(i)+ r*Omega2*s2
            do j=0,N-1
                  Ac(i+1,j+1) = GM*((2*j+1)*Legendre(2*j,x)*Rn**(2*j+1))/(r*a)
            end do
      end do

  !    DO k = 1,N
  !          WRITE (*,'(5(4X, F20.15))') (Ac(k,J), J = 1,N)
  !    end do
  !    WRITE (*,'(5(4X, F18.15))') (B(J), J = 1,N)

      ITASK  =  1
      CALL SGEFS (Ac, LDA, N, B, ITASK, IND, WORK, IWORK, RCOND)
      do i =1,N
         C_Coef(i) = B(i)
      end do
end subroutine


function Calcgravity(latitute,method)
      ! returns normal gravity of katitude
      implicit none
      integer :: method
      real(kind=8) ::  Rphi, x , x2, s , s2
      real(kind=8) ::  Rn , r, Legendre, Calcgravity
      real(kind=8) :: Omega2, Sigma, latitute
      real(kind=8) ,dimension(2) :: ellipsoid
      real(kind=8) ,dimension(5) ::C_Coef, C0, C1

      integer :: j

      real(kind=8), parameter ::  PI = 3.141592653589793238462643D0
      real(kind=8), parameter ::  DEGRAD = PI/180.D0
      real(kind=8), parameter ::  a = 6378137.D0
      real(kind=8), parameter ::  Omega = 0.7292115D-4
      real(kind=8), parameter ::  GM = 0.39860050000D+15

    !  call EXgravity(method,C_Coef)

      data C0 /0.999993551446001216D0 , &
              -0.001084161913952187D0 , &
               0.000004571732568181D0 , &
              -0.000000017043632669D0 , &
              0.000000000081735028D0/

      data C1 /1.000005569571394082D0 , &
                -0.001081298969678650D0 , &
                 0.000000471322575981D0 , &
                -0.000000000922562314D0 , &
                 0.000000000010459033D0 /

      if(method==0) then
            C_Coef = C0
      elseif(method ==1) then
            C_Coef = C1
      end if

      Rphi = latitute*DEGRAD
      x = dsin(Rphi)
      s = dcos(Rphi)
      x2 = x*x
      s2 = s*s
      Omega2 = Omega*Omega

      call ellipsoidal(a,method,latitute,ellipsoid)
      r = ellipsoid(1)
      Rn = ellipsoid(2)
      Sigma = 0.0
      do j=0,4
            Sigma = Sigma + GM*((2*j+1)*C_Coef(j+1)*Legendre(2*j,x)*Rn**(2*j+1))
      end do

     Calcgravity = Sigma/(r*a) - r*Omega2*s2

end function

subroutine FulNormALegenRange(n,x,Pnm)
      implicit none
      integer :: n , l , k , I
      integer(kind = 8) :: factorial
      real(kind=8) :: x , s , Anm , Bnm
      real(kind=8) :: Legendre, ld
      real(kind=8),dimension(0:n+1) :: ALFS,BLFS,Pnm

      s = dsqrt(1.0 - x*x)
      Pnm = 0

      if(n <= 8) then
            call ALegendFnmSeri(n,x,Pnm)
            do k = 1 , n
                  ld = dsqrt(2.D0*(2*n+1)*(factorial(n-k))/factorial(n+k))
            Pnm(k) = ld*Pnm(k)
            end do
       elseif(n > 8) then
            Pnm =  0.D0
            call ALegendFnmSeri(7,x,Pnm)
            do I = 1,7
                  ld = dsqrt(30.D0*real(factorial(7-I),kind=8)/ &
                  real(factorial(7+I),kind=8))
                  BLFS(I) = ld*Pnm(I)
            end do
            Pnm = 0
            call ALegendFnmSeri(8,x,Pnm)
            do I = 1,8
                  ld = dsqrt(34.D0*real(factorial(8-I),kind=8)/ &
                  real(factorial(8+I),kind=8))
                  ALFS(I) = ld*Pnm(I)
            end do
            l = 9
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


      function ALegendreFnm2(n,m,x)
            implicit none
            integer :: n,m,k,l
            real(kind=8) :: x,s
            real(kind=8) :: ALegendreFnm2 , Legendre
            real(kind=8),dimension(0:n) :: ALFS,BLFS,Pnm

            s = dsqrt(1.D0 - x*x)

            select case (n)
                  case (0)
                      ALegendreFnm2  = 1.D0
                  case (1)
                     select case (m)
                        case (0)
                              ALegendreFnm2  = x
                        case (1)
                              ALegendreFnm2  = s
                     end select
                  case (2)
                        select case (m)
                        case (0)
                              ALegendreFnm2  = (3.D0*x*x-1.D0)/2.D0
                        case (1)
                              ALegendreFnm2  = 3.D0*s*x
                        case (2)
                              ALegendreFnm2  = 3.D0*(1.D0-x*x)
                        end select
                  case (3)
                        select case (m)
                        case (0)
                              ALegendreFnm2  = x*(5.D0*x*x-3.D0)/2.D0
                        case (1)
                              ALegendreFnm2  = 3.D0*(5.D0*x*x-1.D0)*s/2.D0
                        case (2)
                              ALegendreFnm2  = 15.D0*x*(1.D0-x*x)
                        case (3)
                              ALegendreFnm2  = 15.D0*s*s*s
                        end select
                  case (4)
                        select case (m)
                        case (0)
                              ALegendreFnm2  = (35.D0*x*x*x*x-30.D0*x*x+3.)/8.D0
                        case (1)
                              ALegendreFnm2  = 5.D0*(7.D0*x*x-3.D0)*x*s/2.D0
                        case (2)
                              ALegendreFnm2  = 15.D0*(7.D0*x*x-1.D0)*(1.D0-x*x)/2.D0
                        case (3)
                              ALegendreFnm2  = 105.D0*s*s*s*x
                        case (4)
                              ALegendreFnm2  = 105.D0*s*s*s*s
                        end select
                  case (5)
                        select case (m)
                        case (0)
                              ALegendreFnm2  = x*(63.D0*x*x*x*x-70.D0*x*x+15.D0)/8.D0
                        case (1)
                              ALegendreFnm2  = 15.D0*s*(21.D0*x*x*x*x-14.D0*x*x+1.D0)/8.D0
                        case (2)
                              ALegendreFnm2  = 105.D0*x*(3.0*x*x-1.D0)*(1.D0-x*x)/2.D0
                        case (3)
                              ALegendreFnm2  = 105.D0*s*s*s*(9.D0*x*x-1.)/2.D0
                        case (4)
                              ALegendreFnm2  = 945.D0*s*s*s*s*x
                        case (5)
                              ALegendreFnm2  = 945.D0*s*s*s*s*s
                        end select
                  case (6)
                        select case (m)

                        case (0)
                              ALegendreFnm2  = (x*x*(x*x*(231.D0*x*x-315.D0)+105.D0)-5.D0)/16.D0
                        case (1)
                              ALegendreFnm2  = 21.D0*x*(x*x*(33.D0*x*x-30.D0)+5.D0)*s/8.D0
                        case (2)
                              ALegendreFnm2  = 105.D0*s*s*(x*x*(33.D0*x*x-18.)+1.D0)/8.D00
                        case (3)
                              ALegendreFnm2  = 315.D0*(11.D0*x*x-3.D0)*x*s*s*s/2.D0
                        case (4)
                              ALegendreFnm2  = 945.D0*s*s*s*s*(11.D0*x*x-1.D0)/2.D0
                        case (5)
                              ALegendreFnm2  = 10395.D0*x*s**5
                        case (6)
                              ALegendreFnm2  = 10395.D0*s**6
                       end select

                  case (7)
                        select case (m)

                        case (0)
                              ALegendreFnm2  = x*(x*x*(429.D0*x**4-693.D0*x*x+315.D0)-35.D0) /16.D0
                        case (1)
                              ALegendreFnm2  = 7.*s*(x*x*(429.D0*x**4-495.D0*x*x+135.D0)-5.D0)/16.D0
                        case (2)
                              ALegendreFnm2  = 63.D0*x*s*s*(x*x*(143.D0*x*x-110.D0)+15.D0)/8.D0
                        case (3)
                              ALegendreFnm2  = 315.D0*s*s*s*(x*x*(143.D0*x*x-66.D0)+3.D0)/8.D0
                        case (4)
                              ALegendreFnm2  = 3465.D0*x*s*s*s*s*(13.D0*x*x-3.D0)/2.D0
                        case (5)
                              ALegendreFnm2  = 10395.D0*(S**5)*(13.D0*x*x-1.D0)/2.D0
                        case (6)
                              ALegendreFnm2  = 135135.D0*x*s**6
                        case (7)
                              ALegendreFnm2  = 135135.D0*s**7
                       end select

                  case default

                       ALFS(0) = (x*x*(x*x*(x*x*(6435.D0*x*x-12012.D0)+6930.D0)-1260.D0)+35.D0)/128.D0
                       ALFS(1) = 9.D0*x*s*(x*x*(x*x*(715.D0*x*x-1001.D0)+385.D0)-35.D0)/16.D0
                       ALFS(2) = 315.D0*s*s*(x*x*(x*x*(143.D0*x*x-143.D0)+33.D0)-1)/16.D0
                       ALFS(3) = 3465.D0*x*s*s*s*(x*x*(39.D0*x*x-26.D0)+3.D0)/8.D0
                       ALFS(4) = 10395.D0*s*s*s*s*(65.D0*x**4-26.D0*x*x+1.D0)/8.D0
                       ALFS(5) = 135135.D0*(x*s**5)*(5.D0*x*x-1.D0)/2.D0
                       ALFS(6) = 135135.D0*(s**6)*(15.D0*x*x-1.D0)/2.D0
                       ALFS(7) = 2027025.D0*x*s*s*s*s*s*s*s
                       ALFS(8) = 2027025.*s*s*s*s*s*s*s*s

                       if(n == 8) then
                              select case (m)
                              case (0)
                                    ALegendreFnm2  = ALFS(0)
                              case (1)
                                    ALegendreFnm2  = ALFS(1)
                              case (2)
                                    ALegendreFnm2  = ALFS(2)
                              case (3)
                                    ALegendreFnm2  = ALFS(3)
                              case (4)
                                    ALegendreFnm2  = ALFS(4)
                              case (5)
                                    ALegendreFnm2  = ALFS(5)
                              case (6)
                                    ALegendreFnm2  = ALFS(6)
                              case (7)
                                    ALegendreFnm2  = ALFS(7)
                              case (8)
                                    ALegendreFnm2  = ALFS(8)
                              end select
                        else
                             if(m == 0) then
                                    ALegendreFnm2  = Legendre(n,x)
                                    return
                              else
                                   BLFS(0) = x*(x*x*(429.D0*x**4-693.D0*x*x+315.D0)-35.D0) /16.D0
                                   BLFS(1) = 7.D0*s*(x*x*(429.D0*x**4-495.D0*x*x+135.D0)-5.D0)/16.D0
                                   BLFS(2) = 63.D0*x*s*s*(x*x*(143.D0*x*x-110.D0)+15.D0)/8.D0
                                   BLFS(3) = 315.D0*s*s*s*(x*x*(143.D0*x*x-66.D0)+3.D0)/8.D0
                                   BLFS(4) = 3465.D0*x*s*s*s*s*(13.D0*x*x-3.D0)/2.D0
                                   BLfS(5) = 10395.D0*(S**5)*(13.D0*x*x-1.D0)/2.D0
                                   BLFS(6) = 135135.D0*x*s**6
                                   BLFS(7) = 135135.D0*s**7

                                    l = 8
                                    do while(l< n)
                                          Pnm(0) = Legendre(l+1,x)
                                          k = 1
                                          do while(k<=m)
                                                if(k<=l-1) then
                                                     Pnm(k) = ((2*l+1)*x*ALFS(k)-(l+k)*BLFS(k))/(l-k+1)
                                                elseif(k>l-1)then
                                                     Pnm(k) = -((l-k+2)*x*Pnm(k-1)-(l+k)*ALFS(k-1))/s
                                                end if
                                                k = k + 1
                                          end do
                                          BLFS = ALFS
                                          ALFS = Pnm
                                          l = l + 1
                                    end do
                                    ALegendreFnm2  = Pnm(m)

                              end if
                       end if
            end select

      end function

      subroutine ALegendFnmSeri(n,x,Pnm)
            implicit none
            integer :: n,k,l
            real(kind=8) :: x,s
            real(kind=8) :: Legendre
            real(kind=8),dimension(0:n) :: ALFS,BLFS,Pnm

            s = dsqrt(1.D0 - x*x)

            select case (n)
                  case (0)
                        Pnm(0) = 1.D0
                  case (1)
                        Pnm(0) = x
                        Pnm(1) = s
                  case (2)
                        Pnm(0) = (3.D0*x*x-1.D0)/2.D0
                        Pnm(1) = 3.D0*s*x
                        Pnm(2) = 3.D0*(1.D0-x*x)
                  case (3)
                        Pnm(0) = x*(5.D0*x*x-3.D0)/2.D0
                        Pnm(1) = 3.D0*(5.D0*x*x-1.D0)*s/2.D0
                        Pnm(2) = 15.0*x*(1.D0-x*x)
                        Pnm(3) = 15.D0*s*s*s
                  case (4)
                        Pnm(0) = (35.D0*x*x*x*x-30.D0*x*x+3.D0)/8.D0
                        Pnm(1) = 5.D0*(7.D0*x*x-3.D0)*x*s/2.D0
                        Pnm(2) = 15.D0*(7.D0*x*x-1.D0)*(1.D0-x*x)/2.D0
                        Pnm(3) = 105.D0*s*s*s*x
                        Pnm(4) = 105.D0*s*s*s*s
                  case (5)
                        Pnm(0) = x*(63.D0*x*x*x*x-70.D0*x*x+15.D0)/8.D0
                        Pnm(1) = 15.D0*s*(21.D0*x*x*x*x-14.D0*x*x+1.D0)/8.D0
                        Pnm(2) = 105.D0*x*(3.D0*x*x-1.D0)*(1.D0-x*x)/2.D0
                        Pnm(3) = 105.D0*s*s*s*(9.D0*x*x-1.D0)/2.D0
                        Pnm(4) = 945.D0*s*s*s*s*x
                        Pnm(5) = 945.D0*s*s*s*s*s
                  case (6)
                       Pnm(0) = (x*x*(x*x*(231.D0*x*x-315.D0)+105.D0)-5.D0)/16.D0
                       Pnm(1) = 21.D0*x*(x*x*(33.D0*x*x-30.D0)+5.D0)*s/8.D0
                       Pnm(2) = 105.D0*s*s*(x*x*(33.D0*x*x-18.D0)+1.D0)/8.D0
                       Pnm(3) = 315.D0*(11.D0*x*x-3.D0)*x*s*s*s/2.D0
                       Pnm(4) = 945.D0*s*s*s*s*(11.D0*x*x-1.D0)/2.D0
                       Pnm(5) = 10395.D0*x*s**5
                       Pnm(6) = 10395.D0*s**6
                  case (7)
                        Pnm(0) = x*(x*x*(429.D0*x**4-693.D0*x*x+315.D0)-35.D0)/16.D0
                        Pnm(1) = 7.D0*s*(x*x*(429.D0*x**4-495.D0*x*x+135.D0)-5.D0)/16.D0
                        Pnm(2) = 63.D0*x*s*s*(x*x*(143.D0*x*x-110.D0)+15.D0)/8.D0
                        Pnm(3) = 315.D0*s*s*s*(x*x*(143.D0*x*x-66.D0)+3.D0)/8.D0
                        Pnm(4) = 3465.D0*x*s*s*s*s*(13.D0*x*x-3.D0)/2.D0
                        Pnm(5) = 10395.D0*(S**5)*(13.D0*x*x-1.D0)/2.D0
                        Pnm(6) = 135135.D0*x*s**6
                        Pnm(7) = 135135.D0*s**7
                  case default

                       ALFS(0) = (x*x*(x*x*(x*x*(6435.D0*x*x-12012.D0)+6930.D0)-1260.D0)+35)/128.D0
                       ALFS(1) = 9.D0*x*s*(x*x*(x*x*(715.D0*x*x-1001.D0)+385.D0)-35.D0)/16.D0
                       ALFS(2) = 315.D0*s*s*(x*x*(x*x*(143.D0*x*x-143.D0)+33.D0)-1.D0)/16.D0
                       ALFS(3) = 3465.D0*x*s*s*s*(x*x*(39.D0*x*x-26.D0)+3.D0)/8.D0
                       ALFS(4) = 10395.D0*s*s*s*s*(65.D0*x**4-26.D0*x*x+1.D0)/8.D0
                       ALFS(5) = 135135.D0*(x*s**5)*(5.D0*x*x-1.D0)/2.D0
                       ALFS(6) = 135135.D0*(s**6)*(15.D0*x*x-1.D0)/2.D0
                       ALFS(7) = 2027025.D0*x*s*s*s*s*s*s*s
                       ALFS(8) = 2027025.D0*s*s*s*s*s*s*s*s

                       if(n == 8) then
                              Pnm = ALFS
                        else
                             BLFS(0) = x*(x*x*(429.D0*x**4-693.D0*x*x+315.D0)-35.D0)/16.D0
                             BLFS(1) = 7.D0*s*(x*x*(429.D0*x**4-495.D0*x*x+135.D0)-5.D0)/16.D0
                             BLFS(2) = 63.D0*x*s*s*(x*x*(143.D0*x*x-110.D0)+15.D0)/8.D0
                             BLFS(3) = 315.D0*s*s*s*(x*x*(143.D0*x*x-66.D0)+3.D0)/8.D0
                             BLFS(4) = 3465.D0*x*s*s*s*s*(13.D0*x*x-3.D0)/2.D0
                             BLfS(5) = 10395.D0*(S**5)*(13.D0*x*x-1.D0)/2.D0
                             BLFS(6) = 135135.D0*x*s**6
                             BLFS(7) = 135135.D0*s**7

                              l = 8
                              do while(l< n)
                                    Pnm(0) = Legendre(l+1,x)
                                    k = 1
                                    do while(k<=l+1)
                                          if(k<=l-1) then
                                                Pnm(k) = ((2*l+1)*x*ALFS(k)-(l+k)*BLFS(k))/(l-k+1)
                                          elseif(k>l-1)then
                                                Pnm(k) = -((l-k+2)*x*Pnm(k-1)-(l+k)*ALFS(k-1))/s
                                          end if
                                          k = k + 1
                                    end do
                                    BLFS = ALFS
                                    ALFS = Pnm
                                    l = l + 1
                              end do
                       end if
            end select
      end subroutine


function FullNormALegendre(n,m,x)
      implicit none
      integer :: n , m , l , k , I
      integer(kind = 8) :: factorial
      real(kind=8) :: x , s , FullNormALegendre, Anm , Bnm
      real(kind=8) :: Legendre,ALegendreFnm2 , ld
       real(kind=8),dimension(0:n) :: ALFS,BLFS,Pnm

      s = dsqrt(1.0 - x*x)
      ld = 0.D0
      if(m == 0) then
            FullNormALegendre = dsqrt(2D0*n+1.D0)*Legendre(n,x)
      elseif(n <= 8) then
            ld = dsqrt(2.D0*((2*n+1))*(factorial(n-m))/(factorial(n+m)))
            FullNormALegendre = ld*ALegendreFnm2 (n,m,x)
      elseif(n > 8) then
            Pnm =  0.D0
            call ALegendFnmSeri(7,x,Pnm)
            do I = 1,7
                  ld = dsqrt(30.D0*factorial(7-I)/factorial(7+I))
                  BLFS(I) = ld*Pnm(I)
            end do
            Pnm = 0
            call ALegendFnmSeri(8,x,Pnm)
            do I = 1,8
                  ld = dsqrt(34.D0*factorial(8-I)/factorial(8+I))
                  ALFS(I) = ld*Pnm(I)
            end do
            l = 9
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
                                    Pnm(k) = dsqrt(2.D0*k+3.D0)*x*ALFS(l-1)
                              elseif(k == l) then
                                    Pnm(k) = s*dsqrt(1.D0+1.D0/(2.D0*k))*ALFS(l-1)
                              endif
                        end if
                        k = k + 1
                  end do
                  BLFS = ALFS
                  ALFS = Pnm
                  l = l + 1
            end do
            FullNormALegendre = Pnm(m)
      end if

end function


subroutine gravitation(phi,lambda,fileName,method,gravity)
      implicit none
      real(kind=8)  :: gravity
      real(kind=8)  :: Rn , rep
      real(kind=8)::  GM, R, C ,S, Phi, Lambda, Rlam
      real(kind=8):: Sum1 , Tsum , x, Coef
      real(kind=8) , allocatable , dimension(:) ::Pn
      integer, parameter :: linelength = 120
      character (len=linelength) :: line,fileName
      integer :: ierror,L,M,N,method
      real(kind=8) , dimension(2) :: ellipsoid
      real(kind=8), parameter ::  PI = 3.141592653589793238462643D0
      real(kind=8), parameter ::  DEGRAD = PI/180.D0

      Rlam = Lambda*DEGRAD

      sum1 = 0.D0
      Tsum = 0.D0

      open (unit=10, file=fileName,action ="read", status="old",iostat=ierror)
      if ( ierror /= 0 ) then
            print*, "Failed to open file!"
            stop
      end if
      readfile : do
            read (unit=10, fmt="(a)", iostat=ierror) line
            if (ierror < 0 ) then
                  print*, "end of file reached"
                  exit readfile
            else if (ierror > 0) then
                  print*, "Error during read"
                  exit readfile
            else

                  if (index (line, "earth_gravity_constant") > 0) then
                        read(line(24:53),*)GM
                        print*,'earth_gravity_constant, GM0= ', GM
                  elseif(index (line, "radius") > 0) then
                        read(line(24:53),*)R
                        print*,'earth radius R= ', R
                  elseif(index (line, "max_degree") > 0) then
                        read(line(24:53),*)n
                         print*,'max_degree n= ', n
                        exit readfile
                  end if
            end if
      end do readfile

      call ellipsoidal(R,method,phi,ellipsoid)
      rep = ellipsoid(1)
      Rn = ellipsoid(2)
      x = dsin(Phi*DEGRAD)
      allocate(Pn(0:n))

   readgfc: do
      read (unit=10, fmt="(a)", iostat=ierror) line
            if (ierror < 0 ) then
                  print*, "end of file reached"
                  exit readgfc
            else if (ierror > 0) then
                  print*, "Error during read"
                  exit readgfc
            else
                  if(index (line, "gfc") > 0) then
                        line = adjustl(line(4:len_trim(line)))
                        read(line(1:scan(line,' ')),*)L
                        line = adjustl(line(scan(line,' '):len_trim(line)))
                        read(line(1:scan(line,' ')),*) M
                        line = adjustl(line(scan(line,' '):len_trim(line)))
                        read(line(1:scan(line,' ')),*) C
                        line = adjustl(line(scan(line,' '):len_trim(line)))
                        read(line(1:scan(line,' ')),*) S
                        if(M==0) then
                            Pn = 0.D0
                            call FulNormALegenSerie(L,x,Pn)
                            Sum1 = Pn(M)*C
                        else
                            Sum1 = Sum1 + Pn(M)*(C*dcos(m*Rlam)+S*dsin(m*Rlam))
                        end if
                        if(M==L) then
                              Tsum = Tsum + Sum1*(L-1)*Rn**L
                        !      print*,L,M,C,S
                        end if
                  end if
            end if

      end do readgfc
      close (unit=10)
      deallocate(Pn)

      Coef = GM/(rep*rep)
      gravity = Tsum*Coef

end subroutine
