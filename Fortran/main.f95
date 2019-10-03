! A fortran95 program for G95
! By Amin Khiabani
! aminkhiabani@outlook.com
program main
  implicit none

  integer :: j
  real(kind=8):: Calcgravity
  real(kind=8):: Phi ,Lambda ,gravity
  real(kind=8):: CFgravity
  real(kind=8):: normalGravity
  real(kind=8):: grav,adaptedGravity,calcg0,calcg1,adapg

  real(kind=8) ,dimension(5) :: C_Coef

  integer, parameter :: linelength = 120
  character (len=linelength) ::  fileName


  print*,
  call EXgravity(0,C_Coef)
  WRITE (*,'(5(4X, F21.18))') (C_Coef(J), J = 1,5)
  print*,

  print*,
  call EXgravity(1,C_Coef)
  WRITE (*,'(5(4X, F21.18))') (C_Coef(J), J = 1,5)
  print*,
   print*,'------------------------------------------------'

  Phi = 0.D0

  do while(Phi<=90.D0)
      calcg0 = Calcgravity(Phi,0)*1D5
      calcg1 = Calcgravity(Phi,1)*1D5
      adapg = normalGravity(Phi)*1D5
      print*,Phi, calcg0,calcg1,adapg,normalGravity(Phi)*1D5
      Phi = Phi + 5.D0
      print*,
  end do
  
  !
  ! Phi = 45.D0
  ! Lambda = 47.D0

  ! fileName = "file directory\dV_ELL_Earth2014.GFC"

  ! call gravitation(phi,Lambda,fileName,1,grav)

  !  write(*,'(A,F8.4,X,A,F8.4,X,A,3X,F16.14)') &
  ! 'Delta Gravity of latitude ',Phi,'Longitude', Lambda ,'is',grav
  !      print*,
  !      print*,'------------------------------------------------'
  !      print*, Phi, Lambda
  !      gravity = adaptedGravity(Phi,0.D0)
  !      print*, gravity, CFgravity(Phi,0.D0)
  !      print*,normalGravity(Phi)+grav,gravity+grav
end
