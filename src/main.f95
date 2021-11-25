
PROGRAM CHEMCPD 

USE RESIDUES
USE MINPACK
USE MODCTE
USE DSKPHY
USE QUADPACK
USE PARTFUN

!Loop variables 
integer :: i

!Physical parameters 
double precision , dimension(p_Nr) :: r
double precision  :: R_c
double precision , dimension(p_Nr) :: cap_lambda
double precision , dimension(p_Nr) :: omegak
double precision , dimension(p_Nr) :: F_vis
double precision , dimension(p_Nr) :: F_acc 
double precision , dimension(p_Nr) :: T_mid
double precision , dimension(p_Nr) :: T_s
double precision , dimension(p_Nr) :: rho_mid
double precision , dimension(p_Nr) :: rho_add
double precision , dimension(p_Nr) :: rho_s
double precision , dimension(p_Nr) :: z_add
double precision , dimension(p_Nr) :: z_s
double precision , dimension(p_Nr) :: sigma
double precision , dimension(p_Nr) :: kappa_p

!Hybr1 subroutine parameters: 
DOUBLE PRECISION , DIMENSION(p_Nr*8) :: X ! array containin all variables (8*p_Nr) 
DOUBLE PRECISION , dimension(p_Nr*8) :: fvec ! Residues array 
double precision :: tol = 1.0d-2 !asked relative error in the resolution 
double precision , dimension(p_Nr*5) :: args ! constant arguments in function to optimize

integer :: info !output code of the solver : 
!    0, improper input parameters.
!    1, algorithm estimates that the relative error between X and the
!       solution is at most TOL.
!    2, number of calls to FCN has reached or exceeded 200*(N+1).
!    3, TOL is too small.  No further improvement in the approximate
!       solution X is possible.
!    4, the iteration is not making good progress.


!!!!!!!!!!!!!!!!!!!!!!!!
!   Initialization     ! 
!!!!!!!!!!!!!!!!!!!!!!!!


!Linear grid 
forall(i = 1:p_Nr) r(i) = (p_R_disk - p_R_p) * float(i) / float(p_Nr) + p_R_p

Write(*,*) "[MAIN] Grid generated "

!Initialize the profiles 
call Init_profiles(p_Nr,r,cap_lambda,R_c,omegak,F_vis,F_acc,T_mid,T_s,rho_mid,rho_add,rho_s,z_add,z_s,sigma,kappa_p)
Write(*,*) "[MAIN] Guesses Initialized "

! Write the initization in a file in table format (for astropy table use)
open(unit=10, file='../Data/initialisation.dat',status='new')

write(10,*) 'r cap_lambda omegak F_vis F_acc T_mid T_s rho_mid rho_add rho_s z_add z_s sigma kappa_p'

do i = 1,p_Nr 
    write(10,*) r(i),cap_lambda(i),omegak(i),F_vis(i),F_acc(i),T_mid(i),T_s(i) &
    &,rho_mid(i),rho_add(i),rho_s(i),z_add(i),z_s(i),sigma(i),kappa_p(i)
end do 
close(unit=10) 
Write(*,*) "[MAIN] Guesses Written "

!!!!!!!!!!!!!!!!!!!!!!!!
!      Resolution      ! 
!!!!!!!!!!!!!!!!!!!!!!!!

!Create the variable to be parsed in the solver subroutine 
x = [sigma,T_mid,T_s,z_s,z_add,rho_mid,rho_add,rho_s]

!Create the argument to be parsed in Equation_system_ms
args = [cap_lambda,omegak,F_vis,F_acc,r]

Write(*,*) "[MAIN] Begining of solving "

!Lauch the solver 
call hybrd1 (Equation_system_ms, p_Nr*8, x, fvec, tol, info , 5*p_Nr , args)
Write(*,*) "[MAIN] End of solving "
!Parse the solution
sigma = x(1 : p_Nr)
T_mid = x(p_Nr+1 : 2*p_Nr)
T_s = x(2*p_Nr+1 : 3*p_Nr)
z_s = x(3*p_Nr+1 : 4*p_Nr) 
z_add = x(4*p_Nr+1 : 5*p_Nr)
rho_mid = x(5*p_Nr+1 : 6*p_Nr)
rho_add = x(6*p_Nr+1 : 7*p_Nr)
rho_s = x(7*p_Nr+1 : 8*p_Nr)

!Write the solution in a file

open(unit=10, file='../Data/Solution.dat',status='new')

write(10,*) 'r cap_lambda omegak F_vis F_acc T_mid T_s rho_mid rho_add rho_s z_add z_s sigma kappa_p'

do i = 1,p_Nr 
    write(10,*) r(i),cap_lambda(i),omegak(i),F_vis(i),F_acc(i),T_mid(i),T_s(i) &
    &,rho_mid(i),rho_add(i),rho_s(i),z_add(i),z_s(i),sigma(i),kappa_p(i)
end do 
close(unit=10) 
Write(*,*) "[MAIN] Solution Written "

END PROGRAM