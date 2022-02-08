!! @author : Antoine Schneeberger 
!! @mail : antoine.schneeberger@lam.fr

PROGRAM CHEMCPD 


USE RESIDUES
USE MODCTE
USE DSKPHY
USE QUADPACK
USE PARTFUN
USE ENV
USE JFNK

IMPLICIT NONE 

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
double precision , dimension(p_Nr*5) :: args ! constant arguments in function to optimize


integer :: info !output code of the solver : 
!    0, improper input parameters.
!    1, algorithm estimates that the relative error between X and the
!       solution is at most TOL.
!    2, number of calls to FCN has reached or exceeded 200*(N+1).
!    3, TOL is too small.  No further improvement in the approximate
!       solution X is possible.
!    4, the iteration is not making good progress.


!!!!!!!!!!!!!!!!
! Env creation !
!!!!!!!!!!!!!!!!

call init_env()

!info=run_test()
!!!!!!!!!!!!!
! LOGS open !
!!!!!!!!!!!!!

!creation of the log file in datapath

open(unit=30,file=trim(env_datapath)//'/logs.log', status='new')

!!!!!!!!!!!!!!!!!!!!!!!!
!   Initialization     ! 
!!!!!!!!!!!!!!!!!!!!!!!!

!Check if asked grid radius is valid
! else it will not continue  further and output an error
if (p_R_grid < 10d0*p_R_p) then 
    write(30,*)'[MAIN] ERROR Asked grid size is to small, is must be at least 10 time the planet radius'
    write(30,*) '   Asked grid size:', p_R_grid, 'Planet Radius:', p_R_p
    stop
end if 

!Grid construction 

if (p_gtype == 1) then  
    !Linear grid

    Write(30,*) '[MAIN] Construction of Linear grid between :' , p_R_p*2.0d0, 'and'&
    &, p_R_grid,'with',p_Nr,'points'

    forall(i = 1:p_Nr) r(i) = (p_R_grid - p_R_p*2.0d0) * float(i) / float(p_Nr) + p_R_p*2.0d0

else if (p_gtype == 2) then
    !Log Grid 

    Write(30,*) '[MAIN] Construction of Log grid between :' , p_R_p*2.0d0/c_R_jup, 'and'&
    &, p_R_grid/c_R_jup,'with',p_Nr,'points'

    forall(i = 1:p_Nr) r(i) = 2.0d0*p_R_p * exp(float(i-1)/float(p_Nr-1) * log(p_R_grid/(2.0d0*p_R_p)) )
else 
    ! if other values, unknown grid type, output error and stop
    write(30,*) '[MAIN] ERROR Unknown gride type, gtype must be either 1 or 2'
    stop
end if 

Write(30,*) "[MAIN] Grid generated "

!Initialize the profiles 
call Init_profiles(p_Nr,r,cap_lambda,R_c,omegak,F_vis,F_acc,T_mid,T_s,rho_mid,rho_add,rho_s,z_add,z_s,sigma,kappa_p)
Write(30,*) "[MAIN] Guesses Initialized "
if (p_verbose) write(30,*) '[GUESSES] Centrigucal Radius Rc computed as:', R_c/c_R_jup , 'R_jup'

! Write the initization in a file in table format (for astropy table use)
open(unit=10, file=Trim(env_datapath)//'/initialisation.dat',status='new')

write(10,*) 'r cap_lambda omegak F_vis F_acc T_mid T_s rho_mid rho_add rho_s z_add z_s sigma kappa_p'

do i = 1,p_Nr 
    write(10,*) r(i),cap_lambda(i),omegak(i),F_vis(i),F_acc(i),T_mid(i),T_s(i) &
    &,rho_mid(i),rho_add(i),rho_s(i),z_add(i),z_s(i),sigma(i),kappa_p(i)
end do 
close(unit=10) 
Write(30,*) "[MAIN] Guesses Written "

!!!!!!!!!!!!!!!!!!!!!!!!
!      Resolution      !
!!!!!!!!!!!!!!!!!!!!!!!!


!Create the variable to be parsed in the solver subroutine 
x = [sigma,T_mid,T_s,z_s,z_add,rho_mid,rho_add,rho_s]


!Create the argument to be parsed in Equation_system_ms
args = [cap_lambda,omegak,F_vis,F_acc,r]


Write(30,*) "[MAIN] Begining of solving "

!Lauch the solver 

x = solve_JFNK(p_Nr*8,Equation_system_ms,x,5*p_Nr,args,1.0d-5,1000)

Write(30,*) "[MAIN] End of solving, exit status :" , info

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

open(unit=10, file=Trim(env_datapath)//'/Solution.dat',status='new')

write(10,*) 'r cap_lambda omegak F_vis F_acc T_mid T_s rho_mid rho_add rho_s z_add z_s sigma kappa_p'

do i = 1,p_Nr 
    write(10,*) r(i),cap_lambda(i),omegak(i),F_vis(i),F_acc(i),T_mid(i),T_s(i) &
    &,rho_mid(i),rho_add(i),rho_s(i),z_add(i),z_s(i),sigma(i),kappa_p(i)
end do 
close(unit=10) 
Write(30,*) "[MAIN] Solution Written "

close(unit=30)
END PROGRAM