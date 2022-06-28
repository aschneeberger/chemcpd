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
DOUBLE PRECISION , DIMENSION(2*p_Nr) :: X, sol, res ! array containin all variables (8*p_Nr) 
double precision , dimension(5*p_Nr) :: args ! constant arguments in function to optimize


!!!!!!!!!!!!!!!!
! Env creation !
!!!!!!!!!!!!!!!!

call init_env()
 

!i=run_test()
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

call write_file('initialisation.csv',&
&[r,cap_lambda,omegak,F_vis,F_acc,T_mid,T_s,rho_mid,rho_add,rho_s,z_add,z_s,sigma,kappa_p],p_Nr,14,&
&'r,cap_lambda,omegak,F_vis,F_acc,T_mid,T_s,rho_mid,rho_add,rho_s,z_add,z_s,sigma,kappa_p,')

! write all used physical constants 
call write_file('physical_constants.csv', &
&[c_au,c_G,c_gamma,c_L_sun,c_M_earth,c_M_jup,c_M_sun,c_R_jup,c_Rg],&
&1,9,'au,G,gamma,L_sun,M_earth,M_jup,M_sun,R_jup,Rg,')

! Write all disk parameters constants
call write_file('disk_parameters.csv',&
&[p_a_p,p_alpha,p_Chi,p_Ks,p_L,p_L_p,p_M_dot,p_M_p,p_mu_gas,p_R_disk,p_R_grid,p_R_hill,p_R_p,p_T_neb,p_Xd],&
&1,15,'a_p,alpha,Chi,Ks,L,L_p,M_dot,M_p,mu_gas,R_disk,R_grid,R_hill,R_p,T_neb,Xd,')

Write(30,*) "[MAIN] Guesses Written "

!!!!!!!!!!!!!!!!!!!!!!!!
!      Resolution      !
!!!!!!!!!!!!!!!!!!!!!!!!


!Create the variable to be parsed in the solver subroutine 

!Create the argument to be parsed in Equation_system_ms


Write(30,*) "[MAIN] Begining of solving "

!Lauch the solver

x = 1d4![T_mid,T_s]
args = [cap_lambda,omegak,F_vis,F_acc,r]

sol = solve_JFNK(p_Nr*2,Heller_eq_sys,boundary_heller_sys,x,5*p_Nr,args,10000,1d-5,1d-5)

close(unit=10)

!Parse the solution
T_mid = sol(1 : p_Nr)
T_s = sol(p_Nr+1 : 2*p_Nr)

res = Heller_eq_sys(2*p_Nr,sol,5*p_Nr,args)



!Write the solution in a file

open(unit=10, file=Trim(env_datapath)//'/Solution.dat',status='new')

write(10,*) 'r cap_lambda omegak F_vis F_acc T_mid T_s res'

do i = 1,p_Nr 
    write(10,*) r(i),cap_lambda(i),omegak(i),F_vis(i),F_acc(i),T_mid(i),T_s(i),res(i)
end do 
close(unit=10) 
Write(30,*) "[MAIN] Solution Written "

close(unit=30)
END PROGRAM CHEMCPD
