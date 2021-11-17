
PROGRAM CHEMCPD 

USE RESIDUES
USE MINPACK
USE MODCTE
USE DSKPHY
USE QUADPACK
USE PARTFUN

integer :: i ! index variable
real(8) , dimension(2) ::  x0, fvec
real(8) :: tol 
integer :: info 

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

!Linear grid 
forall(i = 1:p_Nr) r(i) = (p_R_disk - p_R_p) * float(i) / float(p_Nr) + p_R_p

write(*,*) r

call Init_profiles(p_Nr,r,cap_lambda,R_c,omegak,F_vis,F_acc,T_mid,T_s,rho_mid,rho_add,rho_s,z_add,z_s,sigma)

 
open(unit=10, file='initialisation.dat',status='new')

write(10,*) 'r(i),cap_lambda(i),omegak(i),F_vis(i),F_acc(i),T_mid(i),T_s(i),rho_mid(i),rho_add(i),rho_s(i),z_add(i),z_s(i),sigma(i)'

do i = 1,p_Nr 

write(10,*) r(i),cap_lambda(i),omegak(i),F_vis(i),F_acc(i),T_mid(i),T_s(i),rho_mid(i),rho_add(i),rho_s(i),z_add(i),z_s(i),sigma(i)

end do 
close(unit=10) 

END PROGRAM

