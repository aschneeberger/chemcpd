
PROGRAM CHEMCPD 

USE RESIDUES
USE MINPACK
USE MODCTE
USE DSKPHY

integer , parameter :: n = 5
real(8) , dimension(2) ::  x0, fvec
real(8) :: tol 
integer :: info 
double precision , dimension(n) :: r,b
double precision :: R_c

x0(1) = -1.0 
x0(2) = -3.0 

r = [1.0,2.0,6.0,10.0,30.0] * c_R_jup

R_c = centrifugal_radius()


call hybrd1 (Test_fcn, 2, x0, fvec, tol, info)

write(*,*) "solution ", x0

write(*,*) "residue ", fvec 

write(*,*) c_pi

END PROGRAM