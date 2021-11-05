
PROGRAM CHEMCPD 

USE RESIDUES
USE MINPACK
USE MODCTE

integer :: n = 2
real(8) , dimension(2) ::  x0, fvec
real(8) :: tol 
integer :: info 

x0(1) = -1.0 
x0(2) = -3.0 

call hybrd1 (Test_fcn, n, x0, fvec, tol, info)

write(*,*) "solution ", x0

write(*,*) "residue ", fvec 

write(*,*) p_R_hill/ c_R_jup

END PROGRAM