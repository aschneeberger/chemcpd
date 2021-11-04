
PROGRAM CHEMCPD 
USE RESIDUE
USE SOLVER
USE minipack

integer :: n = 2
real(8) , dimension(2) ::  x0, fvec
real(8) :: tol 
integer :: info 

WRITE(*,*) "HELLO WORLD" 

x0(1) = 1.0 
x0(2) = 1.0 

call hybrd1 (Fcn, n, x0, fvec, tol, info )

write(*,*) x0

END PROGRAM