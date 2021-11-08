
function fcn_int(x)

    double precision :: x,f 

    f = exp(-x*x)

end function 


PROGRAM CHEMCPD 

USE RESIDUES
USE MINPACK
USE MODCTE
USE DSKPHY
USE QUADPACK


integer , parameter :: n = 5
real(8) , dimension(2) ::  x0, fvec
real(8) :: tol 
integer :: info 
double precision , dimension(n) :: r
double precision :: R_c

double precision ::  epsabs = 1.0d-3,epsrel=1.0d-6
integer :: key = 1
double precision :: a = -20.0d0 ,b = +20.0d0

double precision :: results 
double precision ::  abserr
integer  :: neval,ier



call QAG(fcn_int,a,b,epsabs,epsrel,jey,results,abserr,neval,ier)


call hybrd1 (Test_fcn, 2, x0, fvec, tol, info)

write(*,*) 
write(*,*) results

END PROGRAM

