MODULE RESIDUES
 
USE MINPACK
USE PHYCTE
USE MODCTE

IMPLICIT NONE 

CONTAINS 

SUBROUTINE HELLO_RES


WRITE(*,*) "HELLO RES"

END SUBROUTINE

subroutine Test_fcn ( n, x, res ,iflag )
! Test function to test solvers, constructed as stupulated in MINPACK
! Variables:  
! n : integer , IN : number of unknown and equations = 2 here
! x : double precision (n), IN/OUT : estimate of unknowns from previous iteration
! res : double precision (n), OUT : residue of equation system 
! iflag : integer IN/OUT : flag to communicate with solving subroutine
!-----------------
! Test function must solve eq system of type F[X] = 0 

    integer :: n
    real( kind = 8 ) :: res(n)
    integer :: iflag
    real ( kind = 8 ) ::  x(n)

    res(1) = -2 + x(1)**2 - 3*x(2)
    res(2) = 3*x(1)- 4*x(2)
end subroutine

subroutine Equation_system_ms (n, x, res ,iflag)
! Equation system based on Makalkin 1995. It aim to solve the midplane and 
! photosurface temperature profile along with surface temperature. 
! Each residue is named from the original paper equation number, 
! Ex: res_10 : is the residue of equation 10 in makalkin 1995.
! ----------
! Variables:
! ----------
! n : Integer, IN : 
!     Size of equation system (=8) 
!
! x : double precision, IN/OUT: 
!  vector of [sigma, T_mid, T_s, z_s, z_add, rho_mid, rho_add, rho_s]
!
! res : double precision, OUT :
!   vector of equation system residue: [res_10, res_17, res_23, res_24, res_31, res_36, res_37, res_39]
!
! iflag : integer IN/OUT :
!   flag to communicate with solving subroutine api.

    integer :: n , iflag                                     
    double precision ::  res(n) , x(n)

    ! Function variables 
    double precision :: sigma, T_mid, T_s, z_s, z_add, rho_mid, rho_add, rho_s

    ! residue
    double precision :: res_10, res_17, res_23, res_24, res_31, res_36, res_37, res_39



end subroutine 

END MODULE   