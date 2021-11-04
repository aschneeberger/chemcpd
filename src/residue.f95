MODULE RESIDUE 
use minipack
IMPLICIT NONE 
CONTAINS 

SUBROUTINE HELLO_RES


WRITE(*,*) "HELLO RES"

END SUBROUTINE

subroutine Fcn ( n, x, res ,iflag )
  
    integer :: n
    real( kind = 8 ) :: res(n)
    integer :: iflag
    real ( kind = 8 ) ::  x(n)

    res(1) = x(1)**2 + x(1)*x(2)**2 +1
    res(2) = x(1) + x(2)

end subroutine

END MODULE   