! The system to solve is derived from Makalkin et al 1995 and 2014, there are 8 equations to solbve 
! with 8 unknowns: 
!     Σ: the surface density 
!     ρm: mid-plane density 
!     ρs: photosphere density
!     ρa: addiabatique to isotherme frontier density
!     Tm: mid-plane temperature
!     Ts: photosphere temperature
!     zs: photosphere height from the mid-plane
!     za: addiabatique to isotherm frontier height from the mid-plane 

! The equations are equation 10,17,23,24,31-32,36,37 and 39 + Κs opacity table. 
! To solve the problem we use an initial set of guess derived from the model of Anderson et al 2021
! that do take into account the luminosity if the central planet.

module DSKPHY
! Module containing the disqk physics computation, such as mommentum and 
! Energy flux,  mass conservation computation and radiative trnasfere 

USE PHYCTE
USE MODCTE

implicit none 

contains 

function cap_lambda(r, R_c) 
! Function that compute the angular momentum factor Λ. 
! following equation 2 and 3 from Makalkin et al 2014. 
! -------
! Input :
! -------
! 
! r : double precision : radius [m]
! R_c : double precision : Centrifugal radius [m]
!
! -------
! Output :
! -------
! 
! cap_lambda: double precision : angular momentum factor 

    double precision cap_lambda 
    double precision r 
    double precision R_c 

    ! The angular mommentum is devided in two sections, closer and further than R_c 

    ! Closer than R_c, the gas infall from the PSN contribute to the total CPD angular momentum
    if (r .lt. R_c) then 
        cap_lambda =  1 - 1/5 * (r/R_c)**2 - sqrt(p_R_p/r) - 4/5 * sqrt(R_c/p_R_disk) + 4/5 * sqrt((p_R_p*R_c)/(p_R_disk * r)) &
                    & + 0.2 * sqrt(p_R_p/p_R_disk) * (r/R_c)**2 
    
    ! Further than R_c, the gas drift outward and there is no more infall 
    else 
        cap_lambda = 0.8 * sqrt(R_c/r) - sqrt(p_R_p/r) - 0.8*sqrt(R_c/p_R_disk) + sqrt(p_R_p/p_R_disk) 
    end if 

end function 


function spe_ang_mom() 
! Give the CPD specific angular momentum from the planet mass and distance to the sun
! Formula tacken from Heller 2015 
!--------
! Inputs: 
!--------
!
! variables are taken from the constant table, in the future they will evolve with time
!--------
! Output:
!--------
! 
! spe_ang_mom : double precision , specific angular momentum [m2.s-1]

    double precision spe_ang_mom 

    if (p_M_p .lt. c_M_jup) then 
        spe_ang_mom = 7.8d11 * (p_M_p/c_M_jup) * (p_a_p/c_au)**(7/4)
    
    else 
        spe_ang_mom = 9.0d11 * (p_M_p/c_M_jup)**(2/3) * (p_a_p/c_au) **(7/4)
    end if 

end function


end module

!--------------------------------------------------------------------------------------------------------

MODULE RESIDUES
! Module containing the residue computation F[X] = 0
! It contain test function and subroutines to compute the residue. 

USE MINPACK
USE DSKPHY

IMPLICIT NONE 

CONTAINS 

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

end module 