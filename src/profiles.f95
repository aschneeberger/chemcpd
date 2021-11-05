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

function ang_mom_factor(r, R_c, N) 
! Function that compute the angular momentum factor Λ. 
! following equation 2 and 3 from Makalkin et al 2014. 
! -------
! Input :
! -------
! 
! N : Integer : grid size 
! r(N) : double precision : radius [m]
! R_c : double precision : Centrifugal radius [m]
!
! -------
! Output :
! -------
! 
! cap_lambda(N): double precision : angular momentum factor 

    integer :: N
    double precision , dimension(N) ::  cap_lambda 
    double precision , dimension(N) :: r 
    double precision :: R_c 
    double precision , dimension(N) :: ang_mom_factor 

    ! The angular mommentum is devided in two sections, closer and further than R_c 

    ! Closer than R_c, the gas infall from the PSN contribute to the total CPD angular momentum
    where (r .lt. R_c) 
        cap_lambda =  1 - 1/5 * (r/R_c)**2 - sqrt(p_R_p/r) - 4/5 * sqrt(R_c/p_R_disk) + 4/5 * sqrt((p_R_p*R_c)/(p_R_disk * r)) &
                    & + 0.2 * sqrt(p_R_p/p_R_disk) * (r/R_c)**2 
    
    ! Further than R_c, the gas drift outward and there is no more infall 
    elsewhere  
        cap_lambda = 0.8 * sqrt(R_c/r) - sqrt(p_R_p/r) - 0.8*sqrt(R_c/p_R_disk) + sqrt(p_R_p/p_R_disk) 
   
    end where 

end function 


function spe_ang_mom() 
! Give the CPD specific angular momentum from the planet mass and distance to the sun
! Formula tacken from Heller 2015 
!--------
! Inputs: 
!--------
!
! Variables are taken from the constant table, in the future they will evolve with time
!--------
! Output:
!--------
! 
! spe_ang_mom : double precision , specific angular momentum [m2.s-1]

    double precision ::  spe_ang_mom 

    if (p_M_p .lt. c_M_jup) then 
        spe_ang_mom = 7.8d11 * (p_M_p/c_M_jup) * (p_a_p/c_au)**(7/4)
    
    else 
        spe_ang_mom = 9.0d11 * (p_M_p/c_M_jup)**(2/3) * (p_a_p/c_au) **(7/4)
    end if 

end function

function centrifugal_radius() 
! Give the centrifugal radius from the constant table 
! Forumal from Makalkin 2014 
!-------
!Inputs:
!-------
!
! Variables are taken from the constant table, in the future they will evolve with time
!-------
!Output:
!-------
!
! centrifugal_radius : double precision : centrifugal radius (radius at witch Fgrav = Fcentrifuge) [m]  

    double precision  :: centrifugal_radius 
    double precision  :: J ! Special angular momentum 
    
    J = spe_ang_mom()
    centrifugal_radius = J**2.0 /(c_G * p_M_p)  

end function 

function kep_puls(r,N) 
! Keplerian pulsation at radius r 
!-------
!Inputs:
!-------
! N : integer : size of the grid over r 
! r(N) : double precision : array of radii from the disk center [m] 
!-------
!Output:
!-------
!
! kep_puls : double precision : Keplerian pulsation [s-1] 
    integer ::  N
    double precision , dimension(N) :: r
    double precision , dimension(N) :: kep_puls 

    kep_puls = sqrt(c_G * p_M_p / r**3.0)

end function 

function flux_visc(r, N, omegak, cap_lambda) 
!Compute the energy flux at photosphere altitute from viscous heating inside the disk
!-------
!Inputs: 
!-------
!
! N : integer : size of the grid
! r(N) : double precision : array of radii from disk center [m]
! omegak(N) : double precision : keplerian pulsation corresponding to r [s-1]
! cap_lambda(N) : double precision : angular momentum factor corresponding to r 
!
!--------
! Output:
!-------- 
!
! flux_visc(N) : energy flux from viscous heating corresponding to r [W]  

    integer :: N 
    double precision , dimension(N) :: r
    double precision , dimension(N) :: omegak
    double precision , dimension(N) :: cap_lambda
    double precision , dimension(N) :: flux_visc 

    flux_visc = 3.0d0/(8.0d0*c_pi) * cap_lambda/p_L * p_M_dot * omegak**2

end function 

function flux_accretion(r, N, R_c)
! Compute the energy flux from accretion of gas onto the CPD 
! It is cut off at r = R_c as for the angular mommentum factor
!-------
!Inputs:
!-------
!
! N : integer : size of the grid
! r(N) : double precision : array of radii from disk center [m]
! R_c : double precision : Centrifugal radius [m]
!
!-------
!Output: 
!-------
!
! flux_accetion(N) : double precision : energy flux from gas accretion (W)

    integer :: N
    double precision , dimension(N) :: r 
    double precision , dimension(N) :: R_c 
    double precision , dimension(N) :: flux_accretion

    flux_accretion = ((p_Xd * p_Chi * c_G * p_M_p * p_M_dot)/(4.0d0 * c_pi * R_c * R_c * r)) * exp(-r**2/R_c**2)


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


subroutine Equation_system_ms (N, x, res ,iflag)
! Equation system based on Makalkin 1995. It aim to solve the midplane and 
! photosurface temperature profile along with surface temperature. 
! Each residue is named from the original paper equation number, 
! Ex: res_10 : is the residue of equation 10 in makalkin 1995.
! ----------
! Variables:
! ----------
! N : Integer, IN : 
!     Size of equation system (=8) 
!
! x'N) : double precision, IN/OUT: 
!  vector of [sigma, T_mid, T_s, z_s, z_add, rho_mid, rho_add, rho_s]
!
! res(N) : double precision, OUT :
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