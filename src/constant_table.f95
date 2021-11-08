module PHYCTE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Physical constant table, all constant 
! names are preceeded by c_ 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

implicit none
 


! From astropy.constants : 
double precision , parameter :: c_M_earth = 5.972d+24   ! Earth mass [Kg]

double precision , parameter :: c_M_jup = 1.89813d+27   ! Jupiter mass [Kg] 
double precision , parameter :: c_R_jup = 71492000.0    ! Jupiter radius [Kg]

double precision , parameter :: c_M_sun = 1.98841d+30   ! Sun mass [Kg]
double precision , parameter :: c_L_sun = 3.828d+26     ! Sun Limunosity [W]

double precision , parameter :: c_year = 365.0*24.0*3600.0 ! Year [s] 
double precision , parameter :: c_sigma_sb = 5.6703744191844314d-08 ! Steffan-Boltzman constant 
double precision , parameter :: c_au = 149597870700.0 ! Astronomical unit [m]
double precision , parameter :: c_G = 6.6743d-11    ! Gravitational constant [m3.Kg-1.s-1] 

!Other sources 
double precision, parameter :: c_gamma = 1.45  ! Addiabatique compression factor for perfect gas
double precision , parameter :: c_pi = 4.0d0 * atan(1.0d0)

end module 

module MODCTE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Model input parameters, all parameters
! names are preceeded by p_ 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

USE PHYCTE
implicit none 
 

!---------------------
! Physical parameters 
!---------------------

! Dust to gas mass ratio in the PSN at Jupiter's formation location from Heller 2015
double precision , parameter :: p_Xd = 0.006 

! Enchiment factor in dust between the PSN and the cpd chi = Xd(cpd)/Xd(psn). 
! Can be btw 1 and 1.5 (Ruskol 2006)
double precision , parameter :: p_Chi = 1.0 

!mean gas molecular weight from aguichine 2020 [Kg.mol-1] 
double precision , parameter :: p_mu_gas = 2.341d-3 

! CPD radiation absorption factor 
double precision, parameter :: p_Ks = 0.2 

! Alpha parameter for α-turbulence from Shakura & Sunyaev 1994  
double precision, parameter :: p_alpha = 1.0d-3

! PSN temperature in the viscinity of Jupiter forming area [K]
double precision, parameter :: p_T_neb = 100.0 

! Distance of Jupiter from the sun [m]
double precision , parameter :: p_a_p = 5.0 * c_au
!----------------------
! Numerical parameters 
!----------------------

integer, parameter :: p_Nr = 1000 ! number of point in the r grid 


!--------------------
!   Pre-computation 
!--------------------

double precision , parameter :: p_M_p = 300.0 * c_M_earth             ! planet mass  [Kg] 
double precision , parameter :: p_L_p = 1.0d-4 * c_L_sun              ! planet luminosity [W] 
double precision , parameter :: p_R_p = 1.4 * c_R_jup               ! planet radius [m]
double precision , parameter :: p_M_dot = 4.0d-6 * c_M_earth/c_year   ! Accretion rate [Kg.s-1] 
double precision , parameter :: p_R_hill = 5.0 * c_au * ( p_M_p / (3.0 * c_M_sun) )**(1.0/3.0)  ! Hill radius [m]
double precision , parameter :: p_R_disk = p_R_hill / 5.0           ! Disk size [m]
double precision , parameter :: p_L = sqrt(p_R_p/p_R_disk)          ! angular momentum transfert coeficient 

end module 
