module PHYCTE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Physical constant table, all constant 
! names are preceeded by c_ 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

implicit none
 


! From astropy.constants : 
double precision , parameter :: c_M_earth = 5.972d+24   ! earth mass [Kg]
double precision , parameter :: c_M_jup = 1.89813d+27   ! Jupiter mass [Kg] 
double precision , parameter :: c_L_sun = 3.828e+26     ! Sun Limunosity [W]
double precision , parameter :: c_year = 365.0*24.0*3600.0 ! Year [s] 
double precision , parameter :: c_sigma_sb = 5.6703744191844314d-08 ! Steffan-Boltzman constant 

!Other sources 
double precision, parameter :: c_gamma = 1.45  ! Addiabatique compression factor for perfect gas

end module 

module MODCTE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Model input parameters, all parameters
! names are preceeded by p_ 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

implicit none 
 

!---------------------
! Physical parameters 
!---------------------

! Dust to gas mass ratio in the PSN at Jupiter's formation location from Heller 2015
double precision , parameter :: p_Xd = 0.0006 

! Enchiment factor in dust between the PSN and the cpd chi = Xd(cpd)/Xd(psn). 
! Can be btw 1 and 1.5 (Ruskol 2006)
double precision , parameter :: p_Chi = 1.0 

!mean gas molecular weight from aguichine 2020 [Kg.mol-1] 
double precision , parameter :: p_mu_gas = 2.341d-3 

! CPD radiation absorption factor 
double precision, parameter :: p_Ks = 0.2 

! Alpha parameter for α-turbulence from Shakura & Sunyaev 1994  
double precision, parameter :: p_alpha = 1d-3

! PSN temperature in the viscinity of Jupiter forming area [K]
double precision, parameter :: p_T_neb = 100.0 

!----------------------
! Numerical parameters 
!----------------------

integer, parameter :: p_Nr = 1000 ! number of point in the r grid 


end module 
