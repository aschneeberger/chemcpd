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

MODULE OPACITY
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Module computing disk opacity depending on temperature and density  !
! Is isolated in a particular module to allow opcity complexification !
! Without complicting with DSKPHY variables                           !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

USE MODCTE

implicit none

contains
subroutine opacity_table(N,Temp,beta,kappa_0,kappa_p)

    !IN/OUT
    integer , intent(in) :: N 
    double precision , intent(in) , dimension(N) :: Temp 
    double precision , intent(out) , dimension(N) :: beta 
    double precision , intent(out) , dimension(N) :: kappa_0 
    double precision , intent(out) , dimension(N) :: kappa_p

    !INTERNALS
    integer :: i

    forall (i = 1:N) kappa_0(i) = 0.1/10.0d0 
    forall (i = 1:N) beta(i) = 0.5 
    forall (i = 1:N) kappa_p(i) = 0.1 * p_Chi * Temp(i)**0.5 / 10.0d0

end subroutine

END MODULE 

!--------------------------------------------------------------------------------------------------------

module DSKPHY
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Module containing the disk physics computation, such as mommentum and !
! Energy flux,  mass conservation computation and radiative trnasfere   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
USE PHYCTE
USE MODCTE
USE QUADPACK
USE OPACITY

implicit none 

contains 

!---------------------------------------------------!
! Intrasec disk properties (computed only one time) !
!---------------------------------------------------!


function ang_mom_factor(N,r, R_c) 
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
! ang_mom_factor(N): double precision : angular momentum factor 

    !IN/OUT
    integer :: N
    double precision , dimension(N) :: r 
    double precision :: R_c 
    double precision , dimension(N) :: ang_mom_factor 

    ! The angular mommentum is devided in two sections, closer and further than R_c 

    ! Closer than R_c, the gas infall from the PSN contribute to the total CPD angular momentum
    where (r .lt. R_c) 
        ang_mom_factor =  1 - 1/5 * (r/R_c)**2 - sqrt(p_R_p/r) - 4/5 * sqrt(R_c/p_R_disk) + 4/5 * sqrt((p_R_p*R_c)/(p_R_disk * r)) &
                    & + 0.2 * sqrt(p_R_p/p_R_disk) * (r/R_c)**2 
    
    ! Further than R_c, the gas drift outward and there is no more infall 
    elsewhere  
        ang_mom_factor = 0.8 * sqrt(R_c/r) - sqrt(p_R_p/r) - 0.8*sqrt(R_c/p_R_disk) + sqrt(p_R_p/p_R_disk) 
   
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

    !IN/OUT
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

    !IN/OUT
    double precision  :: centrifugal_radius 
    double precision  :: J ! Special angular momentum 
    
    J = spe_ang_mom()
    centrifugal_radius = J**2.0 /(c_G * p_M_p)  

end function 

function kep_puls(N,r) 
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
    
    !IN/OUT
    integer ::  N
    double precision , dimension(N) :: r
    double precision , dimension(N) :: kep_puls 

    kep_puls = sqrt(c_G * p_M_p / r**3.0)

end function 



function flux_visc(N, r, omegak, cap_lambda) 
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

    !IN/OUT
    integer :: N 
    double precision , dimension(N) :: r
    double precision , dimension(N) :: omegak
    double precision , dimension(N) :: cap_lambda
    double precision , dimension(N) :: flux_visc 

    flux_visc = 3.0d0/(8.0d0*c_pi) * cap_lambda/p_L * p_M_dot * omegak**2

end function 

function flux_accretion(N, r, R_c)
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
! flux_accetion(N) : double precision : energy flux from gas accretion [W]

    !IN/OUT
    integer :: N
    double precision , dimension(N) :: r 
    double precision :: R_c 
    double precision , dimension(N) :: flux_accretion

    flux_accretion = ((p_Xd * p_Chi * c_G * p_M_p * p_M_dot)/(4.0d0 * c_pi * R_c * R_c * r)) * exp(-r**2/R_c**2)


end function 



!-----------------------------------------------------------------! 
! Properties depending on temperature used in the equation system !
!-----------------------------------------------------------------!



function flux_planet(r,N,z_s)
! Compute the energy from the planet the disk received, NOT THE ONE 
! THE DISK ABSORB ! In case of debuging code convergence, pay attention 
! to this paramater involving derivative of z_s that depend on T_s
!-------
!Inputs:
!-------
!
! N : integer : size of the grid
! r(N) : double precision : array of radii from disk center [m]
! z_s(N): double precision : array of the photosphere of the disk 
!         (def as the altitude where τ = 2/3  ) [m]
!-------
!Output:
!-------
! 
! flux_planet(N) : double precision : Planet energy flux on the disk. 

    !IN/OUT vars 
    integer :: N
    double precision , dimension(N) :: r  !radius 
    double precision , dimension(N) :: z_s  !altidude 
    double precision , dimension(N) :: flux_planet !flux 

    !Internal vars 
    double precision , dimension(N) :: eps , eta    !Disk curvature parameters 
    double precision , dimension(N) :: dzdr        ! z_s derivative 

    !Central difference derivative of z_s 
    dzdr(1) = 0.0d0
    dzdr(2:N-1) = (z_s(1:N-2) - z_s(3:N)) / (r(1:N-2) - r(3:N)) 
    dzdr(N) = dzdr(N-1)

    !Disk curvature parameters 
    eps = atan( (4.0d0/(3.0d0*c_pi) * p_R_p) / sqrt( r*r + z_s*z_s ))
    eta = atan(dzdr) - atan(z_s/r)

    !Planet energy flux on disk 
    flux_planet = p_L_p * sin(eps * r + eta) / (8.0d0*c_pi*(r*r + z_s*z_s))

end function


function gas_scale_height(N,Temp,omegak)
! Compute the gas scale height above the altitude where the  
! The temperature is taken.
!-------
!Inputs:
!-------
!
! N : integer : size of the grid
! Temp(N) : double precision : Temperature 
! omegak(N) : double precision : keplerian pulsation corresponding to r [s-1]
!
!-------
!Output:
!-------
!
! gas_scale_height(N) : double precision : gas scale height [m] 

    integer :: N 
    double precision , dimension(N) :: Temp
    double precision , dimension(N) :: omegak  
    double precision , dimension(N) :: gas_scale_height

    !internal vars 
    double precision , dimension(N) :: c_s ! sound speed 

    c_s = sqrt( Temp * c_gamma * c_Rg / p_mu_gas )

    gas_scale_height = c_s/omegak 

end function 

!----------------------------------! 
!  Equations of the system itself  !
!----------------------------------!
! The Equation system is derived from makalkin 1995,2014 and Heller 2015 
! They solve the surface and midplane temperature profile of the circum planetary disk 
! in the case of an addiabatique disk, with the photosurface being in an isothermal layer
! The photosphere is defined as the area where the optical depth τ = 2/3. 
! There are 8 equations to solve with 8 unknowns: 
!     Σ: the surface density 
!     ρm: mid-plane density 
!     ρs: photosphere density
!     ρa: addiabatique to isotherme frontier density
!     Tm: mid-plane temperature
!     Ts: photosphere temperature
!     zs: photosphere height from the mid-plane
!     za: addiabatique to isotherm frontier height from the mid-plane 
!
! The equations are equation 10,17,23,24,31-32,36,37 and 39 + Κs opacity table


function temp_surface(N,kappa_p,beta,sigma,f_vis,f_acc,f_planet)
! Compute the disk temperature at the photosphere altitude, defined to at 
! an optical depth τ=2/3. It is assumed that energy flux come from accretion 
! heating, viscous heating, planetary luminosity and sourounding PSN temperature. 
! Taken from Heller 2015 
!
!------
!Input:
!------
!
! N : integer : size of the grid
! kappa_p(N) : double precision : Mean planck opacity [m2.Kg.-1]
! beta(N) : double precision : Mean opacity has a form of κ0 * T^β, it is the exponent 
! sigma(N) : double precision : Gas surface density [Kg.m-2]
! f_vis(N) : double precision : Viscous heating power [W] 
! f_acc(N) : double precision : Accretion heating power [W] 
! f_planet(N) : double precision : Planet heating power [W] 
!
!-------
!Output:
!-------
!
! temp_surface(N) : double precision : Temperature at photosphere surface [K] 

    !IN/OUT
    integer :: N
    double precision , dimension(N) :: kappa_p , beta ! opacity terms 
    double precision , dimension(N) :: sigma ! surface density 
    double precision , dimension(N) :: f_vis, f_acc, f_planet ! energy powers 
    double precision , dimension(N) :: temp_surface  !surface temperature

    !Internal 
    double precision , dimension(N) :: prefactor ! prefactor in equation 

    prefactor = (1.0d0 + (2.0d0 * kappa_p * sigma)**(-1.0d0) ) / c_sigma_sb 

    temp_surface = (prefactor * (f_vis + f_acc + p_Ks* f_planet) + p_T_neb)**(0.25d0)

end function

function temp_mid_equation(N,kappa_p,beta,kappa_0,cap_lambda,q_s,omegak)
! Right hand side of equation 31 (or 32 ), midplane temperature can not be isolated, 
! This function is used in the solving subroutine to find the midplane temperature
! It is taken from makalkin 2014. It assume that the disk is addiabatique until 
! optical depth τ=2/3
! It compute T_m^(5-β) - T_s^(4-β) * T_m
!------
!Input:
!------
!
! N : integer : size of the grid
! kappa_p(N) : double precision : Mean planck opacity [m2.Kg.-1]
! beta(N) : double precision : Mean opacity has a form of κ0 * T^β, it is the exponent 
! kappa_0(N) : double precision : Mean opacity has a form of κ0 * T^β, it is the prefactor [m2.Kg-1]
! cap_lambda(N) : double precision : Angular mommentum factor Λ 
! q_s(N) : Mass coordinate of the disk's photosphere. 
! omegak(N) : double precision : keplerian pulsation corresponding to r [s-1]
!
!-------
!Output:
!-------
!
! temp_mid_equation(N) : double precision : T_m^(5-β) - T_s^(4-β) * T_m value

    !IN/OUT 
    integer :: N
    double precision , dimension(N) :: kappa_p , beta , kappa_0 ! mean opacity vars 
    double precision , dimension(N) :: cap_lambda !angular mommentum factor Λ
    double precision , dimension(N) :: q_s ! mass coordinate of photosphere 
    double precision , dimension(N) :: temp_mid_equation  ! midplane temp residual 
    double precision , dimension(N) :: omegak ! keplerian pulsation

    !Internal vars 
    double precision :: constants ! prefactor composed of all constants 
    constants = 3.0d0 / (512.0d0 * c_pi*c_pi) * p_mu_gas / (c_sigma_sb * c_Rg * c_gamma)  

    !T_m^(5-β) - T_s^(4-β) * T_m
    temp_mid_equation = constants * (4.0d0 - beta) * p_Chi * kappa_0 * (p_M_dot*p_M_dot)/p_alpha * omegak**3.0d0 &
                        & *(cap_lambda/p_L)**2.0d0 * q_s*q_s

end function 

function accretion_rate(N,sigma,T_m,omegak,cap_lambda) 
! Compute the accretion rate from the temperature, function used in the 
! résolution of the equation system, its value must match the parameter p_M_dot
!------
!Input:
!------
!
! N : integer : size of the grid
! sigma(N) : double precision : Gas surface density [Kg.m-2]
! T_m(N) : double precision : midplane temperature [K] 
! omegak(N) : double precision : keplerian pulsation corresponding to r [s-1]
! cap_lambda(N) : double precision : Angular mommentum factor Λ 
!
!-------
!Output:
!-------
!
! accretion_rate(N) : double precision : The accretion rate [Kg.s-1]  

    !IN/OUt 
    integer :: N 
    double precision , dimension(N) :: sigma
    double precision , dimension(N) :: T_m
    double precision , dimension(N) :: omegak 
    double precision , dimension(N) :: cap_lambda 
    double precision , dimension(N) :: accretion_rate

    !Internal vars 
    double precision :: constants ! prefactor with all constants in it 
    
    constants = (3.0 * c_pi * c_gamma * c_Rg * p_alpha * p_L) / p_mu_gas 

    accretion_rate = (constants * sigma * T_m) / (omegak * cap_lambda) 

end function

function rho_add_23(N,rho_s,z_s,z_add,T_s,omegak) 
! In the equation set there are two independent method to compute the 
! gas density at the altitude of transition between addiabatique and 
! isothermal part of the disk. The first is from equation 23 of makalkin 
! 1995, using the pressure in an isothermal fluid
!------
!Input:
!------
! N : integer : size of the grid
! rho_s(N) : double precision : Gas density at photosphere altitude [Kg.m-3] 
! z_s(N) : double precision : altitude of the photosphere [m] 
! z_add(N) : double precision : altitude of the addiabatique to isothermal transition [m]
! T_s(N) : double precision : Temperature at photosphere surface [K] 
! omegak(N) : double precision : keplerian pulsation corresponding to r [s-1]
!
!-------
!Output:
!-------
!
! rho_add_23(N) : double precision : gas density at addiabatique to isothermal transition
!                                   from eq 23 

    !IN/OUT
    integer :: N 
    double precision , dimension(N) :: rho_s 
    double precision , dimension(N) :: z_s  
    double precision , dimension(N) :: z_add  
    double precision , dimension(N) :: rho_add_23
    double precision , dimension(N) :: T_s
    double precision , dimension(N) :: omegak
    
    ! Internal vars 
    double precision , dimension(N) :: h_s 

    h_s = gas_scale_height(N,T_s,omegak)

    rho_add_23 = rho_s * exp( (c_gamma / 2.0d0) * ( z_s*z_s - z_add*z_add ) / (h_s*h_s) )

end function 

function rho_add_36(N,rho_mid,T_s,T_mid) 
! In the equation set there are two independent method to compute the 
! gas density at the altitude of transition between addiabatique and 
! isothermal part of the disk. The second is from equation 36 of makalkin 
! 1995, using the pressure in an addiabatique fluid
!------
!Input:
!------
! N : integer : size of the grid
! rho_mid(N) : double precision : midplane gas density [Kg.m-3] 
! T_s(N) : double precision : photosphere temperature [K] 
! T_mid(N) : double precision : midplane temperature 
!
!-------
!Output:
!-------
!
! rho_add_36(N) : double precision : gas density at addiabatique to isothermal transition
!                                    from eq 36

    !IN/OUT 
    integer :: N
    double precision , dimension(N) :: rho_mid
    double precision , dimension(N) :: T_s
    double precision , dimension(N) :: T_mid 
    double precision , dimension(N) :: rho_add_36 

    rho_add_36 = rho_mid * (T_s / T_mid) ** ( 1.0d0/(c_gamma -1.0d0 ) )

end function 

function addiabatique_height(N,T_mid,T_s,omegak) 
!Compute the altitude of the addiabatique to isothermal zone in the 
! disk. Taken from eq 37 in makalkin 1995. 
!------
!Input:
!------
!
! N : integer : size of the grid
! T_mid(N) : double precision : midplane temperature [K)
! T_s(N) : double precision : photosphere temperature [K] 
! omegak(N) : double precision : keplerian pulsation corresponding to r [s-1]
!
!-------
!Output:
!-------
!
! addiabatique_height(N) : double precision :  altitude of the addiabatique to isothermal zone [m]

    !IN/OUT
    integer :: N
    double precision , dimension(N) :: T_mid  
    double precision , dimension(N) :: T_s  
    double precision , dimension(N) :: omegak  
    double precision , dimension(N) :: addiabatique_height

    !internal vars
    double precision :: constants !  Prefactor made of all the constants 

    constants = sqrt( ( 2.0d0 * c_gamma) / ( c_gamma -1.0d0 ) * c_Rg / p_mu_gas  )

    addiabatique_height = constants * sqrt( T_mid - T_s ) / omegak

end function 

function optical_depth(N,rho_s,kappa_p,z_s,omegak,T_s) 
! Compute the optical depth at the photophere altitude, it must be 
! equal to 2/3, it is taken from eq 24 from makalkin 1995. It assume
! that we are in the isothermal part of the disk
!------
!Input:
!------
!
! N : integer : size of the grid
! rho_s(N) : double precision : Gas density at photosphere altitude [Kg.m-3] 
! kappa_p(N) : double precision : Mean planck opacity [m2.Kg.-1]
! z_s(N) : double precision : altitude of the photosphere [m] 
! omegak(N) : double precision : keplerian pulsation corresponding to r [s-1]
! T_s(N) : double precision : Temperature at photosphere surface [K] 
!
!-------
!Output:
!-------
!
! optical_depth(N) : double precision : Opticla depth [dimensionless]

    !IN/OUT 
    integer :: N 
    double precision , dimension(N) :: rho_s 
    double precision , dimension(N) :: kappa_p 
    double precision , dimension(N) :: z_s
    double precision , dimension(N) :: omegak 
    double precision , dimension(N) :: T_s 
    double precision , dimension(N) :: optical_depth 

    !INTERNAL VARS 
 
    double precision , dimension(N) :: b_s ! dimensionless integration factor  
    double precision , dimension(N) :: sqrt_b_s ! dimensionless integration factor sqrt (to optimize)  
    
    b_s = 0.5d0 * (p_mu_gas / c_Rg ) * (omegak *omegak * z_s *z_s) / T_s  
    sqrt_b_s = sqrt(b_s)

    optical_depth = kappa_p * rho_s *z_s * exp(b_s)/sqrt_b_s * ( 1 - erf(sqrt_b_s) ) ! equation 24 reduced   
end function 

function surface_density(N,z_add,rho_add,z_s,rho_s,T_s,omegak) 
! Compute the surface density of gas from density in the addiabatique zone 
! and density in isothermal zone. It is taken from equation 39 in makalkin
! 1995.
!------
!Input:
!------
!
! N : integer : size of the grid
! z_add(N) : double precision : altitude of the addiabatique to isothermal transition [m]
! rho_add(N) : double precision : gas density at addiabatique to isothermal transition [Kg.m-3] 
! z_s(N) : double precision : altitude of the photosphere [m] 
! rho_s(N) : double precision : Gas density at photosphere altitude [Kg.m-3] 
! T_s(N) : double precision : Temperature at photosphere surface [K] 
! omegak(N) : double precision : keplerian pulsation corresponding to r [s-1]
!
!-------
!Output:
!-------
!
! surface_density(N) : double precision : surface density of gas [Kg.m-2]

    !IN/OUT 
    integer :: N
    double precision , dimension(N) :: z_add 
    double precision , dimension(N) :: rho_add
    double precision , dimension(N) :: z_s 
    double precision , dimension(N) :: rho_s 
    double precision , dimension(N) :: T_s
    double precision , dimension(N) :: omegak
    double precision , dimension(N) :: surface_density

    !INTERNAL VARS 
    !Intermediates computations vars 
    double precision , dimension(N) :: b_a  ! in addiabatique section of integral
    double precision , dimension(N) :: b_s  ! in isothermal section of integral
    double precision , dimension(N) :: integration_add ! addiabatique section integral results array 
    double precision , dimension(N) :: integration_s ! isothermal section integral results array 
    double precision , dimension(N) :: zeta_a       ! isothermal integral lower bound (z_s/z_add)
    integer :: i

    !integration subroutine 
    !inputs 
    double precision :: a = 0.0d0, b=1.0d0 ! integration borders 
    double precision :: epsabs=1.0d-3 , epsrel=1.0e-8 ! absolute and relative deisred integration errors 
    integer :: key = 1 ! api communication with qag, 1 is set for 7 gauss points, 15 Gauss-Kronrod points,
    !outputs 
    double precision :: result ! integration results 
    double precision :: abserr ! integration  absolute error 
    integer :: neval ! number of evaluation of function to integrate 
    integer :: ier ! return code should be 0 for successful integration 
    
    b_a = 0.5 * p_mu_gas/c_Rg * (omegak*omegak * z_add*z_add) / (T_s*T_s)
    b_s = 0.5 * p_mu_gas/c_Rg * (omegak*omegak * z_s*z_s) / (T_s*T_s)

    ! Loop integrating the function addiabatique_density from 0 to 1 on all grid points 
    do i=1 , N 
        ! do integration 
        call qag(addiabatique_density,a,b,epsabs,epsrel,key,result,abserr,neval,ier,b_a(i))
        ! store result in array 
        integration_add(i) = result
    end do

    ! analitical integral results in isotherm
    zeta_a = z_add/z_s 
    !check if z_add > z_s to proceed, avoiding unphysical results 
    if (All(z_add < z_s)) then 
        write(*,*) "WARNING : Unphysical situation, z_addiabatique > z_surface"
    end if 

    integration_s = exp(b_s)/sqrt(b_s) * sqrt(c_pi)/2.0d0 * ( erf( 1/sqrt(b_s) ) - erf(zeta_a / sqrt(b_s)) )

    ! final surface density 
    surface_density = 2.0d0 * z_add * rho_add * integration_add + 2.0d0 * z_s * rho_s * integration_s
end function  

function addiabatique_density(zeta , b_a)
! Function to integrate to compute the gas mass in addiabatique zone
! depends on the parameter b_a that change with radius
!------
!Input:
!------
!
! zeta : double precision : integration variable, provided by quadpack   
! b_a : double precision :  constant parameter throughout the integration 
!
!-------
!Output:
!-------
!
! addiabatique_density : double precision 

    !IN/OUT
    double precision :: zeta 
    double precision :: b_a 
    double precision :: addiabatique_density 

    !Internal 
    double precision :: power ! power factor 
    double precision :: prefac ! prefactor 

    power = 1.0d0 / (c_gamma -1)
    prefac = (c_gamma -1.0d0) / c_gamma * b_a

    addiabatique_density = (1 + prefac*(1-zeta*zeta))**(power)

end function


!---------------------------------------------------------!
!       INITIALIZATION AND INITIAL GUESSES SUBROUTINE     ! 
!---------------------------------------------------------! 
! This part of the module is dedicated to initialize variables and find initial guesses to solve equations systems 
! Different guesses prescription are implemented : 
!   - Anderson et al 2020
!   - Makalkin et al 1995 analytical solution in pure addiabatique 
!   - Aguichine et al 2020 guesses adaptated to a circum planetary disk
! The user is free choose guesses that gives the better results  


subroutine Init_profiles(N,r,cap_lambda,R_c,omegak,F_vis,F_acc,T_mid,T_s,rho_mid,rho_add,rho_s,z_add,z_s,sigma)

    !IN/OUT
    integer , intent(in) :: N 
    double precision , intent(out) :: R_c
    double precision , intent(in) , dimension(N) :: r
    double precision , intent(out) , dimension(N) :: cap_lambda
    double precision , intent(out) , dimension(N) :: omegak
    double precision , intent(out) , dimension(N) :: F_vis
    double precision , intent(out) , dimension(N) :: F_acc 
    double precision , intent(out) , dimension(N) :: T_mid
    double precision , intent(out) , dimension(N) :: T_s
    double precision , intent(out) , dimension(N) :: rho_mid
    double precision , intent(out) , dimension(N) :: rho_add
    double precision , intent(out) , dimension(N) :: rho_s
    double precision , intent(out) , dimension(N) :: z_add
    double precision , intent(out) , dimension(N) :: z_s
    double precision , intent(out) , dimension(N) :: sigma

    !INTERNAL VARS 

    ! Parameters constant with temperature, initialised
    R_c = centrifugal_radius()

    cap_lambda = ang_mom_factor(N,r,R_c)

    omegak = kep_puls(N,r)

    F_vis = flux_visc(N,r,omegak,cap_lambda)

    F_acc = flux_accretion(N,r,R_c)

    ! Real first guesses 



end subroutine

subroutine Guesses_Anderson(N,r,cap_lambda,omegak,F_vis,F_acc,T_mid,T_s,rho_mid,rho_add,rho_s,z_add,z_s,sigma)
    
    !IN/OUT
    integer , intent(in) :: N 
    double precision , intent(in) , dimension(N) :: r
    double precision , intent(in) , dimension(N) :: cap_lambda
    double precision , intent(in) , dimension(N) :: omegak
    double precision , intent(in) , dimension(N) :: F_vis
    double precision , intent(in) , dimension(N) :: F_acc 
    double precision , intent(out) , dimension(N) :: T_mid
    double precision , intent(out) , dimension(N) :: T_s
    double precision , intent(out) , dimension(N) :: rho_mid
    double precision , intent(out) , dimension(N) :: rho_add
    double precision , intent(out) , dimension(N) :: rho_s
    double precision , intent(out) , dimension(N) :: z_add
    double precision , intent(out) , dimension(N) :: z_s
    double precision , intent(out) , dimension(N) :: sigma

    !INTERNAL VARS 

    double precision , dimension(N) :: c_s ! mid plane sound speed 
    double precision , dimension(N) :: h ! gas scale height 
    double precision , dimension(N) :: mu ! dynamical viscosity 
    double precision :: power ! factor in a power 

    !opacity subroutine 
    double precision , dimension(N) :: beta 
    double precision , dimension(N) :: kappa_0 
    double precision , dimension(N) :: kappa_p

    ! From equation 4 of anderson et al 2021 
    T_mid = ((3.0d0*omegak**2.0d0 * p_M_dot * cap_lambda) / (10.0d0 * c_PI * c_sigma_sb)  + p_T_neb**4.0d0)**0.25
    T_s = T_mid / 5.0d0

    c_s = sqrt(c_gamma * c_Rg * T_mid / p_mu_gas)
    h = c_s/omegak
    mu = p_alpha * c_s * c_s / omegak

    !From steady state radial equation 
    sigma = p_M_dot * cap_lambda /(3.0d0 * c_PI * mu) 

    !Approximating the density ρ(r,z) profile on z as a gaussian one
    rho_mid = sqrt( 2.0d0/c_PI) * sigma/(2.0d0*h)

    !One can find ρ_add from eq 36 in makalakin 1995
    power  = 1.0d0 / (c_gamma -1.0d0)
    rho_add = rho_mid * (T_s/T_mid) ** (power)

    !Approximation de z_s from the gas scale height from heller 2015
    call opacity_table(N, T_s, beta, kappa_0, kappa_p) ! use opacity at surface temperature

     

end subroutine 

end module

!--------------------------------------------------------------------------------------------------------

MODULE RESIDUES
! Module containing the residue computation F[X] = 0
! It contain test function and subroutines to compute the residue. 

USE MINPACK
USE DSKPHY

IMPLICIT NONE 

CONTAINS 

function fcn_int(x,param)
! Test function to test integration librairies 
! It gives a normalized gaussian, with inegration 
! over infinit domain must equal 1.0 
!
!------
!Input:
!------
!
! x: double precision 
! param : double precision 
!
!-------
!Output:
!-------
!
!fcn_int : double precision 


    double precision :: x,fcn_int 
    double precision :: param

    fcn_int = param * exp(-x*x) / sqrt(c_pi)

end function 

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