"""
proof of concept script for the makalkin stationary model of cpd on (r,z)  
"""

import numpy as np 
import scipy.optimize as opt 
from scipy.integrate import solve_ivp
import astropy.constants as cte


#------------------------------- Constant table ----------------------------------#

M_earth = cte.M_earth.value
M_jup = cte.M_jup.value
L_sun = cte.L_sun.value
year = 365*24*3600          # year in seconds
sigma_sb = cte.sigma_sb
Chi = 0.0142 # dust to gas mass ratio from makalkin 
mu_gas = 2.34 # mean molecular weight 

Nr = 2000 # Number of points in the grid 

# precomputed model data of Mordasini 2013 taken from the graphs of Heller 2015

M_p_pic_accretion = 200 * M_earth # Jupiter mass at pic accretion runaway 
M_p_high = 300 * M_earth    # Later jupiter mass after gas accretion runaway 

M_p = M_p_high  

L_p_pic_accretion = 2e-3 * L_sun    # Jupiter luminosity at pic accretion runaway 
L_p_high = 1e-3 * L_sun     # later Jupiter luminosity after accreion runaway, at t=1Myr 

L_p = L_p_high 

Mdot_pic_accretion = 1e-5 * M_earth/year #Accretion rate at accretion runaway  
Mdot_high = 4e-6 * M_earth/year     # Later accretion rate after 1Myr

Mdot = Mdot_high 

R_p_pic_accretion = 10 * cte.R_jup.value # Jupiter's radius at pic accretion runaway 
R_p_high = 1.4 * cte.R_jup.value # Jupiter's radius after pic accretion at t = 1Myrs 

R_p = R_p_high

# disk parameters 

R_H = R_p * M_p / (3* cte.M_sun.value)   # Hill radius from makalkin et al 2014 [m]

R_disk = R_H / 5                           # Disk radius from Mordasini et al 2013 [m]


# -------------- Temperature and radius independent quantities -----------------------------------# 

def sepecific_angular_momentum(M_p,a_p = 5) :
    """
    Give the cpd specific angular momentun (heller 2015 eq 1)
    
    inputs :
    M_p : planet mass in kg 
    a_p  ; planet distance to the sun in au 

    return :
    J : specific angular momentum in m².s^{-1} 
    """

    if M_p < M_jup : 
        return 7.8e11 * (M_p/M_jup) * a_p**(7/4)
    else :
        return 9.0e11 * (M_p/M_jup)**(2/3) * a_p **(7/4)

def cap_lambda (r,Rc,Rp,R_disk) :
    """
    Compute the Value capital lambda from makalkin et al 2014 eq 2 and 3

    inputs:
    r : 1D array, computational grid [m]
    Rc : Centrifugal radius [m]
    Rp : planet radius [mp]

    return :
    lambda : Capital lambda factor, angular momentum factor  
    """

    cap_lambda = np.zeros(r.shape)
    
    # At r<rc, the angular momentum of gas infall from the psn is taken into acount
    cap_lambda[r<Rc] = 1 - 1/5 * (r[r<Rc]/Rc)**2 - np.sqrt(Rp/r[r<Rc]) - 4/5 * np.sqrt(Rc/R_disk) + 4/5 * (Rp*Rc)/(R_disk * r[r<Rc]) + 0.2 * np.sqrt(Rp/R_disk) * (r[r<Rc]/Rc)**2 

    #At r>=rc, the angular momentum of gas infall is not considered since at those radii gas is ejected from the cpd and there is no infalls
    cap_lambda[r>=Rc] = 0.8 * np.sqrt(Rc/r[r>=Rc]) - Rp/r[r>=Rc] - 0.8*Rc/R_disk + np.sqrt(Rp/R_disk) 

    return cap_lambda     


L = 1 - np.sqrt(R_p/R_disk)   # momentum transfert coeficient 
J = sepecific_angular_momentum(M_p)  # Specific angular momentum [m².s^-1]
R_c = J^2 / (cte.G * M_p) # Centrifuge radius [m] 

r = np.logspace(R_disk,R_c,Nr) # computational log grid, non linear

lmbd = cap_lambda(r,R_c,R_p,R_disk) # Lambda from equations 2 and 3 from makalkin 2014