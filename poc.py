"""
proof of concept script for the makalkin stationary model of cpd on (r,z)  
"""

import numpy as np 
import scipy.optimize as opt 
from scipy.integrate import solve_ivp
import astropy.constants as cte
from scipy.special import erfinv
###################################################################################
#------------------------------- Constant table ----------------------------------#
###################################################################################

dict_cte = {} # Dictionary initialisation of all contants value in the temperature determination. 

M_earth = cte.M_earth.value
M_jup = cte.M_jup.value
L_sun = cte.L_sun.value
year = 365*24*3600          # year in seconds
sigma_sb = cte.sigma_sb
Xd = 0.0142 # dust to gas mass ratio in the PSN at Jupiter's formation location from makalkin 2014
Chi = 1 # Enchiment factor in dust between the PSN and the cpd chi = Xd(cpd)/Xd(psn). Can be btw 1 and 1.5 (Ruskol 2006)
mu_gas = 2.34 # mean molecular weight 
Ks = 0.2 # CDP radiation absorption factor 

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

R_H = 5*cte.au.value * (M_p / (3* cte.M_sun.value))**(1/3)   # Hill radius from makalkin et al 2014 [m]

R_disk = R_H / 5                           # Disk radius from Mordasini et al 2013 [m]

########################################################################################
# -------------- Temperature independent quantities -----------------------------------# 
########################################################################################

def specific_angular_momentum(M_p,a_p = 5) :
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
    cap_lambda[r<Rc] = 1 - 1/5 * (r[r<Rc]/Rc)**2 - np.sqrt(Rp/r[r<Rc]) - 4/5 * np.sqrt(Rc/R_disk) + 4/5 * np.sqrt((Rp*Rc)/(R_disk * r[r<Rc])) + 0.2 * np.sqrt(Rp/R_disk) * (r[r<Rc]/Rc)**2 

    #At r>=rc, the angular momentum of gas infall is not considered since at those radii gas is ejected from the cpd and there is no infalls
    cap_lambda[r>=Rc] = 0.8 * np.sqrt(Rc/r[r>=Rc]) - np.sqrt(Rp/r[r>=Rc]) - 0.8*np.sqrt(Rc/R_disk) + np.sqrt(Rp/R_disk) 

    return cap_lambda     





##################################################################################################
# ------- Temperature determination along with the surface density by recursive optimisation-----#
##################################################################################################

"""
Since the temperature depends on the surface density and vise-versa, we need to solve this non 
linear set of equations via an iterative methode. The iterative least_square  method is taken from the lib scipy.optimize.
We initialize it with a guess for the midplane temperature profile Tm0.  
"""


#--------------------- Surface density ---------------------------------------------------------#

def sound_speed(gamma,mid_temp,mug):
    return np.sqrt(gamma * 8.314 * mid_temp / mug)

def viscosity(cs,alpha,omegak):
    return alpha * cs*cs / omegak


def gas_scale_height(cs,mid_temp,omgak):
    return cs*mid_temp/omgak

def surface_density(Mdot,mu,cap_lambda,L):
    return Mdot/(3 * np.pi * mu) * cap_lambda/L

#-----------------------Surface temperature computation----------------------------------------#

def F_vis(cap_lambda,Mdot,omegak,L):
    """
    Independent of temperature
    """
    return 3/(8*np.pi) * cap_lambda/L * Mdot * omegak**2

def F_acc(Xd,Chi,M_p,Mdot,Rc,r):
    """
    independent of temperature
    """
    return (Xd * Chi*cte.G.value * M_p * Mdot)/(4 * np.pi * Rc * Rc * r)

def F_planet(L_p, R_p , zs, r):
    
    eps = np.arctan(4/(3*np.pi) * R_p/np.sqrt(r*r + zs*zs))

    eta = np.arctan(np.gradient(zs,r,edge_order=2)) - np.arctan(zs/r)

    return L_p * np.sin(eps * r + eta)/(8 * np.pi *(r*r + zs*zs))

def z_surface(sigmag,kappa_P,h):

    return erfinv(1 -2/3 * 2/(sigmag * kappa_P)) * np.sqrt(2) * h

def mean_planck_opacity(temp_surf,Chi) : 
    
    return 0.1,0.5,Chi*0.1*temp_surf**0.5

def surface_temperature(sigmag,F_vis,F_acc,F_p,ks,temp_neb,kappa_p):
      

      temp_surf = ( ((1 + (2*kappa_p*sigmag)**(-1) )/cte.sigma_sb.value) * (F_acc + F_vis + ks*F_p) + temp_neb**4 )**0.25 

      return temp_surf


#------------------ Midplane temperature and residu résolution -----------------------------# 

def surface_mass_coordinate(kappa_p,sigmag) :
    
    return 1 - 4 / (3 * kappa_p * sigmag)


def midplane_temp_residue(mid_temp_guess,temp_surf,gamma,Chi,kappa0,beta,mug,alpha,qs,cap_lambda) :

    return mid_temp_guess**(5-beta) - temp_surf**(4-beta)*mid_temp_guess - 3/(512*np.pi**2) * mug/(cte.sigma_sb.value * 8.314 * gamma) * (4-beta) * Chi * kappa0 * Mdot**2/alpha * omegak**3 * (cap_lambda/L)**2 * qs**2



def Residue(X,dict_cte,N) :
    """
    Residue computation to determine bpth surface and midplane temperatures, inserted in the optimization
    algorithm as a single array [mid_temp,temp_surf].  
    """
    # Disk properties depending on temperature
    cs = sound_speed(dict_cte['gamma'],X[0:N-1],dict_cte["mu_gas"])
    mu = viscosity(cs,dict_cte['alpha'],dict_cte["omageak"])
    sigmag = surface_density(dict_cte['Mdot'],mu,dict_cte['cap_lambda'],dict_cte['L'])
    h = gas_scale_height(cs,X[0:N-1],dict_cte['omegak'])

    #Surface properties depending on surface temperature
    kappa0,beta,kappa_p = mean_planck_opacity(X[N:],dict_cte['Chi'])
    zs = z_surface(sigmag,kappa_p,h)

    #Heating from planet luminosity depends on surface géometry
    F_p  = F_planet(dict_cte['L_p'],dict_cte['R_p'],zs,dict_cte['r'])

    #Compute of surface temp and resulting residue
    temp_surface= surface_temperature(X[N:],sigmag,dict_cte['F_vis'],dict_cte['F_acc'],F_p,dict_cte['ks'],dict_cte['temp_neb'],dict_cte['Chi'])

    res_temp_surf = temp_surface - X[N:]

    #Compute of midplane temperatureand residue 

    qs = surface_mass_coordinate(kappa_p,sigmag)

    res_mid_temp = midplane_temp_residue(X[0:N-1],temp_surface,dict_cte["gamma"],dict_cte["Chi"],kappa0,beta,dict_cte["mu_gas"],dict_cte['alpha'],qs,dict_cte["cap_lambda"])

    return np.concatenate((res_mid_temp,res_temp_surf))



#######################################################################################
#---------Contant values construction and optimization of temperatures----------------#
#######################################################################################

#-----------------=Construction of the constant value wqith temperature dictionary---------#

dict_cte['L'] = 1 - np.sqrt(R_p/R_disk)   # momentum transfert coeficient 

J = specific_angular_momentum(M_p)  # Specific angular momentum [m².s^-1]
R_c = J**2 / (cte.G.value * M_p) # Centrifuge radius [m] 

dict_cte['r'] = np.logspace(np.log10(2*R_p),np.log10(R_disk),Nr) # computational log grid, non linear

dict_cte['cap_lambda'] = cap_lambda(dict_cte['r'],R_c,R_p,R_disk) # Lambda from equations 2 and 3 from makalkin 2014

dict_cte['omegak'] = np.sqrt(cte.G.value * M_p / dict_cte['r']**3) # keplerian pulsasion