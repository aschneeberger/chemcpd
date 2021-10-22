"""
proof of concept script for the makalkin stationary model of cpd on (r,z)  
"""

from math import inf
import numpy as np
from numpy.typing import _256Bit 
import scipy.optimize as opt 
from scipy.integrate import solve_ivp
import astropy.constants as cte
from scipy.optimize._lsq.least_squares import soft_l1
from scipy.special import erfinv
import matplotlib.pyplot as plt 
import scipy as sc
from autograd import grad 

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
mu_gas = 2.341 # mean molecular weight 
Ks = 0.5 # CDP radiation absorption factor 
alpha = 1e-3 # Alpha tubulence from shakura and sunyaev
gamma = 1.42 # addiabatic compression facteor of Perfect Gas
temp_neb = 100 # PSN temperature at 5 AU [K]

Nr = 1000 # Number of points in the grid 

# precomputed model data of Mordasini 2013 taken from the graphs of Heller 2015

M_p_pic_accretion = 200 * M_earth # Jupiter mass at pic accretion runaway 
M_p_high = 300 * M_earth    # Later jupiter mass after gas accretion runaway 

M_p = M_p_high  

L_p_pic_accretion = 1e-3 * L_sun    # Jupiter luminosity at pic accretion runaway 
L_p_high = 1e-4 * L_sun     # later Jupiter luminosity after accreion runaway, at t=1Myr 

L_p =L_p_high 

Mdot_pic_accretion = 1e-5 * M_earth/year #Accretion rate at accretion runaway  
Mdot_high = 1e-8 * M_jup/year #4e-6 * M_earth/year     # Later accretion rate after 1Myr

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
    """
    Compute the sound speed in the disk midplane 

    inputs:
    gamma: coefficient de chaleur spécifique [no units]
    mid_temp: midplane temperature  [K]
    mug: mean molécular weight [Kg.mol-1]

    output:
    cs: sound speed [m.s-1]
    """
    return np.sqrt(gamma * 8.314 * mid_temp / mug )

def viscosity(cs,alpha,omegak):
    """
    Dynamic viscosity of the gas

    input:
    alpha: alpha viscosity parameter from shakura and sunyaev 
    cs: sound speed [m.s-1]
    omegak: keplerian pulsation [s-1]

    return:
    mu: the dysnamic viscosity [m2.s-1]
    """
    return alpha * cs*cs / omegak


def gas_scale_height(cs,mid_temp,omgak):
    """
    Compute the gas gaussian hydrastatic scale height 

    input:
    cs: the sound speed [m.s-1]
    omegak: keplerian pusation [s-1]


    """
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
    return ((Xd * Chi * cte.G.value * M_p * Mdot)/(4 * np.pi * Rc * Rc * r)) * np.exp(-r**2/Rc**2)

def F_planet(L_p, R_p , zs, r):
    
    eps = np.arctan((4/(3*np.pi) * R_p) / np.sqrt(r*r + zs*zs))

    dzsdr = np.gradient(zs,r,edge_order=2)
    #dzsdr[-1] = dzsdr[-2] # assume the same disk slope at the edge

    eta = np.arctan(dzsdr) - np.arctan(zs/r)

    F_p = L_p * np.sin(eps*r  + eta)/(8 * np.pi *(r*r + zs*zs))
    F_p[F_p<0.0] = 0.0

    return F_p

def z_surface(sigmag,kappa_P,h):

    # value inside erfinv function must be between -1 and 1
    # It imply that sigmag*kappa_p > 2/3
    zs= erfinv(1 -2/3 * 2/(sigmag * kappa_P)) * np.sqrt(2) * h 
    return  zs

def mean_planck_opacity(temp_surf,Chi) : 
    
    return 0.1/10 , 0.5 ,Chi * 0.1 * temp_surf**0.5/10

def surface_temperature(sigmag,F_vis,F_acc,F_p,ks,temp_neb,kappa_p):
      

      temp_surf = ( ((1 + 1/(2*kappa_p*sigmag) )/cte.sigma_sb.value) * (F_acc + F_vis + ks*F_p) + temp_neb**4 )**0.25 
      temp_surf[temp_surf<0.0] = 0.0

      return temp_surf


#------------------ Midplane temperature and residu résolution -----------------------------# 

def surface_mass_coordinate(kappa_p,sigmag) :
    
    return 1 - 4 / (3 * kappa_p * sigmag)


def midplane_temp_residue(mid_temp_guess,temp_surf,gamma,Chi,kappa0,beta,mug,alpha,qs,cap_lambda,omegak,L) :

    return mid_temp_guess**(5-beta) - temp_surf**(4-beta)*mid_temp_guess - 3/(512*np.pi**2) * mug/(cte.sigma_sb.value * 8.314 * gamma) * (4-beta) * Chi * kappa0 * Mdot**2/alpha * omegak**3 * (cap_lambda/L)**2 * qs**2



def Residue(X,dict_cte,N) :
    """
    Residue computation to determine bpth surface and midplane temperatures, inserted in the optimization
    algorithm as a single array [mid_temp,temp_surf].  
    """

    #Values from previous iteration 
    mid_temp_prev = X[:N]
    temp_surf_prev = X[N:2*N]
    sigmag_prev = X[2*N:3*N]
    kappa_p_prev = X[3*N:4*N]
    zs_prev = X[4*N:]

    #Compute new values (5 eq with 5 unknowns)

    # Disk properties depending on temperature
    cs = sound_speed(dict_cte['gamma'],mid_temp_prev,dict_cte["mu_gas"])
    mu = viscosity(cs,dict_cte['alpha'],dict_cte["omegak"])
    sigmag = surface_density(dict_cte['Mdot'],mu,dict_cte['cap_lambda'],dict_cte['L'])
    h = gas_scale_height(cs,mid_temp_prev,dict_cte['omegak'])

    #Surface properties depending on surface temperature
    kappa0,beta,kappa_p = mean_planck_opacity(temp_surf_prev,dict_cte['Chi'])
    zs = z_surface(sigmag_prev,kappa_p_prev,h)

    #Heating from planet luminosity depends on surface géometry
    F_p  = F_planet(dict_cte['L_p'],dict_cte['R_p'],zs_prev,dict_cte['r'])

    #Compute of surface temp and resulting residue
    temp_surface= surface_temperature(sigmag_prev,dict_cte['F_vis'],dict_cte['F_acc'],F_p,dict_cte['Ks'],dict_cte['temp_neb'],kappa_p)

    res_temp_surf = temp_surface - temp_surf_prev

    #Compute of midplane temperatureand residue 

    qs = surface_mass_coordinate(kappa_p,sigmag_prev)

    res_mid_temp = midplane_temp_residue(X[0:N],temp_surface,dict_cte["gamma"],dict_cte["Chi"],kappa0,beta,dict_cte["mu_gas"],dict_cte['alpha'],qs,dict_cte["cap_lambda"],dict_cte['omegak'],dict_cte["L"])

    return np.concatenate((res_mid_temp,res_temp_surf,sigmag-sigmag_prev,kappa_p-kappa_p_prev,zs-zs_prev))



#######################################################################################
#---------Contant values construction and optimization of temperatures----------------#
#######################################################################################

# ------------------Constants from the constant table------------------------------------#

dict_cte['Ks'] = Ks # rayonement absorption factor of cpd 
dict_cte["Chi"] = Chi 
dict_cte["mu_gas"] = mu_gas
dict_cte["alpha"] = alpha
dict_cte["gamma"] = gamma
dict_cte["temp_neb"] = temp_neb
dict_cte["L_p"] = L_p
dict_cte["R_p"] = R_p
dict_cte["Mdot"] = Mdot



#-----------------Construction of the constant value with temperature dictionary---------#

dict_cte['L'] = 1 - np.sqrt(R_p/R_disk)   # momentum transfert coeficient 


J = specific_angular_momentum(M_p)  # Specific angular momentum [m².s^-1]
R_c = J**2 / (cte.G.value * M_p) # Centrifuge radius [m] 

# dict_cte['r'] = np.logspace(np.log10(1.2*R_p),np.log10(100*cte.R_jup.value),Nr) # computational log grid, non linear
dict_cte['r'] = np.logspace(np.log10(1.2*R_p),np.log10(50*cte.R_jup.value),Nr) # computational log grid, non linear

dict_cte['cap_lambda'] = cap_lambda(dict_cte['r'],R_c,R_p,R_disk) # Lambda from equations 2 and 3 from makalkin 2014

dict_cte['omegak'] = np.sqrt(cte.G.value * M_p / dict_cte['r']**3) # keplerian pulsasion

#----------------energy flux independant of temperature------------------------------------#

# Viscosity energy flux to disk surface 
dict_cte['F_vis'] = F_vis(dict_cte['cap_lambda'],dict_cte['Mdot'],dict_cte['omegak'],dict_cte['L'])

#Accretion energy flux to the surface
dict_cte['F_acc'] = F_acc(Xd,Chi,M_p,Mdot,R_c,dict_cte['r'])

#---------------Setting the initial guess from a prescription of Anderson et al 2021------#
"""
The solver need an initial guess close enougth to not diverge. since the tempertaure dependence come the gas surface density
and the heating via planet luminosity. A simple prescription from anderson on an equilibrium between disk black body radiation 
and viscous heating is a close enougth guess.
"""

# From eq 4 in Anderson 2021 guess of surface and mid plane temps 
mid_temp0 = ((3*dict_cte['omegak']**2 * dict_cte['Mdot']*year * dict_cte['cap_lambda']) / (10 * np.pi * cte.sigma_sb.value)  + temp_neb**4)**0.25 

# guess for surface density 
cs = sound_speed(dict_cte['gamma'],mid_temp0,dict_cte["mu_gas"])
mu = viscosity(cs,dict_cte['alpha'],dict_cte["omegak"])
sigmag0 = surface_density(dict_cte['Mdot'],mu,dict_cte['cap_lambda'],dict_cte['L'])

# guess opacity 
kappa0,beta,kappa_p0 = mean_planck_opacity(mid_temp0,dict_cte['Chi'])

#Guess ofr photosphere surface altitude 
h = gas_scale_height(cs,mid_temp0,dict_cte['omegak'])
zs0 = z_surface(sigmag0,kappa_p0,h)


X0 = np.concatenate((mid_temp0/5,mid_temp0,sigmag0,kappa_p0,zs0))
sol = opt.least_squares(Residue,X0,args=(dict_cte,Nr),method='trf',jac='2-point')

print(sol)

plt.plot(dict_cte['r']/cte.R_jup.value,sol.x[:Nr],label='midplane')
plt.plot(dict_cte['r']/cte.R_jup.value,sol.x[Nr:2*Nr],label='surface')
#plt.plot(dict_cte['r']/cte.R_jup.value,X0[Nr:],label='guess')
plt.legend()
plt.xscale('log')
plt.yscale('log')
plt.show()