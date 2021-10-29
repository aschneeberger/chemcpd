"""
Proof of concept python script to solve the 9-equations system of the stationary makalkin cpd model 
The solving method is the MINIPACK HYBRID one to be coherent with the futur F̵̭͈́͠O̷̹̭͌͘RT̶̠͆̚Ṟ̶͎̋AN̵̖͑̿ version of 
this algorithm 
"""

from math import inf
import numpy as np
from numpy.core.numeric import full
from numpy.typing import _256Bit 
import scipy.optimize as opt 
import scipy.integrate as  int
import astropy.constants as cte
from scipy.optimize._lsq.least_squares import soft_l1
from scipy.special import erfinv
import matplotlib.pyplot as plt 
import scipy as sc

###################################################################################
#------------------------------- Constant table ----------------------------------#
###################################################################################

dict_cte = {} # Dictionary initialisation of all contants value in the temperature determination. 

M_earth = cte.M_earth.value
M_jup = cte.M_jup.value
L_sun = cte.L_sun.value
year = 365*24*3600          # year in seconds
sigma_sb = cte.sigma_sb
Xd = 0.006 # dust to gas mass ratio in the PSN at Jupiter's formation location from makalkin 2014
Chi = 1 # Enchiment factor in dust between the PSN and the cpd chi = Xd(cpd)/Xd(psn). Can be btw 1 and 1.5 (Ruskol 2006)
mu_gas = 2.341 # mean molecular weight 
Ks = 0.2 # CDP radiation absorption factor 
alpha = 1e-3 # Alpha tubulence from shakura and sunyaev
gamma = 1.42 # addiabatic compression facteor of Perfect Gas
temp_neb = 100 # PSN temperature at 5 AU [K]

Nr = 100 # Number of points in the grid 

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

"""
Values that will change as jupiter's mass ans accretion rate change but are independant of 
the tempererature and considered as constant in the résolution
"""

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
# ------- Diskparamers determination by recursive optimisation (MINI PACK HYBRID)----------------#
##################################################################################################

"""
The system to solve is derived fomr Makalkin et al 1995 and 2014, there are 8 equations to solbve 
with 8 unknowns: 
    Σ: the surface density 
    ρm: mid-plane density 
    ρs: photosphere density
    ρa: addiabatique to isotherme frontier density
    Tm: mid-plane temperature
    Ts: photosphere temperature
    zs: photosphere height from the mid-plane
    za: addiabatique to isotherm frontier height from the mid-plane 

The equations are equation 10,17,23,24,31-32,36,37 and 39 + Κs opacity table. 
To solve the problem we use an initial set of guess derived from the model of Anderson et al 2021
that do take into account the luminosity if the central planet. 
"""


#--------------------- Planck opcity table ---------------------------------------------------------#

def mean_planck_opacity(temp_surf,Chi) : 
    
    return  0.1/10 , 0.5 ,Chi * 0.1 * temp_surf**0.5/10

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


def surface_temperature(sigmag,F_vis,F_acc,F_p,ks,temp_neb,kappa_p):
      

      temp_surf = (( (1 + (2*kappa_p*sigmag)**(-1) )/cte.sigma_sb.value) * (F_acc + F_vis + ks*F_p) + temp_neb**4 )**0.25 

      return temp_surf


#------------------ Midplane temperature and residu résolution -----------------------------# 

def surface_mass_coordinate(kappa_p,sigmag) :
    
    return 1 - 4 / (3 * kappa_p * sigmag)


def midplane_temp_residue(mid_temp_guess,temp_surf,gamma,Chi,kappa0,beta,mug,alpha,qs,cap_lambda,omegak,L) :

    return mid_temp_guess**(5-beta) - temp_surf**(4-beta)*mid_temp_guess - 3/(512*np.pi**2) * mug/(cte.sigma_sb.value * 8.314 * gamma) * (4-beta) * Chi * kappa0 * Mdot**2/alpha * omegak**3 * (cap_lambda/L)**2 * qs**2

#-----------------------------Residue vector construction----------------------------------# 
"""
The method from scipy only take 1D vectors, ones have to parse the values on the grid from 
a 8 by Nr matrix to a 8*Nr vector. 
The résidue compute the vector function F[X]=R, where R is the residue vectore to minimize to 0 
Such that the method solve the equation F[X]=0.  
"""

def Parse(X,Nr):
    """
    Parser from the Vector form to the matrix form.

    input:
    X:variable vector of size 8*Nr

    output:
    dict_var: dict with variable names as keywords and values as output of size Nr 
    """
    dict_var = {}

    dict_var['sigma'] = X[:Nr]
    dict_var['T_mid'] = X[Nr:2*Nr]
    dict_var['T_s'] = X[2*Nr:3*Nr]
    dict_var['z_s'] = X[3*Nr:4*Nr]
    dict_var['z_add'] = X[4*Nr:5*Nr]
    dict_var['rho_mid'] = X[5*Nr:6*Nr]
    dict_var['rho_add'] = X[6*Nr:7*Nr] 
    dict_var['rho_s'] = X[7*Nr:]

    return dict_var


def Serialize(dict_var):
    """
    Convert the Matrix to the Vector form
    """

    return np.concatenate((dict_var['sigma'], dict_var['T_mid'], dict_var['T_s'], dict_var['z_s'], dict_var['z_add'], dict_var['rho_mid'], dict_var['rho_add'], dict_var['rho_s']))


def Residue(X,dict_cte,N) :
    """
    Compute the equation system vector residues F[X] 
    The equations are equation 10,17,23,24,31-32,36,37 and 39 + Κs opacity table from makalkin et al 1995
    """


    #Values from the previous iteration 
    dict_var = Parse(X,N)

    # Mean planck opacity at photosphere heigh 
    kappa_p0 , beta, kappa_p = mean_planck_opacity(dict_var['T_s'], dict_cte['Chi'])

    #Eq 10: Accretion rate 

    res_10 = dict_cte['Mdot'] - (3*np.pi * dict_cte['gamma'] * 8.314 * dict_cte['alpha'] * dict_var['sigma'] * dict_var['T_mid'] * dict_cte['L']) / (dict_cte['omegak'] * dict_cte['mu_gas'] * dict_cte['cap_lambda']) 

    #Eq 17: Surface temperature with improvemnt from heller 2015 

    # Planet luminosity energy flux
    flux_p = F_planet(dict_cte['L_p'],dict_cte['R_p'],dict_var['z_s'],dict_cte['r'])    

    res_17 = dict_var['T_s'] - surface_temperature(dict_var['sigma'],dict_cte['F_vis'],dict_cte['F_acc'],flux_p,dict_cte['Ks'],dict_cte['temp_neb'],kappa_p)

    #Eq 23:  Addiabatique altitude gas density 
    
    c_s = np.sqrt(dict_cte['gamma'] * 8.314 * dict_var['T_mid'] / dict_cte['mu_gas']) # mid plane sound speed 
    h = c_s / dict_cte['omegak'] # disk scale heghy

    res_23 = dict_var['rho_add'] - dict_var['rho_s'] * np.exp( dict_cte['gamma']/2.0 * (dict_var['z_s'] - dict_var['z_add']) / h )

    #Eq 24: At z=z_s the optical depth is 2/3

    # constant in integration
    b_s = 0.5 * dict_cte['mu_gas'] / 8.314  * dict_cte['omegak']**2 * dict_var['z_s']**2 / dict_var['T_s'] 

    # The int.quad do the integration of the lambda function from 1 to infinity
    res_24 = 2/3 -  kappa_p * dict_var['rho_s'] * dict_var['z_s'] * int.quad_vec(lambda eps : np.exp(b_s*(1-eps**2)), 1, np.inf)[0]

    #Eq 31 (or 32) The residue of the midplane temperature equation

    #surface mass coordinates
    q_s = surface_mass_coordinate(kappa_p,dict_var['sigma'])

    res_31 = midplane_temp_residue(dict_var['T_mid'],dict_var['T_s'],dict_cte['gamma'],dict_cte['Chi'],kappa_p0,beta,dict_cte['mu_gas'],dict_cte['alpha'],q_s,dict_cte['cap_lambda'],dict_cte['omegak'],dict_cte['L'])

    #Eq 36: Addiabatique density as function of mid plane density and temperature and photosphere
    # temperature (assumed to be equal to addiabatique frontier temperature)

    res_36 = dict_var['rho_add'] - dict_var['rho_mid']* (dict_var['T_s'] / dict_var['T_mid'])**(1/(dict_cte['gamma']-1))

    #Eq 37: Addiabatique to isotherm frontier height 

    res_37 = dict_var['z_add'] - np.sqrt( (2 * dict_cte['gamma'] * 8.314) / ( (dict_cte['gamma'] -1 ) * dict_cte['mu_gas'] ) ) * np.sqrt( dict_var['T_mid'] - dict_var['T_s'] ) / dict_cte['omegak']

    #Eq 39 we should be able to retreive the surface density from the density profile rho(r,z)

    # integration parameters
    b_a =  0.5 * dict_cte['mu_gas'] / 8.314  * dict_cte['omegak']**2 * dict_var['z_add']**2 / dict_var['T_s'] 
    zeta_a = dict_var['z_add']/dict_var['z_s']

    # function to integration to get 
    function_a = lambda zeta : (1 + (dict_cte['gamma'] -1)/dict_cte['gamma'] * b_a * (1 - zeta**2) ) ** (1/(dict_cte['gamma']-1))
    function_s = []
    
    # tournaround to integrate a vector with variable boundaries 
    for i in range(Nr) :

        function_s.append(int.quad_vec(lambda zeta : np.exp( b_s[i] * (1- zeta**2)),zeta_a[i], 1)[0])
    function_s = np.array(function_s)

    res_39 = dict_var['sigma'] - 2 * dict_var['z_add'] * dict_var['rho_add'] * int.quad_vec(function_a,0,1)[0] + 2 *dict_var['z_s'] * dict_var['rho_s'] * function_s

    # return all residues as single vector of size 8*Nr

    res_vec = np.concatenate( (res_10/dict_cte['Mdot'],res_17,res_23,res_24,res_31,res_36,res_37,res_39) )

    # plt.loglog(dict_cte['r'],dict_var['T_s'])
    # plt.loglog(dict_cte['r'],res_17)
    # plt.pause(0.1)
    # plt.cla()

    return res_vec



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



#-----------------Dictionary of values constant withh respect to temperature construction---------#

dict_cte['L'] = 1 - np.sqrt(R_p/R_disk)   # momentum transfert coeficient 


J = specific_angular_momentum(M_p)  # Specific angular momentum [m².s^-1]
R_c = J**2 / (cte.G.value * M_p) # Centrifuge radius [m] 

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
The solver need an initial guess close enougth to not diverge.
We use a prescription of temperature from Anderson et al 2021 and a gaussian gas density profile on z
"""

dict_var = {}

# From eq 4 in Anderson 2021 guess of surface and mid plane temps 
dict_var['T_mid'] = ((3*dict_cte['omegak']**2 * dict_cte['Mdot']*year * dict_cte['cap_lambda']) / (10 * np.pi * cte.sigma_sb.value)  + temp_neb**4)**0.25 
dict_var['T_s'] = dict_var['T_mid']/5.0

c_s = np.sqrt(dict_cte['gamma'] * 8.314 * dict_var['T_mid'] / dict_cte['mu_gas'])
h = c_s / dict_cte['omegak']

# approximation de z_s from the gas scale height from heller 2015
dict_var['z_s'] = 0.03 * h 

# We find sigma from the steady state radial equation
mu = dict_cte['alpha'] * c_s**2 / dict_cte["omegak"]
dict_var['sigma'] = dict_cte['Mdot'] * dict_cte['cap_lambda'] / (3*np.pi * mu) + 0.1

#Approximating the density ρ(,rz) profile on z as a gaussian one
dict_var['rho_mid'] = np.sqrt(2/np.pi) * dict_var['sigma'] / (2*h) +0.1


#One can find ρ_add from eq 36 in makalakin 1995
dict_var['rho_add'] = dict_var['rho_mid'] * (dict_var['T_s']/dict_var['T_mid'])**(1/(dict_cte['gamma']-1)) +0.1

#And z_add from eq 37 
dict_var['z_add'] = np.sqrt( (2 * dict_cte['gamma'] * 8.314) / ( (dict_cte['gamma'] -1 ) * dict_cte['mu_gas'] ) ) * np.sqrt( dict_var['T_mid'] - dict_var['T_s'] ) / dict_cte['omegak'] + 0.1

#We find z_s by solving equation 23 for z_s, all other variables are known
dict_var['rho_s'] = dict_var['rho_add'] * np.exp(- dict_cte['gamma']/2.0 * (dict_var['z_s']**2 - dict_var['z_add']**2)/h**2) +0.1

X0 = Serialize(dict_var)



sol = opt.least_squares(Residue,X0,args=(dict_cte,Nr),method='trf')

dict_sol = Parse(sol.x,Nr)

print(sol.message)


plt.plot(dict_cte['r']/cte.R_jup.value,dict_sol['T_mid'],label='midplane')
plt.plot(dict_cte['r']/cte.R_jup.value,dict_sol['T_s'],label='surface')
#plt.plot(dict_cte['r']/cte.R_jup.value,X0[Nr:],label='guess')
plt.legend()
plt.xscale('log')
plt.yscale('log')
plt.show()