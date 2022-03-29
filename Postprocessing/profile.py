#!/Users/aschneeberger/miniforge3/bin/python3.9

# @author : Antoine Schneeberger
# @mail : antoine.schneeberger@lam.fr

# Post processing script to compute the 2D temperature and  
# density profiles from the surface and midplane temperatures
# outputted by the model. 
# 
# USAGE
# -----
#
# python3 profile.py [PATH_TO_FILE] 
# 
# OUTPUT
# -----

import matplotlib.pyplot as plt  # Plot library 
import numpy as np               # Array manipulation  
import pandas as pd              # CSV manipulation 
import sys                       # Sys calls 
import os                        # File management
import astropy.constants as cte  # Astrophysics constants  
import scipy.special as spe 

PATH = sys.argv[1]   # get file path 
FILE = sys.argv[2]
THRESHOLD = 1e-2     # Residue threshold to consider that a point converged  


data_table = pd.read_csv(PATH+'/'+FILE,sep=',')                     # Open the data csv file 
init_table = pd.read_csv(PATH+'/initialisation.csv',sep=',')        # Open the init csv file for grid and constants disk param 
disk_table = pd.read_csv(PATH+'/disk_parameters.csv',sep=',')       # Open csv where all disk parameters are stored 
cte_table  = pd.read_csv(PATH+'/physical_constants.csv',sep=',')     # Open csv file with all used physics constants



# Create a mask where both surface and midplane temperature  residuals are superior to the threshold 
mask = (np.abs(data_table['res_Tm'].values) > THRESHOLD) * (np.abs(data_table['res_Ts'].values) > THRESHOLD)

# Get the radial coordinates 
r = init_table["r"].values

rmin =r.min() 
rmax =r.max() 


# fileter all quantities 
Tm = data_table['Tm'].values
Ts = data_table['Ts'].values
sigma = data_table["sigma"].values
scale_height = data_table["scale_height"].values


# Create the grid in z 
zmin=0
zmax = cte_table["R_jup"].values[0]

z = np.linspace(zmin,zmax,1000)


# Create the 2D mesh 
R,Z = np.meshgrid(r,z)

rho = np.zeros((len(z),len(r)))
q = np.zeros((len(z),len(r)))
T2D = np.zeros((len(z),len(r)))
kappa = disk_table['Chi'].values *   1.6e-5 * Ts**2.1

# Compute all 1D quantities 

rho0 = np.sqrt(2/np.pi) * sigma / (2*scale_height) 
for j in range(len(z)) :
    for i in range(len(r)):
        rho[j,i] = rho0[i] * np.exp(-z[j]**2 / (2*scale_height[i]**2)) 
        q[j,i] = 2*rho0[i]/sigma[i] * ( np.sqrt(np.pi/2) * scale_height[i] * spe.erf(z[j]/(np.sqrt(2) * scale_height[i]) ) )

        T2D[j,i] = (1 + 3/64 * 1.9 * init_table['F_vis'].values[i]/(cte.sigma_sb.value * Ts[i]**4) * kappa[i] * sigma[i] * (1-q[j,i]**2) ) ** (1/1.9) * Ts[i]

rho[:,mask] = 0
T2D[:,mask] = 0

xpos = np.linspace(0,1,5)
xval = np.around(np.geomspace(rmin/cte_table['R_jup'][0],rmax/cte_table['R_jup'][0],5),2)

plt.figure()
plt.title("gas density")
plt.imshow(rho[-1::-1,:],cmap="hot" , extent=(0,1,zmin/cte_table["R_jup"].values[0],zmax/cte_table["R_jup"].values[0]))
plt.xlabel('r in $R_{jup}$')
plt.ylabel('z in $R_{jup}$')
plt.xticks(xpos,xval)
#plt.clim(-10,-4)
plt.colorbar()
plt.savefig(PATH + '/Density.pdf',dpi=500)

plt.figure()
plt.title("disk temperature")
plt.imshow(T2D[-1::-1,:],cmap='hot',extent=(0,1,zmin/cte_table["R_jup"].values[0],zmax/cte_table["R_jup"].values[0]))
plt.colorbar()
plt.xlabel('r in $R_{jup}$')
plt.ylabel('z in $R_{jup}$')
plt.xticks(xpos,xval)
plt.savefig(PATH + '/Temperature.pdf',dpi=500)