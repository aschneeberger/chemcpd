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

PATH = sys.argv[1]   # get file path 
FILE = sys.argv[2]
THRESHOLD = 1     # Residue threshold to consider that a point converged  


data_table = pd.read_csv(PATH+'/'+FILE,sep=',')                     # Open the data csv file 
init_table = pd.read_csv(PATH+'/initialisation.csv',sep=',')        # Open the init csv file for grid and constants disk param 
disk_table = pd.read_csv(PATH+'/disk_parameters.csv',sep=',')       # Open csv where all disk parameters are stored 
cte_table  = pd.read_csv(PATH+'/physical_constants.csv',sep=',')     # Open csv file with all used physics constants


# Create a mask where both surface and midplane temperature  residuals are under the threshold 
mask = (np.abs(data_table['res_Tm'].values) < THRESHOLD) * (np.abs(data_table['res_Ts'].values) < THRESHOLD)

# Get the radial coordinates 
r = init_table["r"].values

# fileter all quantities 
Tm = data_table['Tm'][mask].values
Ts = data_table['Ts'][mask].values
sigma = data_table["sigma"][mask].values
scale_height = data_table["scale_height"][mask].values
r = r[mask]

# Create the grid in z 
z = np.logspace(1,np.log10(cte_table["R_jup"]),1000)

# Create the 2D mesh 
R,Z = np.meshgrid(r,z)

rho = np.zeros((len(z),len(r)))


# Compute all 1D quantities 

rho0 = np.sqrt(2/np.pi) * sigma / (2*scale_height) 

for j in range(len(z)) :
    for i in range(len(r)):
        rho[j,i] = rho0[i] * np.exp(-z[j]**2 / (2*scale_height[i]))



plt.contourf(R,Z,np.log10(rho),cmap="hot")
plt.clim(-10,-2)
plt.colorbar()
plt.show()