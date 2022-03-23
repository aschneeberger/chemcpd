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

import matplotlib.pyplot as plt 
import numpy as np 
import pandas as pd 
import sys 
import os 

FILE = sys.argv[1] 
THRESHOLD = 1e-2

data_table = pd.read_csv(FILE,sep=',')

mask = (np.abs(data_table['res_Tm']) < THRESHOLD) * (np.abs(data_table['res_Ts']) < THRESHOLD)

r = np.loadtxt('Data_test/initialisation.dat',skiprows=1)[:,0]

Tm = data_table['Tm'][mask]
Ts = data_table['Ts'][mask]
r = r[mask]

plt.plot(r,Tm,'+')
plt.plot(r,Ts,'+')
plt.show()