"""
Lib to parse output files from the model into astropy table. 
@author : Antoine Schneeberger 
@mail : antoine.schneeberger@lam.fr
"""

import numpy as np 
import matplotlib.pyplot as plt 
from astropy.table import Table
import os 
import astropy.constants as cte

DIRECTORY  = '../Data' 
FILENAME = 'initialisation.dat'

def Parse(dir,fname) :
    """
    Read a data file in a table format en return it as an astropy table
    """
    data_table = Table.read(dir+'/'+fname, format = 'ascii')

    print('Columns found in datafile : ' , data_table.colnames)

    return data_table

def control_plot(data_table,key,xlog=False,ylog=False) :
    """
    Plot a control graph and save it in ag rpah directory as a subdirectory of 
    the data directory 
    """

    #Check if the asked column exist
    if not (key in data_table.colnames) :
        raise os.error("Asked Column not present in table, change key variable to an existing table entry ")

    # Do the plot 
    plt.plot(data_table['r']/cte.R_jup.value,data_table[key])

    # change to log scale if needed 
    if xlog :
        plt.xscale('log')
    if ylog :
        plt.yscale('log')

    plt.xlabel('radius')
    plt.ylabel(key)

    # Verify is there is already a graph directory in data dir
    # if does not exist creat it 
    list_dir = os.listdir(DIRECTORY) 
    if not ('Graph' in list_dir) :
        os.mkdir(DIRECTORY+'/'+'Graph')

    # Save the figure in the dedicated directory 
    plt.savefig(DIRECTORY + '/' +'Graph/control_'+ key +'.png', dpi=300)
    plt.close()


if __name__ == '__main__':
    fnames = np.array(os.listdir(DIRECTORY),dtype='str')
    for elem in fnames :
        if not '.dat' in elem :
            fnames = np.delete(fnames,np.where(fnames==elem))
    print(fnames)
    table = Parse(DIRECTORY,FILENAME)
    for key in table.colnames :
        control_plot(table,key,xlog=True)
