"""
Lib to parse output files from the model into astropy table. 
@author : Antoine Schneeberger 
@mail : antoine.schneeberger@lam.fr
"""

import sys 
import matplotlib.pyplot as plt 
from astropy.table import Table
import os 
import astropy.constants as cte

DIRECTORY  = '../Archives/' + sys.argv[1]
FILENAME = sys.argv[2] 

def Parse(dir,fname) :
    """
    Read a data file in a table format en return it as an astropy table
    """
    data_table = Table.read(dir+'/'+fname, format = 'ascii')

    print('Columns found in datafile : ' , data_table.colnames)

    return data_table

def control_plot(data_table,key,xlog=False,ylog=False,r=None) :
    """
    Plot a control graph and save it in ag rpah directory as a subdirectory of 
    the data directory 
    """

    #Check if the asked column exist
    if not (key in data_table.colnames) :
        raise os.error("Asked Column not present in table, change key variable to an existing table entry ")


    # If no external r given, take the one in the table 
    if r is None:
        #Verify if there is a column r 
        if not('r' in data_table.colnames):
            raise os.error("There is no radius column in the table and radius array was not provided")
        
        r = data_table['r']
   
    # Do the plot 
    plt.plot(r/cte.R_jup.value,data_table[key])

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
        os.mkdir(DIRECTORY+'/Graph')
    if not (FILENAME[:-4] in os.listdir(DIRECTORY+'/Graph')):
        os.mkdir(DIRECTORY+'/Graph/' + FILENAME[:-4])
            

    # Save the figure in the dedicated directory 
    plt.savefig(DIRECTORY + '/Graph/'+ FILENAME[:-4]+'/control_'+ key +'.png', dpi=300)
    plt.close()


if __name__ == '__main__':

    table = Parse(DIRECTORY,FILENAME)
    #If there is r in the colnames of the file, use its radius 
    if 'r' in table.colnames :
        for key in table.colnames :
            control_plot(table,key,xlog=True)
    #For other files such as debbug files or res files, use the radius from initialisation
    else :
        r_init = Parse(DIRECTORY,'initialisation.dat')['r']
        for key in table.colnames :
            control_plot(table,key,xlog=True,r=r_init)
