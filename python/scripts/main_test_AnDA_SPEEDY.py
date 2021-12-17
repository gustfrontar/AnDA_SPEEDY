#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  9 12:33:56 2018
@author: jruiz
"""
# analog data assimilation
import AnDA_SPEEDY_regional_cme_mod as AnDA
import os
import copy
import pickle as pkl
import numpy as np
import matplotlib.pyplot as plt

current_path = os.getcwd()  #Path were the script is being excecuted.

#=====================================================================================
# Define a configuration for the experiment
#=====================================================================================

configuration=dict()

configuration['NEns']=30                    #Number of ensemble members in AnDA
configuration['data_path']=current_path + '/../../DATA/SPEEDY_catalogs/'  #Location of the binary files storing the catalog.
configuration['output_path'] = current_path + '/../../DATA/pkl/'          #Location of output files
configuration['nx_max']=96                  #Maximum size of SPEEDY grid in x
configuration['ny_max']=48                  #Maximum size of SPEEDY grid in y
configuration['nz_max']=7                   #Maximum size of SPEEDY grid in z
configuration['NAnalogs']=250               #Number of analogs used in forecasting step.
configuration['alpha']=-1.2                 #Multiplicative inflation factor (negative values activates adaptive inflation) Miyoshi et al. 2011 (Mon. Wea. Rev.)
configuration['model']='SPEEDY'             #Name of the model
configuration['ctl_true']='1969010100.ctl'    #ctl file corresponding to the true catalog (to generate the true state and the noisy observations)
configuration['ctl_catalog']='1969010100.ctl' #ctl file corresponding to the catalog to be used in AnDA (catalog to be used in the analog forecasting step in AnDA)
configuration['center_point']=[0,0,0]       #Center of the box (grid index)
configuration['grid_span']=[1,1,1]          #How many grid points will the local domain span (in each direction). This is the local domain for the implementation of AnDA.
configuration['regression']='local_linear'  # chosen regression ('locally_constant', 'increment', 'local_linear'), see Lguensat (2017). 
configuration['sampling']='gaussian'        # chosen sampler ('gaussian', 'multinomial'), see Lguensat (2017).
configuration['method']='AnEnKF'            # chosen method ('AnEnKF', 'AnETKF', 'AnFor')
configuration['expname_true']='RHBL011'     # Number to be used to identify the catalog.
configuration['time_limits_catalog']=['1970010100','1994113018'] #Starting date and end date of the long-range simulation used as catalog (yyyymmddhh)
configuration['Q']=None                     #Model for model error (currently None)
configuration['std_data']=False             #Wheter the data will be standarized (i.e. remove seasonal cycle and standarize the variability) before runing AnDA.
configuration['dt_states'] = 4              #Assimilation frequency (in number of 6 hour periods, e.g. 4 -> 24 hours)
configuration['dt_obs']    = 1              #Frequency of available observations (in number of 6 hour periods, e.g. 4 -> 24 hours)
configuration['perturb_catalog_sigma']  = None  #Sigma of simulated model errors
configuration['perturb_catalog_offset'] = None  #Mean of simulated model errors
configuration['vert_lev']=3                 #Vertical levels in the local domain used to perform AnDA.
configuration['kernel'] = 'gaussian'        #Analog kernel ('tricube','gaussian')
configuration['random_seed_obs'] = 1        #Random seed for observation generation.
configuration['var_list'] = ['T']           #Variables to be used in the analog search and DA.
configuration['vert_lev'] = 3               #Vertical level at which the local AnDA domains will be centered.
configuration['expname_catalog']='RHBL011'  #The experiment will be repeated using the catalogs indicated by the elements of this list.
configuration['time_limits_true']=['1996010100','1996033100']  #Starting data and end date of the AnDA experiment (yyyymmddhh)
configuration['domain_name'] = 'Test_domain'

grid_point_start_we = 80    #Start and end grid point of the domain in which AnDA will be applied (West-East direction)
grid_point_end_we   = 81
grid_point_start_sn = 24    #Start and end grid point of the domain in which AnDA will be applied (South-North direction)
grid_point_end_sn   = 25
configuration['domain_limits'] = [grid_point_start_we,grid_point_end_we,grid_point_start_sn,grid_point_end_sn]
configuration['exp_name'] = configuration['time_limits_true'][0] + '_' + configuration['method'] 

var_list_fn = ''
for my_var in configuration['var_list'] :
   var_list_fn = var_list_fn + my_var

configuration['output_file_name'] = configuration['output_path'] + 'output_par.ens.' + str(configuration['NEns']) + 'cat.' + configuration['expname_catalog'] + 'true.' + configuration['expname_true'] + 'vars.' + var_list_fn  + 'lev.' + str(configuration['vert_lev']) + 'reg.' + configuration['domain_name'] + '.' + configuration['exp_name'] + '.pkl'

#=====================================================================================
#Run SPEEDY AnDA
#=====================================================================================

AnDA.run_sanda( configuration )

#=====================================================================================
# Make a simple plot of the time series of the true, the observations and the analysis
#=====================================================================================

start_time = 100
end_time   = 200
exp_data = pkl.load( open( configuration['output_file_name'] , "rb" ) )

plt.figure()

plt.plot( exp_data['xt'][0,0,0,start_time:end_time] , '-b' , label='True')
plt.plot( exp_data['xa'][0,0,0,start_time:end_time] , '-r' , label='Analysis')
plt.plot( exp_data['yo'][0,0,0,start_time:end_time] , 'ok' , label='Obs.')
plt.xlabel('Time (analysis cycles)')
plt.ylabel('Temperature (K)')
plt.title('AnDA Example Time series')
plt.legend()

plt.savefig('./Plot_test.png')

#plt.show()




   





