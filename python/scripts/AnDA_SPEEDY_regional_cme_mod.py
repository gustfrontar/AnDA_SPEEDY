#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  9 12:33:56 2018

@author: jruiz
"""
import sys
sys.path.append('../common_functions/')
sys.path.append('../common_modules/')
sys.path.append('../AnDA/')
sys.path.append('../')
                

import numpy as np
import matplotlib.pyplot as plt
# analog data assimilation
from AnDA_codes.AnDA_analog_forecasting import AnDA_analog_forecasting
from AnDA_codes.AnDA_data_assimilation  import AnDA_data_assimilation
from AnDA_codes.AnDA_generate_data_SPEEDY import AnDA_generate_data_SPEEDY_grid_based
from scipy.signal import savgol_filter
import pickle
#from tqdm import tqdm
import gc


def run_sanda(  configuration  )  :

   output = dict()
   output['configuration'] = configuration   #Store experiment configuration into output

   domain_limits       = configuration['domain_limits'] + [ configuration['vert_lev'] , configuration['vert_lev'] ]    #Add vertical levels to domain description.

   var_list_fn = ''
   for my_var in configuration['var_list'] :
      var_list_fn = var_list_fn + my_var

   domain_nx = domain_limits[1] - domain_limits[0] + 1
   domain_ny = domain_limits[3] - domain_limits[2] + 1
   domain_nz = domain_limits[5] - domain_limits[4] + 1

   #output_path= configuration['data_path'] + '/pkl/'

   if configuration['std_data'] :
       tmp_str='_STD'
   else        :
       tmp_str=''

   output_file_name = configuration['output_file_name']    
   #output_file_name = output_path + 'output_par.ens.' + str(configuration['NEns']) + 'cat.' + configuration['expname_catalog'] + 'true.' + configuration['expname_true'] + 'vars.' + var_list_fn  + tmp_str + 'lev.' + str(configuration['vert_lev']) + 'reg.' + configuration['domain_name'] + '.' + configuration['exp_name'] + '.pkl'


   ### GET DATA FROM THE SPEEDY MODEL

   class GD:
    model = configuration['model']
    mainpath_true=configuration['data_path']
    expname_true=configuration['expname_true']
    mainpath_catalog=configuration['data_path']
    expname_catalog=configuration['expname_catalog']
    ctl_true=configuration['ctl_true']
    ctl_catalog=configuration['ctl_catalog']
    center_point=configuration['center_point']   #Center of the box (grid index)
    grid_span=configuration['grid_span']      #How many grid points will the box span (in each direction)
    
    #grid_limits=np.array([315.0,325.0,-5.0,5.0,0.510,0.510])
    time_limits_catalog= configuration['time_limits_catalog'] #Time range for catalog.
    time_limits_true   = configuration['time_limits_true']    #Time range for true simulation
    var_list=configuration['var_list']               #Variable list for analog identification.
    nvar = 0                                         # Number of model variables (This will be computed by AnDA_generate_data_SPEEDY)      
    dt_integration = 3600*6                          # delta T in between outputs
    dt_states = configuration['dt_states']           # number of integeration times between consecutive states (for xt and catalog)
    dt_obs = configuration['dt_obs']                 # number of integration times between consecutive observations (for yo)
    nt_true = 0                                      # Number of times in the true integration
    nt_catalog = 0                                   # Number of times in the catalog
    obs_density = 1.0                                # Obs density in percentaje (1 means full observation)
    obs_locations = 0                                # Indices to obs location (to be computed by AnDA_generate_data_SPEEDY)
    sigma2_state = 0.0                               # variance of the model error to generate the catalog
    sigma2_obs = 0.5                                 # variance of the observation error to generate observation
    
    random_seed_obs = configuration['random_seed_obs']        #Random seed that will be combined with the grid point index for the generation of the
                                                              #noise in the observations. Change this seed to generate a new realization of the observation error.

    perturb_catalog_sigma = configuration['perturb_catalog_sigma']
    perturb_catalog_offset= configuration['perturb_catalog_offset']

    #Data transformation / normalization.
    standarize_data = configuration['std_data']    #Standarize data tacking into account the seasonal cycle.
    standarization_window = 30.0  #Number of days (+/- around a center day to define the time dependent seasonal mean)

   # ANALOG DATA ASSIMILATION (dynamical model given by the catalog)

   # parameters of the analog forecasting method
   class AF:
    k = configuration['NAnalogs']      # number of analogs
    neighborhood = 0  # global analogs
    catalog = 0       # catalog with analogs and successors
    regression = configuration['regression'] # chosen regression ('locally_constant', 'increment', 'local_linear')
    sampling = configuration['sampling'] # chosen sampler ('gaussian', 'multinomial')
    cov_est=True
    kdt=None
    initialized=False
    global_linear=None
    kernel = configuration['kernel']
    
   class DA:
    method = configuration['method'] # chosen method ('AnEnKF', 'AnEnKS', 'AnPF', 'An2PF_T')
    N = configuration['NEns']  # number of members (AnEnKF/AnEnKS) or particles (AnPF)
    Ns = 10 # number of smoothing trajectories AnPS
    xb = 0
    B = 0
    Xcond = []
    H = 0
    R = 0
    alpha = configuration['alpha']   #Multiplicative inflation factor.
    minalpha = 1.05                   #Minimum allowed multiplicative inflation in online inflation estimation
    Q     = configuration['Q']
    dt_states = configuration['dt_states'] 

    @staticmethod
    def m(x):
        return AnDA_analog_forecasting(x,AF)

    
   ii=0
   for igrid_iter in  range( domain_limits[0] , domain_limits[1] + 1)   :
       #Apply cyclic boundary conditions in X 
       igrid = igrid_iter
       if( igrid_iter < 0 )  :
           igrid = configuration['nx_max'] + igrid_iter
       if( igrid_iter > configuration['nx_max'] ) :
           igrid = igrid_iter - nx_max 

       jj = 0 
       for jgrid in range( domain_limits[2] , domain_limits[3] + 1)   :
           kk = 0 
           for kgrid in range( domain_limits[4] , domain_limits[5] + 1)   :
       
              #
              GD.center_point = [ igrid , jgrid , kgrid ]
          
              #Observational error for a particular grid point will be the same upon
              #repeated calls to this function.
              AF.catalog, xt, yo , GD , xgrid , ygrid , zgrid , R = AnDA_generate_data_SPEEDY_grid_based(GD) 
              AF.neighborhood = np.ones([xt.values.shape[1],xt.values.shape[1]])
              DA.xb = xt.values[0,:]
              if configuration['method'] == 'An3DVAR' :
                 DA.B = configuration['B3DVAR']
              else                                    :
                 DA.B = 0.1*np.eye(xt.values.shape[1])
              DA.H = np.eye(xt.values.shape[1])
              DA.R = R
              #DA.R = GD.sigma2_obs*np.eye(xt.values.shape[1])

              ## run the analog data assimilation
              x_hat_EnKF = AnDA_data_assimilation(yo, DA , xtrue=xt )
              AF.initialized = False
        
              ##Find the column corresponding to the domain center (this may change from box to box)
              for ivar in range( np.size( xgrid ) )  :
                  if( xgrid[ivar] == GD.center_point[0] and ygrid[ivar] == GD.center_point[1] and zgrid[ivar] == GD.center_point[2] ) :
                      center_domain_var_ind = ivar 
                 
              if kk == 0 and jj == 0 and ii == 0  :
                  output['loglike'] = np.zeros( (domain_nx , domain_ny , domain_nz , GD.nt_true) ) 
                  output['rmsea']    = np.zeros( (domain_nx , domain_ny , domain_nz , GD.nt_true) )
                  output['biasa']    = np.zeros( (domain_nx , domain_ny , domain_nz , GD.nt_true) )
                  output['rmsef']    = np.zeros( (domain_nx , domain_ny , domain_nz , GD.nt_true) )
                  output['biasf']    = np.zeros( (domain_nx , domain_ny , domain_nz , GD.nt_true) )
                  output['rmseao']    = np.zeros( (domain_nx , domain_ny , domain_nz , GD.nt_true) )
                  output['biasao']    = np.zeros( (domain_nx , domain_ny , domain_nz , GD.nt_true) )
                  output['rmsefo']    = np.zeros( (domain_nx , domain_ny , domain_nz , GD.nt_true) )
                  output['biasfo']    = np.zeros( (domain_nx , domain_ny , domain_nz , GD.nt_true) )

                  output['xt']    = np.zeros( (domain_nx , domain_ny , domain_nz , GD.nt_true) )
                  output['xa']    = np.zeros( (domain_nx , domain_ny , domain_nz , GD.nt_true) )
                  output['xf']    = np.zeros( (domain_nx , domain_ny , domain_nz , GD.nt_true) )
                  output['yo']    = np.zeros( (domain_nx , domain_ny , domain_nz , GD.nt_true) )
                  output['xc']    = np.zeros( (domain_nx , domain_ny , domain_nz , AF.catalog.analogs.shape[0] ) )
                  output['spread_f'] = np.zeros( (domain_nx , domain_ny , domain_nz , GD.nt_true) )
                  output['spread_a'] = np.zeros( (domain_nx , domain_ny , domain_nz , GD.nt_true) )
                  output['multinf']  = np.zeros( (domain_nx , domain_ny , domain_nz , GD.nt_true) )
    
              output['loglike'][ii,jj,kk,:] = x_hat_EnKF.loglik[:]
              output['rmsea'][ii,jj,kk,:] = np.sqrt( np.mean( np.power( x_hat_EnKF.xamean - xt.values , 2 ) , 1 ) )
              output['rmsef'][ii,jj,kk,:] = np.sqrt( np.mean( np.power( x_hat_EnKF.xfmean - xt.values , 2 ) , 1 ) )
              output['biasa'][ii,jj,kk,:] = np.mean( x_hat_EnKF.xamean - xt.values , 1 )
              output['biasf'][ii,jj,kk,:] = np.mean( x_hat_EnKF.xfmean - xt.values , 1 )

              output['rmseao'][ii,jj,kk,:] = np.sqrt( np.nanmean( np.power( x_hat_EnKF.xamean - yo.values , 2 ) , 1 ) )
              output['rmsefo'][ii,jj,kk,:] = np.sqrt( np.nanmean( np.power( x_hat_EnKF.xfmean - yo.values , 2 ) , 1 ) )
              output['biasao'][ii,jj,kk,:] = np.nanmean( x_hat_EnKF.xamean - xt.values , 1 )
              output['biasfo'][ii,jj,kk,:] = np.nanmean( x_hat_EnKF.xfmean - xt.values , 1 )

              output['xt'][ii,jj,kk,:] = xt.values[:,center_domain_var_ind]
              output['xa'][ii,jj,kk,:] = x_hat_EnKF.xamean[:,center_domain_var_ind]
              output['xf'][ii,jj,kk,:] = x_hat_EnKF.xfmean[:,center_domain_var_ind]
              output['yo'][ii,jj,kk,:] = yo.values[:,center_domain_var_ind]
              output['xc'][ii,jj,kk,:] = AF.catalog.analogs[:,center_domain_var_ind]

              output['spread_f'][ii,jj,kk,:] = np.std( x_hat_EnKF.xfens[:,:,center_domain_var_ind] , 1 )
              output['spread_a'][ii,jj,kk,:] = np.std( x_hat_EnKF.xaens[:,:,center_domain_var_ind] , 1 )

              output['multinf'][ii,jj,kk,:] = x_hat_EnKF.alpha[:] 
          
              kk = kk + 1

              del x_hat_EnKF , AF.catalog , xt , yo , AF.kdt

              gc.collect()
          
           jj = jj + 1
       
       ii = ii + 1

   with open(output_file_name, 'wb') as my_file  :
       pickle.dump( output , my_file , protocol=pickle.HIGHEST_PROTOCOL)

