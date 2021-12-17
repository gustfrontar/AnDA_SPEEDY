
import sys
sys.path.append('../common_functions/')
sys.path.append('./')

import numpy as np
#import matplotlib.pyplot as plt
import datetime as dt
import ctl_reader as ctlr
import os
from scipy import linalg as la


def reformat_speedy_catalog( MAINPATH , EXPNAME , INI_DATE , END_DATE , CTL='yyyymmddhh.ctl' )   :

   #This version usually won't work due to memory issues. Need to update this function
   #to be more memory efficient.
 
   #INI_DATE='1995010100' 
   #END_DATE='1995011000'

   #MAINPATH='/home/jruiz/Dropbox/DATA/analog_data_assimilation_SPEEDY/TEST_DATA/'
   #EXPNAME='TEST_CATALOG'

   OUTDIR=MAINPATH + EXPNAME + '/single_point_data/'


   TIME_DELTA=dt.timedelta( hours = 6 )

   itime=dt.datetime.strptime( INI_DATE , '%Y%m%d%H%M' )
   etime=dt.datetime.strptime( END_DATE , '%Y%m%d%H%M' )

   total_times = int( np.ceil( ( etime - itime ).total_seconds() / ( 6.0 * 3600 ) ) ) + 1

   #=========================================================
   #  READ CTL FILE
   #=========================================================

   CTL= MAINPATH + '/' + EXPNAME + '/' + CTL
   ctl_dict = ctlr.read_ctl( CTL )

   #Estimate number of bins in histogram. 
   nx=np.int(ctl_dict['nx'])
   ny=np.int(ctl_dict['ny'])
   nz=np.array(ctl_dict['end_record']).max() + 1 #Total number of records in binary file.
   nlev=ctl_dict['nz']                 #Number of vertical levels for 3D variables.

   undef=np.float32( ctl_dict['undef'] )

   if  ctl_dict['big_endian']   :
      dtypein = '>f4'
      endian='big_endian'
   else                         :
      dtypein = 'f4'
      endian='little_endian'
   if  ctl_dict['sequential']   :
      access='sequential'
   else                         :
      access='direct'
      sequential=ctl_dict['sequential']
 
   #All the fields will be read so..
   n_selected_fields = nz
   selected_fields = np.arange(0,nz) + 1

   #[ lon , lat ] = np.meshgrid( ctl['xlevels'] , ctl['ylevels'] )

   #=========================================================
   #  ALLOCATE MEMORY
   #=========================================================

   speedy_data = dict()

   for ivar,my_var in enumerate(ctl_dict['var_list']) :
      my_var_size = int( ctl_dict['var_size'][ivar] )
      if my_var_size > 0 :
          speedy_data[my_var] = np.zeros( ( nx , ny , my_var_size , total_times ) ).astype('float32')
      else                         :
          speedy_data[my_var] = np.zeros( ( nx , ny , 1 , total_times ) ).astype('float32')


   #=========================================================
   #  READ DATA
   #=========================================================

   #Write the dates to a text file. 
   my_file = OUTDIR + '/dates.txt'

   os.makedirs( OUTDIR , exist_ok=True )

   fascii=open( my_file , 'w+')


   ctime=itime

   it = 0 

   while ( ctime <= etime )  :

      my_file = MAINPATH + '/' + EXPNAME + '/' + ctime.strftime("%Y%m%d%H") + '.grd'

      print( ' Reading file ', my_file )
   
      my_data = ctlr.read_data_grads(my_file,ctl_dict,use_nt_from_ctl=False,nt=1,masked=False)

      #Transfer data into structure.

      for my_var in my_data :
         if my_data[my_var].shape[2] == 1 :
            speedy_data[my_var][:,:,0,it] = np.squeeze( my_data[my_var] )
         else                             :
            speedy_data[my_var][:,:,:,it] = np.squeeze( my_data[my_var] )

   
      fascii.write(ctime.strftime("%Y%m%d%H") + "\n")
      it = it + 1
      ctime = ctime + TIME_DELTA


   #=========================================================
   #  SPLIT DATA INTO SINGLE LOCATION - SINGLE VARIABLE FILES
   #=========================================================

   OUTDIR=MAINPATH + EXPNAME + '/single_point_data/'

   for ivar , my_var in enumerate( speedy_data ) :

      for ilon , my_lon in enumerate(ctl_dict['xlevels'])  :

         for ilat , my_lat in enumerate(ctl_dict['ylevels'])  :

            if ctl_dict['var_size'][ivar] == 0  :
               nlev=1
            else                                :
               nlev=int( ctl_dict['var_size'][ivar] )

            for ilev  in range( nlev )  :

               tmp_data=speedy_data[my_var][ilon,ilat,ilev,:]

               my_file = OUTDIR + '/var_' + my_var + '_lon_' + str(ilon) + '_lat_' + str(ilat) + '_lev_' + str(ilev) + '.bin' 

               tmp_data.astype('float32').tofile(my_file)


def read_speedy_data( MAINPATH , EXPNAME , CTL , GRID_LIMITS , TIME_LIMITS , VAR_LIST )  :

   #MAINPATH is the mainpath of the data
   #EXPNAME is the experiment name (last folder in the path)
   #GRID_LIMITS describe the domain limits ( min_lon , max_lon , min_lat , max_lat , min_lev , max_lev )
   #TIME_LIMITS describe the time range limits ( min_date , max_date ) * character list (date format yyyymmddhh)
   #VAR_LIST = list of variables that will be incorporated.

   #=========================================================
   #  READ CTL FILE
   #=========================================================
   MY_DIR = MAINPATH + '/' + EXPNAME + '/' + 'single_point_data/'
   my_data = dict()

   my_data['var_list'] = VAR_LIST

   CTL= MAINPATH + '/' + EXPNAME + '/' + CTL

   ctl_dict = ctlr.read_ctl( CTL )

   #Estimate number of bins in histogram. 
   #nx=np.int(ctl_dict['nx'])
   #ny=np.int(ctl_dict['ny'])
   #nz=np.array(ctl_dict['end_record']).max() + 1 #Total number of records in binary file.
   #nlev=ctl_dict['nz']                 #Number of vertical levels for 3D variables.

   #undef=np.float32( ctl_dict['undef'] )

   #READ THE LIST OF TIMES AND ADJUST THE DATE RANGE ACCORDINGLY.

   itime_req=dt.datetime.strptime( TIME_LIMITS[0] , '%Y%m%d%H%M' )
   etime_req=dt.datetime.strptime( TIME_LIMITS[1] , '%Y%m%d%H%M' )   

   
   tmp_dates_str=[]
   my_file = MY_DIR + '/dates.txt'
   fascii=open( my_file , 'r')
   for line in fascii.readlines() :
     tmp_dates_str.append( line[0:10] )
  
   itime_input=dt.datetime.strptime( tmp_dates_str[0] , '%Y%m%d%H%M' )     
   etime_input=dt.datetime.strptime( tmp_dates_str[-1] , '%Y%m%d%H%M' )

   total_times_input = int( np.ceil( ( etime_input - itime_input ).total_seconds() / ( 6.0 * 3600 ) ) ) + 1

   if itime_req >= itime_input :
      itime = itime_req
   else                        :
      itime = itime_input
      it_ini = 0
      print('Warning: time range updated new ini time is :', itime_input)
   if etime_req <= etime_input :
      etime = etime_req 
   else                        :
      print('Warning: time range updated new end time is :', etime_input)
      etime = etime_input

   it_ini = int( ( itime - itime_input ).total_seconds() / ( 3600.0 * 6.0 ) )
   it_end = int( ( etime - itime_input ).total_seconds() / ( 3600.0 * 6.0 ) ) + 1

   total_times = ( it_end - it_ini )

   if GRID_LIMITS[0] > GRID_LIMITS[1]   :  #Then we assume that Greenwich is in 
                                           #the middle of the domin.
      #Use np.logical_or instead of np.logical_and.                                     
      xind=np.where( np.logical_or( ctl_dict['xlevels'] >= GRID_LIMITS[0] , ctl_dict['xlevels'] <= GRID_LIMITS[1] ) )
   else                                 :
      xind=np.where( np.logical_and( ctl_dict['xlevels'] >= GRID_LIMITS[0] , ctl_dict['xlevels'] <= GRID_LIMITS[1] ) )
   yind=np.where( np.logical_and( ctl_dict['ylevels'] >= GRID_LIMITS[2] , ctl_dict['ylevels'] <= GRID_LIMITS[3] ) )
   zind=np.where( np.logical_and( ctl_dict['vlevels'] >= GRID_LIMITS[4] , ctl_dict['vlevels'] <= GRID_LIMITS[5] ) )

   nvar= len( VAR_LIST )   

   ctime = itime
   my_data['dates_str']=[]
   my_data['dates']=[] 
   while ctime <= etime    :

      my_data['dates'].append(ctime)
      my_data['dates_str'].append( ctime.strftime("%Y%m%d%H") )

      ctime = ctime + dt.timedelta( hours = 6 )

   #Allocate variable:

   my_data['data'] = np.zeros( ( np.size(xind)*np.size(yind)*np.size(zind)*nvar , total_times ) )
      
   my_data['xind']=[]
   my_data['yind']=[]
   my_data['zind']=[]
   my_data['var'] =[]
   
   print('==================================================================')
   print('NX=',np.size(xind))
   print('NY=',np.size(yind))
   print('NZ=',np.size(zind))
   print('XRANGE=',ctl_dict['xlevels'][xind])
   print('YRANGE=',ctl_dict['ylevels'][yind])
   print('ZRANGE=',ctl_dict['vlevels'][zind])   
   print('==================================================================')
   
   my_ind = 0

   for ivar , my_var in enumerate( VAR_LIST )  :
      for ix in np.nditer(xind) :
         for iy in np.nditer(yind) :
            for iz in np.nditer(zind) :
                               
               my_file = MY_DIR + '/var_' + my_var + '_lon_' + str(ix) + '_lat_' + str(iy) + '_lev_' + str(iz) + '.bin'

               print('I m reading file ', my_file )

               my_data['data'][my_ind,:] = np.fromfile(my_file,dtype='float32',count=total_times_input)[it_ini:it_end]

               my_data['xind'].append( ix )
               my_data['yind'].append( iy )
               my_data['zind'].append( iz )
               my_data['var'].append( my_var )

               my_ind = my_ind + 1

   return my_data 


def read_speedy_data_grid_based( MAINPATH , EXPNAME , CTL , CENTER_POINT , GRID_SPAN , TIME_LIMITS , VAR_LIST , STANDARIZATION = False , STANDARIZATION_WINDOW = 30 )  :

   #MAINPATH is the mainpath of the data
   #EXPNAME is the experiment name (last folder in the path)
   #CENTER_POINT (i,j,k) of the box center
   #GRID_SPAN (delta_x , delta_y , delta_k ) that will be taken 
   #TIME_LIMITS describe the time range limits ( min_date , max_date ) * character list (date format yyyymmddhh)
   #VAR_LIST = list of variables that will be incorporated.

   #=========================================================
   #  READ CTL FILE
   #=========================================================
   MY_DIR = MAINPATH + '/' + EXPNAME + '/' + 'single_point_data/'
   my_data = dict()       #We will store only the requested data.

   my_raw_data = dict()   #We will store some intermediate data (full catalog length)
                          #this will be useful for reading the data and for standarizing the data.

   my_data['var_list'] = VAR_LIST

   CTL= MAINPATH + '/' + EXPNAME + '/' + CTL

   ctl_dict = ctlr.read_ctl( CTL )

   #Estimate number of bins in histogram. 
   #nx=np.int(ctl_dict['nx'])
   #ny=np.int(ctl_dict['ny'])
   #nz=np.array(ctl_dict['end_record']).max() + 1 #Total number of records in binary file.
   #nlev=ctl_dict['nz']                 #Number of vertical levels for 3D variables.

   #undef=np.float32( ctl_dict['undef'] )

   #READ THE LIST OF TIMES AND ADJUST THE DATE RANGE ACCORDINGLY.

   itime_req=dt.datetime.strptime( TIME_LIMITS[0] , '%Y%m%d%H%M' )
   etime_req=dt.datetime.strptime( TIME_LIMITS[1] , '%Y%m%d%H%M' )   

   
   tmp_dates_str=[]
   my_raw_data['dates']=[]
   my_raw_data['dates_str']=[]
   my_file = MY_DIR + '/dates.txt'
   fascii=open( my_file , 'r')
   for line in fascii.readlines() :
     #tmp_dates_str.append( line[0:10] )
     my_raw_data['dates'].append( dt.datetime.strptime(line[0:10] ,"%Y%m%d%H" ) )
     my_raw_data['dates_str'].append(line[0:10]  )

  
   itime_input=dt.datetime.strptime( my_raw_data['dates_str'][0] , '%Y%m%d%H%M' )     
   etime_input=dt.datetime.strptime( my_raw_data['dates_str'][-1] , '%Y%m%d%H%M' )

   total_times_input = int( np.ceil( ( etime_input - itime_input ).total_seconds() / ( 6.0 * 3600 ) ) ) + 1

   if itime_req >= itime_input :
      itime = itime_req
   else                        :
      itime = itime_input
      it_ini = 0
      print('Warning: time range updated new ini time is :', itime_input)
   if etime_req <= etime_input :
      etime = etime_req 
   else                        :
      print('Warning: time range updated new end time is :', etime_input)
      etime = etime_input

   it_ini = int( ( itime - itime_input ).total_seconds() / ( 3600.0 * 6.0 ) )
   it_end = int( ( etime - itime_input ).total_seconds() / ( 3600.0 * 6.0 ) ) + 1

   total_times = ( it_end - it_ini )

   nx=np.size( ctl_dict['xlevels'])
   ny=np.size( ctl_dict['ylevels'])
   nz=np.size( ctl_dict['vlevels'])

   xind = np.arange( CENTER_POINT[0] - GRID_SPAN[0] , CENTER_POINT[0] + GRID_SPAN[0] + 1 )
   yind = np.arange( CENTER_POINT[1] - GRID_SPAN[1] , CENTER_POINT[1] + GRID_SPAN[1] + 1 )
   zind = np.arange( CENTER_POINT[2] - GRID_SPAN[2] , CENTER_POINT[2] + GRID_SPAN[2] + 1 ) 
   
   #Apply boundary conditions in X
   for ii in range( np.size( xind ) ) :
       if xind[ii] < 0  :
           xind[ii] = nx + xind[ii]
       if xind[ii] > nx-1 :
           xind[ii] = xind[ii] - nx 

   xind = xind[np.logical_and( xind >= 0 , xind <= nx -1 )]           
   yind = yind[np.logical_and( yind >= 0 , yind <= ny -1 )]
   zind = zind[np.logical_and( zind >= 0 , zind <= nz -1 )]
   
   nvar= len( VAR_LIST )   

   ctime = itime
   my_data['dates_str']=[]
   my_data['dates']=[] 

   while ctime <= etime    :

      my_data['dates'].append(ctime)
      my_data['dates_str'].append( ctime.strftime("%Y%m%d%H") )

      ctime = ctime + dt.timedelta( hours = 6 )

   #Allocate variable:

   my_data['data'] = np.zeros( ( np.size(xind)*np.size(yind)*np.size(zind)*nvar , total_times ) )
      
   my_data['xind']=[]
   my_data['yind']=[]
   my_data['zind']=[]
   my_data['var'] =[]
   
   print('==================================================================')
   print('NX=',np.size(xind))
   print('NY=',np.size(yind))
   print('NZ=',np.size(zind))
   print('XRANGE=',ctl_dict['xlevels'][xind])
   print('YRANGE=',ctl_dict['ylevels'][yind])
   print('ZRANGE=',ctl_dict['vlevels'][zind])   
   print('==================================================================')
   
   my_ind = 0

   year_length = 366
   my_data['seasonal_mean'] = np.zeros( ( np.size(xind)*np.size(yind)*np.size(zind)*nvar , year_length ) )
   my_data['seasonal_std'] = np.zeros( ( np.size(xind)*np.size(yind)*np.size(zind)*nvar , year_length ) )
   my_data['data_count'] = np.zeros( ( np.size(xind)*np.size(yind)*np.size(zind)*nvar , year_length ) )

   for ivar , my_var in enumerate( VAR_LIST )  :
      for ix in np.nditer(xind) :
         for iy in np.nditer(yind) :
            for iz in np.nditer(zind) :
                               
               my_file = MY_DIR + '/var_' + my_var + '_lon_' + str(ix) + '_lat_' + str(iy) + '_lev_' + str(iz) + '.bin'

               #print('I m reading file ', my_file )

               my_raw_data['data'] = np.fromfile(my_file,dtype='float32',count=total_times_input)

               if STANDARIZATION :
                   my_raw_data  = normalize_data_seasonal( my_raw_data , STANDARIZATION_WINDOW )
                   #print( my_raw_data['data_count'] )
                   my_data['data'][my_ind,:]  = my_raw_data['data'][it_ini:it_end]
                   #Save seasonal mean, standard deviation and data count
                   my_data['seasonal_mean'][my_ind,:]=my_raw_data['seasonal_mean']
                   my_data['seasonal_std'][my_ind,:]=my_raw_data['seasonal_std']
                   my_data['data_count'][my_ind,:]=my_raw_data['data_count']

               else  :
                   my_data['data'][my_ind,:] = my_raw_data['data'][it_ini:it_end]

               my_data['xind'].append( ix )
               my_data['yind'].append( iy )
               my_data['zind'].append( iz )
               my_data['var'].append( my_var )

               my_ind = my_ind + 1

      #print('max data count', np.max( my_data['data_count'] ) )
      #quit()


   return my_data  


def normalize_data( my_data )  :

   nvar = np.shape( my_data['data'] )[0]
   nt   = np.shape( my_data['data'] )[1]

   my_data['mean']=np.mean( my_data['data'] , 1 )
   my_data['std'] =np.std( my_data['data']  , 1 )

   my_data['data_std']=np.zeros( (nvar , nt) )

   for ii in range( nvar )   :

       my_data['data_std'][ii,:] = ( my_data['data'][ii,:] - my_data['mean'][ii] ) / my_data['std'][ii] 

   return my_data


def normalize_data_seasonal( my_data , normalize_window ) :

    #Normalize data tacking into account the seasonal cycle.
    #my_data a numpy array with the following dimensions [times]
    #my_dates a numpy array with as many elements as columns in my_data (times) in the ctime internal numerical format.
    #normalize_window the number of days to consider for computing local mean and standard deviation. 
    #data_freq time frequency in hours.

    #Two new entries are added to my_data:
    #seasonal_mean = containing the climatological mean for each day of the year
    #seasonal_std  = containing the climatological standard deviation for each day of the year.
    #data_count    = containing the number of intances used to compute the mean and standar deviation for each day of the year.

    #WARNING: This routine is designed for SPEEDY data which has no diurnal cycle. 
    #The diurnal cycle is not removed from the data and the climatology is computed only as a function of the day of the year.

    #nvar = np.shape( my_data['data'] )[0]

    year_length = 366

    #Compute the seasonal mean and standard deviation.
    my_data['seasonal_mean'] = np.zeros( ( year_length  ) )
    my_data['seasonal_std']  = np.zeros( ( year_length  ) )

    my_data['data_count'] = np.zeros( ( year_length  ) )
     
    for ii , date in enumerate( my_data['dates'] ) :

       doy = date.timetuple().tm_yday  #Get day of year.

       max_index = int( doy + normalize_window - 1 )
       min_index = int( doy - normalize_window - 1 )

       if max_index > year_length - 1 :
          max_index = int( max_index - year_length ) 
       if min_index < 0 :
          min_index = int( year_length + min_index )  

       if max_index < min_index :

          my_data['seasonal_mean'][0:max_index+1] = my_data['seasonal_mean'][0:max_index+1] + my_data['data'][ii]
          my_data['seasonal_std'][0:max_index+1] = my_data['seasonal_std'][0:max_index+1] + my_data['data'][ii] ** 2
          my_data['data_count'][0:max_index+1] = my_data['data_count'][0:max_index+1] + 1

          my_data['seasonal_mean'][min_index:] = my_data['seasonal_mean'][min_index:] + my_data['data'][ii]
          my_data['seasonal_std'][min_index:] = my_data['seasonal_std'][min_index:] + my_data['data'][ii] ** 2
          my_data['data_count'][min_index:] = my_data['data_count'][min_index:] + 1

       else                     :

          my_data['seasonal_mean'][min_index:max_index+1] = my_data['seasonal_mean'][min_index:max_index+1] + my_data['data'][ii]
          my_data['seasonal_std'][min_index:max_index+1] = my_data['seasonal_std'][min_index:max_index+1] + my_data['data'][ii] ** 2
          my_data['data_count'][min_index:max_index+1] = my_data['data_count'][min_index:max_index+1] + 1

    my_data['seasonal_mean'] = np.where( my_data['data_count'] > 0 , my_data['seasonal_mean'] / my_data['data_count'] , np.nan + np.zeros( ( year_length  ) ) )   

    my_data['seasonal_std']  = np.where( my_data['data_count'] > 0 , ( my_data['seasonal_std'] / my_data['data_count'] - ( my_data['seasonal_mean'] ** 2 ) ) ** 0.5 , np.nan + np.zeros( ( year_length  ) ) )

   
    #Standarize the data using the seasonal mean and standard deviation.

    for ii , date in enumerate( my_data['dates'] ) :

        doy = date.timetuple().tm_yday  #Get day of year.

        #print( ii , my_data['data'][ii] , my_data['seasonal_mean'][doy-1] , my_data['seasonal_std'][doy-1] )

        my_data['data'][ii] = ( my_data['data'][ii] - my_data['seasonal_mean'][doy-1] ) / my_data['seasonal_std'][doy-1]

    return my_data 

def pca_data( my_data , variable_wise = False , variance_threshold = 1.0 )  :

    #Perform PCA analysis on the input data (speedy format)
    #If variable_wise == True then perform PCA analysis for each variable
    #separatedely. 

    nvar , nt = np.shape( my_data['data'] )

    if not 'data_std' in my_data  :

       my_data = normalize_data( my_data )

    if not( variable_wise )  :

       #Consider all the variables , levels and grid locations at the same time.

       tmp_data_trans , eivals , eivecs = pca( np.transpose( my_data['data_std'] ) ) 

       my_data['data_trans'] = np.transpose( tmp_data_trans )
       my_data['eivecs'] = np.transpose( eivecs )
       my_data['eivals'] = eivals

       #Perform data compression. 

       my_data['eivals_norm'] = my_data['eivals'] / np.sum( my_data['eivals'] )

       acum_variance = np.cumsum( my_data['eivals_norm'] )

       my_mask = acum_variance <= variance_threshold 

       my_data['data_trans_compress'] = my_data['data_trans'][my_mask,:]

       my_data['eival_compress'] = my_data['eivals_norm'][my_mask]
    
       my_data['eivec_compress'] = my_data['eivecs'][my_mask,:]


    else                     :

       my_data['data_trans'] = np.zeros( np.shape( my_data['data'] ) )
       my_data['eivecs']     = np.zeros( nvar )
       my_data['eivals']     = np.zeros( nvar )

       compress_mask = np.zeros( nvar ).astype('bool')
      
       for ivar , my_var in enumerate( my_data['var_list'] )  :

          var_mask = np.zeros( nvar ).astype('bool')

          for tmp_ivar , tmp_my_var in enumerate( my_data['var'] ) :
             if tmp_my_var == my_var  :
                var_mask = True

          tmp_data_trans , tmp_eivals , tmp_eivecs = pca( np.transpose( my_data['data_std'][var_mask,:] ) )

          my_data['data_trans'][var_mask,:] = np.transpose( tmp_data_trans )
          my_data['eivecs'][var_mask,:]       = np.transpose( tmp_eivecs )
          my_data['eivals'][var_mask]       = tmp_eivals

          #Perform data compression. 

          tmp_eivals_norm = tmp_eivals / np.sum( tmp_eivals )

          my_data['eivals_norm'][var_mask] = tmp_eivals_norm

          acum_variance = np.cumsum( tmp_eivals_norm )

          my_mask = acum_variance <= variance_threshold 
 
          compress_mask[var_mask][my_mask] = True

       my_data['data_trans_compress'] = my_data['data_trans'][compress_mask,:]

       my_data['eival_compress'] = my_data['eivals_norm'][compress_mask]

       my_data['eivec_compress'] = my_data['eivecs'][compress_mask,:]

    return my_data           


def pca(data) :
    """
    returns: data transformed in 2 dims/columns + regenerated original data
    pass in: data as 2D NumPy array
    """
    m, n = data.shape
    # mean center the data
    data -= data.mean(axis=0)
    # calculate the covariance matrix
    R = np.cov(data, rowvar=False)
    # calculate eigenvectors & eigenvalues of the covariance matrix
    # use 'eigh' rather than 'eig' since R is symmetric, 
    # the performance gain is substantial
    evals, evecs = la.eigh(R)
    # sort eigenvalue in decreasing order
    idx = np.argsort(evals)[::-1]
    evecs = evecs[:,idx]
    # sort eigenvectors according to same index
    evals = evals[idx]
    # select the first n eigenvectors (n is desired dimension
    # of rescaled data array, or dims_rescaled_data)
    evecs = evecs[:,:]
    # carry out the transformation on the data using eigenvectors
    # and return the re-scaled data, eigenvalues, and eigenvectors
    return np.dot(evecs.T, data.T).T, evals, evecs

 
