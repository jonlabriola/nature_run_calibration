import os,argparse
import numpy as np
from netCDF4 import Dataset
import scipy.stats
import matplotlib
import matplotlib.pyplot as plt
from multiprocessing import Pool
import scipy.ndimage as ndimage
matplotlib.rcParams['axes.linewidth'] = 1.5

parser = argparse.ArgumentParser()
parser.add_argument("time",type=int,help='File Path')
parser.add_argument("fcst_var",type=str,help = 'Var to plot')
parser.add_argument("ob_var",type=str,help = 'Var to plot')
parser.add_argument("lev",type=int,help = 'Var to plot')
parser.add_argument("--vertical_localization",action='store_true',help='Look at the vertical localization')
parser.add_argument("--storm_obs",action='store_true',help='Observations Taken From Storms')
parser.add_argument("--storm_fcst",action='store_true',help='Observations Taken From Forecasts')
arguments = parser.parse_args()


fcst_var = arguments.fcst_var  #--- The forecast variable
ob_var = arguments.ob_var      #--- The observation variable
lev = arguments.lev            #--- The vertical model level to look at
nproc = 20                     #--- The number of processors if you need to run multi processor
nen = 10                       #--- Ensemble Size
thin_incr = 25                 #--- The number of grid points to thin observations
if arguments.vertical_localization:
   dincr = 0.25                      #--- The correlation radius increment 
   eval_dis = np.arange(dincr,6.0,dincr) #--- The range of correlation radii to look at
   thin_incr = 10
else:
   dincr = 2                      #--- The correlation radius increment 
   eval_dis = np.arange(dincr,70,dincr) #--- The range of correlation radii to look at


ens_cor = np.zeros((nen,len(eval_dis)))

#--- Collect information from each ensemble member
for mindex,mem in enumerate(range(1,nen+1)):
   filepath = '/work/jonathan.labriola/TurbPBLExp/all/LES_%02d/cm1out_%06d.nc'%(mem,arguments.time)
   print('Opening ... %s'%filepath)
   dumpfile = Dataset(filepath,"r",fortmat="NETCDF4")


   #--- Step 1: Grab Observations over subdomain, thin every few grid points
   if arguments.storm_obs:
      xmin = 100
      xmax = 300
      ymin = 0
      ymax = 800
   else:
      xmin = 500
      xmax = 800
      ymin = 0
      ymax = 800

   #--- Collect observation info
   #--- The observation array is flattened so you can loop over obs
   obs = np.squeeze(dumpfile.variables[ob_var][0,lev,ymin:ymax,xmin:xmax])
   obs = ndimage.gaussian_filter(obs,sigma=8.0)[::thin_incr,::thin_incr]
   [ny,nx] = obs.shape 
   obx = np.tile(dumpfile.variables['xh'][xmin:xmax:thin_incr],(ny,1)).flatten(order='F')
   oby = np.transpose(np.tile(dumpfile.variables['yh'][ymin:ymax:thin_incr],(nx,1))).flatten(order='F')
   obs = obs.flatten(order='F')
   nob = obs.shape[0]
   obz = np.tile(dumpfile.variables['zh'][lev],(nob)) 

 
   #--- Step 2: Forecast Section
   if arguments.storm_fcst:
      xmin = 100
      xmax = 300
      ymin = 0
      ymax = 800
   else:
      xmin = 500
      xmax = 800
      ymin = 0
      ymax = 800

   #--- Vertical Localization
   #--- Preserve all forecast gridpoints above and below the selected observations
   if arguments.vertical_localization:
      fcst = np.squeeze(dumpfile.variables[fcst_var][0,:,ymin:ymax,xmin:xmax])
      [nz,ny,nx] = fcst.shape
      for zindex in range(0,nz):
         #fcst[zindex] = fcst[zindex] - np.mean(fcst[zindex])
         fcst[zindex] =  ndimage.gaussian_filter(fcst[zindex],sigma=8.0)
      #--- make forecast array [nob,nz]
      fcst = np.transpose(np.reshape(fcst[:,::thin_incr,::thin_incr],(nz,nob),order='F'))
      zh_fcst = dumpfile.variables['zh'][:] #--- zh is the same at all vertical grid points   
      
   #--- Horizontal Localization 
   #--- Preserve all forecast gridpoints on the same horizontal plane as observations
   else:
      fcst = np.squeeze(dumpfile.variables[fcst_var][0,lev,ymin:ymax,xmin:xmax])
      fcst = ndimage.gaussian_filter(fcst,sigma=8.0)
      [ny,nx] = fcst.shape
      xh_fcst = np.tile(dumpfile.variables['xh'][xmin:xmax],(ny,1))
      yh_fcst = np.transpose(np.tile(dumpfile.variables['yh'][ymin:ymax],(nx,1)))

 


   #--- Calculate vertical localization
   #--- reform arrays into [nob,nz]
   if arguments.vertical_localization:
     if mindex == 0:
        [nz] = zh_fcst.shape
        zh_fcst = np.tile(zh_fcst,(nob,1))
        obz2d = np.transpose(np.tile(obz,(nz,1)))
        tot_dis = np.absolute(zh_fcst-obz2d)
        #tot_dis = -1 *(zh_fcst-obz2d)
        #tot_dis = (zh_fcst-obz2d)

     for dindex,dis in enumerate(eval_dis):
        [obindex,zindex] = np.where((tot_dis >= (dis-float(dincr)/2.)) & (tot_dis < (dis + float(dincr)/2.)))
        forecast = fcst[obindex,zindex]
        observation = obs[obindex]
        if observation.shape[0] > 0:
           ens_cor[mindex,dindex] = np.absolute(scipy.stats.pearsonr(observation, forecast)[0])
        else:
           ens_cor[mindex,dindex] = np.nan
  
   #--- Calculate horizontal localization over a "small" number of obs
   #--- reform arrays in [nob,ny,nx]
   elif nob < 10000:  
      if mindex == 0:
         [ny,nx] = xh_fcst.shape
         xh_fcst = np.tile(xh_fcst,(nob,1,1))
         yh_fcst = np.tile(yh_fcst,(nob,1,1))
         obx3d = np.transpose(np.tile(obx,(nx,ny,1)))
         oby3d = np.transpose(np.tile(oby,(nx,ny,1)))
         tot_dis = ((yh_fcst - oby3d)**2 + (xh_fcst - obx3d)**2)**0.5

      for dindex,dis in enumerate(eval_dis):
         #[obindex,yindex,xindex] = np.where((tot_dis >= (dis-0.25)) & (tot_dis < (dis + 0.25)))
         [obindex,yindex,xindex] = np.where((tot_dis >= (dis-float(dincr)/2.)) & (tot_dis < (dis + float(dincr)/2.)))
         forecast = fcst[yindex,xindex]
         observation = obs[obindex]
         if observation.shape[0] > 0:
            print(observation.shape)
            # JDL  - Come back here, this is where correlations are being messed up
            ens_cor[mindex,dindex] = np.absolute(scipy.stats.pearsonr(observation, forecast)[0])
         else:
            ens_cor[mindex,dindex] = np.nan     

   #--- Multi-Processor Function if calculating a large number of obs
   #--- Memory friendly, but much slower
   else:       
      def obloop(oindex):
         tot_dis = ((yh_fcst - oby[oindex])**2 + (xh_fcst - obx[oindex])**2)**0.5
         [yindex,xindex] = np.where((tot_dis >= (dis-0.25)) & (tot_dis < (dis + 0.25)))
         combined_array = np.zeros((2,len(yindex)))
         combined_array[0,:] = obs[oindex]
         combined_array[1,:] = fcst[yindex,xindex]
         return combined_array 


      distance  = [] 
      cor = []
      for dindex,dis in enumerate(eval_dis):
         print("Evolatulating Distance = ",dis)
         #--- Multi Processing Section
         if nproc > 1:
            p = Pool(nproc)
            cor_info = p.map(obloop,np.arange(0,nob))
            p.close()
            for oindex in np.arange(0,nob):
               if oindex == 0:
                  observation = cor_info[oindex][0]
                  forecast = cor_info[oindex][1]
               else: 
                  observation = np.concatenate([observation,cor_info[oindex][0]])
                  forecast = np.concatenate([forecast,cor_info[oindex][1]])
         else:
            for oindex in range(0,nob):
               cor_info = obloop(oindex)
               if oindex == 0:
                  observation = cor_info[0]
                  forecast = cor_info[1]
               else:
                  observation = np.concatenate([observation,cor_info[0]])
                  forecast = np.concatenate([forecast,cor_info[1]])

         ens_cor[mindex,dindex] = np.absolute(scipy.stats.pearsonr(observation, forecast)[0])

      #print('The Correlation is = ',A)
      #print('P value =',p)
      #distance.append(dis)
      #cor.append(np.absolute(A))
for mindex in range(0,nen):
   plt.plot(eval_dis,ens_cor[mindex],color='slategrey',linewidth=2,alpha=0.25)
   
plt.plot(eval_dis,np.percentile(ens_cor,50,axis=0),color='crimson',linewidth=4,label='Ensemble Median')
plt.xlabel('Distance (km)')
plt.ylabel('Absolute Correlation')
plt.ylim([0,1])#np.nanmax(ens_cor)])
plt.xlim([0,np.nanmax(eval_dis)])
plt.legend()
if arguments.vertical_localization:
   outpath = '%s_%s_t%d_lev%03d_nen%03d_vert.png'%(arguments.fcst_var,arguments.ob_var,arguments.time,arguments.lev,nen)
else:
   outpath = '%s_%s_t%d_lev%03d_nen%03d.png'%(arguments.fcst_var,arguments.ob_var,arguments.time,arguments.lev,nen)

plt.title('Correlation between %s (obs) and %s (fcst)'%(arguments.ob_var,arguments.fcst_var))
plt.savefig(outpath)
#plt.show()
 
