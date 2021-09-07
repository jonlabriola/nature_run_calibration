import os,argparse
import numpy as np
from netCDF4 import Dataset
import scipy.stats
import matplotlib
import matplotlib.pyplot as plt
from multiprocessing import Pool
import jl_modules
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


fcst_var = arguments.fcst_var
ob_var = arguments.ob_var
lev = arguments.lev
nproc = 20
nen = 36
dincr = 2
eval_dis = np.arange(dincr,70,dincr)

ens_cor = np.zeros((nen,len(eval_dis)))
for mindex,mem in enumerate(range(1,nen+1)):
   filepath = '/work/jonathan.labriola/TurbPBLExp/all/LES_%02d/cm1out_%06d.nc'%(mem,arguments.time)
   print('Opening ... %s'%filepath)
   dumpfile = Dataset(filepath,"r",fortmat="NETCDF4")

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



   if fcst_var in ['wmax']:
      fcst = np.amax(np.squeeze(dumpfile.variables['winterp'][0,:,ymin:ymax,xmin:xmax]),axis=0)
   else:
      fcst = np.squeeze(dumpfile.variables[fcst_var][0,lev,ymin:ymax,xmin:xmax])
     
   fcst = ndimage.gaussian_filter(fcst,sigma=8.0)

   [ny,nx] = fcst.shape
   xh_fcst = np.tile(dumpfile.variables['xh'][xmin:xmax],(ny,1))
   yh_fcst = np.transpose(np.tile(dumpfile.variables['yh'][ymin:ymax],(nx,1)))

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


   obs = np.squeeze(dumpfile.variables[ob_var][0,lev,ymin:ymax,xmin:xmax])
   obs = ndimage.gaussian_filter(obs,sigma=8.0)
   [ny,nx] = obs.shape
   xh_obs = np.tile(dumpfile.variables['xh'][xmin:xmax],(ny,1))
   yh_obs = np.transpose(np.tile(dumpfile.variables['yh'][ymin:ymax],(nx,1)))


   #--- Step A: Collect Pseudo Obs
   increment = 25
   obs = obs[::increment,::increment].flatten()
   obx = xh_obs[::increment,::increment].flatten()
   oby = yh_obs[::increment,::increment].flatten()
   nob = obs.shape[0]


   def obloop(oindex):
      tot_dis = ((yh_fcst - oby[oindex])**2 + (xh_fcst - obx[oindex])**2)**0.5
      [yindex,xindex] = np.where((tot_dis >= (dis-0.25)) & (tot_dis < (dis + 0.25)))
      combined_array = np.zeros((2,len(yindex)))
      combined_array[0,:] = obs[oindex]
      combined_array[1,:] = fcst[yindex,xindex]
      return combined_array 



   if nob < 10000:
      [ny_fcst,nx_fcst] = xh_fcst.shape
      if mindex == 0:
         xh_fcst = np.tile(xh_fcst,(nob,1,1))
         yh_fcst = np.tile(yh_fcst,(nob,1,1))
    
         tmp = np.tile(obx,(nx_fcst,ny_fcst,1))
         obx3d = np.transpose(np.tile(obx,(nx_fcst,ny_fcst,1)))
         oby3d = np.transpose(np.tile(oby,(nx_fcst,ny_fcst,1)))
         tot_dis = ((yh_fcst - oby3d)**2 + (xh_fcst - obx3d)**2)**0.5

      for dindex,dis in enumerate(eval_dis):
         #[zindex,yindex,xindex] = np.where((tot_dis >= (dis-0.25)) & (tot_dis < (dis + 0.25)))
         [zindex,yindex,xindex] = np.where((tot_dis >= (dis-float(dincr)/2.)) & (tot_dis < (dis + float(dincr)/2.)))
         gridpts = fcst[yindex,xindex]
         observation = obs[zindex]
         ens_cor[mindex,dindex] = np.absolute(scipy.stats.pearsonr(observation, gridpts)[0])
        
   else:       
      #def obloop(oindex):
      #   tot_dis = ((yh_fcst - oby[oindex])**2 + (xh_fcst - obx[oindex])**2)**0.5
      #   [yindex,xindex] = np.where((tot_dis >= (dis-0.25)) & (tot_dis < (dis + 0.25)))
      #   combined_array = np.zeros((2,len(yindex)))
      #   combined_array[0,:] = obs[oindex]
      #   combined_array[1,:] = fcst[yindex,xindex]
      #   return combined_array 


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
                  gridpts = cor_info[oindex][1]
               else: 
                  observation = np.concatenate([observation,cor_info[oindex][0]])
                  gridpts = np.concatenate([gridpts,cor_info[oindex][1]])
         else:
            for oindex in range(0,nob):
               cor_info = obloop(oindex)
               if oindex == 0:
                  observation = cor_info[0]
                  gridpts = cor_info[1]
               else:
                  observation = np.concatenate([observation,cor_info[0]])
                  gridpts = np.concatenate([gridpts,cor_info[1]])

         ens_cor[mindex,dindex] = np.absolute(scipy.stats.pearsonr(observation, gridpts)[0])

      #print('The Correlation is = ',A)
      #print('P value =',p)
      #distance.append(dis)
      #cor.append(np.absolute(A))
for mindex in range(0,nen):
   plt.plot(eval_dis,ens_cor[mindex],color='slategrey',linewidth=2,alpha=0.25)
   
plt.plot(eval_dis,np.percentile(ens_cor,50,axis=0),color='crimson',linewidth=4,label='Ensemble Median')
plt.xlabel('Distance (km)')
plt.ylabel('Absolute Correlation')
plt.ylim([0,np.amax(ens_cor)])
plt.xlim([0,np.amax(eval_dis)])
plt.legend()
outpath = '%s_%s_t%d_lev%03d_nen%03d.png'%(arguments.fcst_var,arguments.ob_var,arguments.time,arguments.lev,nen)
plt.title('Correlation between %s (obs) and %s (fcst)'%(arguments.ob_var,arguments.fcst_var))
plt.savefig(outpath)
#plt.show()
 
