from netCDF4 import Dataset
#import readvar
import post.io 
import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser()
#parser.add_argument("mem",type=int,help = 'Ensemble member to upscale')
parser.add_argument("nen",type=int,help = 'N_Ensemble Members')
parser.add_argument("res",type=int,default=2000,help = 'Forecast Resolution (m)')
parser.add_argument("var",type=str,help = 'The analyzed variable')
parser.add_argument("time",type=int,help='Model Time')
parser.add_argument("--filter",action='store_true',help='Only evaluate Filter Forecasts')
parser.add_argument("--random",action='store_true',help='Only evaluate Random Forecasts')
arguments = parser.parse_args()
eval_hgts = [1.,2.,3.]

if arguments.filter:  #--- Only compare filtered observations against real observations
   if arguments.res == 2000:
      epsilon = 25000
      nskip = 5 #8 
      xmin = 400
      xmax = 800
      ymin = 0 
      ymax = 800

else:
   if arguments.res == 2000:
      xmin = 50 # km
      xmax = 100 # km
      ymin = 0 # km
      ymax = 100 # km
      nskip = 5
   elif arguments.res == 4000:
      xmin = 25
      xmax = 50
      ymin = 0
      ymax = 50
      nskip = 1  
   elif arguments.res == 1000:
      xmin = 100
      xmax = 200
      ymin = 0
      ymax = 200
      nskip = 1
   elif arguments.res == 500:
      xmin = 200
      xmax = 400
      ymin = 0
      ymax = 400
      nskip = 1



#nprofiles = (xmax-xmin) * (ymax-ymin)
nprofiles = int(np.ceil(float(xmax-xmin)/float(nskip)) * np.ceil(float(ymax-ymin)/float(nskip)))
print('The number of profiles = ',str(nprofiles))
nlevs = len(eval_hgts)
sim_profiles = np.zeros((nlevs,arguments.nen,nprofiles))
obs_profiles = np.zeros((nlevs,arguments.nen,nprofiles))
for mem in range(1,arguments.nen+1):
   print('Working on Member %s'%str(mem))
   #--- observation is considered the "high-res" run
   if arguments.random:
      obs_path = '/work/jonathan.labriola/TurbPBLExp/rand/LES_%02d/cm1out_%06d.nc'%(mem,arguments.time)
   else:
      obs_path = '/work/jonathan.labriola/TurbPBLExp/all/LES_%02d/cm1out_%06d.nc'%(mem,arguments.time)
   obsfile = Dataset(obs_path,"r",fortmat="NETCDF4")
   obsgrid = post.io.grdbas.readGrdbas(obsfile)
   for hindex, hgt in enumerate(eval_hgts):
      if hindex == 0:
         tmp = post.io.cm1.readvar(obsfile,arguments.var,interp_hgt=hgt)
         obs = np.zeros((len(eval_hgts),tmp.shape[0],tmp.shape[1]))
         obs[hindex] = tmp
      else:
         obs[hindex] = post.io.cm1.readvar(obsfile,arguments.var,interp_hgt=hgt)

   #--- Simulations are on the "Coarse Grids"
   if arguments.filter:
      if arguments.random:
         sim_path = '/work/jonathan.labriola/TurbPBLExp/rand/LES_%02d/cm1out_%06d.nc_%sepsilon'%(mem,arguments.time,str(epsilon))
      else:
         sim_path = '/work/jonathan.labriola/TurbPBLExp/all/LES_%02d/cm1out_%06d.nc_%sepsilon'%(mem,arguments.time,str(epsilon))
   else:
      if arguments.random:
         sim_path = '/scratch/jonathan.labriola/idealized/250m/LES_Injection_Turbulent/rand_%sm/LES_%02d/cm1out_%06d.nc'%(str(arguments.res),mem,arguments.time)
      else:
         #sim_path = '/scratch/jonathan.labriola/idealized/250m/LES_Injection_Turbulent/all_%sm/LES_%02d/cm1out_%06d.nc'%(str(arguments.res),mem,arguments.time)
         sim_path = '/work/jonathan.labriola/TurbPBLExp/2km_runs/LES_%02d/cm1out_%06d.nc'%(mem,arguments.time)
   simfile = Dataset(sim_path,"r",fortmat="NETCDF4")
   simgrid = post.io.grdbas.readGrdbas(simfile)
   for hindex, hgt in enumerate(eval_hgts):
      if hindex == 0:
         tmp = post.io.cm1.readvar(simfile,arguments.var,interp_hgt=hgt)
         sim = np.zeros((len(eval_hgts),tmp.shape[0],tmp.shape[1]))
         sim[hindex] = tmp
      else:
         sim[hindex] = post.io.cm1.readvar(simfile,arguments.var,interp_hgt=hgt)

   for hindex, hgt in enumerate(eval_hgts):
      #--- Wind Speed Biases
      f = interpolate.RectBivariateSpline(obsgrid.yh, obsgrid.xh, obs[hindex])
      obs_interp = f(simgrid.yh[ymin:ymax:nskip],simgrid.xh[xmin:xmax:nskip])
      sim_profiles[hindex,mem-1] = sim[hindex,ymin:ymax:nskip,xmin:xmax:nskip].flatten()
      obs_profiles[hindex,mem-1] = obs_interp.flatten()

sim_profiles_collapse = np.reshape(sim_profiles,(nlevs,arguments.nen*nprofiles))
obs_profiles_collapse = np.reshape(obs_profiles,(nlevs,arguments.nen*nprofiles))
obs_diff = sim_profiles_collapse - obs_profiles_collapse

#--- Plotting Errors on a Histogram
if arguments.var in ['wspd','u','ua','v','va']:
   bins = np.arange(-4.,4.2,0.2)
   xmin_plot = -4.
   xmax_plot = 4.
   ymin_plot = 0
   ymax_plot = 2
elif arguments.var == 'qv':
   bins = np.arange(-1.5,1.53,0.075)
   xmin_plot = -1.5
   xmax_plot = 1.5
   ymin_plot = 0
   ymax_plot = 6.5
else:
   bins = np.arange(-1.5,1.53,0.075)
   xmin_plot = -1.5
   xmax_plot = 1.5
   ymin_plot = 0
   ymax_plot = 6.5


for hindex, hgt in enumerate(eval_hgts):
   if hindex == 0: color = 'goldenrod'
   elif hindex == 1: color = 'forestgreen'
   elif hindex == 2: color = 'royalblue'
   plt.hist(obs_diff[hindex],bins=bins,color=color,histtype='step',density='True',linewidth=3.0,label='%1.2f km (STD = %1.2f)'%(hgt,np.std(obs_diff[hindex])))
plt.xlim([xmin_plot,xmax_plot])
plt.ylim([ymin_plot,ymax_plot])
plt.legend()


title_base = 'profile_differences_nen%02d_%03d_%s_%06dm'%(arguments.nen,arguments.time,str(arguments.var),arguments.res)
if arguments.filter:
   outpath='%s_filter'%(title_base)
else:
   outpath='%s_fcst'%(title_base)
if arguments.random:
   outpath = outpath+'_rand'
else:
   outpath = outpath+'_all'

plt.savefig('./%s.png'%(outpath), figsize = (13, 13), dpi=300)
#plt.semilogy()
#plt.savefig('./%s_log.png'%(outpath), figsize = (13, 13), dpi=300)





