import argparse
import os
import numpy as np
from netCDF4 import Dataset
# Personal Libraries
import datetime
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('Agg')
from scipy import interpolate
from datetime import datetime

parser = argparse.ArgumentParser()
parser.add_argument("mem",type=int,help = 'Ensemble member to upscale')
parser.add_argument("res",type=int,help='Upscaled resolution (in m)')
parser.add_argument("--filter",type=int,default=-25000,help='Upscaled resolution (in m)')
parser.add_argument("--random",action='store_true',help='Random Perturbations')
arguments = parser.parse_args()

#--- Get the directories that cm1 output is stored
if arguments.filter > 0:
    if arguments.random:
       infile = '/work/jonathan.labriola/TurbPBLExp/restart_files_turbulent_PBL/rand/LES_%02d/cm1rst_000006.nc_%depsilon'%(arguments.mem,arguments.filter)
    else:
      infile = '/work/jonathan.labriola/TurbPBLExp/restart_files_turbulent_PBL/all/LES_%02d/cm1rst_000006.nc_%depsilon'%(arguments.mem,arguments.filter)
else:
   if arguments.random:
      infile = '/work/jonathan.labriola/TurbPBLExp/restart_files_turbulent_PBL/rand/LES_%02d/cm1rst_000006.nc'%(arguments.mem)
   else:
      infile = '/work/jonathan.labriola/TurbPBLExp/restart_files_turbulent_PBL/all/LES_%02d/cm1rst_000006.nc'%(arguments.mem)

#--- Create the directories used to output the upscaled information
if arguments.random:
   os.system('mkdir /work/jonathan.labriola/TurbPBLExp/obs_error/rand_%sm/LES_%02d/'%(str(arguments.res),arguments.mem))
   outfile = '/work/jonathan.labriola/TurbPBLExp/obs_error/rand_%sm/LES_%02d/cm1rst_000006.nc'%(str(arguments.res),arguments.mem)
else:
   os.system('mkdir /work/jonathan.labriola/TurbPBLExp/obs_error/all_%sm/LES_%02d/'%(str(arguments.res),arguments.mem)) 
   outfile = '/work/jonathan.labriola/TurbPBLExp/obs_error/all_%sm/LES_%02d/cm1rst_000006.nc'%(str(arguments.res),arguments.mem)

os.system('rm %s'%outfile)
dumpfile = Dataset(infile,"r",fortmat="format=")
wrtfile =  Dataset(outfile,"w",fortmat="NETCDF4")


nx = 800 
ny = 800 
dx = .25


dxnew = arguments.res/1000. #2.
nxnew = int(nx * (dx/dxnew))
nynew = int(ny * (dx/dxnew))


new_values={
'ni':  nxnew,
'nj':  nynew,
'nk': 120,
'nip1': nxnew+1,
'njp1': nynew+1,
'nkp1': 121,
'time': None,
'nbudget': 10,
'numq': 16,
}


def interpolate_domain(var,dims):
   """
   A Function used to interpolate data

   """
   now = datetime.now()
   current_time = now.strftime("%H:%M:%S")
   print("Current Time =", current_time)
   #xh Want to be between integer values i.e., -49.5 - 49.5
   #xhp1 Want to be on integer values i.e., -50, 50

   xh = np.arange(-((float(nx)*dx)/2.)+(dx/2.),((float(nx)*dx)/2.)+(dx/2.),dx)
   xh_new = np.arange(-((float(nxnew)*dxnew)/2.)+(dxnew/2.),((float(nxnew)*dxnew)/2.)+(dxnew/2.),dxnew)
   yh = np.arange(-((float(ny)*dx)/2.)+(dx/2.),((float(ny)*dx)/2.)+(dx/2.),dx)
   yh_new = np.arange(-((float(nynew)*dxnew)/2.)+(dxnew/2.),((float(nynew)*dxnew)/2.)+(dxnew/2.),dxnew)



   xhp1 = np.arange(-((float(nx)*dx)/2.),((float(nx+1.)*dx)/2.)+(dx/2.),dx)
   xhp1_new = np.arange(-((float(nxnew)*dxnew)/2.),((float(nxnew+1.)*dxnew)/2.)+(dxnew/2.),dxnew) 
   yhp1 = np.arange(-((float(ny)*dx)/2.),((float(ny+1.)*dx)/2.)+(dx/2.),dx)
   yhp1_new = np.arange(-((float(nynew)*dxnew)/2.),((float(nynew+1.)*dxnew)/2.)+(dxnew/2.),dxnew)  

   if 'time' in dims:
      var = np.squeeze(var[0])
   print('dims',dims)

   if ('ni' in dims or 'nip1' in dims) and ('nj' in dims  or 'njp1' in dims):
     
      if 'ni' in dims and 'nj' in dims:
        yy = yh
        xx = xh
        yynew = yh_new
        xxnew = xh_new

      elif 'ni' in dims and 'njp1'in dims:
        yy = yhp1
        xx = xh
        yynew = yhp1_new
        xxnew = xh_new

      elif 'nip1' in dims and 'nj' in dims:
        yy = yh
        xx = xhp1
        yynew = yh_new
        xxnew = xhp1_new

      if 'nk' in dims or 'nkp1' in dims:
         if 'nk' in dims:
             nz = new_values['nk']
         else:
             nz = new_values['nkp1']
         var_updt = np.zeros((nz,len(yynew),len(xxnew)))
         for index in range(0,nz): #--- Time to perform 2-D interpolation over different levels
            #f = interpolate.interp2d(yy, xx, var[index], kind='linear') 
            f = interpolate.RectBivariateSpline(yy, xx, var[index])
            var_updt[index] = f(yynew,xxnew)
      else: #--- Time to perform 2-D interpolation at 1 level
          f = interpolate.RectBivariateSpline(yy, xx, var)
          # JDL RectBivariateSpline is optimal, don't need mesh grid just need x and y -components
          var_updt = f(yynew,xxnew)


   elif 'ni'in dims or 'nip1' in dims or 'nj'in dims or 'njp1' in dims: #--- 2D and 1D interpolation
      print('HEY JDL dims = ',dims)
      if 'ni' in dims:
        xx = xh
        xxnew = xh_new
      elif 'nip1' in dims:
         xx = xhp1
         xxnew = xhp1_new
      elif 'nj' in dims: 
        xx = yh
        xxnew = yh_new
      elif 'njp1' in dims:
        xx = yhp1
        xxnew = yhp1_new 


      if 'nk' in dims or 'nkp1' in dims: #--- 1D interpolation layer by layer
         if 'nk' in dims:
           nz = new_values['nk']
         else:
           nz = new_values['nkp1']
         var_updt = np.zeros((nz,len(xxnew)))
         for index in range(0,nz):
           f = interpolate.InterpolatedUnivariateSpline(xx, var[index])
           var_updt[index] = f(xxnew)
      else: #--- 1D interpolation on one layer:
         f = interpolate.InterpolatedUnivariateSpline(xx, var)
         var_updt = f(xxnew)
   else:
      var_updt = var

  
   return var_updt

#---- Creating the dimensions you want
for value in new_values:
   wrtfile.createDimension(value, new_values[value])

#--- Creating nc attributes
attrs = dumpfile.ncattrs()
for nx_attr in attrs:
   print('%s = %s'%(nx_attr,dumpfile.getncattr(nx_attr)))
   if nx_attr == 'nx':
      wrtfile.setncattr(nx_attr,np.int32(new_values['ni']))
   elif nx_attr == 'ny':
      wrtfile.setncattr(nx_attr,np.int32(new_values['nj'])) #--- attribute is int32
   else:
      wrtfile.setncattr(nx_attr,dumpfile.getncattr(nx_attr))
      
#--- Creating the Variables
variables = dumpfile.variables
for var in variables:
   dims = []
   for dim in dumpfile.variables[var].dimensions:
      dims.append(dim)

   if len(dims) == 1:
      wrtfile.createVariable(var,dumpfile.variables[var].datatype,(dims[0]))
   elif len(dims) == 2:
      wrtfile.createVariable(var,dumpfile.variables[var].datatype,(dims[0],dims[1],))
   elif len(dims) == 3:
      wrtfile.createVariable(var,dumpfile.variables[var].datatype,(dims[0],dims[1],dims[2],))
   elif len(dims) == 4:
      wrtfile.createVariable(var,dumpfile.variables[var].datatype,(dims[0],dims[1],dims[2],dims[3],))
  
   try:
      wrtfile.variables[var].long_name = dumpfile.variables[var].long_name
      wrtfile.variables[var].units = dumpfile.variables[var].units
   except:
      pass

   #if 'time' in dims:
   #   var_tmp = np.squeeze(dumpfile.variables[var][0])
   #else: 
   var_tmp = dumpfile.variables[var][:]
   var_tmp = interpolate_domain(var_tmp,dims)

   #if var == 'umove':
   #   var_tmp = 17.0
   #if var == 'vmove':
   #   var_tmp = 22.0

   #for grid_dim in ['ni','nip1','nj','njp1']:
   #   if grid_dim  in dims: 
   #     index = dims.index(grid_dim) 
   #     print()
   #     print(index)
   #     print(grid_dim)
   #     print(dims)
   #     print('Before transform shape = ',var_tmp.shape) 
   #     var_tmp = contract_domain(var_tmp,grid_dim,index)      
   #     print('After transform shape = ',var_tmp.shape)  
   
   if dims[0] == 'time':
      try:
         wrtfile.variables[var][0] = var_tmp[:] #np.squeeze(dumpfile.variables[var][0])
      except:
         wrtfile.variables[var][0] = var_tmp
   else: 
     wrtfile.variables[var][:] = var_tmp #dumpfile.variables[var]
   print(wrtfile.variables[var].shape)
#--- Looking to see how shifted gridpoints work 
#xh = np.squeeze(dumpfile.variables['xh'][:])
#xf = np.squeeze(dumpfile.variables['xf'][:])
#print('xf = ',xf[0:3],xf[-3:])
#print('xh = ',xh[0:3],xh[-3:])
wrtfile.close()
