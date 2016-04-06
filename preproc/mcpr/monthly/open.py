import numpy as np
from netCDF4 import Dataset
import glob 
import sys

funcs = ['get_pattresid']
__all__ = funcs

def get_pattresid(var, rcp, pattern, residual, patt_dir, resid_dir, xmap, ymap):
    ''' Get pattern and residual '''

    # OPEN PATTERN
 
    # Note: No gridded projections avilable for NORESM1-ME for Tmin/Tmax,
    # so just use average temperature instead.

    if var == "tasmin" or var == "tasmax" and pattern == "noresm1-me" or pattern == "access1-3":
      name = 'tas/seasonal_pattern_tas_bcsd_'+rcp+'_'+pattern+'.nc' # use avg. temperature
    else:
      name = var+'/seasonal_pattern_'+var+'_bcsd_'+rcp+'_'+pattern+'.nc'

    try:
        inf_patt = glob.glob(patt_dir+'/'+name)[0]
        try:
            print 'Opening %s' %inf_patt
            rootgrp = Dataset(inf_patt, 'r', format='NETCDF4_CLASSIC')
        except:
            sys.exit('Sorry, I cannot open: %s' %inf_patt)
        slope_grd = rootgrp.variables['slope'][:,:,:] # (season, lat, lon)

        nstn = len(ymap)
        slope_stn = np.zeros((4,nstn), dtype=np.float)
        slope_stn.fill(-999.99)
        for i in range(nstn):
         slope_stn[:,i] = slope_grd[:,ymap[i],xmap[i]]
    except IndexError:
        print 'I cannot open: '+patt_dir+'/'+var+'/'+name
        sys.exit(0)

    if var == "tasmin" or var == "tasmax" and residual == "noresm1-me" or residual == "access1-3":
      name = 'tas/tas_mon_residual_bcsd_'+rcp+'_'+residual+'_*19*-2*.nc' # use avg. temperature
    else:
      name = var+'/'+var+'_mon_residual_bcsd_'+rcp+'_'+residual+'_*19*-2*.nc'

    try:
        inf_resid = glob.glob(resid_dir+'/'+name)[0]
        try:
            print 'Opening %s' %inf_resid
            rootgrp = Dataset(inf_resid, 'r', format='NETCDF4_CLASSIC')
        except:
            sys.exit('Sorry, I cannot open: %s' %inf_resid)
        resid_grd = rootgrp.variables['residual'][:,:,:] # (time, lat, lon)

        nmon = np.shape(resid_grd)[0]
        resid_stn = np.zeros((nmon,nstn), dtype=np.float)
        resid_stn.fill(-999.99)
        for i in range(nstn):
         resid_stn[:,i] = resid_grd[:,ymap[i],xmap[i]]

        rootgrp.close()

        yrs = inf_resid.split('_')[-1].strip('.nc').split('-')
        rYr = []
        rYr.append(int(yrs[0][0:4])) # strip start and end month
        rYr.append(int(yrs[1][0:4]))
    except IndexError:
        print 'I cannot open: '+resid_dir+'/'+var+'/'+name
        sys.exit(0)

    return slope_stn, resid_stn, rYr

