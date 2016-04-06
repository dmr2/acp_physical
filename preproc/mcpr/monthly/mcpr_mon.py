# mcpr_mon.py

# rasmussen; last updated: Fri 27 Nov 2015 09:14:36 AM PST

# generates GHCN station-level monthly data using the MCPR method

import sys
import glob
import pdb
import numpy as np
import os
from open import *
from temporal import *
from writeNetCDF import write_netcdf

# Use only these forcing scenarios
rcp_list = ['rcp26','rcp45','rcp60','rcp85']

# These variables only! 
var_list = ['tas','pr','tasmax','tasmin']

monthly = True
checkMonthly = True

# for what range of years will use 
yrStart = 1981 # must be .le. to 1981
yrEnd = 2200

# obs years start end
obYrStr = 1981
obYrEnd = 2010
ObYrs = []; ObYrs.append(obYrStr); ObYrs.append(obYrEnd)
 
if monthly:
  # generate yyyymm time record
  yyyymm = monthDates(yrStart, yrEnd)

# location of patterns and residuals
patt_dir = '/home/dmr/jnk/patterns'
resid_dir = '/home/dmr/jnk/residuals'


# Extract pattern and residual values at each station

for rcp in rcp_list:

    mcpr_fil = "/home/dmr/acp_code/preproc/magicc/modelweights/"+rcp+"_MCPR.tsv"
    print "Opening list of pattern weights: %s\n" % mcpr_fil 
    try:
      data = np.genfromtxt(mcpr_fil, delimiter='\t', dtype="|S", autostrip=True, skip_header=0)
    except:
      sys.exit("Sorry, I cannot find %s" %mcpr_fil)
    
    istr = (yrStart - 1950) + 3 
    iend = (yrEnd - 1950) + 1 + 3

    quant = np.array(data[1:,0], dtype=np.float)
    pattern_list = np.array(data[1:,1]).tolist()
    residual_list = np.array(data[1:,2]).tolist()
    years = np.array(data[0,istr:iend]).tolist()
    global_tas = np.array(data[1:,istr:iend], dtype=float)
    wxStrt = (2200 - 1950) + 1 + 3
    wxYrs = np.array(data[1:,wxStrt:], dtype=int)
    
    print 'Going to create daily data from %s to %s \n' %(years[0], years[-1])
   
    for var in var_list:
        # Map GHCN stations to GCM grid cells
        meta = []
        
        if var == 'tas' or var == 'tasmax' or var == 'tasmin':
           map_fil = 'tas_conus_bcsd_xy_map.csv'
        else: 
           map_fil = 'pr_conus_bcsd_xy_map.csv'
        
        print "Opening list of CONUS stations: %s\n" % map_fil
        data = np.genfromtxt(map_fil, delimiter=',', dtype="|S", skip_header=1)
        print np.shape(data)
        meta.append(np.array(data[:,0]).tolist()) # code
        meta.append(np.array(data[:,1]).tolist()) # lat
        meta.append(np.array(data[:,2]).tolist()) # lon
        meta.append(np.array(data[:,4]).tolist()) # state
        meta.append(np.array(data[:,5]).tolist()) # name
        xmap = np.array(data[:,6], dtype=int) # longitude
        ymap = np.array(data[:,7], dtype=int) # latitude

        for i,pattern in enumerate(pattern_list):
            if monthly and checkMonthly:
              outDir = '/home/dmr/jnk/mcpr/001/' + rcp + '/'+str(i+1).zfill(3)+'/' + var 
              testFil = outDir +'/'+var+'_mcpr-'+str(i+1).zfill(3)+'_mon_ghcn_'+rcp+'_'+str(yrStart)+'01-'+str(yrEnd)+'12.nc'
              if os.path.isfile(testFil):
                  continue
  
            print 'Working on set: %s' %(i + 1)

            # get everything we need to generate monthlies
            slope_stn, resid_stn, rYrs = get_pattresid(var, rcp, \
                                                       pattern, residual_list[i], \
                                                       patt_dir, resid_dir, xmap, ymap)
            # generate monthly anomalies through pattern scaling
            yr = []
            yr.append(int(yrStart)); yr.append(int(yrEnd))
            modMonthly = get_monthly(yr, rYrs, slope_stn, resid_stn,\
                                          global_tas[i,:])
           
            if monthly:
                outDir = '/home/dmr/jnk/mcpr/001/' + rcp + '/'+str(i+1).zfill(3)+'/' + var 
                testFil = outDir +'/'+var+'_mcpr-'+str(i+1).zfill(3)+'_mon_ghcn_'+rcp+'_'+str(yrStart)+'01-'+str(yrEnd)+'12.nc'
                if not os.path.exists(outDir):
                    os.makedirs(outDir)

                name = var+'_mcpr-'+str(i+1).zfill(3)+'_mon_ghcn_'+rcp+'_'+str(yrStart)+'01-'+str(yrEnd)+'12.nc'
                
                outFil = outDir +'/'+ name
                args = [outFil, var, rcp, pattern, residual_list[i], quant[i], i + 1]
                write_netcdf(args, meta, yyyymm, modMonthly)
