import numpy as np
import sys
import pdb
from time import strftime
from netCDF4 import Dataset, stringtoarr

__all__ = ['write_netcdf']

def write_netcdf(xargs, meta, dates, xvar):
    ''' Write a netCDF file ''' 
    
    # unpack args
    outFil = xargs[0]
    var = xargs[1]
    rcp = xargs[2]
    pattern = xargs[3]
    residual = xargs[4]
    quant = xargs[5]
    set = str(xargs[6])

    if len(dates[0]) == 7:
      units = { 'tas': 'K', 'tasmin': 'K','tasmax': 'K', 'pr': 'mm'}
      long_name = { 'tas': 'daily avg. temperature', 'tasmin': 'daily minimum temperature',
                    'tasmax': 'daily maximum temperature', 'pr': 'daily precipitation total'}
    elif len(dates[0]) == 6:
      units = { 'tas': 'K', 'tasmin': 'K','tasmax': 'K', 'pr': 'mm/d'}
      long_name = { 'tas': 'monthly avg. temperature', 'tasmin': 'monthly avg. daily minimum temperature',
                    'tasmax': 'monthly avg. daily maximum temperature', 'pr': 'monthly avg. daily precip. rate'}

    rootgrp = Dataset(outFil, 'w')

    # some global attributes
    rootgrp.history = 'Created '+strftime("%Y-%m-%d %H:%M:%S")
    rootgrp.author = 'Processed by DJ Rasmussen for Rhodium Group, LLC; email: '\
                                            'd.m.rasmussen.jr@gmail.com'
    rootgrp.address = 'Rhodium Group, LLC; 312 Clay St.; Oakland, Calif. 94607'
    rootgrp.copyright = 'Commercial use of this data with out the written consent of' \
                         ' Rhodium Group, LLC is prohibited; (c) Rhodium Group, LLC, 2014'
    rootgrp.experiment = rcp
    rootgrp.pattern = pattern
    rootgrp.residual = residual
    
    if len(dates[0]) == 7:
      rootgrp.frequency = 'daily'
      rootgrp.reference = 'Daily data from monthly using method from Wood et al. (2002) ' \
                                           '(section 2.3.2) from JGR-atmospheres'
    elif len(dates[0]) == 6:
      rootgrp.frequency = 'monthly'

    rootgrp.quantile = str(quant)
    rootgrp.set = set.zfill(3)
    rootgrp.contents = 'Anomalies at GHCN station level'
    
    # time dimension (UNLIMITED)
    rootgrp.createDimension('time', None)
    time = rootgrp.createVariable('time', np.dtype('int32'), dimensions=['time'], fill_value=-999)
    time.calendar = 'noleap'
    if len(dates[0]) == 7:
      time.units = 'yyyyddd'
    elif len(dates[0]) == 6:
      time.units = 'yyyymm'
 
    time[:] = dates
    del dates, time
    
    # county dimension
    nstn = xvar.shape[1]
    rootgrp.createDimension('station', nstn)
    
    # lat
    lat = rootgrp.createVariable('lat', np.dtype('float32'), dimensions=['station'], fill_value=-999.99)
    lat.units = 'degrees_north'
    lat.long_name = 'latitude of station'
    lat[:] = np.array(meta[1][:], dtype='float')
    del lat
    
    # lon
    lon = rootgrp.createVariable('lon', np.dtype('float32'), dimensions=['station'], fill_value=-999.99)
    lon.units = 'degrees_east'
    lon.long_name = ' longitude of station'
    lon[:] = np.array(meta[2][:], dtype='float')
    del lon

    # code
    code = rootgrp.createVariable('code', np.dtype('|S'), dimensions=['station'])
    code.long_name = ' station code'
    code[:] = np.array(meta[0][:], dtype='|S')
    del code

    # state
    state = rootgrp.createVariable('state', np.dtype('|S'), dimensions=['station'])
    state.long_name = ' station state'
    state[:] = np.array(meta[3][:], dtype='|S')
    del state

    # name
    name = rootgrp.createVariable('name', np.dtype('|S'), dimensions=['station'])
    name.long_name = ' station name'
    name[:] = np.array(meta[4][:], dtype='|S')
    del name

    xout = rootgrp.createVariable(var, np.dtype('float32'), dimensions=['time','station'], fill_value=-999.99)
    xout.actual_range = [np.min(xvar), np.max(xvar)]
    xout.domain_avg = np.mean(xvar)
    
    try:
       xout.long_name = long_name[var]
       xout.units = units[var]
    except KeyError:
       xout.long_name = 'missing'
       xout.units = 'missing'

    xout[:,:] = xvar
   # for (i,j), value in np.ndenumerate(xvar):
   #   xout[i,j] = xvar[i,j]

    del xout, xvar

    rootgrp.close()
    print '*** Success! Wrote %s ****' %outFil
