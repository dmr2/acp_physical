#!/usr/local/bin/python
"""
write_cdf_table.py 

@author rasmussen; last updated: Mon 01 Feb 2016 06:29:35 PM PST

DESCRIPTION: Extract every quantile from CDF for every county in US, state, and NCA reg. 

"""

import sys
import numpy as np
import acp_tools


def find_nearest(array,value):
    return int((array <= value).argmin())


# USER EDITS START

kind = 'CMIP5'  # SMME | CMIP5 | MCPR

# VARIABLE TO PROCESS
var = "tasmin" # pr | tas | tasmin | tasmax | wmax

# AVERAGING WEIGHT
weight = "POP" # AREA | POP

# VARIABLE UNIT
unit = "C" # C | K | mm | in | F

# EXCEEDANCE LABEL
label = ["lt0"]
#label = ["annual"]
#label = ["DJF","MAM","JJA","SON"]

# NUMBER OF COUNTIES TO PROCESS
ncnty = 3109

want_anomaly = True # for temperature variables only
want_percent = True # for precipitation only
want22nd = False # processing 22nd century results

# ROOT OUTPUT DIRECTORY
root = "/home/dmr/jnk"
path_in = root + "/cnty_txt"

# SMME MODEL PROBABILITY WEIGHTS DIRECTORY
mod_wdir = "/home/dmr/acp_code/preproc/magicc/modelweights"

# USER EDITS END


rcp_list = ["rcp26","rcp45","rcp60","rcp85"]

if want22nd: 
  yr1 = [2020,2040,2080,2100,2120,2140,2180]
  yr2 = [2039,2059,2099,2119,2139,2159,2199]
else:
  yr1 = [2020,2040,2080]
  yr2 = [2039,2059,2099]

if label[0] == 'annual':
  want_annual = True # for annual averages (doesn't work for exceedances)
else:
  want_annual = False

if kind == 'SMME':
  # Do we want the pattern scaled projections?
  want_pattern = True
else:
  want_pattern = False

quants = [i*.01 for i in xrange(1,100)]

season = label

reg_list = ["CONUS","NEA","SEA","UGP","MWE","NWE","SWE",
               "MTN","SCL","CAL"]

states = ["AL", "AZ", "AR", "CA", "CO", "CT", "DC", "DE", "FL", "GA", 
          "ID", "IL", "IN", "IA", "KS", "KY", "LA", "ME", "MD",  
          "MA", "MI", "MN", "MS", "MO", "MT", "NE", "NV", "NH", "NJ", 
          "NM", "NY", "NC", "ND", "OH", "OK", "OR", "PA", "RI", "SC", 
          "SD", "TN", "TX", "UT", "VT", "VA", "WA", "WV", "WI", "WY"]


rnorm = np.zeros((len(season), len(reg_list)), dtype=float)
snorm = np.zeros((len(season), len(states)), dtype=float)

qcnty = np.zeros((len(season), len(rcp_list), ncnty, len(quants), len(yr1)), dtype=float)
qcnty.fill(-999.999)
qstate = np.zeros((len(season), len(rcp_list), len(states), len(quants), len(yr1)), dtype=float)
qstate.fill(-999.999)
qreg = np.zeros((len(season), len(rcp_list), len(reg_list), len(quants), len(yr1)), dtype=float)
qreg.fill(-999.999)

for ircp,rcp in enumerate(rcp_list):

    # LOAD MODEL PROBABILITY WEIGHTS AND LIST OF MODELS/ SURROGATES
    prob_weight, models = acp_tools.get_mod_weights(rcp, kind, mod_wdir)

    for iy in xrange(0,len(yr1)):

        cmod = np.zeros((len(season), ncnty, len(models)), dtype=float)
        cmod.fill(-999.99)

        norm = np.zeros((len(season), ncnty), dtype=float)
        norm.fill(-999.99)

        imod = 0
        for imod,mod in enumerate(models):
            if mod[0:5] =='patte':
              if "DJF" in season or "annual" in season:
               fil = path_in+'/smme_county_seasonal_'+var+'_'+rcp+'_'+mod.rsplit('_')[0]+ \
                                       '_'+str(yr1[iy])+'-'+str(yr2[iy])+'.csv'
              else:
               fil = path_in+'/smme_county_'+var+'_num_day_'+season[0]+'_'+mod.rsplit('_')[0]+ \
                                       '_'+rcp+'_'+str(yr1[iy])+'-'+str(yr2[iy])+'.csv'
            else:
              if kind == 'MCPR':
               if "DJF" in season or "annual" in season:
                fil = path_in+'/mcpr_county_seasonal_'+var+'_'+rcp+'_'+mod.rsplit('_')[0]+ \
                                       '_'+str(yr1[iy])+'-'+str(yr2[iy])+'.csv'
               else:
                fil = path_in+'/mcpr_county_'+var+'_num_day_'+season[0]+'_'+mod+ \
                                     '_'+rcp+'_'+str(yr1[iy])+'-'+str(yr2[iy])+'.csv'
	      else:
               if "DJF" in season or "annual" in season:
                fil = path_in+'/smme_county_seasonal_'+var+'_'+rcp+'_'+mod.rsplit('_')[0]+ \
                                       '_'+str(yr1[iy])+'-'+str(yr2[iy])+'.csv'
               else:
                fil = path_in+'/smme_county_'+var+'_num_day_'+season[0]+'_'+mod+ \
                                       '_'+rcp+'_'+str(yr1[iy])+'-'+str(yr2[iy])+'.csv'

            print "Opening file: %s " %fil
            data = np.genfromtxt(fil,delimiter=',',dtype="|S",autostrip=True,skip_header=1)

            fips = np.array(data[:,0], dtype=int) # FIPS
            counties = np.array(data[:,1], dtype="|S") # County Name
            cstate = np.array(data[:,2], dtype="|S") # County State
            lat = np.array(data[:,3], dtype=float) # County Lat
            lon = np.array(data[:,4], dtype=float) # County Lon
            area = np.array(data[:,5], dtype=float) # County Area
            ncnty = area.size
            
            if want_annual and "num_day" not in fil:
               for iseason in xrange(0, len(season)):
                if var == "pr":
                  cmod[iseason,:,imod] = np.sum(np.array(data[:,11:16], dtype=float),axis=1) # Model
                  norm[iseason,:] = np.sum(np.array(data[:,7:11], dtype=float),axis=1) # Normals
                else:
                  cmod[iseason,:,imod] = np.mean(np.array(data[:,11:16], dtype=float),axis=1) # Model
                  norm[iseason,:] = np.mean(np.array(data[:,7:11], dtype=float),axis=1) # Normals
            else:
              if "num_day" in fil:
               for iseason in xrange(0, len(season)):
                   cmod[iseason,:,imod] = np.array(data[:,8+iseason], dtype=float) # Model
                   norm[iseason,:] = np.array(data[:,7+iseason], dtype=float) # Normals
              else:
               for iseason in xrange(0, len(season)):
                   cmod[iseason,:,imod] = np.array(data[:,11+iseason], dtype=float) # Model
                   norm[iseason,:] = np.array(data[:,7+iseason], dtype=float) # Normals

        if "num_day" not in fil:
          if unit == "F":
              cmod *= 9/5.
              norm = norm*9/5. - 459.67
          elif unit == "C":
              norm = norm - 273.15
          elif unit == "in":
              cmod /= 25.4 # mm to inches
              norm /= 25.4 # mm to inches
        else:
          unit = "days"

        # AGGREGATE TO REGIONAL LEVEL
        rmod = np.zeros((len(season), len(reg_list), len(models)), dtype=float)
        rmod.fill(-999.99)

        rnorm = np.zeros((len(season), len(reg_list)), dtype=float)
        rnorm.fill(-999.99)

        for ire,reg in enumerate(reg_list):

            reg_dict = acp_tools.region_defs(reg)
            weights = np.array(reg_dict[weight], dtype=float)
            regfips = map(int,reg_dict['FIPS'])

            temp = np.zeros((len(season), len(regfips), len(models)), dtype=float)
            temp.fill(-999.99)

            ntemp = np.zeros((len(season), len(regfips)), dtype=float)
            ntemp.fill(-999.99)

            print "Pulling out the counties that are in %s..." %reg
            for i,ifip in enumerate(regfips):
                for imod in xrange(0, len(models)):
                    indx = np.where(ifip == fips)[0][0]
                    for iseason in xrange(0, len(season)):
                        temp[iseason,i,imod] = cmod[iseason,indx,imod]

                for iseason in xrange(0, len(season)):
                    ntemp[iseason,i] = norm[iseason,indx]

            # WEIGHT BY POP OR AREA
            for imod in xrange(0, len(models)):
                for iseason in xrange(0, len(season)):
                    rmod[iseason,ire,imod] = np.average(temp[iseason,:,imod], weights=weights)

            for iseason in xrange(0, len(season)):
                rnorm[iseason,ire] = np.average(ntemp[iseason,:], weights=weights)

        # AGGREGATE TO STATE LEVEL
        smod = np.zeros((len(season),len(states), len(models)), dtype=float)
        smod.fill(-999.99)
        
        snorm = np.zeros((len(season),len(states)), dtype=float)
        snorm.fill(-999.99)

        reg_dict = acp_tools.region_defs('CONUS')
        weights = np.array(reg_dict[weight], dtype=float)

        for istate,state in enumerate(states):
            print "Aggregating county level data for %s..." % state

            indx = np.where(state == cstate)[0][:]
            temp = np.zeros((len(season), len(indx), len(models)), dtype=float)
            temp.fill(-999.99)
           
            ntemp = np.zeros((len(season), len(indx)), dtype=float)
            ntemp.fill(-999.99)

            for imod in xrange(0, len(models)):
                for iseason in xrange(0, len(season)):
                    temp[iseason,:,imod] = cmod[iseason,indx,imod]

            for ii,i in enumerate(indx):
                for iseason in xrange(0, len(season)):
                    ntemp[iseason,ii] = norm[iseason,i]

            # WEIGHT BY POP OR AREA
            for imod in xrange(0, len(models)):
                for iseason in xrange(0, len(season)):
                    smod[iseason,istate,imod] = np.average(temp[iseason,:,imod], weights=weights[indx])

            for iseason in xrange(0, len(season)):
                snorm[iseason,istate] = np.average(ntemp[iseason,:], weights=weights[indx])

        # DETERMINE QUANTILES

        # REGION
        for ire,reg in enumerate(reg_list):
            print "Determining the quants we want for %s" % reg

            for iq,iquant in enumerate(quants):
                for iseason in xrange(0, len(season)):
                    # sort from lowest to highest and return permutation vector only
                    ind_sort = np.argsort(rmod[iseason,ire,:])

                    sortmod = []
                    for i in range(len(models)): 
                        sortmod.append(models[ind_sort[i]]) # sort models 

                    weightsort = prob_weight[ind_sort]
                    weightsort = weightsort/sum(weightsort) # re-normalize

                    wq = np.cumsum(weightsort, dtype=float)  # calculate probability weight factor
                    imod = models.index(sortmod[find_nearest(wq, iquant)]) #ordered to unordered
                    if var != 'pr':
                      if want_anomaly:
                        qreg[iseason,ircp,ire,iq,iy] = rmod[iseason,ire,imod]
                      else:
                        qreg[iseason,ircp,ire,iq,iy] = rmod[iseason,ire,imod] + \
                                                        rnorm[iseason,ire]
                    elif var == 'pr':
                      if not want_percent:
                          qreg[iseason,ircp,ire,iq,iy] = rmod[iseason,ire,imod]
                      else:
                          qreg[iseason,ircp,ire,iq,iy] = (rmod[iseason,ire,imod] - \
                                                            rnorm[iseason,ire])/(rnorm[iseason,ire]) * 100.
        # STATE
        for istate,state in enumerate(states):
            print "Determining the quants we want for %s" % state

            for iq,iquant in enumerate(quants):
                for iseason in xrange(0, len(season)):
                    ind_sort = np.argsort(smod[iseason,istate,:])

                    sortmod = [ ]
                    for i in range(len(models)): 
                        sortmod.append(models[ind_sort[i]])

                    weightsort = prob_weight[ind_sort]
                    weightsort = weightsort/sum(weightsort)

                    wq = np.cumsum(weightsort, dtype=float)
                    imod = models.index(sortmod[find_nearest(wq, iquant)])

                    if var != 'pr':
                      if want_anomaly:
                        qstate[iseason,ircp,istate,iq,iy] = smod[iseason,istate,imod]
                      else:
                        qstate[iseason,ircp,istate,iq,iy] = smod[iseason,istate,imod] + snorm[iseason,istate]
                    elif var == 'pr':
                      if not want_percent:
                          qstate[iseason,ircp,istate,iq,iy] = smod[iseason,istate,imod]
                      else:
                          qstate[iseason,ircp,istate,iq,iy] = (smod[iseason,istate,imod] - \
                                                            snorm[iseason,istate])/ snorm[iseason,istate] * 100.

     
        # COUNTY
        for icnty,county in enumerate(counties):
            #print "Determining the quants we want for counties in state: %s" %cstate[icnty]
            sys.stdout.write("\rDetermining the quants we want for CONUS counties: %d/%d" %(icnty+1,ncnty))
            sys.stdout.flush()
            for iq,iquant in enumerate(quants):
                for iseason in xrange(0, len(season)):
                    ind_sort = np.argsort(cmod[iseason,icnty,:])

                    sortmod = []
                    for i in range(len(models)): 
                        sortmod.append(models[ind_sort[i]])

                    weightsort = prob_weight[ind_sort]
                    weightsort = weightsort/sum(weightsort)

                    wq = np.cumsum(weightsort, dtype=float)
                    imod = models.index(sortmod[find_nearest(wq, iquant)])
                    if var != 'pr':
                      if want_anomaly:
                        qcnty[iseason,ircp,icnty,iq,iy] = cmod[iseason,icnty,imod]
                      else:
                        qcnty[iseason,ircp,icnty,iq,iy] = cmod[iseason,icnty,imod] + norm[iseason,icnty]
                    elif var == 'pr':
                      if not want_percent:
                        qcnty[iseason,ircp,icnty,iq,iy] = cmod[iseason,icnty,imod]
                      else:
                        qcnty[iseason,ircp,icnty,iq,iy] = (cmod[iseason,icnty,imod] - \
                                                              norm[iseason,icnty])/ norm[iseason,icnty] * 100.

    # WRITE TO TEXT FILES

    cdfwriter = acp_tools.CDF_write(root,var=var,kind=kind,rcp=rcp,\
                unit=unit, percent=want_percent, season=season, \
                quants=quants, weight=weight)

    cdfwriter.region(reg_list, rnorm, qreg)
    cdfwriter.state(states, snorm, qstate)
    cdfwriter.county(counties, cstate, fips, lat, lon, norm, qcnty)
