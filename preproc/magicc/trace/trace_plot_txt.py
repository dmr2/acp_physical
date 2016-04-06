#!/usr/bin/python
'''
@author rasmussen
last updated Sun Mar 13 11:54:29 2016

Generate global temperature percentiles as a function of MAGICC model runs

'''

import pandas as pd
import numpy as np
import csv
import sys

class prettyfloat(float):
    def __repr__(self):
        return "%0.2f" % self

fil_list = ["IPCCAR5climsens_rcp85_DAT_SURFACE_TEMP_BO_15Nov2013_185227.OUT", \
            "IPCCAR5climsens_rcp3pd_DAT_SURFACE_TEMP_BO_16Nov2013_085858.OUT",\
            "IPCCAR5climsens_rcp45_DAT_SURFACE_TEMP_BO_16Nov2013_070923.OUT", \
            "IPCCAR5climsens_rcp6_DAT_SURFACE_TEMP_BO_16Nov2013_064508.OUT"]

diri = "MAGICC_CMIP"

head = ["YEAR"] + [str(run+1).zfill(3) for run in range(600)]

year_list = [[1981,2010],[2020,2039],[2040,2059],[2080,2099],[2140,2159],[2180,2199]]

perc_list = [5,17,50,83,95,99,99.9]

vals = np.zeros((len(head)-1,len(year_list)),dtype=float)
p = np.zeros((len(head)-1,len(year_list),len(perc_list)),dtype=float)

# EACH RCP SCENARIO
for fil in fil_list:

  df = pd.read_csv(diri+"/"+fil,skiprows=23,sep=r"\s*")
  df.columns = head

  # EACH MAGICC RUN
  for i,run in enumerate(head[1:len(head)]):
    print "Crunching MAGICC run ... %s" %(run)
     
    base_line = []
    # EACH AVERAGING PERIOD
    for j,pair in enumerate(year_list):
      tmp = []
      for yr in xrange(pair[0],pair[1]+1):
       tmp.append(df[df["YEAR"] == yr][run].values[0])

      if j == 0:
       base_line =  np.mean(tmp)
      else:
       vals[i,j] =  np.mean(tmp) - base_line

      for k,perc in enumerate(perc_list):
        p[i,j,k] = np.percentile(vals[0:i+1,j],perc)


  # WRITE OUT AS TEXT FILE
  for j,pair in enumerate(year_list):
    with open(fil+"_"+str(pair[0])+"-"+str(pair[1])+"_perc", 'wb') as outf:
      writer = csv.writer(outf,delimiter="\t")
      writer.writerow(["MAGICC runs"] + ["q"+str(perc) for perc in perc_list])
      for i,run in enumerate(head[1:len(head)]):
       writer.writerow([i+1]+map(prettyfloat,p[i,j,:]))

