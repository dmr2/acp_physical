import numpy as np
import sys

__all__ = ['get_monthly','date','monthDates']

def get_monthly(yr, Ryr, patt, resid, globT):
    ''' Get monthly averages from patts and resids '''

    nyr = ((yr[1] - yr[0]) + 1) # all years
    nryr = ((Ryr[1] - Ryr[0]) + 1) # residual years
    nstn = patt.shape[1]
    resid = np.reshape(resid, (nryr, 12, nstn), order='C') # (year, month, station)
 
    stn_monthly = np.zeros((nyr*12, nstn), dtype=np.float)
    stn_monthly.fill(-999.99)

    month_list = range(1, 12+1)*nyr
  
    # locate start year, month in the residual array
    k = ((yr[0] - Ryr[0]) + 1) - 1

    if k < 0:
      print 'Residual start year is greater than desired start year. Exiting.'
      sys.exit(0)

    iy = 0
    print 'Generating monthly projections... \n'
    for i,month in enumerate(month_list):
        if month == 12 or month == 1 or month == 2: 
            slope = patt[0,:]
        elif month >= 3 and month <= 5:
            slope = patt[1,:]
        elif month >= 6 and month <= 8:
            slope = patt[2,:]
        elif month >= 9 and month <= 11:
            slope = patt[3,:]

        stn_monthly[i,:] = slope[:]*globT[iy] + resid[k,month-1,:]

        if (i + 1) % 12 == 0:
          iy += 1
          # for resids
          if iy+yr[0] > 2099:
            k -= 1 # reverse 
          else:
            k += 1

    return stn_monthly


def date(yrStart, yrEnd):
   import datetime
   start_dt = datetime.datetime.strptime(str(yrStart)+'001', "%Y%j")
   end_dt = datetime.datetime.strptime(str(yrEnd)+'365', "%Y%j")
   while start_dt <= end_dt:
       yield start_dt.strftime("%Y%j")
       start_dt += datetime.timedelta(days=1)

def monthDates(yrStart, yrEnd):
    yyyymm=[]
    for year in xrange(yrStart, yrEnd + 1):
      for month in xrange(1, 12 + 1):
        yyyymm.append(str(year) + str(month).zfill(2))
    return yyyymm
