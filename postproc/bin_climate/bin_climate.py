'''

 @author: Rasmussen

 For each model, county and year, extract and bin daily temperature and precipitation

 Write to text file

'''
from netCDF4 import Dataset
import numpy as np
import glob
import os
import csv

rcp_list = ['26','45','60','85']

var_list = ['TAS','PR','TASMIN','TASMAX']
_dict = {"TAS":"average","TASMIN":"minimum","TASMAX":"maximum"}

# Where to write to
outdir_root = '/shares/gcp/climate/ACP/SMME/BIN'

#location of netCDFs
indir_root = '/shares/gcp/climate/ACP/SMME'

for RCP in rcp_list:

    for VAR in var_list:

        idir = indir_root + '/RCP' +RCP+ '/' +VAR+ '_DAILY'
        file_list = glob.glob( idir + '/*19810101*' )
        
        for _file in file_list:
        
            model = _file.split('_')[-3].upper()

            # get first and last year
            no_ext = _file.split(".")[0]
            yrs = no_ext.split("_")[-1]

            yr1 = int( yrs.split("-")[0][0:4] )
            yr2 = int( yrs.split("-")[1][0:4] )
            year_list = range(yr1,yr2+1)

            print "Working with MODEL: %s" %model.lower()

            # Extract to a temporary file that we will delete later
            base = os.path.dirname( _file )

            print "Unzipping %s ..." %_file
            os.system("gunzip -c " + _file + "> " + base + "/temp.nc")

            print "Opening %s ..." %(base + "/temp.nc")

            ncfile = Dataset(base + "/temp.nc", "r", format="NETCDF4")
            data = ncfile.variables[VAR.lower()][:,:] # time,county
            county_list = ncfile.variables["name"][:]
            state_list = ncfile.variables["state"][:]
            fips_list = ncfile.variables["fips"][:]

            if VAR[0:2] == 'PR':
               bins = np.linspace(0, 500, 251)
               header = "Description: binned total daily liquid precipitation (millimeters; rain, snow). Bin width is 2 mm. All but the last (righthand-most) bin is half-open. In other words, if bins is: [0, 2, 4, 6] then the first bin is [0, 2) (including 0, but excluding 2) and the second [2, 4). The last bin, however, is [4, 6], which includes 6."
            else:
               header = "Description: binned daily "+_dict[ VAR ]+ " temperature (degrees Celsius). Bin width is 1 C. All but the last (righthand-most) bin is half-open. In other words, if bins is: [1, 2, 3, 4] then the first bin is [1, 2) (including 1, but excluding 2) and the second [2, 3). The last bin, however, is [3, 4], which includes 4."
               data = data - 273.15
               bins = np.linspace(-40, 60, 101)

            # delete temporary file
            print "Removing %s ..." %(base + "/temp.nc")
            os.remove(base + "/temp.nc")
        
            xarr = np.zeros((len(year_list),len(bins)-1),np.int)
            for j,county in enumerate(county_list):
        
              county_lab = str("_".join([county_list[j].replace(" ","_"),state_list[j]]))
              odir = outdir_root + '/RCP' +RCP+ '/' +VAR + "/"+model+"/"

              if not os.path.exists(odir):
                 os.makedirs(odir)

              outf = odir + VAR.lower()+"_smme_daily_"+county_lab+"_"+str(fips_list[j])+"_"+model+"_rcp"+RCP+"_"+str(yr1)+"-"+str(yr2)+".tsv"

              print "Writing file: %s" %outf
              with open(outf,"wb") as out:

                  writer = csv.writer(out,delimiter="\t")
                  writer.writerow([header])
                  writer.writerow( ["Year"] + map(str, map(int,bins[0:(len(bins)-1)])) )

                  for i,year in enumerate(year_list):

                    i1 = i*365
                    i2 = (i+1)*365
  
                    writer.writerow([year] + map(int,np.histogram(data[i1:i2,j], bins)[0]))
