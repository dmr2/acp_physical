#!/bin/bash

# This script downloads the CONUS BCSD data set from the Bureau of Reclamation


# Select ensemble to download
ensemble=r1i1p1

# Local directory where downloaded BCSD data will be placed
root_dir=/home/dmr/jnk/

   for (( count=1; count<=39; count++))
   do
     if [ "$count" -eq 1 ]; then
       model=ACCESS1-0
       elif [ "$count" -eq 2 ]; then
       model=ACCESS1-3
       elif [ "$count" -eq 3 ]; then
       model=bcc-csm1-1-m
       elif [ "$count" -eq 4 ]; then
       model=bcc-csm1-1
       elif [ "$count" -eq 5 ]; then
       model=BNU-ESM
       elif [ "$count" -eq 6 ]; then
       model=CanESM2
       elif [ "$count" -eq 7 ]; then
       model=CCSM4
       elif [ "$count" -eq 8 ]; then
       model=CESM1-BGC
       elif [ "$count" -eq 9 ]; then
       model=CESM1-CAM5
       elif [ "$count" -eq 10 ]; then
       model=CMCC-CM
       elif [ "$count" -eq 11 ]; then
       model=CNRM-CM5
       elif [ "$count" -eq 12 ]; then
       model=CSIRO-Mk3-6-0
       elif [ "$count" -eq 13 ]; then
       model=EC-EARTH
       elif [ "$count" -eq 14 ]; then
       model=FGOALS-g2
       elif [ "$count" -eq 15 ]; then
       model=FGOALS-s2
       elif [ "$count" -eq 16 ]; then
       model=FIO-ESM
       elif [ "$count" -eq 17 ]; then
       model=GFDL-CM3
       elif [ "$count" -eq 18 ]; then
       model=GFDL-ESM2G
       elif [ "$count" -eq 19 ]; then
       model=GFDL-ESM2M
       elif [ "$count" -eq 20 ]; then
       model=GISS-E2-H-CC
       elif [ "$count" -eq 21 ]; then
       model=GISS-E2-R
       elif [ "$count" -eq 22 ]; then
       model=giss-E2-R-CC
       elif [ "$count" -eq 23 ]; then
       model=HadCM3
       elif [ "$count" -eq 24 ]; then
       model=HadGEM2-AO
       elif [ "$count" -eq 25 ]; then
       model=HadGEM2-CC
       elif [ "$count" -eq 26 ]; then
       model=HadGEM2-ES
       elif [ "$count" -eq 27 ]; then
       model=INMCM4
       elif [ "$count" -eq 28 ]; then
       model=IPSL-CM5A-LR
       elif [ "$count" -eq 29 ]; then
       model=IPSL-CM5A-MR
       elif [ "$count" -eq 30 ]; then
       model=IPSL-CM5B-LR
       elif [ "$count" -eq 31 ]; then
       model=MIROC-ESM
       elif [ "$count" -eq 32 ]; then
       model=MIROC-ESM-CHEM
       elif [ "$count" -eq 33 ]; then
       model=MIROC4h
       elif [ "$count" -eq 34 ]; then
       model=MIROC5
       elif [ "$count" -eq 35 ]; then
       model=MPI-ESM-LR
       elif [ "$count" -eq 36 ]; then
       model=MPI-ESM-MR
       elif [ "$count" -eq 37 ]; then
       model=MRI-CGCM3
       elif [ "$count" -eq 38 ]; then
       model=NorESM1-M
       elif [ "$count" -eq 39 ]; then
       model=NorESM1-ME
     else
        echo "ERROR"
        ehco "marker1"
     fi

     modlow=`echo $model | awk '{print tolower($0)}'`

     for RCP in rcp60 rcp85 rcp45 rcp26 historical
     do
       for var in pr tas tasmin tasmax
       do

# check if the file exists on the FTP...
        url=ftp://gdo-dcp.ucllnl.org/pub/dcp/archive/cmip5/bcsd/BCSD/${modlow}/${RCP}/
        wget -O/dev/null -q $url && exist=true || exist=false
      
        if $exist ; then

          dir=${root_dir}/${modlow}/${RCP}/mon/${var}/

          if [ ! -d ${dir} ]; then
           echo "Creating directory: ${dir}" 
           mkdir -p ${dir}
          fi


          filcount=`ls -1 ${dir}*.nc 2>/dev/null | wc -l` 
          if [ $filcount -eq 0 ]; then
            size=0
          else
            fil=`ls -1 ${dir}*.nc 2>/dev/null`
            size=$(wc -c < ${fil[0]})
          fi


          if [ $size -ge 200000000 ]; then
             echo "File already downloaded: ${fil[0]}"
          else
             url=ftp://gdo-dcp.ucllnl.org/pub/dcp/archive/cmip5/bcsd/BCSD/${modlow}/${RCP}/mon/${ensemble}/${var}/
             wget -nH --no-parent --directory-prefix=${dir} ${url}/*.nc
  
             find $dir -depth -exec rename 's/(.*)\/([^\/]*)/$1\/\L$2/' {} \;

          fi
        fi
       done # var
     done # rcp
   done # each model

# convert all characters to lower case

echo "COMPLETE!"  
