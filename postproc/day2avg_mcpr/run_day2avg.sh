#!/bin/bash

# Create monthly averages from daily 

# location of NCL
ncl=/usr/local/ncl-6.3.0/bin/ncl

# location of daily county-level climate projections of Tmin, Tavg, Tmax
root=/home/dmr/RUCORE

SET=001
METHOD=mcpr

OUTDIR=$root/

if [ ! -e $OOUTDIR ]; then
   mkdir -p $OUTDIR
fi

for RCP in "rcp26" "rcp45" "rcp60" "rcp85"; do
for VAR in "pr" "tas" "tasmin" "tasmax"; do

fils=`ls $root/${METHOD^^}/${RCP^^}/${VAR^^}_DAILY/${METHOD}_county_daily*.nc 2>/dev/null`

if [ ${#fils[@]} -eq 0 ]; then continue; fi

for FIL in ${fils[@]}; do

MOD=`echo $FIL | rev | cut -d'_' -f3 | rev`

echo "Working with [ RCP: $RCP ] [ Model: $MOD ] [Variable: $VAR ]..."

cat > day2avg.ncl.tmp << EOF

; rasmussen; last updated: Sun 27 Mar 2016 07:01:39 PM PDT

; Averages daily to monthly

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
 begin

 print("Opening file: ${FIL}")
 ifil = addfile("${FIL}","r")

 xtime = ifil->time
 yyyymmdd = yyyyddd_to_yyyymmdd(xtime)

 yr1 = toint(floor(xtime(0)/1000))
 yr2 = toint(floor(xtime(dimsizes(xtime)-1)/1000))

 xvar = ifil->${VAR}

 cnty_fips = ifil->fips
 cnty_state = ifil->state
 cnty_name = ifil->name
 cnty_lat = ifil->lat
 cnty_lon = ifil->lon
 cnty_area = ifil->area
 cnty_pop = ifil->pop

 dims = dimsizes(xvar)
 ntim = dims(0)
 ncnty = dims(1)

 ; Convert time to seconds/hours/minutes since...
 xtime2 = new((/ntim/),"double")
 
 option = 1
 option@calendar = "noleap"

 do i=0, ntim-1
   yyyy = yyyymmdd(i)/10000
   mmdd = yyyymmdd(i) - (yyyy*10000)
   mm   = mmdd/100
   dd   = mmdd - (mm*100)
   xtime2(i) = cd_inv_calendar(yyyy,mm,dd,12,0,0,"days since 1950-01-01 00:00:00",option)
 end do

 delete(xvar&time)
 xvar&time = xtime2

 ; Average daily data to monthly
 print("")
 print("Calculating monthly means from daily...")
 if "${VAR}" .eq. "pr" then
   xmon = calculate_monthly_values(xvar,"sum",0,False) ; sum daily precipitation
 else
   xmon = calculate_monthly_values(xvar,"avg",0,False)
 end if

 xmon!0 = "time"
 xmon!1 = "county"

 yyyymm = yyyymm_time(yr1, yr2, "integer")
 yyyymm!0 = "time"

 ; Average monthly to annual
 if "${VAR}" .eq. "pr" then
   xann = month_to_annual(xmon,0) ; sum monthly precipitation
 else
   xann = month_to_annual(xmon,1) 
 end if

 xann!0 = "time"
 xann!1 = "county"

 yyyy_time = new((/yr2-yr1+1/),"integer")
 yyyy_time!0 = "time"

 do i=0, (yr2-yr1+1)-1
   yyyy_time(i) = yr1 + i
 end do

 xann&time = yyyy_time

 ; write to netCDF
 ncDir = "$root/${METHOD^^}/${RCP^^}/${VAR^^}_MONTH/"
 system("if ! test -d " + ncDir +" ; then mkdir -p " + ncDir + " ; fi")
 ncFil = "${METHOD}_county_mon_${VAR}_${MOD}_${RCP}_"+yr1+"01-"+yr2+"12.nc"
 ncFile1 = ncDir + "/" +ncFil

 print("Writing:" +ncFile1)
 print("")
 setfileoption("nc","Format","NetCDF4") ; explicitly start file definition mode.
 ncdf = addfile(ncFile1,"c")

 globAtt = True
 globAtt@creation_date = systemfunc ( "date" ) ; update creation time
 fileattdef(ncdf, globAtt) ; update attributes
 setfileoption(ncdf,"DefineMode",True)

 ; Check for values missing in arrays
 if any(ismissing(xmon)) then
  xarr_1d = ndtooned(xmon)
  dsizes = dimsizes(xmon)
  indices  = ind_resolve(ind(ismissing(xarr_1d)),dsizes)
  ncells = dimsizes(indices)
  do ii=0, ncells(0)-1
    print(" FOUND MISSING! x: " +indices(ii,0)+" y: "+indices(ii,1))
  end do
  printVarSummary(xmon)
  exit
 end if

 ; coordinate variables
 dimNames = (/"time","county"/)
 dimSizes = (/-1,ncnty/)
 dimUnlim = (/True,False/)
 filedimdef( ncdf,dimNames,dimSizes,dimUnlim) ; make time UNLIMITED dimension

 filevardef(ncdf,"time",typeof(yyyymm),getvardims(yyyymm))
 filevardef(ncdf,"fips",typeof(cnty_fips),getvardims(cnty_fips))
 filevardef(ncdf,"state",typeof(cnty_state),getvardims(cnty_state))
 filevardef(ncdf,"name",typeof(cnty_name),getvardims(cnty_name))
 filevardef(ncdf,"lat",typeof(cnty_lat),getvardims(cnty_lat))
 filevardef(ncdf,"lon",typeof(cnty_lon),getvardims(cnty_lon))
 filevardef(ncdf,"area",typeof(cnty_area),getvardims(cnty_area))
 filevardef(ncdf,"pop",typeof(cnty_pop),getvardims(cnty_pop))
 filevardef(ncdf,"${VAR}",typeof(xmon),getvardims(xmon))

 filevarattdef(ncdf,"time",yyyymm)
 filevarattdef(ncdf,"fips",cnty_fips)
 filevarattdef(ncdf,"state",cnty_state)
 filevarattdef(ncdf,"name",cnty_name)
 filevarattdef(ncdf,"lat",cnty_lat)
 filevarattdef(ncdf,"lon",cnty_lon)
 filevarattdef(ncdf,"area",cnty_area)
 filevarattdef(ncdf,"pop",cnty_pop)
 filevarattdef(ncdf,"${VAR}",xmon)
  
 varAtt = True
 varAtt@contents = "County level from GHCN stations"
 if "${VAR}" .eq. "pr" then
   varAtt@comments = "Monthly total from daily"
   varAtt@units = "mm"
 else
   varAtt@comments = "Monthly average from daily"
 end if
 varAtt@model = "${MOD}"
 varAtt@experiment = "${RCP}"
 varAtt@history = "Processed by DJ Rasmussen for Rhodium Group, LLC; email: d.m.rasmussen.jr@gmail.com"
 varAtt@frequency ="monthly"
 varAtt@actual_range = (/ min(xmon), max(xmon) /)
 varAtt@time = yyyymm(0)
 filevarattdef(ncdf,"${VAR}",varAtt)

 ncdf->time = (/yyyymm/)
 ncdf->fips = (/cnty_fips/)
 ncdf->state = (/cnty_state/)
 ncdf->name = (/cnty_name/)
 ncdf->area = (/cnty_area/)
 ncdf->pop = (/cnty_pop/)
 ncdf->lat = (/cnty_lat/)
 ncdf->lon = (/cnty_lon/)
 ncdf->${VAR}  =  (/xmon/)

 setfileoption(ncdf,"DefineMode",False) ; explicitly exit file definition mode.

 ; annual
 ncDir = "$root/${METHOD^^}/${RCP^^}/${VAR^^}_ANNUAL/"
 system("if ! test -d " + ncDir +" ; then mkdir -p " + ncDir + " ; fi")
 ncFil = "${METHOD}_county_ann_${VAR}_${MOD}_${RCP}_"+yr1+"-"+yr2+".nc"
 ncFile1 = ncDir + "/" +ncFil

 print("Writing:" +ncFile1)
 print("")
 setfileoption("nc","Format","NetCDF4") ; explicitly start file definition mode.
 ncdf = addfile(ncFile1,"c")

 globAtt = True
 globAtt@creation_date = systemfunc ( "date" ) ; update creation time
 fileattdef(ncdf, globAtt) ; update attributes
 setfileoption(ncdf,"DefineMode",True)

 ; Check for values missing in arrays
 if any(ismissing(xann)) then
  xarr_1d = ndtooned(xann)
  dsizes = dimsizes(xann)
  indices  = ind_resolve(ind(ismissing(xarr_1d)),dsizes)
  ncells = dimsizes(indices)
  do ii=0, ncells(0)-1
    print(" FOUND MISSING! x: " +indices(ii,0)+" y: "+indices(ii,1))
  end do
  printVarSummary(xann)
  exit
 end if

 ; coordinate variables
 dimNames = (/"time","county"/)
 dimSizes = (/-1,ncnty/)
 dimUnlim = (/True,False/)
 filedimdef( ncdf,dimNames,dimSizes,dimUnlim) ; make time UNLIMITED dimension

 filevardef(ncdf,"time",typeof(yyyy_time),getvardims(yyyy_time))
 filevardef(ncdf,"fips",typeof(cnty_fips),getvardims(cnty_fips))
 filevardef(ncdf,"state",typeof(cnty_state),getvardims(cnty_state))
 filevardef(ncdf,"name",typeof(cnty_name),getvardims(cnty_name))
 filevardef(ncdf,"lat",typeof(cnty_lat),getvardims(cnty_lat))
 filevardef(ncdf,"lon",typeof(cnty_lon),getvardims(cnty_lon))
 filevardef(ncdf,"area",typeof(cnty_area),getvardims(cnty_area))
 filevardef(ncdf,"pop",typeof(cnty_pop),getvardims(cnty_pop))
 filevardef(ncdf,"${VAR}",typeof(xann),getvardims(xann))

 filevarattdef(ncdf,"time",yyyy_time)
 filevarattdef(ncdf,"fips",cnty_fips)
 filevarattdef(ncdf,"state",cnty_state)
 filevarattdef(ncdf,"name",cnty_name)
 filevarattdef(ncdf,"lat",cnty_lat)
 filevarattdef(ncdf,"lon",cnty_lon)
 filevarattdef(ncdf,"area",cnty_area)
 filevarattdef(ncdf,"pop",cnty_pop)
 filevarattdef(ncdf,"${VAR}",xann)
  
 varAtt = True
 varAtt@contents = "County level from GHCN stations"
 if "${VAR}" .eq. "pr" then
   varAtt@comments = "Annual total from daily"
   varAtt@units = "mm"
 else
   varAtt@comments = "Annual average from daily"
 end if
 varAtt@model = "${MOD}"
 varAtt@experiment = "${RCP}"
 varAtt@history = "Processed by DJ Rasmussen for Rhodium Group, LLC; email: d.m.rasmussen.jr@gmail.com"
 varAtt@frequency ="annual"
 varAtt@actual_range = (/ min(xann), max(xann) /)
 varAtt@time = yyyy_time(0)
 filevarattdef(ncdf,"${VAR}",varAtt)

 ncdf->time = (/yyyy_time/)
 ncdf->fips = (/cnty_fips/)
 ncdf->state = (/cnty_state/)
 ncdf->name = (/cnty_name/)
 ncdf->area = (/cnty_area/)
 ncdf->pop = (/cnty_pop/)
 ncdf->lat = (/cnty_lat/)
 ncdf->lon = (/cnty_lon/)
 ncdf->${VAR}  =  (/xann/)

 setfileoption(ncdf,"DefineMode",False) ; explicitly exit file definition mode.
 end
EOF

rm ncl.log
$ncl day2avg.ncl.tmp | tee ncl.log
rm day2avg.ncl.tmp

done # each model
done # each variable
done # each scenario
