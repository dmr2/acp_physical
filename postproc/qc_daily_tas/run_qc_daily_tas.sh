#!/bin/bash

# Quality Control # 

# Constrain daily Tmin < Tavg (Tas) < Tmax

# location of NCL
ncl=/usr/local/ncl-6.3.0/bin/ncl

# location of daily county-level climate projections of Tmin, Tavg, Tmax
root=/home/dmr/RUCORE

SET=001
METHOD=smme

OUTDIR=$root/

if [ ! -e $OOUTDIR ]; then
   mkdir -p $OUTDIR
fi

for RCP in "rcp26" "rcp45" "rcp60" "rcp85"; do

fils=`ls $root/${METHOD^^}/${RCP^^}/TAS_DAILY/${METHOD}_county_daily*.nc 2>/dev/null`

if [ ${#fils[@]} -eq 0 ]; then continue; fi

for FIL in ${fils[@]}; do

MOD=`echo $FIL | rev | cut -d'_' -f3 | rev`

echo "Working with [ RCP: $RCP ] [ Model: $MOD ] [Variables: tasmin, tas, tasmax]..."

fil1=`ls $root/${METHOD^^}/${RCP^^}/TAS_DAILY/${METHOD}_county_daily*_${MOD}_*.nc 2>/dev/null`
if [ ${#fil1[@]} -eq 0 ]; then continue; fi
fil2=`ls $root/${METHOD^^}/${RCP^^}/TASMIN_DAILY/${METHOD}_county_daily*_${MOD}_*.nc 2>/dev/null`
if [ ${#fil2[@]} -eq 0 ]; then continue; fi
fil3=`ls $root/${METHOD^^}/${RCP^^}/TASMAX_DAILY/${METHOD}_county_daily*_${MOD}_*.nc 2>/dev/null`
if [ ${#fil3[@]} -eq 0 ]; then continue; fi

TEST=`ls $root/${METHOD^^}/${RCP^^}/TASMAX_DAILY/ADJ/${METHOD}_county_daily_001_tasmax_${MOD}_${RCP}_19810101-*.nc 2>/dev/null`
echo  $TEST
if [ ! -z $TEST ]; then continue; fi

cat > qc_daily_tas.ncl.tmp << EOF

; rasmussen; last updated: Tue Feb 18 11:29:38 PST 2014

; assures that Tmin < Tavg < Tmax

; Where Tmin < Tavg < Tmax does not hold, calculates Tmin as Tavg-2.5 
; and Tmax as Tavg+2.5 where Tmin < Tavg < Tmax is not true

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
external QAQC "f90/qa_tmax_tmin.so"
 begin

 print("Opening file: ${fil1[0]}")
 ifil_tas = addfile("${fil1[0]}","r")

 print("Opening file: ${fil2[0]}")
 ifil_tmin = addfile("${fil2[0]}","r")

 print("Opening file: ${fil3[0]}")
 ifil_tmax = addfile("${fil3[0]}","r")

 xtime = ifil_tas->time

 yr1 = toint(floor(xtime(0)/1000))
 yr2 = toint(floor(xtime(dimsizes(xtime)-1)/1000))

 xtas = ifil_tas->tas
 xtmin = ifil_tmin->tasmin
 xtmax = ifil_tmax->tasmax

 cnty_fips = ifil_tmax->fips
 cnty_state = ifil_tmax->state
 cnty_name = ifil_tmax->name
 cnty_lat = ifil_tmax->lat
 cnty_lon = ifil_tmax->lon
 cnty_area = ifil_tmax->area
 cnty_pop = ifil_tmax->pop

 dims = dimsizes(xtas)
 ntime = dims(0)
 ncnty = dims(1)

 ; check if Tmin < Tas < Tmax
 print("")
 print("Constraining Tmin < Tas < Tmax")
 icnt = 0
 QAQC::qa_minmax(ncnty,ntime,icnt,xtmin(county|:,time|:),xtmax(county|:,time|:),xtas(county|:,time|:))

 print("Adjusted Tmin and Tmax for "+icnt+" days.")

 ; write to netCDF
 
 ; Tmax
 ncDir = "$root/${METHOD^^}/${RCP^^}/TASMAX_DAILY/ADJ"
 system("if ! test -d " + ncDir +" ; then mkdir -p " + ncDir + " ; fi")
 ncFil = "smme_county_daily_${SET}_tasmax_${MOD}_${RCP}_"+yr1+"0101-"+yr2+"1231_tmp.nc"
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
 if any(ismissing(xtmax)) then
  xarr_1d = ndtooned(xtmax)
  dsizes = dimsizes(xtmax)
  indices  = ind_resolve(ind(ismissing(xarr_1d)),dsizes)
  ncells = dimsizes(indices)
  do ii=0, ncells(0)-1
    print(" FOUND MISSING! x: " +indices(ii,0)+" y: "+indices(ii,1))
  end do
  printVarSummary(xtmax)
  exit
 end if

 ; coordinate variables
 dimNames = (/"time","county"/)
 dimSizes = (/-1,ncnty/)
 dimUnlim = (/True,False/)
 filedimdef( ncdf,dimNames,dimSizes,dimUnlim) ; make time UNLIMITED dimension

 filevardef(ncdf,"time",typeof(xtime),getvardims(xtime))
 filevardef(ncdf,"fips",typeof(cnty_fips),getvardims(cnty_fips))
 filevardef(ncdf,"state",typeof(cnty_state),getvardims(cnty_state))
 filevardef(ncdf,"name",typeof(cnty_name),getvardims(cnty_name))
 filevardef(ncdf,"lat",typeof(cnty_lat),getvardims(cnty_lat))
 filevardef(ncdf,"lon",typeof(cnty_lon),getvardims(cnty_lon))
 filevardef(ncdf,"area",typeof(cnty_area),getvardims(cnty_area))
 filevardef(ncdf,"pop",typeof(cnty_pop),getvardims(cnty_pop))
 filevardef(ncdf,"tasmax",typeof(xtmax),getvardims(xtmax))

 filevarattdef(ncdf,"time",xtime)
 filevarattdef(ncdf,"fips",cnty_fips)
 filevarattdef(ncdf,"state",cnty_state)
 filevarattdef(ncdf,"name",cnty_name)
 filevarattdef(ncdf,"lat",cnty_lat)
 filevarattdef(ncdf,"lon",cnty_lon)
 filevarattdef(ncdf,"area",cnty_area)
 filevarattdef(ncdf,"pop",cnty_pop)
 filevarattdef(ncdf,"tasmax",xtmax)
  
 varAtt = True
 varAtt@contents = "County level from GHCND stations"
 varAtt@comments = "Adjusted for Tmin < Tavg < Tmax"
 varAtt@model = "${MOD}"
 varAtt@experiment = "${RCP}"
 varAtt@history = "Processed by DJ Rasmussen for Rhodium Group, LLC; email: d.m.rasmussen.jr@gmail.com"
 varAtt@frequency ="daily"
 varAtt@actual_range = (/ min(xtmax), max(xtmax) /)
 varAtt@time = xtime(0)
 filevarattdef(ncdf,"tasmax",varAtt)

 ncdf->time = (/xtime/)
 ncdf->fips = (/cnty_fips/)
 ncdf->state = (/cnty_state/)
 ncdf->name = (/cnty_name/)
 ncdf->area = (/cnty_area/)
 ncdf->pop = (/cnty_pop/)
 ncdf->lat = (/cnty_lat/)
 ncdf->lon = (/cnty_lon/)
 ncdf->tasmax  =  (/xtmax/)

 setfileoption(ncdf,"DefineMode",False) ; explicitly exit file definition mode.

 ; replace original file
 ncFil = "smme_county_daily_${SET}_tasmax_${MOD}_${RCP}_"+yr1+"0101-"+yr2+"1231.nc"
 ncFile2 = ncDir + "/" +ncFil
 system("mv " + ncFile1 + " " +ncFile2)

 ; Tmin
 ncDir = "$root/${METHOD^^}/${RCP^^}/TASMIN_DAILY/ADJ"
 system("if ! test -d " + ncDir +" ; then mkdir -p " + ncDir + " ; fi")
 ncFil = "smme_county_daily_${SET}_tasmin_${MOD}_${RCP}_"+yr1+"0101-"+yr2+"1231_tmp.nc"
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
 if any(ismissing(xtmin)) then
  xarr_1d = ndtooned(xtmin)
  dsizes = dimsizes(xtmin)
  indices  = ind_resolve(ind(ismissing(xarr_1d)),dsizes)
  ncells = dimsizes(indices)
  do ii=0, ncells(0)-1
    print(" FOUND MISSING! x: " +indices(ii,0)+" y: "+indices(ii,1))
  end do
  printVarSummary(xtmin)
  exit
 end if

 ; coordinate variables
 dimNames = (/"time","county"/)
 dimSizes = (/-1,ncnty/)
 dimUnlim = (/True,False/)
 filedimdef( ncdf,dimNames,dimSizes,dimUnlim) ; make time UNLIMITED dimension

 filevardef(ncdf,"time",typeof(xtime),getvardims(xtime))
 filevardef(ncdf,"fips",typeof(cnty_fips),getvardims(cnty_fips))
 filevardef(ncdf,"state",typeof(cnty_state),getvardims(cnty_state))
 filevardef(ncdf,"name",typeof(cnty_name),getvardims(cnty_name))
 filevardef(ncdf,"lat",typeof(cnty_lat),getvardims(cnty_lat))
 filevardef(ncdf,"lon",typeof(cnty_lon),getvardims(cnty_lon))
 filevardef(ncdf,"area",typeof(cnty_area),getvardims(cnty_area))
 filevardef(ncdf,"pop",typeof(cnty_pop),getvardims(cnty_pop))
 filevardef(ncdf,"tasmin",typeof(xtmin),getvardims(xtmin))

 filevarattdef(ncdf,"time",xtime)
 filevarattdef(ncdf,"fips",cnty_fips)
 filevarattdef(ncdf,"state",cnty_state)
 filevarattdef(ncdf,"name",cnty_name)
 filevarattdef(ncdf,"lat",cnty_lat)
 filevarattdef(ncdf,"lon",cnty_lon)
 filevarattdef(ncdf,"area",cnty_area)
 filevarattdef(ncdf,"pop",cnty_pop)
 filevarattdef(ncdf,"tasmin",xtmin)
  
 varAtt = True
 varAtt@contents = "County level from GHCND stations"
 varAtt@comments = "Adjusted for Tmin < Tavg < Tmax"
 varAtt@model = "${MOD}"
 varAtt@experiment = "${RCP}"
 varAtt@history = "Processed by DJ Rasmussen for Rhodium Group, LLC; email: d.m.rasmussen.jr@gmail.com"
 varAtt@frequency ="daily"
 varAtt@actual_range = (/ min(xtmin), max(xtmin) /)
 varAtt@time = xtime(0)
 filevarattdef(ncdf,"tasmin",varAtt)

 ncdf->time = (/xtime/)
 ncdf->fips = (/cnty_fips/)
 ncdf->state = (/cnty_state/)
 ncdf->name = (/cnty_name/)
 ncdf->area = (/cnty_area/)
 ncdf->pop = (/cnty_pop/)
 ncdf->lat = (/cnty_lat/)
 ncdf->lon = (/cnty_lon/)
 ncdf->tasmin  =  (/xtmin/)

 setfileoption(ncdf,"DefineMode",False) ; explicitly exit file definition mode.

 ; replace original file
 ncFil = "smme_county_daily_${SET}_tasmin_${MOD}_${RCP}_"+yr1+"0101-"+yr2+"1231.nc"
 ncFile2 = ncDir + "/" +ncFil
 system("mv " + ncFile1 + " " +ncFile2)

 end
EOF

rm ncl.log
$ncl qc_daily_tas.ncl.tmp | tee ncl.log
rm qc_daily_tas.ncl.tmp

done # each model
done # each scenario
