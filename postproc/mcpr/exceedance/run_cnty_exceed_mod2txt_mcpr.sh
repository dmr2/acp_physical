#!/bin/bash

# netCDF to text

# location of NCL
ncl=/usr/local/ncl-6.3.0/bin/ncl

# location of concatenated historical+projected BCSD files
root=/home/dmr/jnk

station=/home/dmr/acp_code/preproc/downscale/station_data # location of daily station data

SET=001 # weather realization set

# Set desired variable
VAR=tasmin # wmax | temperature (tas, tasmin, tasmax)

# Set desired threshold for variable
THRESH=0 # degrees C
ABOVE=False # count number of days above threshold (False for below threshold)

# Probabilistic model ensemble method
METHOD=mcpr # SMME | MCPR 

OUTDIR=$root/cnty_txt

if [ ! -e $OOUTDIR ]; then
   mkdir -p $OUTDIR
fi

# Multi-year averages between these years
YR1=(2020 2040 2080)
YR2=(2039 2059 2099)

for RCP in "rcp26" "rcp45" "rcp60" "rcp85"; do

fils=`ls $root/${METHOD}/${SET}/${RCP}/*/${VAR}/${METHOD}_county_daily*.nc 2>/dev/null`

if [ ${#fils[@]} -eq 0 ]; then continue; fi

for FIL in ${fils[@]}; do

MOD=`echo $FIL | rev | cut -d'/' -f3 | rev`

echo "Working with [ RCP: $RCP ] [ Model: $MOD ] [ Variable: $VAR ] ..."


for i in `seq 0 2` ; do # each year pair


cat > cnty_exceed_mod2txt_mcpr.ncl.tmp << EOF

; Uses daily model projections at the county level

; 1. Calculates the average number of days over a threshold 
;    in a multiyear period for both observational and model record

; 2. Writes to a text file


 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
 external COUNT_DAYS "f90/count.so" ; Count number of days
 begin

 above = ${ABOVE}

 ; Open daily observational-reanalysis combined historical station data 
 daily_filname = "${station}/ghcnd_station_${VAR}_daily_19810101-20101231.nc"
 print("Opening daily observational-reanalysis combined historical station data: "+daily_filname)
 dfil = addfile(daily_filname,"r")
 stn_daily_all = dfil->${VAR} - 273.15 ; Kelvin to Celsius
 stn_time = dfil->time ; daily

 stn_lat = dfil->lat
 stn_lon = dfil->lon
 stn_code = dfil->code

 nstn = dimsizes(stn_code)

 nyear = 2010 - 1981 + 1
 obs_count_stn = new((/nstn/),"float",-999.99) ; store threshold counts
 obs_count_stn = 0.

 COUNT_DAYS::count(nstn, dimsizes(stn_time), stn_daily_all, ${THRESH}, obs_count_stn, above)
 obs_count_stn = obs_count_stn * 1/tofloat(nyear)

 ; Open county-level model projections 
 print("Opening: ${FIL}")
 fil = addfile("${FIL}","r")
 model_time = fil->time

 ; find start and end year
 fstr = stringtoint(tostring(${YR1[$i]})+sprinti("%0.3i",1))
 fend = stringtoint(tostring(${YR2[$i]})+tostring(day_of_year(2001,12,31))); force leap year
 istrt =ind(model_time.eq.fstr)
 iend = ind(model_time.eq.fend)

 model_cnty_all = fil->${VAR}(istrt:iend,:) - 273.15 ; Kelvin to Celsius
 dims = dimsizes(model_cnty_all)
 ntime = dims(0)

 ; County-specific meta data
 cnty_fips = fil->fips
 cnty_name = fil->name
 cnty_state = fil->state
 cnty_lat = fil->lat
 cnty_lon = fil->lon
 cnty_area = fil->area
 cnty_pop = fil->pop
 ncnty = dimsizes(cnty_fips)

 ; Map station observations to counties
 print("Mapping GHCN stations to county lat/lon centroids...")
 print("")
 isite = new((/ncnty/),"integer",-999)
 ii = 0
 do i=0, ncnty-1
  if ( cnty_lon(i)+360 .gt. 235 ) then
    r_polar = 6370000.
    r_equi = 6370000.
    lat_dist = 2*(3.141532)*r_polar*(abs(stn_lat-cnty_lat(i))/360.)
    lon_dist = 2*(3.141532)*r_equi*(abs((stn_lon+360)-(cnty_lon(i)+360))/360.)
    dist = sqrt((lat_dist)^2 + (lon_dist)^2)
    indx = minind(dist)
    print("Mapping GHCN site: "+stn_code(indx)+ \
                     " to county/ equivalent: "+cnty_name(i)+", "+cnty_state(i))
    isite(ii) = indx
    ii = ii + 1
  end if
 end do

 obs_count = new((/ncnty/),"float",-999.99)
 obs_count!0 = "county"

 do icnty=0, ncnty-1
  obs_count(icnty) = (/obs_count_stn(isite(icnty))/)
 end do


 ; Calculate number of model exceedances
 nyear = ${YR2[$i]} - ${YR1[$i]} + 1
 model_count = new((/ncnty/),"float",-999.99) ; store threshold counts
 model_count = 0.0
 
 COUNT_DAYS::count(ncnty, ntime, model_cnty_all, ${THRESH}, model_count, above)
 model_count = model_count * 1/tofloat(nyear)

 if ( above ) then
  outFil = "${OUTDIR}/${METHOD}_county_${VAR}_num_day_gt${THRESH}_${MOD}_${RCP}_${YR1[$i]}-${YR2[$i]}.csv"
 else
  outFil = "${OUTDIR}/${METHOD}_county_${VAR}_num_day_lt${THRESH}_${MOD}_${RCP}_${YR1[$i]}-${YR2[$i]}.csv"
 end if
 
 print("Writing: "+outFil)
 asciiwrite(outFil,"FIPS,county,state,lat,lon,area,population,nday_obs,nday_model")
 write_table(outFil,"a",[/cnty_fips,cnty_name,cnty_state,cnty_lat,cnty_lon, \
            cnty_area,cnty_pop,obs_count,model_count/], \
            "%4i,%s,%s,%7.3f,%7.3f,%5e,%5e,%7.3f,%7.3f")
 end

EOF

rm ncl.log
$ncl cnty_exceed_mod2txt_mcpr.ncl.tmp | tee ncl.log
rm cnty_exceed_mod2txt_mcpr.ncl.tmp

done # each year pair 
done # each model
done # each scenario
