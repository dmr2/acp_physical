#!/bin/bash

# netCDF to text

# location of NCL
ncl=/usr/local/ncl-6.3.0/bin/ncl

# location of concatenated historical+projected BCSD files
root=/home/dmr/jnk

path_in=root+"/cnty_txt"

# location of monthly climate normals at GHCN stations
normals=/home/dmr/acp_code/preproc/downscale/normals

OUTDIR=$root/cnty_txt

if [ ! -e $OOUTDIR ]; then
   mkdir -p $OUTDIR
fi

# Multi-year averages between these years
YR1=(2020 2040 2080)
YR2=(2039 2059 2099)

for RCP in "rcp26" "rcp45" "rcp60" "rcp85"; do
for VAR in "pr" "tas"; do

fils=`ls $root/mcpr/*/${RCP}/*/${VAR}/*ghcn_*198101-*.nc 2>/dev/null`


if [ ${#fils[@]} -eq 0 ]; then continue; fi

for FIL in ${fils[@]}; do

MOD=`echo $FIL | rev | cut -d'/' -f3 | rev`


echo "Working with [ RCP: $RCP ] [ Model: $MOD ] [ Variable: $VAR ] ..."

for i in `seq 0 2` ; do # each year pair

cat > cnty_season_mod2txt_mcpr.ncl.tmp << EOF

; rasmussen; last updated: Thu 17 Mar 2016 03:37:49 PM PDT

; Extracts multi-year averages of seasonal temperature change and 
; precipitation totals, maps to the county level and writes to a text file

; 1.  Calculates county-level seasonal normals from GHCN monthly station data
; 2.  Extracts station-level monthly data from netCDF files and calculates
;     seasonal averages as anomalies
; 3.  Map GHCN stations model and normals to CONUS counties 

; Assumes input are anomalies at the station level

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "utils.ncl"

begin

   ; open county meta data
   meta_fil = "conus_county_info.csv"
   print("Opening county meta data: "+meta_fil)
   lines = asciiread(meta_fil,-1,"string")
   nlines = dimsizes(lines)-1   ; First line is a header
   delim = ","
   cnty_fips = toint(str_get_field(lines(1:),1,delim))
   cnty_fips!0 = "county"
   cnty_fips@long_name = "county FIPS"
   cnty_fips@_FillValue = -999
   cnty_state = str_get_field(lines(1:),2,delim)
   cnty_state!0 = "county"
   cnty_state@long_name = "county state"
   cnty_state@_FillValue = "missing"
   cnty_name  = str_get_field(lines(1:),3,delim)
   cnty_name!0 = "county"
   cnty_name@long_name = "county name"
   cnty_name@_FillValue = "missing"
   cnty_lat = tofloat(str_get_field(lines(1:),4,delim))
   cnty_lat!0 = "county"
   cnty_lat@long_name = "county centroid latitude"
   cnty_lat@units = "degrees_north"
   cnty_lat@_FillValue = -999.99
   cnty_lon = tofloat(str_get_field(lines(1:),5,delim))
   cnty_lon!0 = "county"
   cnty_lon@long_name = "county centroid longitude"
   cnty_lon@units = "degrees_east"
   cnty_lon@_FillValue = -999.99
   cnty_area = tofloat(str_get_field(lines(1:),6,delim))
   cnty_area!0 = "county"
   cnty_area@long_name = "county land area"
   cnty_area@units = "km**2"
   cnty_area@_FillValue = -999.99
   cnty_pop = tofloat(str_get_field(lines(1:),7,delim))
   cnty_pop!0 = "county"
   cnty_pop@long_name = "county population"
   cnty_pop@units = "number"
   cnty_pop@_FillValue = -999.99
   ncnty = dimsizes(lines)-1
   delete(lines)

   ; get list of normals and GHCN sites to process (no Alaska or Hawaii)
   norm_fil = "${normals}/ghcnd_${VAR}_month_normals_198101-201012.csv"
   print("Opening average monthly station normals: "+norm_fil)
   lines = asciiread(norm_fil,-1,"string")

   stn_code_all = str_get_field(lines,1,delim)
   stn_code_all!0 = "station"
   stn_code_all@_FillValue = "missing"
   stn_code_all@long_name = "GHCN site codes"
   nstn = dimsizes(stn_code_all)
 
   stn_state = str_get_field(lines,5,delim)
   stn_state!0 = "station"
   stn_state@_FillValue = "missing"
   stn_state@long_name = "GHCN site name"

   stn_name = str_get_field(lines,6,delim)
   stn_name!0 = "station"
   stn_name@_FillValue = "missing"
   stn_name@long_name = "GHCN site state"

   stn_lat_all = tofloat(str_get_field(lines,2,delim))
   stn_lat_all!0 = "station"
   stn_lat_all@_FillValue = -999.99
   stn_lat_all@units = "degrees_north"
   stn_lat_all@long_name = "lat"
   stn_lon_all = tofloat(str_get_field(lines,3,delim))
   stn_lon_all!0 = "station"
   stn_lon_all@_FillValue = -999.99
   stn_lon_all@units = "degrees_east"
   stn_lon_all@long_name = "lon"

   ; convert GHCN normal units
   stn_norm = new((/12,nstn/),"float",-999.999)
   do im=0, 11
    if "${VAR}" .eq. "pr" then
      ; Precip projections are avg. daily precip totals while GHCN normals are monthly totals, so
      ; we need to divide the normal monthly precip totals by the number of days in each month
      stn_norm(im,:) = tofloat(str_get_field(lines(:),7+im,delim))*25.4; inches to mm
    else
      stn_norm(im,:) = tofloat(str_get_field(lines(:),7+im,delim)) + 273.15
    end if
   end do
   delete(lines)

    
   ; calculate seasonal normals from monthly
   stn_norm_season = new((/4,nstn/),"float",-999.999)
   do is=0,3
     if "${VAR}" .eq. "pr" then ; total seasonal precip
        if is .eq. 0 then
          stn_norm_season(is,:) = stn_norm(11,:) + stn_norm(0,:) + stn_norm(1,:)
        else if is .eq. 1 then
          stn_norm_season(is,:) = stn_norm(2,:) + stn_norm(3,:) + stn_norm(4,:)
        else if is .eq. 2 then
          stn_norm_season(is,:) = stn_norm(5,:) + stn_norm(6,:) + stn_norm(7,:)
        else 
          stn_norm_season(is,:) = stn_norm(8,:) + stn_norm(9,:) + stn_norm(10,:)
        end if
        end if
        end if
     else 
        if is .eq. 0 then
          stn_norm_season(is,:) = ( stn_norm(11,:) + stn_norm(0,:) + stn_norm(1,:) )/3.
        else if is .eq. 1 then
          stn_norm_season(is,:) = ( stn_norm(2,:) + stn_norm(3,:) + stn_norm(4,:) )/3.
        else if is .eq. 2 then
          stn_norm_season(is,:) = ( stn_norm(5,:) + stn_norm(6,:) + stn_norm(7,:) )/3.
        else 
          stn_norm_season(is,:) = ( stn_norm(8,:) + stn_norm(9,:) + stn_norm(10,:) )/3.
        end if
        end if
        end if
     end if
   end do

   print("Opening anomalized monthly model projections: ${FIL}")
   fil = addfile("${FIL}","r")
   xtime = fil->time
    
   fstart = tointeger(tostring(${YR1[$i]}+"01"))
   istart = closest_val(fstart,xtime)

   fend = tointeger(tostring(${YR2[$i]}+"12"))
   iend = closest_val(fend,xtime)

   ; Model projections must be anomalized
   var_stn = fil->${VAR}(istart:iend,:) ; Read in monthly mean anomalies (time,station)
   stn_lat = fil->lat
   stn_lon = fil->lon
   stn_code = fil->code
   nstn = dimsizes(stn_lat)

   stn_state_tmp = new((/nstn/),"string")
   stn_name_tmp = new((/nstn/),"string")
   stn_code_tmp = new((/nstn/),"string")


   do i=1,dimsizes(stn_lat)-1
     ii = ind(stn_code_all.eq.stn_code(i))
     stn_state_tmp(i) = stn_state(ii)
     stn_name_tmp(i) = stn_name(ii)
     stn_code_tmp(i) = stn_code_all(ii)
   end do

   var_stn&time = xtime(istart:iend)
   xtime2 = xtime(istart:iend)
   ntim = dimsizes(xtime2)

   ; Convert precip to monthly totals
   if "${VAR}" .eq. "pr" then
     print("Converting avg. daily precipitation to monthly totals")
     do i=0, ntim-1
      yyyy = tointeger(floor(xtime2(i)/100))
      mm = floor(((xtime2(i)/100.) - tointeger(floor(xtime2(i)/100)))*100)
      var_stn(i,:) = var_stn(i,:)*days_in_month(yyyy,tointeger(mm)) ; assume non-leap year
     end do
   end if 


   cnty_state_temp = new((/ncnty/),"string")
   cnty_state_temp!0 = "county"
   cnty_name_temp = new((/ncnty/),"string")
   cnty_name_temp!0 = "county"
   cnty_fips_temp = new((/ncnty/),"string")
   cnty_fips_temp!0 = "county"
   cnty_area_temp = new((/ncnty/),"string")
   cnty_area_temp!0 = "county"
   cnty_lat_temp = new((/ncnty/),"float",-999.99)
   cnty_lat_temp!0 = "county"
   cnty_lon_temp = new((/ncnty/),"float",-999.99)
   cnty_lon_temp!0 = "county"

   model_month_tmp = new((/(iend-istart)+1,ncnty/),"float",-999.99)
   model_month_tmp!0 = "time"
   model_month_tmp!1 = "county"

   cnty_obs_season = new((/4,ncnty/),"float",-999.99)
   cnty_obs_season!0 = "time"
   cnty_obs_season!1 = "county"

   if "${VAR}" .eq. "pr" then
    model_month_tmp@units = "mm"
    model_month_tmp@long_name = "monthly avg. daily precipitation (liquid equivalent)"
   else
    model_month_tmp@units = "degrees C"
    model_month_tmp@long_name = "daily average temperature"
   end if
   model_month_tmp@_FillValue = -999.999
   model_month_tmp@comment = "monthly model projection anomaly at county centroids"

   ; Map county centroid to nearest GHCN station
   r_polar = 6370000.
   r_equi = 6370000.
   ii = 0
   do i=0, ncnty-1
    if cnty_lon(i)+360 .gt. 235 then; CONUS only
     lat_dist = 2*(3.141532)*r_polar*(abs(stn_lat-cnty_lat(i))/360.)
     lon_dist = 2*(3.141532)*r_equi*(abs((stn_lon+360)-(cnty_lon(i)+360))/360.)
     dist = sqrt((lat_dist)^2 + (lon_dist)^2)
     indx = minind(dist)
     print("Mapping GHCN site: "+stn_name_tmp(indx)+", "+stn_state_tmp(indx)+" to county/ equivalent: "+cnty_name(i)+", "+cnty_state(i))
     model_month_tmp(:,ii) = (/var_stn(:,indx)/)
     cnty_obs_season(:,ii) = (/stn_norm_season(:,indx)/)
     ii = ii + 1
    end if
   end do

   ; Calculate seasonal average anomalies and totals
   xvar_out = new((/4,ncnty/),"float",-999.99)
   season = (/"DJF","MAM","JJA","SON"/)

   if ( "${VAR}" .eq. "pr" ) then ; sum monthly precip totals over each season, leave as totals
     print("Calculating model seasonal precipitation total...")
     do iseason=0,3 
      if iseason .eq. 0 then ; Exclude last DJF in average due to no data for following year
       DJF = month2season_pr_total(model_month_tmp, season(iseason))
       xvar_out(iseason,:) = dim_avg_n_Wrap(DJF(0:(ntim/12)-2,:),0) + cnty_obs_season(iseason,:)
      else
       xvar_out(iseason,:) = dim_avg_n_Wrap(month2season_pr_total(model_month_tmp, season(iseason)),0) + cnty_obs_season(iseason,:)
      end if
     end do
   else
     print("Calculating model seasonal average anomalies...")
     do iseason=0,3
      xvar_out(iseason,:) = dim_avg_n_Wrap(month2season(model_month_tmp, season(iseason)),0) ;+ cnty_obs_season(iseason,:)
     end do
   end if
    
   txtFil = "${OUTDIR}/mcpr_county_seasonal_${VAR}_${RCP}_${MOD}_${YR1[$i]}-${YR2[$i]}.csv"

   print("Writing text file: "+txtFil)
   asciiwrite(txtFil,\
   "FIPS,county,state,lat,lon,area,population,DJF_obs,MAM_obs,JJA_obs,SON_obs,DJF_model,MAM_model,JJA_model,SON_model")
   write_table(txtFil,"a",[/cnty_fips,cnty_name,cnty_state,cnty_lat,cnty_lon,cnty_area,cnty_pop, \
                         cnty_obs_season(0,:),cnty_obs_season(1,:),cnty_obs_season(2,:),cnty_obs_season(3,:), \
                         xvar_out(0,:),xvar_out(1,:),xvar_out(2,:),xvar_out(3,:)/], \
                         "%4i,%s,%s,%7.3f,%7.3f,%5e,%5e,%7.3f,%7.3f,%7.3f,%7.3f,%7.3f,%7.3f,%7.3f,%7.3f")
    

 end
EOF

rm ncl.log
$ncl cnty_season_mod2txt_mcpr.ncl.tmp | tee ncl.log
rm cnty_season_mod2txt_mcpr.ncl.tmp

done # each year pair 
done # each model
done # each variable
done # each scenario
