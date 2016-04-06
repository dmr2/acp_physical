#!/bin/bash

# netCDF to text

# location of NCL
ncl=/usr/local/ncl-6.3.0/bin/ncl

# location of concatenated historical+projected BCSD files
root=/home/dmr/jnk

path_in=root+"/cnty_txt"

normals=/home/dmr/acp_code/preproc/downscale/normals # location of monthly climate normals

OUTDIR=$root/cnty_txt

if [ ! -e $OOUTDIR ]; then
   mkdir -p $OUTDIR
fi

# Multi-year averages between these years
YR1=(2020 2040 2080)
YR2=(2039 2059 2099)

for RCP in "rcp26" "rcp45" "rcp60" "rcp85"; do
for VAR in "pr" "tas"; do

fils=`ls $root/merged/*/${RCP}/${VAR}/*195001-*.nc 2>/dev/null`

if [ ${#fils[@]} -eq 0 ]; then continue; fi

for FIL in ${fils[@]}; do

MOD=`echo $FIL | rev | cut -d'/' -f4 | rev`

echo "Working with [ RCP: $RCP ] [ Model: $MOD ] [ Variable: $VAR ] ..."

for i in `seq 0 2` ; do # each year pair

cat > cnty_season_mod2txt.ncl.tmp << EOF

; rasmussen; last updated: Sat 23 Jan 2016 11:45:26 AM PST

; Extracts multi-year averages of seasonal temperature change and 
; precipitation totals at county level and writes to a text file

; 1.  Calculates county-level seasonal normals from GHCND monthly station data
; 2.  Extracts county-level monthly data from netCDF files and calculates
;     seasonal averages as anomalies
; 3.  Map GHCND stations model and normals to CONUS counties 

; Assumes inputted gridded climate data are anomalies

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

   ; get list of normals and GHCND sites to process (no Alaska or Hawaii)
   norm_fil = "${normals}/ghcnd_${VAR}_month_normals_198101-201012.csv"
   print("Opening average monthly station normals: "+norm_fil)
   lines = asciiread(norm_fil,-1,"string")

   stn_code = str_get_field(lines,1,delim)
   stn_code!0 = "station"
   stn_code@_FillValue = "missing"
   stn_code@long_name = "GHCND site codes"
   nstn = dimsizes(stn_code)
 
   stn_state = str_get_field(lines,5,delim)
   stn_state!0 = "station"
   stn_state@_FillValue = "missing"
   stn_state@long_name = "GHCND site name"

   stn_name = str_get_field(lines,6,delim)
   stn_name!0 = "station"
   stn_name@_FillValue = "missing"
   stn_name@long_name = "GHCND site state"

   stn_lat = tofloat(str_get_field(lines,2,delim))
   stn_lat!0 = "station"
   stn_lat@_FillValue = -999.99
   stn_lat@units = "degrees_north"
   stn_lat@long_name = "lat"
   stn_lon = tofloat(str_get_field(lines,3,delim))
   stn_lon!0 = "station"
   stn_lon@_FillValue = -999.99
   stn_lon@units = "degrees_east"
   stn_lon@long_name = "lon"

   ; convert GHCND normal units
   stn_norm = new((/12,nstn/),"float",-999.999)
   do im=0, 11
    if "${VAR}" .eq. "pr" then
      ; Precip projections are avg. daily precip totals while GHCND normals are monthly totals, so
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
    
   ;Check time dimension
   utc_date = cd_calendar(xtime(dimsizes(xtime)-1), 0)
   dd = tointeger(utc_date(:,2))
   hh = tointeger(utc_date(:,3))
   EndCheck = tointeger(utc_date(:,0))
   if ( EndCheck .lt. ${YR2[$i]} ) then
     print( EndCheck+" ${YR2[$i]}" )
     print( "Model last year is less than desired end year...Skipping..." )
     exit
   end if
    
   ; find start and end year
   hh = tointeger(utc_date(:,3))
   option = 1
   option@calendar=xtime@calendar
   fstart = cd_inv_calendar(${YR1[$i]},1,dd,hh,0,0,xtime@units,option)
   istart = closest_val(fstart,xtime)
   fend = cd_inv_calendar(${YR2[$i]},12,dd,hh,0,0,xtime@units,option)
   iend = closest_val(fend,xtime)

   ; Model projections must be anomalized
   var_grid = fil->${VAR}(istart:iend,:,:) ; Read in monthly mean anomalies (time,lat,lon)
   var_grid&time = xtime(istart:iend)
   xtime2 = xtime(istart:iend)
   ntim = dimsizes(xtime2)

   ; Convert precip to monthly totals
   if "${VAR}" .eq. "pr" then
     print("Converting avg. daily precipitation to monthly totals")
     do i=0, ntim-1
      utc_date = cd_calendar(xtime2(i),0)
      yyyy = tointeger(utc_date(:,0))
      mm = tointeger(utc_date(:,1))
      var_grid(i,:,:) = var_grid(i,:,:)*days_in_month(yyyy,mm) ; assume non-leap year
     end do
   end if 

   lat_grid = fil->lat
   lon_grid = fil->lon

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

   ; Map gridcells to GHCND stations to 
   ; 2d lat/lon for GCM-GHCND Station Mapping

   lat2d = new((/dimsizes(lat_grid),dimsizes(lon_grid)/),"float",-999.999)
   lon2d = new((/dimsizes(lat_grid),dimsizes(lon_grid)/),"float",-999.999)

   do i=0,dimsizes(lon_grid)-1
    lat2d(:,i) = lat_grid
   end do
   do i=0,dimsizes(lat_grid)-1
    lon2d(i,:) = lon_grid
   end do
    
   ; Map GHCND stations to counties and GCM grid cells to county centroids

   ; Before mapping, check to see if BCSD grid cell mapping file already exists
   if fileexists("./${VAR}_conus_bcsd_xy_map.csv") then
    lines = asciiread("./${VAR}_conus_bcsd_xy_map.csv",-1,"string")
    x = toint(str_get_field(lines(1:),3,","))
    y = toint(str_get_field(lines(1:),4,","))
    do i=0, ncnty-1
      model_month_tmp(:,i) = (/var_grid(:,y(i),x(i))/)
    end do
    ncnty_no_AK_HI = i
   else
     print("Mapping GHCND stations and model grid cells to county lat/lon centroids...")
     print("")
     x = new((/ncnty/),"integer",-999) ; GCM mapping (x,y) to GHCND station
     y = new((/ncnty/),"integer",-999)
     ii=0
     do i=0, ncnty-1
      if cnty_lon(i)+360 .gt. 235 then; CONUS only
        ij = getind_closest_latlon2d(lat2d,lon2d,cnty_lat(i),cnty_lon(i)+360,var_grid(0,:,:),cnty_name(i),cnty_state(i)) 
        print((ii+1)+" Mapping GCM grid cell: ("+ij(0)+","+ij(1)+") to county centroid: "+ \
                                              cnty_name(i)+", "+cnty_state(i))
        x(i) = ij(1)
        y(i) = ij(0)
        model_month_tmp(:,ii) = (/var_grid(:,ij(0),ij(1))/)
        ii = ii + 1
      end if
     end do
     ; write table to text file
     txtFil = "${VAR}_conus_bcsd_xy_map.csv"
     print("Writing text file: "+txtFil)
     asciiwrite(txtFil,"x,y")
     write_table(txtFil,"a",[/cnty_name,cnty_state,x,y/], \
         "%s,%s,%3i,%3i")
     print("")
     ncnty_no_AK_HI = ii
   end if

   ; Map county centroid to nearest GHCND station
   r_polar = 6370000.
   r_equi = 6370000.
   ii = 0
   do i=0, ncnty-1
    if cnty_lon(i)+360 .gt. 235 then; CONUS only
     lat_dist = 2*(3.141532)*r_polar*(abs(stn_lat-cnty_lat(i))/360.)
     lon_dist = 2*(3.141532)*r_equi*(abs((stn_lon+360)-(cnty_lon(i)+360))/360.)
     dist = sqrt((lat_dist)^2 + (lon_dist)^2)
     indx = minind(dist)
     print("Mapping GHCND site: "+stn_name(indx)+", "+stn_state(indx)+" to county/ equivalent: "+cnty_name(i)+", "+cnty_state(i))
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
    
   txtFil = "${OUTDIR}/smme_county_seasonal_${VAR}_${RCP}_${MOD}_${YR1[$i]}-${YR2[$i]}.csv"

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
$ncl cnty_season_mod2txt.ncl.tmp | tee ncl.log
rm cnty_season_mod2txt.ncl.tmp

done # each year pair 
done # each model
done # each variable
done # each scenario
