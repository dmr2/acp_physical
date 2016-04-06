#!/bin/bash


# location of NCL
ncl=/usr/local/ncl-6.3.0/bin/ncl

# location of concatenated historical+projected BCSD files
root=/home/dmr/jnk

normals=../downscale/normals # location of monthly climate normals
daily=../downscale/station_data # location of daily hybrid observation/reanalysis dataset

rm ncl.log

SET=001 # number of weather realizations (must be <= 20)

METHOD=smme

for RCP in "rcp26" "rcp45" "rcp60" "rcp85"; do

for MOD in "access1-3" "noresm1-me"; do

#Location of gridded 'tas' model output
FIL=`ls ${root}/merged/${MOD}/${RCP}/tas/tas_mon_bcsd_${RCP}_${MOD}_19*-2*.nc 2>/dev/null`

if [ -z ${FIL[@]} ]; then continue; fi

echo "Working with [ RCP: $RCP ] [ Model: $MOD ] ..."

cat > gen_missing.ncl.tmp << EOF

; last updated: Wed 27 Jan 2016 04:52:36 PM PST

; The Bureau of Reclaimation does not provide Tmin and/or Tmax output 
; for some BCSD models (see list below). This script calculates 
; monthly Tmax and Tmin from observed relationships between Tmin, 
; Tmax, and Tavg at GHCND stations.

; This is intended to produce monthly Tmin and Tmax projections 
; for models 1) NorESM1-ME 2) ACCESS1-3 where monthly Tmin and Tmax 
; projections are not provided, but models have monthly Tavg projections.
; Note: Expects to input Tavg as an anomaly; output are monthly anomalies.

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "utils.ncl"

 begin

  month_abbr = (/"","Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep", \
                    "Oct","Nov","Dec"/)

; produce lat/lon data (for BCSD ONLY!!)
  lat_grid = fspan(25.1875,52.8125,222)
  lon_grid = fspan(235.3125,292.9375,462)

  lat2d = new((/dimsizes(lat_grid),dimsizes(lon_grid)/),"float",-999.999)
  lon2d = new((/dimsizes(lat_grid),dimsizes(lon_grid)/),"float",-999.999)

  do i=0, dimsizes(lon_grid) - 1
   lat2d(:,i) = lat_grid
  end do
  do i=0, dimsizes(lat_grid) - 1
   lon2d(i,:) = lon_grid
  end do

  ; Open GHCND monthly normals for only site-related meta data
  ; get list of normals and GHCND sites to process (no Alaska or Hawaii)
  norm_fil = "${normals}/ghcnd_tasmin_month_normals_198101-201012.csv"
  print("Opening average monthly station normals: "+norm_fil)
  lines = asciiread(norm_fil,-1,"string")

  delim = ","
  stn_code_all = str_get_field(lines,1,delim)
  stn_code_all!0 = "station"
  stn_code_all@_FillValue = "missing"
  stn_code_all@long_name = "GHCND site codes"

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

  stn_state_all = str_get_field(lines,5,delim)
  stn_state_all!0 = "station"
  stn_state_all@_FillValue = "missing"
  stn_state_all@long_name = "GHCND site name"

  stn_name_all = str_get_field(lines,6,delim)
  stn_name_all!0 = "station"
  stn_name_all@_FillValue = "missing"
  stn_name_all@long_name = "GHCND site state"
  nstation = dimsizes(stn_code_all)

  ; convert GHCND normal units
  norm_tmin = new((/12,dimsizes(lines)/),"float",-999.999)
  do im=0, 11
     norm_tmin(im,:) = tofloat(str_get_field(lines(:),7+im,delim)) + 273.15
  end do
  delete(lines)

  norm_fil = "${normals}/ghcnd_tasmax_month_normals_198101-201012.csv"
  print("Opening average monthly station normals: "+norm_fil)
  lines = asciiread(norm_fil,-1,"string")

  ; convert GHCND normal units
  norm_tmax = new((/12,dimsizes(lines)/),"float",-999.999)
  do im=0, 11
     norm_tmax(im,:) = tofloat(str_get_field(lines(:),7+im,delim)) + 273.15
  end do
  delete(lines)

  norm_fil = "${normals}/ghcnd_tas_month_normals_198101-201012.csv"
  print("Opening average monthly station normals: "+norm_fil)
  lines = asciiread(norm_fil,-1,"string")

  ; convert GHCND normal units
  norm_tavg = new((/12,dimsizes(lines)/),"float",-999.999)
  do im=0, 11
     norm_tavg(im,:) = tofloat(str_get_field(lines(:),7+im,delim)) + 273.15
  end do
  delete(lines)

  ; Open daily observational-reanalysis combined historical station data 
  daily_filname = "${daily}/ghcnd_station_tasmax_daily_19810101-20101231.nc"
  print("Opening daily observational-reanalysis combined historical station data: "+daily_filname)
  dfil = addfile(daily_filname,"r")
  obsDayTmax = dfil->tasmax ; (time,station)
  stn_time = dfil->time ; daily

  daily_filname = "${daily}/ghcnd_station_tasmin_daily_19810101-20101231.nc"
  print("Opening daily observational-reanalysis combined historical station data: "+daily_filname)
  dfil = addfile(daily_filname,"r")
  obsDayTmin = dfil->tasmin ; (time,station)
  stn_time = dfil->time ; daily

  ; calculate monthly means from daily
  print("Calculating monthly means from daily...")
  obsMonTmin = calculate_monthly_values(obsDayTmin, "avg", 0, False); includes leap days
  obsMonTmax = calculate_monthly_values(obsDayTmax, "avg", 0, False)
  dims = dimsizes(obsMonTmax)
  nmon = dims(0)


; for each month, randomly select a historical year from the daily climatology 
  low = 1981 ; first obs year
  high = 2010 ; last obs year
  randyrs = floattointeger(random_uniform(low,high,(/151*12/))); assume 120 year record from CMIP5 archive

  ; Open this model's monthly average projections for Tavg

  print("Opening projected monthly avg. anomalies: ${FIL[0]}")
  conusfil = addfile("${FIL[0]}","r")
  grid_tas = conusfil->tas(:,:,:); (time, lat, lon)
  xtime = conusfil->time(:)

  utc_date = cd_calendar(xtime(0), 0)
  yr1 = tointeger(utc_date(:,0))
  utc_date = cd_calendar(xtime(dimsizes(xtime)-1), 0)
  yr2 = tointeger(utc_date(:,0))

  ; Map GHCND stations to grid cells
  stn_daily_temp = new((/dimsizes(stn_time),nstation/),"float",-999.99)
  stn_daily_temp!0 = "time"
  stn_daily_temp!1 = "station"
  stn_daily_temp&time = stn_time

  stn_code_temp = new((/nstation/),"string")
  stn_code_temp!0 = "station"
  stn_state_temp = new((/nstation/),"string")
  stn_state_temp!0 = "station"
  stn_name_temp = new((/nstation/),"string")
  stn_name_temp!0 = "station"

  lat_temp = new((/nstation/),"float",-999.99)
  lat_temp!0 = "station"
  lon_temp = new((/nstation/),"float",-999.99)
  lon_temp!0 = "station"

  x_tmp = new((/nstation/),"integer",-999)
  y_tmp = new((/nstation/),"integer",-999)

  obs_mon_tmin_tmp = new((/nmon,nstation/),"float",-999.99)
  obs_mon_tmax_tmp = new((/nmon,nstation/),"float",-999.99)

  ii = 0
  do istn=0, nstation-1
   if stn_lon(istn)+360 .gt. 235 then ; no AK or HI

      ; determine model lat and lon indices of GHCND station
      ij = getind_closest_latlon2d(lat2d,lon2d,stn_lat(istn),stn_lon(istn)+360, \
                           grid_tas(0,:,:),stn_name_all(istn),stn_state_all(istn)) 
      print((ii+1)+" Mapping GCM grid cell: ("+ij(0)+","+ij(1)+") to GHCND station: "\ 
                                   + stn_name_all(istn)+", "+stn_state_all(istn))

      x_tmp(ii) = ij(0)
      y_tmp(ii) = ij(1)

      ; select only non-Alaska and Hawaii stations
      lat_temp(ii) = stn_lat(istn)
      lon_temp(ii) = stn_lon(istn)
      stn_code_temp(ii) = stn_code_all(istn)
      stn_name_temp(ii) = stn_name_all(istn)
      stn_state_temp(ii) = stn_state_all(istn)

      obs_mon_tmin_tmp(:,ii) = (/obsMonTmin(:,istn)/)
      obs_mon_tmax_tmp(:,ii) = (/obsMonTmax(:,istn)/)

     ii = ii + 1
   end if ; if site is not Alaska or Hawaii
  end do ;  all stations
  print("")
  nstn_no_AK_HI = ii

  ; reshape arrays
  stn_daily = stn_daily_temp(:,0:nstn_no_AK_HI-1)
  stn_code = stn_code_temp(0:nstn_no_AK_HI-1)
  stn_name = stn_name_temp(0:nstn_no_AK_HI-1)
  stn_state = stn_state_temp(0:nstn_no_AK_HI-1)
  obs_mon_tmin = obs_mon_tmin_tmp(:,0:nstn_no_AK_HI-1)
  obs_mon_tmax = obs_mon_tmax_tmp(:,0:nstn_no_AK_HI-1)
  lat = lat_temp(0:nstn_no_AK_HI-1)
  lon = lon_temp(0:nstn_no_AK_HI-1)
  x = x_tmp(0:nstn_no_AK_HI-1)
  y = y_tmp(0:nstn_no_AK_HI-1)

  dTmin = new((/nstn_no_AK_HI/),"float",-999.999)
  dTmax = new((/nstn_no_AK_HI/),"float",-999.999)

  xTmin = new((/dimsizes(xtime),nstn_no_AK_HI/),"float",-999.999)
  xTmin!0 = "time"
  xTmin!1 = "station"
  xTmin@long_name = "monthly average minimum daily temperature"
  xTmin@units = "K"

  xTmax = new((/dimsizes(xtime),nstn_no_AK_HI/),"float",-999.999)
  xTmax!0 = "time"
  xTmax!1 = "station"
  xTmax@long_name = "monthly average maximum daily temperature"
  xTmax@units = "K"

; loop over all months, years

 
 yyyymm = yyyymm_time(1950,2100,"integer")
 months = yyyymm - (yyyymm/100)*100

 do imon=0,dimsizes(xtime)-1
    rndyr = randyrs(imon)
    utc_date = cd_calendar(xtime(imon), 0)
    yr = tointeger(utc_date(:,0))
    mon = tointeger(utc_date(:,1))
    print("Working on year, month: "+yr+", "+month_abbr(mon))
  
    ; Select random month's Tmin, Tavg, Tmax relationship
    idx = ((12*(rndyr-1981))+mon-1) - 1
    dTmin = (obs_mon_tmin(((12*(rndyr-1981))+mon)-1,:) + obs_mon_tmax(((12*(rndyr-1981))+mon)-1,:))/2. \
                       - obs_mon_tmin(((12*(rndyr-1981))+mon)-1,:)
    dTmax = obs_mon_tmax(((12*(rndyr-1981))+mon)-1,:) - (obs_mon_tmin(((12*(rndyr-1981))+mon)-1,:) \ 
                       + obs_mon_tmax(((12*(rndyr-1981))+mon)-1,:))/2.

    do isite=0, nstn_no_AK_HI-1
       xTmin(imon,isite) = (/grid_tas(imon,x(isite),y(isite))/) - dTmin(isite) \
                            + norm_tavg(months(imon)-1,isite) - norm_tmin(months(imon)-1,isite)
       xTmax(imon,isite) = (/grid_tas(imon,x(isite),y(isite))/) + dTmax(isite) \
                            + norm_tavg(months(imon)-1,isite) - norm_tmax(months(imon)-1,isite)
    end do; site
 end do ; each month

; write to netCDF
; QA before we write to netCDF
   if any(ismissing(xTmin)) .or. any(ismissing(xTmax))then
     print("Error! Found a missing value in the out array ")
     print("Exiting...")
     exit
   end if

 ; write out netCDF file for Tmax
   ncDir = "${root}/smme/${SET}/${RCP}/${MOD}/tasmax"
   system("if ! test -d " + ncDir +" ; then mkdir -p " + ncDir + " ; fi")
   ncFil = "smme_ghcnd_mon_${SET}_tasmax_${MOD}_${RCP}_"+yr1+"01-"+yr2+"12.nc"
   NCFILE = ncDir + "/" +ncFil
   
   print("")
   print("Writing: " +NCFILE)
   setfileoption("nc","Format","NetCDF4") ; explicitly start file definition mode.
   ncdf = addfile(NCFILE,"c")
   
   globAtt = True
   globAtt@Conventions = "None"
   globAtt@frequency = "monthly"
   globAtt@creation_date = systemfunc ( "date" )
   fileattdef( ncdf, globAtt) ; update attributes
   
   ;   predefine the coordinate variables
   dimNames = (/"time","station"/)
   dimSizes = (/-1,nstn_no_AK_HI/)
   dimUnlim = (/True,False/)
   filedimdef( ncdf,dimNames,dimSizes,dimUnlim) ; make time UNLIMITED dimension
   
   filevardef(ncdf,"time",typeof(xtime),getvardims(xtime))
   filevardef(ncdf,"code",typeof(stn_code),getvardims(stn_code))
   filevardef(ncdf,"name",typeof(stn_name),getvardims(stn_name))
   filevardef(ncdf,"state",typeof(stn_state),getvardims(stn_state))
   filevardef(ncdf,"lat",typeof(lat),getvardims(lat))
   filevardef(ncdf,"lon",typeof(lon),getvardims(lon))
   filevardef(ncdf,"tasmax",typeof(xTmax),getvardims(xTmax))

   filevarattdef(ncdf,"time",xtime)
   filevarattdef(ncdf,"code",stn_code_all)
   filevarattdef(ncdf,"name",stn_name_all)
   filevarattdef(ncdf,"state",stn_state_all)
   filevarattdef(ncdf,"lat",stn_lat)
   filevarattdef(ncdf,"lon",stn_lon)
   filevarattdef(ncdf,"tasmax",xTmax)

   varAtt = True
   varAtt@contents = "projections at GHCND stations"
   varAtt@model = "${MOD}"
   varAtt@experiment = "${RCP}"
   varAtt@set = "${SET}"
   varAtt@history = "Processed by DJ Rasmussen (Rhodium Group, LLC); email: d.m.rasmussen.jr@gmail.com"
   varAtt@frequency ="monthly"
   varAtt@reference = "Method from Wood et al. (2002) (section 2.3.2) from JGR-atmospheres"
   varAtt@actual_range = (/ min(xTmax), max(xTmax) /)
   varAtt@time = xtime(0)
   filevarattdef(ncdf,"tasmax",varAtt)
   
   ncdf->time = (/xtime/)
   ncdf->code = (/stn_code/)
   ncdf->name = (/stn_name/)
   ncdf->state = (/stn_state/)
   ncdf->lat = (/lat/)
   ncdf->lon = (/lon/)
   ncdf->tasmax  =  (/xTmax/)
    
   setfileoption(ncdf,"DefineMode",False) ; explicitly exit file definition mode.


 ; write out netCDF file for Tmin
   ncDir = "${root}/smme/${SET}/${RCP}/${MOD}/tasmin"
   system("if ! test -d " + ncDir +" ; then mkdir -p " + ncDir + " ; fi")
   ncFil = "smme_ghcnd_mon_${SET}_tasmin_${MOD}_${RCP}_"+yr1+"01-"+yr2+"12.nc"
   NCFILE = ncDir + "/" +ncFil
   
   print("")
   print("Writing: " +NCFILE)
   setfileoption("nc","Format","NetCDF4") ; explicitly start file definition mode.
   ncdf = addfile(NCFILE,"c")
   
   globAtt = True
   globAtt@Conventions = "None"
   globAtt@frequency = "monthly"
   globAtt@creation_date = systemfunc ( "date" )
   fileattdef( ncdf, globAtt) ; update attributes
   
   ;   predefine the coordinate variables
   dimNames = (/"time","station"/)
   dimSizes = (/-1,nstn_no_AK_HI/)
   dimUnlim = (/True,False/)
   filedimdef( ncdf,dimNames,dimSizes,dimUnlim) ; make time UNLIMITED dimension
   
   filevardef(ncdf,"time",typeof(xtime),getvardims(xtime))
   filevardef(ncdf,"code",typeof(stn_code),getvardims(stn_code))
   filevardef(ncdf,"name",typeof(stn_name),getvardims(stn_name))
   filevardef(ncdf,"state",typeof(stn_state),getvardims(stn_state))
   filevardef(ncdf,"lat",typeof(lat),getvardims(lat))
   filevardef(ncdf,"lon",typeof(lon),getvardims(lon))
   filevardef(ncdf,"tasmin",typeof(xTmin),getvardims(xTmin))

   filevarattdef(ncdf,"time",xtime)
   filevarattdef(ncdf,"code",stn_code_all)
   filevarattdef(ncdf,"name",stn_name_all)
   filevarattdef(ncdf,"state",stn_state_all)
   filevarattdef(ncdf,"lat",stn_lat)
   filevarattdef(ncdf,"lon",stn_lon)
   filevarattdef(ncdf,"tasmin",xTmin)

   varAtt = True
   varAtt@contents = "projections at GHCND stations"
   varAtt@model = "${MOD}"
   varAtt@experiment = "${RCP}"
   varAtt@set = "${SET}"
   varAtt@history = "Processed by DJ Rasmussen (Rhodium Group, LLC); email: d.m.rasmussen.jr@gmail.com"
   varAtt@frequency ="monthly"
   varAtt@reference = "Method from Wood et al. (2002) (section 2.3.2) from JGR-atmospheres"
   varAtt@actual_range = (/ min(xTmin), max(xTmin) /)
   varAtt@time = xtime(0)
   filevarattdef(ncdf,"tasmin",varAtt)
   
   ncdf->time = (/xtime/)
   ncdf->code = (/stn_code/)
   ncdf->name = (/stn_name/)
   ncdf->state = (/stn_state/)
   ncdf->lat = (/lat/)
   ncdf->lon = (/lon/)
   ncdf->tasmin  =  (/xTmin/)
    
   setfileoption(ncdf,"DefineMode",False) ; explicitly exit file definition mode.

end
EOF

$ncl gen_missing.ncl.tmp | tee ncl.log
rm gen_missing.ncl.tmp

done # each model
done # each scenario
