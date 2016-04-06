#!/bin/bash

ncl=/usr/local/ncl-6.3.0/bin/ncl

rm ncl.log2

for var in 'pr' 'tas' 'tasmin' 'tasmax' ; do
for rcp in 'rcp26' 'rcp45' 'rcp60' 'rcp85'; do

root="/home/dmr/jnk"

# Output directory for seasonal means
seasonDir=$root"/merged/seasonal/"

fils=`ls $root/merged/*/${rcp}/${var}/${var}_mon_bcsd_${rcp}_*195001*.nc 2>/dev/null`
if [ -z "${fils[@]}" ]; then continue; fi

for fili in ${fils[@]}; do

model=`echo $fili | rev | cut -d'_' -f2 | rev`


echo -e "Generating patterns and residuals for [Model: $model] [Variable: $var] [Scenario: $rcp] \n"

cat > pattern_gen.tmp.ncl << EOF
; pattern_gen.ncl

; Written by DJ Rasmussen, last updated: Sun Sep  6 09:48:36 PDT 2015
; email: d-dot-m-dot-rasmussen-dot-jr-AT-gmail-dot-com
 
; In this script, a 30-year running average of seasonal temperatures
; are regressed against a 30-year running average of annual global mean 
; temperatures

; This scripts generates a pattern for pattern scaling using
; seasonal and annual data from 30-yr running averages

; Projected values are anomalized to the historical period of 1981-2010

; localT_anom(t2,i,j) = globalT_anom(t1)*pattern(t3,i,j) + residual(t2,i,j)

; Where t1 = year; t2 = month; t3 = season; i = latitude; j = longitude
; and localT_anom is with respect to each model's 1981-2010 period and the
; pattern and residual are from the same model (the models used are 
; defined below). The regression intercept is not included as it is 
; assumed that there is no local change under no global change.

; References:

; Mitchell, T., Pattern Scaling: An Examination of the Accuracy of the 
; Technique for Describing Future Climates. Climatic Change. 2003

 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

 begin

 month_abbr = (/"","Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep", \
                   "Oct","Nov","Dec"/)

 season = (/"DJF","MAM","JJA","SON"/)

 ; open monthly BCSD data
 print("")
 print("Opening: $fili ...")
 print("")
 fil = addfile("$fili", "r")

 ; Start by calculating the anomalies in each grid cell wrt 1981-2010
 xlocal_anom = fil->${var}(:,:,:)
 xtime = fil->time
 lat = fil->lat
 lon = fil->lon

 ntim = dimsizes(xtime)
 nlat = dimsizes(lat)
 nlon = dimsizes(lon)

 utc_date = cd_calendar( xtime(0), 0 )
 yyyy = tointeger( utc_date(:,0) )
 dd = tointeger( utc_date(:,2) )
 hh = tointeger( utc_date(:,3) )
 option = 1
 option@calendar = xtime@calendar
 yrStart = tointeger(yyyy)

 utc_date = cd_calendar( xtime(ntim-1), 0 )
 yyyy = tointeger( utc_date(:,0) )
 yrLast = tointeger(yyyy)

; fstart = cd_inv_calendar( 1981, 1, dd, hh, 0, 0, xtime@units, option )
; fend = cd_inv_calendar( 2010, 12, dd, hh, 0, 0, xtime@units, option )
; ist = closest_val( fstart, xtime )
; ied = closest_val( fend, xtime )
; xbcsd_norm = clmMonTLL( xbcsd(ist:ied,:,:) ) ; monthly averages
;
; xlocal_anom = new( (/(((yrLast-yrStart)+1)*12), nlat, nlon/), "float", -999.999 )
; print("")
; print("Calculating local monthly anomalies...")
; print("")
;
; xlocal_anom = calcMonAnomTLL( xbcsd, xbcsd_norm ) ; monthly avg anomalies
; xlocal_anom!0 = "time"
; xlocal_anom!1 = "lat"
; xlocal_anom!2 = "lon"
; xlocal_anom@units = xbcsd@units
; xlocal_anom@long_name = xbcsd@long_name
; xlocal_anom@_FillValue = -999.999
; xlocal_anom@actual_range = (/ min(xlocal_anom), max(xlocal_anom) /)
; xlocal_anom@comment = "local monthly anomaly wrt 1981-2010 monthly average"
       
 xtime_local = new( (/(((yrLast-yrStart)+1)*12)/), "double" )
 xtime_local!0 = "time"
 xtime_local@comment = "coordinate time for local anomalies"

 icount = 0
 do iy = yrStart,yrLast
  do im = 1, 12
   xtime_local(icount) = cd_inv_calendar( iy, im, dd, hh, 0, 0, xtime@units, option )
   icount = icount + 1
  end do
 end do
 xlocal_anom&time = xtime_local

 ; Calculate 30-yr run avg of local variable
 istart_month = 1
 istart_year = 1985 ; always start Jan 1985
 iend_month = 12
 iend_year = 2115 ; always end Dec 2115
 
 ; last year and month in the BCSD data
 utc_date = cd_calendar( xtime(ntim-1), 0 )
 yyyy_end = tointeger( utc_date(:,0) )
 mm_end = tointeger( utc_date(:,1) )
   
 ; create arrays to store 1985 to 2115 record
 nmonth = ((iend_year-istart_year) + 1)*12
 nyear = nmonth/12
 xwork = new( (/nmonth,nlat,nlon/), "float", 1.E+20 )
 xwork!0 = "time"
 xwork!1 = "lat"
 xwork!2 = "lon"
 
 xwork_time = new( (/nmonth/), "double", -999.999 )
 xwork_time!0 = "time"
 
 option = 1
 option@calendar = xtime@calendar
 fstart = cd_inv_calendar( 1985, 1, dd, hh, 0, 0, xtime@units, option )
 ist = closest_val( fstart, xtime )

 xwork(0:(((yrLast-1985)+1)*12)-1,:,:) = xlocal_anom(ist:ntim-1,:,:)

 it = 0
 do iy = 1985, 2115
  do im = 1, 12
    xwork_time(it) = cd_inv_calendar( iy, im, dd, hh, 0, 0, xtime@units, option )
    it = it + 1
  end do
 end do

 xseason = new( (/4,nyear,nlat,nlon/), "float", 1.E+20 )
 xseason!0 = "season"
 xseason!1 = "time"
 xseason!2 = "lat"
 xseason!3 = "lon"

 xDJF = new( (/nyear,nlat,nlon/), "float", 1.E+20)
 xDJF!0 = "time"
 xDJF!1 = "lat"
 xDJF!2 = "lon"
 xDJFtime = new( (/nyear/), "double", 1.E+20)
 xDJFtime!0 = "time"

 xMAM = new( (/nyear,nlat,nlon/), "float", 1.E+20)
 xMAM!0 = "time"
 xMAM!1 = "lat"
 xMAM!2 = "lon"
 xMAMtime = new( (/nyear/), "double", 1.E+20)
 xMAMtime!0 = "time"

 xJJA = new( (/nyear,nlat,nlon/), "float", 1.E+20)
 xJJA!0 = "time"
 xJJA!1 = "lat"
 xJJA!2 = "lon"
 xJJAtime = new( (/nyear/), "double", 1.E+20)
 xJJAtime!0 = "time"

 xSON = new( (/nyear,nlat,nlon/), "float", 1.E+20)
 xSON!0 = "time"
 xSON!1 = "lat"
 xSON!2 = "lon"
 xSONtime = new( (/nyear/), "double", 1.E+20)
 xSONtime!0 = "time"
 
 utc_date = cd_calendar( xwork_time(0), 0 )
 mm_str = tointeger( utc_date(:,1) )
 
 print("Calculating seasonal averages...")
 print("")
 if (mm_str .eq. 12) then 
  do i=0, nyear-1
    xDJF(i,:,:) = dim_avg_n( (/xwork(0+i*12,:,:),xwork(1+i*12,:,:),xwork(2+i*12,:,:)/), 0)
    xDJFtime(i) = xwork_time(0+i*12)
    xMAM(i,:,:) = dim_avg_n( (/xwork(3+i*12,:,:),xwork(4+i*12,:,:),xwork(5+i*12,:,:)/), 0)
    xMAMtime(i) = xwork_time(3+i*12)
    xJJA(i,:,:) = dim_avg_n( (/xwork(6+i*12,:,:),xwork(7+i*12,:,:),xwork(8+i*12,:,:)/), 0)
    xJJAtime(i) = xwork_time(6+i*12)
    xSON(i,:,:) = dim_avg_n( (/xwork(9+i*12,:,:),xwork(10+i*12,:,:),xwork(11+i*12,:,:)/), 0)
    xSONtime(i) = xwork_time(9+i*12)
  end do
 else if (mm_str .eq. 1) then
  do i=0, nyear-1
    xDJF(i,:,:) = dim_avg_n( (/xwork(11+i*12,:,:),xwork(0+i*12,:,:),xwork(1+i*12,:,:)/), 0)
    xDJFtime(i) = xwork_time(0+i*12)
    xMAM(i,:,:) = dim_avg_n( (/xwork(2+i*12,:,:),xwork(3+i*12,:,:),xwork(4+i*12,:,:)/), 0)
    xMAMtime(i) = xwork_time(2+i*12)
    xJJA(i,:,:) = dim_avg_n( (/xwork(5+i*12,:,:),xwork(6+i*12,:,:),xwork(7+i*12,:,:)/), 0)
    xJJAtime(i) = xwork_time(5+i*12)
    xSON(i,:,:) = dim_avg_n( (/xwork(8+i*12,:,:),xwork(9+i*12,:,:),xwork(10+i*12,:,:)/), 0)
    xSONtime(i) = xwork_time(8+i*12)
  end do
 else 
  print("Error. Model month must be 1 or 12. Exiting.")
  exit
 end if
 end if
  
 ; Historical seasonal averages
 bcsd_HISTseason = new( (/4,nlat,nlon/), "float", 1E+20)

 xtmp = new( (/(2010-1981)+1,nlat,nlon/), "float", 1E+20)
 do iy=1981,2010
    ftarg = cd_inv_calendar( iy, 1, dd, hh, 0, 0, xwork_time@units, option)
    itarg = closest_val( ftarg, xDJFtime )
    xtmp(iy-1981,:,:) = xDJF(itarg,:,:)
 end do
 bcsd_HISTseason(0,:,:) = dim_avg_n_Wrap( xtmp, 0)
 delete(xtmp)

 xtmp = new( (/(2010-1981)+1,nlat,nlon/), "float", 1E+20)
 do iy=1981,2010
    ftarg = cd_inv_calendar( iy, 3, dd, hh, 0, 0, xwork_time@units, option)
    itarg = closest_val( ftarg, xMAMtime )
    xtmp(iy-1981,:,:) = xMAM(itarg,:,:)
 end do
 bcsd_HISTseason(1,:,:) = dim_avg_n_Wrap( xtmp, 0)
 delete(xtmp)

 xtmp = new( (/(2010-1981)+1,nlat,nlon/), "float", 1E+20)
 do iy=1981,2010
    ftarg = cd_inv_calendar( iy, 6, dd, hh, 0, 0, xwork_time@units, option)
    itarg = closest_val( ftarg, xJJAtime )
    xtmp(iy-1981,:,:) = xJJA(itarg,:,:)
 end do
 bcsd_HISTseason(2,:,:) = dim_avg_n_Wrap( xtmp,0)
 delete(xtmp)

 xtmp = new( (/(2010-1981)+1,nlat,nlon/), "float", 1E+20)
 do iy=1981,2010
    ftarg = cd_inv_calendar( iy, 9, dd, hh, 0, 0, xwork_time@units, option)
    itarg = closest_val( ftarg, xSONtime )
    xtmp(iy-1981,:,:) = xSON(itarg,:,:)
 end do
 bcsd_HISTseason(3,:,:) = dim_avg_n_Wrap( xtmp, 0)
 delete(xtmp)

; linearly extrapolate to 2115 using 2069/70-2099/2100 rate

  xyears = ispan( yyyy_end-29, yyyy_end, 1); interpolate over these years

; DJF
  print("Linearly extrapolating DJF from "+yyyy_end+" to "+iend_year+"...")
  sm = cd_inv_calendar( yyyy_end-29, mm_str, dd, hh, 0, 0, xwork_time@units, option)
  em = cd_inv_calendar( yyyy_end, mm_str, dd, hh, 0, 0, xwork_time@units, option)
  is = closest_val( sm, xDJFtime)
  ie = closest_val( em, xDJFtime)
 
  rate = new( (/nlat,nlon/), "float", 1.E+20)
  yintercept = new( (/nlat,nlon/), "float", 1.E+20)
  yintercept!0 = "lat"
  yintercept!1 = "lon"
  yintercept&lat = lat
  yintercept&lon = lon
  
  rate(:,:) = regCoef( xyears, xDJF(lat|:,lon|:,time|is:ie) )
  rate!0 = "lat" 
  rate!1 = "lon"
  rate&lat = lat
  rate&lon = lon
  yintercept = onedtond( rate@yintercept, (/nlat,nlon/) )
  delete([/rate@nptxy,rate@rstd,rate@tval,rate@yintercept/]) ; delete unnecessary attributes

  do iy=(yyyy_end+1),iend_year
   iloc = (yyyy_end-1985) + (iy-(yyyy_end))
   xDJF(iloc,:,:) = (rate*iy) + yintercept
   xDJFtime(iloc) = cd_inv_calendar( iy, mm_str, dd, hh, 0, 0, xwork_time@units, option)
  end do
  delete([/rate,yintercept/])

; MAM
  print("Linearly extrapolating MAM from "+yyyy_end+" to "+iend_year+"...")
  sm = cd_inv_calendar( (yyyy_end-29), 3, dd, hh,0, 0, xwork_time@units, option)
  em = cd_inv_calendar( yyyy_end, 3, dd, hh, 0, 0, xwork_time@units, option)
  is = closest_val( sm, xMAMtime)
  ie = closest_val( em, xMAMtime)

  rate = new( (/nlat,nlon/), "float", 1.E+20)
  yintercept = new( (/nlat,nlon/), "float", 1.E+20)
  yintercept!0 = "lat"
  yintercept!1 = "lon"
  yintercept&lat = lat
  yintercept&lon = lon
  
  rate(:,:) = regCoef( xyears, xMAM(lat|:,lon|:,time|is:ie) )
  rate!0 = "lat"
  rate!1 = "lon"
  rate&lat = lat
  rate&lon = lon
  yintercept = onedtond( rate@yintercept, (/nlat,nlon/) )
  delete([/rate@nptxy,rate@rstd,rate@tval,rate@yintercept/]) ; delete unnecessary attributes

  do iy=(yyyy_end+1),iend_year
   iloc = (yyyy_end-1985) + (iy-(yyyy_end))
   xMAM(iloc,:,:) = (rate*iy) + yintercept
   xMAMtime(iloc) = cd_inv_calendar( iy, 3, dd, hh, 0, 0, xwork_time@units, option)
  end do
  delete([/rate,yintercept/])
 
; JJA
  print( "Linearly extrapolating JJA from "+yyyy_end+" to "+iend_year+"..." )
  sm = cd_inv_calendar( (yyyy_end-29), 6, dd, hh, 0, 0, xwork_time@units, option)
  em = cd_inv_calendar( (yyyy_end), 6, dd, hh, 0, 0, xwork_time@units, option)
  is = closest_val( sm, xJJAtime)
  ie = closest_val( em, xJJAtime)

  rate = new( (/nlat,nlon/), "float", 1.E+20)
  yintercept = new( (/nlat,nlon/), "float", 1.E+20)
  yintercept!0 = "lat"
  yintercept!1 = "lon"
  yintercept&lat = lat
  yintercept&lon = lon
  
  rate(:,:) = regCoef( xyears, xJJA(lat|:,lon|:,time|is:ie) )
  rate!0 = "lat" 
  rate!1 = "lon"
  rate&lat = lat
  rate&lon = lon
  yintercept = onedtond( rate@yintercept, (/nlat,nlon/) )
  delete([/rate@nptxy,rate@rstd,rate@tval,rate@yintercept/]) ; delete unnecessary attributes

  do iy=(yyyy_end+1),iend_year
   iloc = (yyyy_end-1985) + (iy-(yyyy_end))
   xJJA(iloc,:,:) = (rate*iy) + yintercept
   xJJAtime(iloc) = cd_inv_calendar( iy, 6, dd, hh, 0, 0, xwork_time@units, option)
  end do
  delete([/rate,yintercept/])

; SON
  print( "Linearly extrapolating SON from "+yyyy_end+" to "+iend_year+"..." )
  sm = cd_inv_calendar( (yyyy_end-29), 9, dd, hh, 0, 0, xwork_time@units, option)
  em = cd_inv_calendar( (yyyy_end),9, dd, hh, 0, 0, xwork_time@units, option)
  is = closest_val( sm, xSONtime)
  ie = closest_val( em, xSONtime)

  rate = new( (/nlat,nlon/), "float", 1.E+20)
  yintercept = new( (/nlat,nlon/), "float", 1.E+20)
  yintercept!0 = "lat"
  yintercept!1 = "lon"
  yintercept&lat = lat
  yintercept&lon = lon
  
  rate(:,:) = regCoef( xyears, xSON(lat|:,lon|:,time|is:ie) )
  rate!0 = "lat" 
  rate!1 = "lon"
  rate&lat = lat
  rate&lon = lon
  yintercept = onedtond( rate@yintercept, (/nlat,nlon/) )
  delete([/rate@nptxy,rate@rstd,rate@tval,rate@yintercept/]) ; delete unnecessary attributes

  do iy=(yyyy_end+1),iend_year
   iloc = (yyyy_end-1985) + (iy-(yyyy_end))
   xSON(iloc,:,:) = (rate*iy) + yintercept
   xSONtime(iloc) = cd_inv_calendar( iy, 9, dd, hh, 0, 0, xwork_time@units, option)
  end do
  delete([/rate,yintercept/])

  print("")
  print("Calculating 30-yr running averages...")
  x30DJF = runave_n_Wrap( xDJF, 30, 0, 0 )
  x30MAM = runave_n_Wrap( xMAM, 30, 0, 0 )
  x30JJA = runave_n_Wrap( xJJA, 30, 0, 0 )
  x30SON = runave_n_Wrap( xSON, 30, 0, 0 )

  ; each season: anomalize the 30-year running seasonal averages with the observed 1981-2010 period
  print("")
  print("Anomalizing BCSD data with 1981-2010 historical period...")
  print("")
  xworkanom = new( (/(2100-2000)+1,4,nlat,nlon/), "float", 1.E+20)
  xworkanom!0 = "time"
  xworkanom!1 = "season"
  xworkanom!2 = "lat"
  xworkanom!3 = "lon"

  ; DJF  
  if (mm_str .eq. 12) then ; if first month is December
    sm = cd_inv_calendar( 2000, 12, dd, hh, 0, 0, xwork_time@units, option )
    em = cd_inv_calendar( 2100, 12, dd, hh, 0, 0, xwork_time@units, option )
    is = closest_val( sm, xDJFtime )
    ie = closest_val( em, xDJFtime )
  else if (mm_str .eq. 1) then ; if first month is January
    sm = cd_inv_calendar( 2000, 1, dd, hh, 0, 0, xwork_time@units, option )
    em = cd_inv_calendar( 2100, 1, dd, hh, 0, 0, xwork_time@units, option )
    is = closest_val( sm, xDJFtime )
    ie = closest_val( em, xDJFtime )
  else 
    print("DJF month neither 12 or 1...Exiting...")
    exit
  end if
  end if

  icount = 0
  ; if ( "${var}" .eq. "pr" ) then ; get % change for precip
  ;   do ii=is,ie ; each year
  ;     xanomBCSD(icount,0,:,:)=( (bcsd_HISTseason(0,:,:)-x30DJF(ii,:,:))/x30DJF(ii,:,:) )*100. ;don't want the end years (no values)
  ;     icount=icount+1
  ;   end do
  ; else ; calculate absolute change for all other variables
     do ii=is,ie ; each year
       xworkanom(icount,0,:,:) = x30DJF(ii,:,:) - bcsd_HISTseason(0,:,:) ;don't want the end years (no values)
       icount = icount + 1
     end do
  ; end if

; MAM
   sm = cd_inv_calendar( 2000, 3, dd, hh, 0, 0, xwork_time@units, option )
   em = cd_inv_calendar( 2100, 3, dd, hh, 0, 0, xwork_time@units, option )
   is = closest_val( sm, xMAMtime )
   ie = closest_val( em, xMAMtime )

   icount = 0
  ; if ( "${var}" .eq. "pr" ) then ; get % change for precip
  ;   do ii=is,ie ; each year
  ;     xanomBCSD(icount,0,:,:)=( (bcsd_HISTseason(0,:,:)-x30MAM(ii,:,:))/x30MAM(ii,:,:) )*100. ;don't want the end years (no values)
  ;     icount=icount+1
  ;   end do
  ; else ; calculate absolute change for all other variables
     do ii=is,ie ; each year
       xworkanom(icount,1,:,:) = x30MAM(ii,:,:) - bcsd_HISTseason(1,:,:) ;don't want the end years (no values)
       icount = icount + 1
     end do
  ; end if

; JJA

   sm = cd_inv_calendar( 2000, 6, dd, hh, 0, 0, xwork_time@units, option )
   em = cd_inv_calendar( 2100, 6, dd, hh, 0, 0, xwork_time@units, option )
   is = closest_val( sm, xJJAtime)
   ie = closest_val( em, xJJAtime)

   icount = 0
 ;  if ( "${var}" .eq. "pr" ) then ; get % change for precip
 ;    do ii=is,ie ; each year
 ;      xanomBCSD(icount,0,:,:)=( (bcsd_HISTseason(0,:,:)-x30JJA(ii,:,:))/x30JJA(ii,:,:) )*100. ;don't want the end years (no values)
 ;      icount=icount+1
 ;    end do
 ;  else ; calculate absolute change for all other variables
     do ii=is,ie ; each year
       xworkanom(icount,2,:,:) = x30JJA(ii,:,:) - bcsd_HISTseason(2,:,:) ;don't want the end years (no values)
       icount = icount + 1
     end do
 ;  end if

;SON
   sm = cd_inv_calendar( 2000, 9, dd, hh, 0, 0, xwork_time@units, option )
   em = cd_inv_calendar( 2100, 9, dd, hh, 0, 0, xwork_time@units, option )
   is = closest_val( sm, xSONtime )
   ie = closest_val( em, xSONtime )

   icount = 0
 ;  if ( "${var}" .eq. "pr" ) then ; get % change for precip
 ;    do ii=is,ie ; each year
 ;      xanomBCSD(icount,0,:,:)=( (bcsd_HISTseason(0,:,:)-x30SON(ii,:,:))/x30SON(ii,:,:) )*100. ;don't want the end years (no values)
 ;      icount = icount + 1
 ;    end do
 ;  else ; calculate absolute change for all other variables
     do ii=is,ie ; each year
       xworkanom(icount,3,:,:) = x30SON(ii,:,:) - bcsd_HISTseason(3,:,:) ;don't want the end years (no values)
       icount = icount + 1
     end do
 ;  end if


     print( "Writing seasonal time series to disk..." )

     ; DJF
     ncDir = "$seasonDir"
     system("if ! test -d " + ncDir +" ; then mkdir -p " + ncDir + " ; fi")
     ncName = "${var}_DJF_bcsd_${model}_${rcp}_1981-2100.nc"
     ncFil = ncDir + ncName
     system("if test -s " + ncFil +" ; then rm " + ncFil + " ; fi")

     ncdf = addfile( ncFil, "c")
     setfileoption( ncdf, "DefineMode", True)
     globAtt = True
     globAtt@Conventions = "DJF average"
     globAtt@creation_date = systemfunc ( "date" )
     fileattdef( ncdf, globAtt )
     
     dimNames = (/"time","lat","lon"/)
     dimSizes = (/-1,nlat,nlon/)
     dimUnlim = (/True,False,False/)
     filedimdef( ncdf, dimNames, dimSizes, dimUnlim)
     filevardef( ncdf, "time", typeof(xDJFtime), getvardims(xDJFtime) )
     filevardef( ncdf, "lon", typeof(lon), getvardims(lon) )
     filevardef( ncdf, "lat", typeof(lat), getvardims(lat) )
     filevardef( ncdf, "${var}_anom", typeof(xlocal_anom), getvardims(xlocal_anom) )

     filevarattdef( ncdf, "time", xDJFtime )
     filevarattdef( ncdf, "lon", lon )
     filevarattdef( ncdf, "lat", lat )
     
     if (mm_str .eq. 12) then 
       sm = cd_inv_calendar( 1981, 12, dd, hh, 0, 0, xwork_time@units, option )
       em = cd_inv_calendar( 2100, 12, dd, hh, 0, 0, xwork_time@units, option )
       is = closest_val( sm, xDJFtime )
       ie = closest_val( em, xDJFtime )
     else if (mm_str .eq. 1) then 
       sm = cd_inv_calendar( 1981, 1, dd, hh, 0, 0, xwork_time@units, option )
       em = cd_inv_calendar( 2100, 1, dd, hh, 0, 0, xwork_time@units, option )
       is = closest_val( sm, xDJFtime )
       ie = closest_val( em, xDJFtime )
     else 
       print("DJF month neither 12 or 1...Exiting...")
       exit
     end if
     end if
     
     varAtt = True
     varAtt@info = "anomalized to 1981-2010"
     varAtt@history = "DJF seasonal average anomalies"
     varAtt@actual_range = (/ min(xDJF), max(xDJF) /)
     filevarattdef( ncdf, "${var}_anom", varAtt )

     print("Writing: " +ncFil)
     ncdf->${var}_anom = (/xDJF(is:ie,:,:)/)
     ncdf->time = (/xDJFtime(is:ie)/)
     ncdf->lat = (/lat/)
     ncdf->lon = (/lon/)
     setfileoption( ncdf, "DefineMode", False )

     ; MAM
     ncName = "${var}_MAM_bcsd_${model}_${rcp}_1981-2100.nc"
     ncFil = ncDir + ncName
     system("if test -s " + ncFil +" ; then rm " + ncFil + " ; fi")

     ncdf = addfile( ncFil, "c" )
     setfileoption( ncdf, "DefineMode", True )
     globAtt = True
     globAtt@Conventions = "MAM average"
     globAtt@creation_date = systemfunc ( "date" )
     fileattdef( ncdf, globAtt )
     
     dimNames = (/"time","lat","lon"/)
     dimSizes = (/-1,nlat,nlon/)
     dimUnlim = (/True,False,False/)
     filedimdef( ncdf, dimNames, dimSizes, dimUnlim )
     filevardef( ncdf, "time", typeof(xMAMtime), getvardims(xMAMtime) )
     filevardef( ncdf, "lon", typeof(lon), getvardims(lon) )
     filevardef( ncdf, "lat", typeof(lat), getvardims(lat) )
     filevardef( ncdf, "${var}_anom", typeof(xlocal_anom), getvardims(xlocal_anom) )

     filevarattdef( ncdf, "time", xMAMtime )
     filevarattdef( ncdf, "lon", lon )
     filevarattdef( ncdf, "lat", lat )
     
     sm = cd_inv_calendar( 1981, 3, 16, hh, 0, 0, xMAMtime@units, option )
     em = cd_inv_calendar( 2100, 3, 16, hh, 0, 0, xMAMtime@units, option )
     is = closest_val( sm, xMAMtime )
     ie = closest_val( em, xMAMtime )

     varAtt = True
     varAtt@info = "anomalized to 1981-2010"
     varAtt@history = "MAM seasonal average anomalies"
     varAtt@actual_range = (/ min(xMAM), max(xMAM) /)
     filevarattdef( ncdf,"${var}_anom", varAtt )

     print( "Writing: " +ncFil )
     ncdf->${var}_anom = (/xMAM(is:ie,:,:)/)
     ncdf->time = (/xMAMtime(is:ie)/)
     ncdf->lat = (/lat/)
     ncdf->lon = (/lon/)
     setfileoption( ncdf, "DefineMode", False )

     ; JJA
     ncName = "${var}_JJA_bcsd_${model}_${rcp}_1981-2100.nc"
     ncFil = ncDir + ncName
     system("if test -s " + ncFil +" ; then rm " + ncFil + " ; fi")

     ncdf = addfile( ncFil, "c" )
     setfileoption( ncdf, "DefineMode", True )
     globAtt = True
     globAtt@Conventions = "JJA average"
     globAtt@creation_date = systemfunc ( "date" )
     fileattdef( ncdf, globAtt )
     
     dimNames = (/"time","lat","lon"/)
     dimSizes = (/-1,nlat,nlon/)
     dimUnlim = (/True,False,False/)
     filedimdef( ncdf, dimNames, dimSizes, dimUnlim )
     filevardef( ncdf, "time", typeof(xJJAtime), getvardims(xMAMtime) )
     filevardef( ncdf, "lon", typeof(lon), getvardims(lon) )
     filevardef( ncdf, "lat", typeof(lat), getvardims(lat) )
     filevardef( ncdf, "${var}_anom", typeof(xlocal_anom), getvardims(xlocal_anom) )

     filevarattdef( ncdf, "time", xJJAtime )
     filevarattdef( ncdf, "lon", lon )
     filevarattdef( ncdf, "lat", lat )

     sm = cd_inv_calendar( 1981, 6, 16, hh, 0, 0, xJJAtime@units, option )
     em = cd_inv_calendar( 2100, 6, 16, hh, 0, 0, xJJAtime@units, option )
     is = closest_val( sm, xJJAtime )
     ie = closest_val( em, xJJAtime )
     
     varAtt = True
     varAtt@info = "anomalized to 1981-2010"
     varAtt@history = "JJA seasonal average anomalies"
     varAtt@actual_range = (/ min(xJJA), max(xJJA) /)
     filevarattdef( ncdf, "${var}_anom", varAtt )

     print( "Writing: " +ncFil )
     ncdf->${var}_anom = (/xJJA(is:ie,:,:)/)
     ncdf->time = (/xJJAtime(is:ie)/)
     ncdf->lat = (/lat/)
     ncdf->lon = (/lon/)
     setfileoption( ncdf, "DefineMode", False )

     ; SON
     ncName = "${var}_SON_bcsd_${model}_${rcp}_1981-2100.nc"
     ncFil = ncDir + ncName
     system( "if test -s " + ncFil +" ; then rm " + ncFil + " ; fi" )

     ncdf = addfile( ncFil, "c" )
     setfileoption( ncdf, "DefineMode", True )
     globAtt = True
     globAtt@Conventions = "SON average"
     globAtt@creation_date = systemfunc ( "date" )
     fileattdef( ncdf, globAtt )
     
     dimNames = (/"time","lat","lon"/)
     dimSizes = (/-1,nlat,nlon/)
     dimUnlim = (/True,False,False/)
     filedimdef( ncdf, dimNames, dimSizes, dimUnlim )
     filevardef( ncdf, "time", typeof(xSONtime), getvardims(xSONtime) )
     filevardef( ncdf, "lon", typeof(lon), getvardims(lon) )
     filevardef( ncdf, "lat", typeof(lat), getvardims(lat) )
     filevardef( ncdf, "${var}_anom", typeof(xlocal_anom), getvardims(xlocal_anom) )

     filevarattdef( ncdf, "time", xSONtime )
     filevarattdef( ncdf, "lon", lon )
     filevarattdef( ncdf, "lat", lat )

     sm = cd_inv_calendar( 1981, 9, 16, hh, 0, 0, xSONtime@units, option )
     em = cd_inv_calendar( 2100, 9, 16, hh, 0, 0, xSONtime@units, option )
     is = closest_val( sm, xSONtime )
     ie = closest_val( em, xSONtime )

     varAtt = True
     varAtt@info = "anomalized to 1981-2010"
     varAtt@history = "SON seasonal average anomalies"
     varAtt@actual_range = (/ min(xSON), max(xSON) /)
     filevarattdef( ncdf, "${var}_anom", varAtt )

     print( "Writing: " +ncFil )
     ncdf->${var}_anom = (/xSON(is:ie,:,:)/)
     ncdf->time = (/xSONtime(is:ie)/)
     ncdf->lat = (/lat/)
     ncdf->lon = (/lon/)
     setfileoption( ncdf, "DefineMode", False )

     ; generate 30-yr running avg. of global mean temperature record 
     print("")
     print( "Processing of BCSD data for "+str_upper("${model}")+" complete!" )
     print("")

     ; open global mean temperature
     txtDir = "global_tas/${rcp}/"+str_lower("${model}")+"/global_mean_tas/"
     txtFilName = "global_tas_aann_"+str_lower("${model}")+"_${rcp}_r1i1p1_1875-*.txt"
     txtFil = systemfunc( "ls "+txtDir+txtFilName )
     print( "Opening 30-yr running average of global mean temperatures...  "+txtFil )

     if ( ismissing(txtFil) ) then
       print( "Global temperatures not found: "+txtDir+txtFilName )
       print( "Skipping ${model}..." )
       continue
     end if
    
     lines = asciiread( txtFil, -1, "string" )
     delim = ","
     glob_anom_yr = toint( str_get_field(lines(1:), 1, delim) ) ; skip header
     glob_anom_tas = tofloat( str_get_field(lines(1:), 3, delim) ) ; anomalized to 1981-2010
     print( glob_anom_yr+" "+glob_anom_tas )
     delete( lines )
     
     ; generating patterns
     istrt = closest_val( 2000, glob_anom_yr )
     iend = closest_val( 2100, glob_anom_yr )
     xtime_glob = (/glob_anom_yr(istrt:iend)/) ; year (time)
     xtime_glob!0 = "year"
     xglob_runave30 = glob_anom_tas(istrt:iend)
     xglob_runave30!0 = "year"

     ; for residuals (maybe beyond 2000-2100)
     istrt = closest_val( yrStart, glob_anom_yr )
     iend = closest_val( yrLast, glob_anom_yr )
     xglob2_runave30 = glob_anom_tas(istrt:iend)
     xglob2_runave30!0 = "year"

     ; regress 30-yr running averages of seasonal BCSD data onto global mean temperature data
     slope = new( (/4,nlat,nlon/), "float", 1.E+20)
     intercept = new( (/4,nlat,nlon/), "float", 1.E+20)
     intercept!0 = "season" 
     intercept!1 = "lat" 
     intercept!2 = "lon" 
     intercept&lat = lat 
     intercept&lon = lon 
   
     ; slope standard errors
     rstd = new((/4,nlat,nlon/), "float", 1.E+20)
     rstd!0 = "season" 
     rstd!1 = "lat" 
     rstd!2 = "lon" 
     rstd&lat = lat 
     rstd&lon = lon 
    
     do iseason=0, 3
      slope(iseason,:,:) = regCoef( xglob_runave30, xworkanom(season|iseason,lat|:,lon|:,time|:) )
      rstd(iseason,:,:) = onedtond( slope@rstd, (/nlat,nlon/) )
      intercept(iseason,:,:) = onedtond( slope@yintercept, (/nlat,nlon/) )
     end do
    
     if ( "${var}" .eq."tas" ) then
      slope@units = "degree change per degree change in the global mean temperature"
     else if ( "${var}" .eq. "pr" ) then
      slope@units = "mm/d change per degree change in the global mean temperature"
     else
      slope@units = "per degree change in the global mean temperature"
     end if
     end if
   
     slope!0 = "season" 
     slope!1 = "lat" 
     slope!2 = "lon" 
     slope&lat = lat
     slope&lon = lon

     ; delete these attributes because they attach too much info and bloat filesize
     delete( [/slope@tval,slope@nptxy,slope@yintercept,slope@rstd/] )

     ; calculate the local residuals
     residual = new( (/(((yrLast-yrStart)+1)*12),nlat,nlon/), "float", -999.999 )
     residual!0 = "time"
     residual!1 = "lat"
     residual!2 = "lon"

     imonth = 0 ; count months
     do iyear=0,((yrLast-yrStart)+1)-1
       do im=1, 12
        if (im .eq. 12 .or. im .le. 2) then
         iseason = 0
        else if (im .ge. 3 .and. im .le. 5) then
         iseason = 1
        else if (im .ge. 6 .and. im .le. 8) then
         iseason = 2
        else if (im .ge. 9 .and. im .le. 11) then
         iseason = 3
        else
         print("Error. im is "+im)
         exit
        end if
        end if
        end if
        end if
       
       residual(imonth,:,:) = xlocal_anom(imonth,:,:) - xglob2_runave30(iyear)*slope(iseason,:,:)
       imonth = imonth + 1
       end do
     end do

     ; write everything we need to pattern scale to a netCDF file
     ncDir = "${root}/residuals/${var}/"
     system( "if ! test -d " + ncDir +" ; then mkdir -p " + ncDir + " ; fi" )
     ncName = "${var}_mon_residual_bcsd_${rcp}_${model}_"+yrStart+"-"+yrLast+".nc"
     ncFil = ncDir + ncName
     system( "if test -s " + ncFil +" ; then rm " + ncFil + " ; fi" )
     print("Going to write residuals to netCDF file: "+ncFil+" ...")

     ncdf = addfile( ncFil, "c" )
     setfileoption( ncdf, "DefineMode", True ) ; explicitly start file definition mode.
     globAtt = True
     globAtt@Conventions = "monthly residuals derived from DJF, MAM, JJA, SON patterns"
     globAtt@creation_date = systemfunc ( "date" )
     fileattdef( ncdf, globAtt ) ; update attributes
    
     dimNames = (/"time","lat","lon"/)
     dimSizes = (/-1,nlat,nlon/)
     dimUnlim = (/True,False,False/)
     filedimdef( ncdf, dimNames, dimSizes, dimUnlim ) ; make time UNLIMITED dimension
     delete( [/dimNames,dimSizes,dimUnlim/] )
    
     filevardef( ncdf, "time", typeof(xtime_local), getvardims(xtime_local) )
     filevardef( ncdf, "lon", typeof(lon), getvardims(lon) )
     filevardef( ncdf, "lat", typeof(lat), getvardims(lat) )
     filevardef( ncdf, "residual", typeof(residual), getvardims(residual) )
    
     filevarattdef( ncdf, "time", xtime_local )
     filevarattdef( ncdf, "lon", lon )
     filevarattdef( ncdf, "lat", lat )
     filevarattdef( ncdf, "residual", residual )
    
     varAtt = True
     varAtt@history = "monthly residuals from seasonal model patterns produced by DJ Rasmussen; email: d.m.rasmussen.jr@gmail.com"
     varAtt@info = "From model: ${model}"
     varAtt@actual_range = (/min(residual), max(residual)/)
     varAtt@time = xtime_local(0)
     filevarattdef( ncdf, "residual", varAtt )

     print( "Writing: " +ncFil )
     ncdf->time = (/xtime_local/)
     ncdf->lat = (/lat/)
     ncdf->lon = (/lon/)
     ncdf->residual = (/residual/)
     setfileoption( ncdf, "DefineMode", False )

     ncDir = "${root}/patterns/${var}/"
     system( "if ! test -d " + ncDir +" ; then mkdir -p " + ncDir + " ; fi" )
     ncName = "seasonal_pattern_${var}_bcsd_${rcp}_${model}.nc"
     ncFil = ncDir + ncName
     system( "if test -s " + ncFil +" ; then rm " + ncFil + " ; fi" )
     print( "Going to write patterns to netCDF file: "+ncFil+" ..." )

     xseason_time = new(4, "double", -999.999)
     xseason_time = (/1,2,3,4/)
     xseason_time!0 = "season"
    
     ncdf = addfile( ncFil, "c" )
     setfileoption( ncdf, "DefineMode", True )
     globAtt = True
     globAtt@Conventions = "seasons are: DJF, MAM, JJA, SON"
     globAtt@creation_date = systemfunc ( "date" )
     fileattdef( ncdf, globAtt )
    
     dimNames = (/"time","season","year","lat","lon"/)
     dimSizes = (/-1,4,dimsizes(xtime_glob),nlat,nlon/)
     dimUnlim = (/True,False,False,False,False/)
     filedimdef( ncdf, dimNames, dimSizes, dimUnlim )
     delete([/dimNames,dimSizes,dimUnlim/])
    
     filevardef( ncdf, "season", typeof(xseason_time), getvardims(xseason_time) )
     filevardef( ncdf, "year", typeof(xtime_glob), getvardims(xtime_glob) )
     filevardef( ncdf, "time", typeof(xtime_local), getvardims(xtime_local) )
     filevardef( ncdf, "lon", typeof(lon), getvardims(lon) )
     filevardef( ncdf, "lat", typeof(lat), getvardims(lat) )
     filevardef( ncdf, "slope", typeof(slope), getvardims(slope) )
     filevardef( ncdf, "standard_error", typeof(rstd), getvardims(rstd) )
     filevardef( ncdf, "intercept", typeof(intercept), getvardims(intercept) )
     filevardef( ncdf, "local_anom", typeof(xlocal_anom), getvardims(xlocal_anom) )
     filevardef( ncdf, "global_temp", typeof(xglob_runave30), getvardims(xglob_runave30) )
    
     filevarattdef( ncdf, "season", xseason_time )
     filevarattdef( ncdf, "year", xtime_glob )
     filevarattdef( ncdf, "time", xtime_local )
     filevarattdef( ncdf, "lon", lon )
     filevarattdef( ncdf, "lat", lat )
     filevarattdef( ncdf, "slope", slope )
     filevarattdef( ncdf, "local_anom", xlocal_anom )
     filevarattdef( ncdf, "global_temp", xglob_runave30 )
    
     varAtt = True
     varAtt@history = "seasonal patterns produced by DJ Rasmussen; email: d.m.rasmussen.jr@gmail.com"
     varAtt@actual_range = (/ min(slope), max(slope) /)
     varAtt@time = xseason_time(0)
     filevarattdef( ncdf, "slope", varAtt )
    
     varAtt = True
     varAtt@info = "30-yr running average of global mean temperature (degrees Celsius); anomalized to 1981-2010"
     varAtt@actual_range = (/ min(xglob_runave30), max(xglob_runave30) /)
     filevarattdef( ncdf, "global_temp", varAtt )

     print( "Writing: " +ncFil )
     ncdf->season = (/xseason_time/)
     ncdf->lat = (/lat/)
     ncdf->lon = (/lon/)
     ncdf->slope = (/slope/)
     ncdf->standard_error = (/rstd/)
     ncdf->intercept = (/intercept/)
     ncdf->year = (/xtime_glob/)
     ncdf->time = (/xtime_local/)
     ncdf->local_anom = (/xlocal_anom/)
     ncdf->global_temp = (/xglob_runave30/)

     setfileoption( ncdf, "DefineMode", False )

 end

EOF

$ncl pattern_gen.tmp.ncl >> ncl.log2
rm pattern_gen.tmp.ncl

done

done
done
