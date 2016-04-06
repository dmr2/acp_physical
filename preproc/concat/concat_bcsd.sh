#!/bin/bash

# Concatenate historical and future BCSD files

# location of NCL
ncl=/usr/local/ncl-6.3.0/bin/ncl

# location of downloaded BCSD files
root=/home/dmr/jnk

rm ncl.log

for RCP in "rcp26" "rcp45" "rcp60" "rcp85"; do
for VAR in "pr" "tas" "tasmin" "tasmax"; do

# Download directory
fils=`ls ${root}/*/historical/mon/${VAR}/*.nc 2>/dev/null`

if [ ${#fils[@]} -eq 0 ]; then continue; fi

for FIL in ${fils[@]}; do

MOD=`echo $FIL | rev | cut -d'/' -f5 | rev`

ncDir=$root/merged/${MOD}/${RCP}/${VAR}

if [ ! -e $ncDir ]; then
   mkdir -p $ncDir
fi

echo "Working with $RCP $MOD $VAR ..."

cat > concat_bcsd.ncl.tmp << EOF
; concat_bcsd.ncl

; Written by DJ Rasmussen, last modified Mon Dec 23 15:14:01 PST 2013
; email: d-dot-m-dot-rasmussen-dot-jr-AT-gmail-dot-com

; 1. Concatenates fragmented BCSD output (i.e. historical and projected records)
; 2. Anomalize all model output to 1981-2010 period from their own historical record

; This is necessary because we want the time dimension to have a consistent calendar
; between historical and projected periods and also one file is easier to work with

 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
 load "utils.ncl"

 begin

 ; open historical 
 path = "${root}/${MOD}/historical/mon/${VAR}"
 hfils = systemfunc("ls "+path+"/bcsd_0.125deg_${VAR}_amon_${MOD}_historical_r1i1p1_*.nc 2>/dev/null")

 ; open projected
 path = "${root}/${MOD}/${RCP}/mon/${VAR}"
 pfils = systemfunc("ls "+path+"/bcsd_0.125deg_${VAR}_amon_${MOD}_${RCP}_r1i1p1_*.nc 2>/dev/null")
   
 if (.not.ismissing(pfils(0)) .and. .not.ismissing(hfils(0))); historical and present must be available
    
 ; concatenate historical files first

 print("Opening files: "+hfils)
 fils = addfiles(hfils,"r")
 ListSetType(fils,"cat")
 historicaltime = fils[:]->time
 xhistorical = fils[:]->$VAR
 xhistorical!0 = "lat"
 xhistorical!0 = "lon"
 delete(fils)
   
 ; for lat and lon

 fil = addfile(hfils(0),"r")
 lat = fil->latitude
 lon = fil->longitude
 lat!0 = "lat"
 lon!0 = "lon"
 nlat = dimsizes(lat)
 nlon = dimsizes(lon)
 delete(fil)
     
 nhtim = dimsizes(historicaltime)
 option = 1
 option@calendar = historicaltime@calendar
 utc_date = cd_calendar(historicaltime(nhtim-1), 0)
 yyyy  = tointeger(utc_date(:,0))
 mm  = tointeger(utc_date(:,1))
 dd  = tointeger(utc_date(:,2))
 hh  = tointeger(utc_date(:,3))
   
 ; determine the next month and year 
 mm = mm + 1
 if (mm .gt. 12 ) then
   yyyy = yyyy + 1 
  mm = 1
 end if
 
 print("Opening files: "+pfils) 

 nfils = dimsizes(pfils)
 fils = addfiles(pfils,"r")
 ListSetType(fils,"cat")
 projectedtime = fils[:]->time
 nprojtim = dimsizes(projectedtime)
 xprojected = fils[:]->$VAR
 xprojected!1 = "lat"
 xprojected!2 = "lon"
 delete(fils)
 
 if ( nfils .gt. 1) then

   delete( projectedtime )
   projectedtime = new( (/nprojtim/), "double")
   fil = addfile(pfils(0),"r")
   xtp =  fil->time
   projTimUnit = xtp@units
   projCal = xtp@calendar
   option = 1
   option@calendar = projCal
   tstart = xtp(0) 
   delete([/fil, xtp/]) 
   fil = addfile(pfils(nfils-1),"r")
   xtp = fil->time
   ntim = dimsizes(xtp)
   tend = xtp(ntim-1)
   delete([/fil, xtp/]) 
 
 ; files have multiple time units...
 ; ...re-generate the time record to have consistent units    
   utc_date = cd_calendar(tstart, 0)
   ystrt  = tointeger(utc_date(:,0))
   mstrt  = tointeger(utc_date(:,1))
   dd  = tointeger(utc_date(:,2))
   hh  = tointeger(utc_date(:,3))
   
   utc_date = cd_calendar(tend, 0)
   yend  = tointeger(utc_date(:,0))
   mend  = tointeger(utc_date(:,1))
   ;print(mend)
   delete(utc_date)
   
   ii = 0
   imend = 12
   do iy = ystrt, yend
    if (iy .eq. yend) then
      imend = mend
    end if
    do im = mstrt, imend
     projectedtime(ii) = cd_inv_calendar(iy,im,dd,hh,0,0,projTimUnit,option)
     if(im .eq. 12) then
       mstrt = 1
     end if
     ii = ii + 1
    end do
   end do

 end if ; if multiple files
  
 start = cd_inv_calendar(yyyy,mm,dd,hh,0,0,projectedtime@units,option)
 is = closest_val(start,projectedtime)
 nptim = dimsizes(projectedtime(is::))
 ntottim = nptim + nhtim
 
 ; determine last month of projected
 utc_date = cd_calendar(projectedtime(dimsizes(projectedtime)-1), 0)
 endmon = utc_date(:,1)
 
 xall = new((/ntottim,nlat,nlon/),"float",-999.999)
 xall!0 = "time"
 xall!1 = "lat"
 xall!2 = "lon"
 
 xalltime = new((/ntottim/),"double")
 xalltime!0 = "time"
 
 ; put historical and projected in to "xall"
 xall(0:(nhtim-1),:,:) = (/xhistorical(:,:,:)/)
 xall(nhtim::,:,:) = (/xprojected(is::,:,:)/)
   
 ; re-generate the entire time record
 
 utc_date = cd_calendar( historicaltime(0), 0)
 ystrt  = tointeger( utc_date(:,0) )
 mstrt  = tointeger( utc_date(:,1) )
 
 utc_date = cd_calendar( projectedtime( dimsizes(projectedtime) - 1 ), 0)
 yend  = tointeger( utc_date(:,0) )
 mend  = tointeger( utc_date(:,1) )
 
 option@calendar = "standard"
 ii = 0
 imend = 12
 do iy = ystrt, yend
  if (iy .eq. yend) then
    imend = mend
  end if
  do im = mstrt, imend
   xalltime(ii) = cd_inv_calendar( iy, im, 15, 12, 0, 0, "days since 1850-01-01 00:00:00", option )
   utc_date = cd_calendar( xalltime(ii), 0 )
   print( "Re-working time... Year: "+utc_date(:,0)+"  Month: "+utc_date(:,1) )
   delete( utc_date )
   if( im .eq. 12 ) then
     mstrt = 1
   end if
   ii = ii + 1
  end do
 end do

 ; In some cases, HADGEM2-ES is missing December 2099--extrapolate 2069-2098 trend to fill in this month
 if ( "${MOD}" .eq. "hadgem2-es" .and. mend .eq. 12 ) then

   ; test to see if there is data for this month
   ; find Topeka, KS ; lat = 39.055, lon=-95.689
     i = closest_val_AnyOrder(39.055,lat)
     j = closest_val_AnyOrder((360-95.689),lon)
    if ismissing(xall(dimsizes(xall(:,0,0))-1,i,j)) then ; delete empty months
      ; delete last month in record
      xtemp = new((/dimsizes(xall(:,0,0))-1,nlat,nlon/),"float",-999.999)
      xtemp = xall(0:dimsizes(xall(:,0,0))-2,:,:)

      delete(xall)
      xall = xtemp
      delete(xtemp)
      xall@_FillValue = -999.999
      ; delete last time record
      xtemp_time = new((/dimsizes(xalltime)-1/),"double")
      xtemp_time = xalltime(0:dimsizes(xalltime)-2)
      delete(xalltime)
      xalltime = new((/dimsizes(xtemp_time)/),"double")
      xalltime!0 = "time"
      xalltime = xtemp_time
      delete(xtemp_time)
      mend = 11
    else
      mend = 12
    end if
  end if

  if mend .eq. 11 then ; model end in november
       print("Going to extrapolate to fill in missing December...")
       mend = 12
       xdec = new((/30,nlat,nlon/),"float",-999.999)
       rate = new((/nlat,nlon/),"float",-999.999)
       it = 0
       xyears = ispan(2069,2098,1) ; get rate of Dec increase over these years
       do iy=2069,2098
         ftime = cd_inv_calendar(iy,12,dd,hh,0,0,xalltime@units,option)
         iloc = closest_val(ftime,xalltime)
         xdec(it,:,:) = xall(iloc,:,:) ; aggregate decembers
         it = it + 1
       end do
       xdec!0 = "time"
       xdec!1 = "lat"
       xdec!2 = "lon"
 
       rate(:,:) = regCoef(xyears,xdec(lat|:,lon|:,time|:))
       rate!0 = "lat"
       rate!1 = "lon"
       rate&lat = lat
       rate&lon = lon

       yintercept = onedtond( rate@yintercept, (/nlat,nlon/) )
 
       xtemp = new((/dimsizes(xalltime)+1,nlat,nlon/),"float",-999.999) ; add 1 for december 2099
       xtemp_time = new((/dimsizes(xalltime)+1/),"double")

       xtemp_time(0:dimsizes(xalltime)-1) = xalltime
       xtemp_time(dimsizes(xalltime)) = cd_inv_calendar(2099,12,dd,hh,0,0,xalltime@units,option)
    
       delete(xalltime)
       xalltime = new((/dimsizes(xtemp_time)/),"double")
       xalltime!0 = "time"
       xalltime = xtemp_time
 
       delete(xtemp_time)
       xtemp(0:dimsizes(xalltime)-2,:,:) = xall
       xtemp(dimsizes(xalltime)-1,:,:) = (rate*2099) + yintercept

       delete(xall)
       xall = xtemp
       xall!0 = "time"
       xall!1 = "lat"
       xall!2 = "lon"
       xall&time = xalltime

     end if

 ; Anomalize all output to 1981-2010

 fstart = cd_inv_calendar( 1981, 1, dd, hh, 0, 0, xalltime@units, option )
 fend = cd_inv_calendar( 2010, 12, dd, hh, 0, 0, xalltime@units, option )
 ist = closest_val( fstart, xalltime )
 ied = closest_val( fend, xalltime )
 xall_norm = clmMonTLL( xall(ist:ied,:,:) ) ; model 1981-2010 monthly average

 print("")
 print("Calculating monthly anomalies...")
 print("")
 xall_anom = calcMonAnomTLL( xall, xall_norm ) ; monthly avg anomalies
   
 ; write to a netCDF file

 ncDir = "$ncDir"
 system( "if ! test -d " + ncDir +" ; then mkdir -p " + ncDir + " ; fi" )

 ncName = "${VAR}_mon_bcsd_${RCP}_${MOD}_"+ystrt+sprinti("%0.2i",mstrt)+"-"+yend+sprinti("%0.2i",mend)+".nc"
 ncFil = ncDir +"/"+ ncName
 
 if ( isfilepresent(ncFil) ) then
     system( "rm "+ncFil)
 end if
 print( "Going to write: "+ncFil )
 
 ncdf = addfile( ncFil, "c" )
 setfileoption (ncdf, "DefineMode", True ) ; explicitly start file definition mode.
 fileattdef ( ncdf, xprojected ) ; global attributes
 globAtt = True
 strP0p = str_insert( pfils, " ", 0 )
 strP0h = str_insert( hfils, " ", 0 )
 globAtt@sourcefiles_historical = str_concat( strP0p )
 globAtt@sourcefiles_future = str_concat( strP0h )
 delete( [/strP0p,strP0h/] )
 globAtt@Conventions = "None"
 globAtt@frequency = "monthly average"
 globAtt@creation_date = systemfunc ( "date" )
 fileattdef( ncdf, globAtt) ; update attributes
 
 dimNames = (/"time","lat","lon"/)
 dimSizes = (/-1,nlat,nlon/)
 dimUnlim = (/True,False,False/)
 filedimdef( ncdf, dimNames, dimSizes, dimUnlim ) ; make time UNLIMITED dimension
 
 filevardef( ncdf, "time", typeof(xalltime), getvardims(xalltime) )
 filevardef( ncdf, "lon", typeof(lon), getvardims(lon) )
 filevardef( ncdf, "lat", typeof(lat), getvardims(lat) )
 filevardef( ncdf, "$VAR", typeof(xprojected), getvardims(xall) )
 
 filevarattdef( ncdf, "time", xalltime )
 filevarattdef( ncdf, "lon", lon )
 filevarattdef( ncdf, "lat", lat )
 filevarattdef( ncdf, "$VAR", xprojected )
 
 varAtt = True
 varAtt@history = "Processed by DJ Rasmussen; email: d.m.rasmussen.jr@gmail.com"
 varAtt@actual_range = (/min(xall_anom), max(xall_anom)/)
 varAtt@time = xalltime(0)
 varAtt@comment = "Monthly anomaly with respect to 1981-2010 monthly average"
 filevarattdef( ncdf, "$VAR", varAtt )
 
 print("Writing:" +ncFil)
 ncdf->time  =   (/xalltime/)
 ncdf->lat = (/lat/)
 ncdf->lon = (/lon/)
 ncdf->$VAR  =  (/xall_anom/)
 
 setfileoption( ncdf, "DefineMode", False ) ; explicitly exit file definition mode.

 end if ; if historical and projected exist
 end

EOF

$ncl concat_bcsd.ncl.tmp >> ncl.log
rm concat_bcsd.ncl.tmp

done # each model
done # each variable
done # each scenario



