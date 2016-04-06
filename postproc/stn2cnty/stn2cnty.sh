#!/bin/bash

# Map downscaled dailies at GHCN stations to CONUS county
# centroids for all variables

# location of NCL
ncl=/usr/local/ncl-6.3.0/bin/ncl

# location of concatenated historical+projected BCSD files
root=/home/dmr/jnk

rm ncl.log

NSET=1 # number of weather realizations (must be <= 20)

for RCP in "rcp26" "rcp45" "rcp60" "rcp85"; do
for VAR in "tas" "pr" "tasmin" "tasmax"; do

fils=`ls $root/smme/*/${RCP}/*/${VAR}/*ghcnd*19810101-*.nc 2>/dev/null`


if [ ${#fils[@]} -eq 0 ]; then continue; fi

for FIL in ${fils[@]}; do

MOD=`echo $FIL | rev | cut -d'/' -f3 | rev`

#if [ ! -z `ls /home/dmr/jnk/smme/001/${RCP}/${MOD}/${VAR}/smme_county_daily_001_${VAR}_${MOD}_${RCP}_19810101-*.nc 2>/dev/null` ]; then continue; fi

echo "Working with [ RCP: $RCP ] [ Model: $MOD ] [ Variable: $VAR ] ..."

cat > stn2cnty.ncl.tmp << EOF
; stn2cnty.ncl

; rasmussen; last updated: Tue 29 Dec 2015 05:40:44 AM PST
; Maps county centroid to nearest GHCN station county lat/lon

 begin

 yr1 = 1981

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

 ; Random weather realizations (must be less than 20)
 do iset = 1, $NSET
      
     print("Opening file: "+"${FIL}")
     xfil = addfile("${FIL}","r")
     xtime = xfil->time

     ; determine last year in time record
     utc_date = cd_calendar(xtime(dimsizes(xtime)-1), 0)
     yr2 = tointeger(utc_date(:,0))

     xvar = xfil->${VAR}
     stn_name = xfil->name
     stn_state = xfil->state
     stn_lat = xfil->lat
     stn_lon = xfil->lon

     ; generate time record in the form of yyyyddd
     outtime = new((/((yr2-yr1)+1)*365/),"integer",-999)
     outtime!0 = "time"
     outtime@calendar = "noleap"
     outtime@format = "yyyyddd"
     
     print("Generating time record...") ; yyyyddd format
     it = 0
     do iy=yr1,yr2
      ddd = 1
      do im=1,12
       do id=1,days_in_month(2001,im) ; force non-leap year
        outtime(it) = toint(tostring(iy)+sprinti("%0.3i", ddd))
        it = it + 1
        ddd = ddd + 1
       end do
      end do
     end do

     print("Mapping GHCN stations to county lat/lon centroids...")
     print("")
     isite = new((/ncnty/),"integer",-999)
     ii = 0
     do i=0, ncnty-1
      if cnty_lon(i)+360 .gt. 235 then
        r_polar = 6370000.
        r_equi = 6370000.
        lat_dist = 2*(3.141532)*r_polar*(abs(stn_lat-cnty_lat(i))/360.)
        lon_dist = 2*(3.141532)*r_equi*(abs((stn_lon+360)-(cnty_lon(i)+360))/360.)
        dist = sqrt((lat_dist)^2 + (lon_dist)^2)
        indx = minind(dist)
        print("Mapping GHCN site: "+stn_name(indx)+", "+stn_state(indx)+" to county/ equivalent: "+cnty_name(i)+", "+cnty_state(i))
        isite(ii) = indx
        ii = ii + 1
      end if
     end do
     ncnty_no_AK_HI = ii

     xout = new((/((yr2-yr1)+1)*365,ncnty_no_AK_HI/),"float",-999.99)
     xout!0 = "time"
     xout!1 = "county"
     do icnty=0, ncnty_no_AK_HI-1
      xout(:,icnty) = (/xvar(:,isite(icnty))/)
     end do

     ; write to netCDF
     ncDir = "${root}/smme/"+sprinti("%0.3i", iset)+"/${RCP}/${MOD}/${VAR}"
     system("if ! test -d " + ncDir +" ; then mkdir -p " + ncDir + " ; fi")
     ncFil = "smme_county_daily_"+sprinti("%0.3i", iset)+"_${VAR}_${MOD}_${RCP}_"+yr1+"0101-"+yr2+"1231.nc"
     ncFile = ncDir + "/" +ncFil

     print("Writing:" +ncFile)
     print("")
     setfileoption("nc","Format","NetCDF4") ; explicitly start file definition mode.
     ncdf = addfile(ncFile,"c")

     globAtt = True
     globAtt@creation_date = systemfunc ( "date" ) ; update creation time
     fileattdef(ncdf, globAtt) ; update attributes
     setfileoption(ncdf,"DefineMode",True)

     ; Check for values missing in arrays
     if any(ismissing(xout)) then
      xarr_1d = ndtooned(xout)
      dsizes = dimsizes(xout)
      indices  = ind_resolve(ind(ismissing(xarr_1d)),dsizes)
      ncells = dimsizes(indices)
      do ii=0, ncells(0)-1
        print(" FOUND MISSING! x: " +indices(ii,0)+" y: "+indices(ii,1))
      end do
      printVarSummary(xout)
      exit
     end if

     ; coordinate variables
     dimNames = (/"time","county"/)
     dimSizes = (/-1,ncnty/)
     dimUnlim = (/True,False/)
     filedimdef( ncdf,dimNames,dimSizes,dimUnlim) ; make time UNLIMITED dimension

     filevardef(ncdf,"time",typeof(outtime),getvardims(outtime))
     filevardef(ncdf,"fips",typeof(cnty_fips),getvardims(cnty_fips))
     filevardef(ncdf,"state",typeof(cnty_state),getvardims(cnty_state))
     filevardef(ncdf,"name",typeof(cnty_name),getvardims(cnty_name))
     filevardef(ncdf,"lat",typeof(cnty_lat),getvardims(cnty_lat))
     filevardef(ncdf,"lon",typeof(cnty_lon),getvardims(cnty_lon))
     filevardef(ncdf,"area",typeof(cnty_area),getvardims(cnty_area))
     filevardef(ncdf,"pop",typeof(cnty_pop),getvardims(cnty_pop))
     filevardef(ncdf,"${VAR}",typeof(xout),getvardims(xout))

     filevarattdef(ncdf,"time",outtime)
     filevarattdef(ncdf,"fips",cnty_fips)
     filevarattdef(ncdf,"state",cnty_state)
     filevarattdef(ncdf,"name",cnty_name)
     filevarattdef(ncdf,"lat",cnty_lat)
     filevarattdef(ncdf,"lon",cnty_lon)
     filevarattdef(ncdf,"area",cnty_area)
     filevarattdef(ncdf,"pop",cnty_pop)
     filevarattdef(ncdf,"${VAR}",xout)
      
     varAtt = True
     varAtt@contents = "County level from GHCN stations"
     varAtt@model = "${MOD}"
     varAtt@experiment = "${RCP}"
     varAtt@history = "Processed by DJ Rasmussen for Rhodium Group, LLC; email: d.m.rasmussen.jr@gmail.com"
     varAtt@frequency ="daily"
     varAtt@actual_range = (/ min(xout), max(xout) /)
     varAtt@time = outtime(0)
     filevarattdef(ncdf,"${VAR}",varAtt)

     ncdf->time = (/outtime/)
     ncdf->fips = (/cnty_fips/)
     ncdf->state = (/cnty_state/)
     ncdf->name = (/cnty_name/)
     ncdf->area = (/cnty_area/)
     ncdf->pop = (/cnty_pop/)
     ncdf->lat = (/cnty_lat/)
     ncdf->lon = (/cnty_lon/)
     ncdf->${VAR}  =  (/xout/)
 end do

end
EOF

$ncl stn2cnty.ncl.tmp | tee ncl.log
rm stn2cnty.ncl.tmp

done # each model
done # each variable
done # each scenario

