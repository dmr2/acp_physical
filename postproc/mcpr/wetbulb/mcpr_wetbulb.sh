#!/bin/bash

ncl=/usr/local/ncl-6.3.0/bin/ncl

rm ncl.log

YR1=1981
YR2=2200

SET=001

ROOT="/home/dmr/jnk/mcpr"

WGT_DIR="/home/dmr/acp_code/preproc/magicc/modelweights" # MCPR global temperature pathway
PARM_DIR="/home/dmr/acp_code/postproc/wetbulb/param" # regression wetbulb parameters
HIST_DIR="/home/dmr/acp_code/postproc/wetbulb/narr" # historical wetbulb from NARR
PATT_DIR="/home/dmr/jnk/patterns" # seasonal patterns

for RCP in 'rcp26' 'rcp45' 'rcp60' 'rcp85'; do

OUT=$ROOT"/${RCP}/wetbulb" # output directory

#FILS=`ls /home/dmr/RUCORE/MCPR/${RCP^^}/TAS_DAILY/mcpr_county_daily_001_tas_*_${RCP}_19810101-2*1231.nc 2>/dev/null`

#echo /home/dmr/RUCORE/MCPR/${RCP^^}/TAS_DAILY/mcpr_county_daily_001_tas_*_${RCP}_19810101-2*1231.nc

#FILS=`ls ROOT/${RCP}/${VAR}/${VAR}_mcpr-*_daily_county_${RCP}_*198101*.nc 2>/dev/null`


FILS=`ls /home/dmr/RUCORE/MCPR/${RCP^^}/TAS_DAILY/mcpr_county_daily_001_tas_*_${RCP}_19810101-2*1231.nc 2>/dev/null`


#if [ -z ${FILS[@]} ]; then continue; fi

for FILI in ${FILS[@]}; do


MOD=`echo $FILI | rev | cut -d'_' -f3 | rev`

if [ -e ~/RUCORE/MCPR/${RCP^^}/WETBULB_DAILY/wetbulb_mcpr-${MOD}_daily_county_${RCP}_19810601-22000831.nc ]; then continue; fi

echo -e "Generating daily max. wetbulb temperature for [Model: $MOD] [Scenario: $RCP] \n"


cat > mcpr_wetbulb.tmp.ncl << EOF
; rasmussen; last updated: Fri Mar 21 11:02:41 PDT 2016

; estimates maximum daily wetbulb temperature at the county level for 
; June, July, August only

external wet "f90/calc_wetbulb.so"

begin

          ; MCPR global temperature pathway and model-pattern pairing
          print("Opening global MAGICC global mean temperature pathways: ${WGT_DIR}/${RCP}_MCPR.tsv")
          rlines = asciiread("${WGT_DIR}/${RCP}_MCPR.tsv",-1,"string")
          delim = "	" ; tab delimited
          patterns = str_get_field(rlines(1:),2,delim) ; skip first line
          pattern_model = patterns(tointeger(${MOD})-1)

          ;tavgf = "${ROOT}/${RCP}/tas/mcpr_county_daily_${SET}_tas_${MOD}_${RCP}_19810101-22001231.nc"
          ;tavgf = "/home/dmr/RUCORE/MCPR/${RCP^^}/TAS_DAILY/mcpr_county_daily_${SET}_tas_${MOD}_${RCP}_19810101-22001231.nc"
          tavgf = "$FILI"

          ; get temperature pathway associated with this model
           delim = "	"
           nyear = (${YR2}-${YR1})+1
           tas_path = new((/nyear/),"float",-999.999)
           i = 0
           do iy=${YR1},${YR2}
            tas_path(i) = tofloat(str_get_field(rlines(tointeger(${MOD})),4+(iy-1950),delim))
            i = i + 1
           end do
           delete(rlines)

           print("Working with Tavg [  RCP: ${RCP} model: ${MOD} pattern: "+pattern_model+" ]")

           diri = "${PATT_DIR}/tas/"
    
           ; BCSD grid to county centroid mapping
           lines = asciiread("bcsd_grid2cnty.csv",-1,"string")
           y = toint(str_get_field(lines(1:),3,","))
           x = toint(str_get_field(lines(1:),4,","))
           delete(lines)

           pattf = systemfunc("ls "+diri+"/seasonal_pattern_tas_bcsd_${RCP}_"+pattern_model+".nc")
           
           print("Opening file: "+tavgf)
           ifil_tas = addfile(tavgf,"r")
           print("Opening file: "+pattf)
           ifil_pat = addfile(pattf,"r")
                  
           ; daily Tavg
           xtas_all = ifil_tas->tas
           xtime = ifil_tas->time
           cnty_fips = ifil_tas->fips
           cnty_state = ifil_tas->state
           cnty_name = ifil_tas->name
           cnty_lat = ifil_tas->lat
           cnty_lon = ifil_tas->lon

           ncnty =  dimsizes(cnty_fips)
           ntim = dimsizes(xtime)
           ;yr2 = toint(floor(xtime(dimsizes(xtime)-1)/1000))
           nyr = ((${YR2}-${YR1})+1)
           xtas = new((/nyr*92,ncnty/),"float",-999.99)
           xtas!0 = "time"
           xtas!1 = "county"

           outtime = new((/nyr*92/),"float",-999.99)
           outtime!0 = "time"

           it = 0
           ; We only want JJA days from Tavg record
           ;June 1 is 152, August 31 is 243
           do i=0,ntim-1
             ddd = xtime(i) - toint(xtime(i)/1000)*1000
             if (ddd .ge. 152 .and. ddd .le. 243) then
                xtas(it,:) = xtas_all(i,:)
                outtime(it) = xtime(i)
                it = it + 1
             end if
           end do

           ; seasonal patterns
           slope = ifil_pat->slope(2,:,:)  ; (season, lat, lon) JJA only

           ntim = dimsizes(xtas(:,0))
           xout = new((/ntim,ncnty/),"float",-999.99)
           xout!0 = "time"
           xout!1 = "county"
           xout@unit = "K"
           xout@long_name = "daily maximum wet bulb temperature"

           do icnty=0, ncnty-1

             ; open linear model parameters
              name = str_concat(str_split(cnty_name(icnty), " "))
              lines = asciiread("${PARM_DIR}/wetbulb_tmax_JJA_"+name+"_"+cnty_state(icnty)+"_scatter_change_pt.txt",-1,"string")
              nparm = 12
              params = new(nparm,"float",-999.99)
              do i = 0, 11
                params(i) = tofloat(str_get_field(lines,i+4,","))
              end do
              delete(lines)

              ; get historical minimum max daily wetbulb temperature
              lines = asciiread("${HIST_DIR}/wetbulb_tavg_JJA_"+cnty_fips(icnty)+"_" \
                            +name+"_"+cnty_state(icnty)+"_daily_1981-2010.txt",-1,"string")
              hmin = min(tofloat(str_get_field(lines(1:),3,",")))

              xwet = new((/ntim/),"float",-999.99)

              nexceed = 0
              wet::estimwet(nparm,params,ntim,nexceed,nyear,slope(y(icnty),x(icnty)),tas_path,xwet,xtas(:,icnty),hmin)
              xout(:,icnty) = (/xwet/)

              delete([/lines,params/])
           end do
          
           ; write to netcdf
           ncDir = "${OUT}"
           ncFileName = "wetbulb_mcpr-${MOD}_daily_county_${RCP}_${YR1}0601-${YR2}0831.nc"

           system("if ! test -d " + ncDir +" ; then mkdir -p " + ncDir + " ; fi")
           ncFile = ncDir+"/"+ncFileName
          
           print("Writing:" +ncFile)
           setfileoption("nc","Format","NetCDF4") ; explicitly start file definition mode.
           ncdf = addfile(ncFile,"c")
          
           globAtt = True
           globAtt@creation_date = systemfunc ( "date" ) ; update creation time
           fileattdef(ncdf, globAtt) ; update attributes
           setfileoption(ncdf,"DefineMode",True)
          
           
           ; QA output before writing
            if any(ismissing(xout)) then
             print("Error! Found missing values in out array!")
             print("Exiting...")
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
            filevardef(ncdf,"wetbulb",typeof(xout),getvardims(xout))
            ;
            filevarattdef(ncdf,"time",outtime)
            filevarattdef(ncdf,"fips",cnty_fips)
            filevarattdef(ncdf,"state",cnty_state)
            filevarattdef(ncdf,"name",cnty_name)
            filevarattdef(ncdf,"lat",cnty_lat)
            filevarattdef(ncdf,"lon",cnty_lon)
            filevarattdef(ncdf,"wetbulb",xout)
            ;
            varAtt = True
            varAtt@contents = "County level from NCDC sites"
            varAtt@model = "${MOD}"
            varAtt@experiment = "${RCP}"
            varAtt@history = "Processed by DJ Rasmussen for Rhodium Group, LLC; email: d.m.rasmussen.jr@gmail.com"
            varAtt@frequency ="daily"
            varAtt@reference = "Daily data from monthly using method from Wood et al. (2002) (section 2.3.2) from JGR-atmospheres"
            varAtt@actual_range = (/ min(xout), max(xout) /)
            varAtt@time = outtime(0)
            filevarattdef(ncdf,"wetbulb",varAtt)
            
            ncdf->time = (/outtime/)
            ncdf->fips = (/cnty_fips/)
            ncdf->state = (/cnty_state/)
            ncdf->name = (/cnty_name/)
            ncdf->lat = (/cnty_lat/)
            ncdf->lon = (/cnty_lon/)
            ncdf->wetbulb  =  (/xout/)

            setfileoption(ncdf,"DefineMode",False) ; explicitly exit file definition mode.

end 
EOF



$ncl mcpr_wetbulb.tmp.ncl >> ncl.log
rm mcpr_wetbulb.tmp.ncl



done # each model
done # each RCP
