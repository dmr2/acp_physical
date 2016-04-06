#!/bin/bash


# NOTE: Requires GCM patters of Forced Change

ncl=/usr/local/ncl-6.3.0/bin/ncl

rm ncl.log

# Start Year
YR1=1981

SET=001

ROOT="/home/dmr/jnk"

WGT_DIR="/home/dmr/acp_code/preproc/magicc/modelweights" # SMME global temperature pathway
PARM_DIR="/home/dmr/acp_code/postproc/wetbulb/param" # regression wetbulb parameters
HIST_DIR="/home/dmr/acp_code/postproc/wetbulb/narr" # historical wetbulb from NARR
PATT_DIR="/home/dmr/jnk/patterns" # seasonal patterns


for RCP in 'rcp26' 'rcp45' 'rcp60' 'rcp85'; do

case $RCP in
  rcp26)
  # Models used for SMME pattern scaling
  MODELS="(/"'"gfdl-esm2g"'","'"giss-e2-r"'","'"fgoals-g2"'","'"mpi-esm-lr"'", \
           "'"miroc-esm-chem"'","'"hadgem2-es"'","'"miroc-esm-chem"'"/)"
  # Names for pattern scaled climate model output
  PATTNAME="(/"'"pattern2"'","'"pattern4"'","'"pattern5"'","'"pattern7"'", \
              "'"pattern24"'","'"pattern28"'","'"pattern29"'"/)"
  ;;
  rcp45)
  MODELS="(/"'"gfdl-esm2g"'","'"gfdl-esm2g"'","'"miroc-esm-chem"'","'"hadgem2-ao"'", \
              "'"miroc-esm-chem"'","'"hadgem2-ao"'","'"miroc-esm-chem"'","'"hadgem2-ao"'"/)"
  PATTNAME="(/"'"pattern1"'","'"pattern4"'","'"pattern38"'","'"pattern39"'", \
              "'"pattern40"'","'"pattern41"'","'"pattern42"'","'"pattern43"'"/)"
  ;;
  rcp60)
  MODELS="(/"'"gfdl-esm2m"'","'"fio-esm"'","'"gfdl-esm2m"'","'"fio-esm"'","'"miroc-esm-chem"'",  \
           "'"hadgem2-es"'","'"miroc-esm-chem"'","'"hadgem2-es"'","'"miroc-esm-chem"'",  \
              "'"hadgem2-es"'"/)"
  PATTNAME="(/"'"pattern1"'","'"pattern2"'","'"pattern3"'","'"pattern4"'",  \
               "'"pattern23"'","'"pattern24"'",  \
              "'"pattern25"'","'"pattern26"'","'"pattern27"'","'"pattern28"'"/)"
  ;;
  rcp85)
  MODELS="(/"'"giss-e2-r"'","'"inmcm4"'","'"giss-e2-r"'","'"inmcm4"'","'"miroc-esm-chem"'",  \
           "'"gfdl-cm3"'","'"miroc-esm-chem"'","'"gfdl-cm3"'","'"miroc-esm-chem"'", \
           "'"gfdl-cm3"'","'"miroc-esm-chem"'"/)"
  PATTNAME="(/"'"pattern1"'","'"pattern2"'","'"pattern3"'","'"pattern4"'", \
            "'"pattern38"'","'"pattern39"'","'"pattern40"'",  \
            "'"pattern41"'","'"pattern42"'","'"pattern43"'",  \
            "'"pattern44"'"/)"
  ;;
      *) echo "Not a valid RCP: $RCP" ;exit
esac


OUT=$ROOT"/${RCP}/wetbulb" # output directory

# daily temperature files
FILS=`ls /home/dmr/RUCORE/SMME/${RCP^^}/TAS_DAILY/smme_county_daily_001_tas_pattern*_${RCP}_19810101-2*1231.nc 2>/dev/null`
#echo /home/dmr/RUCORE/SMME/${RCP^^}/TAS_DAILY/smme_county_daily_001_tas_*_${RCP}_19810101-2*1231.nc

#FILS=`ls $ROOT/smme/${SET}/${RCP}/*/tas/smme_county_daily_${SET}_tas_*_${RCP}_*198101*.nc 2>/dev/null`
#echo $ROOT/smme/${SET}/${RCP}/*/tas/smme_county_daily_${SET}_tas_*_${RCP}_*198101*.nc

#FILS=`ls $ROOT/smme/${SET}/${RCP}/*/tas/smme_county_daily_${SET}_tas_*_${RCP}_*198101*.nc 2>/dev/null`
#echo $ROOT/smme/${SET}/${RCP}/*/tas/smme_county_daily_${SET}_tas_*_${RCP}_*198101*.nc

#/home/dmr/jnk/smme/001/rcp85/pattern38/tas/smme_county_daily_001_tas_pattern38_rcp85_19810101-21001231.nc


for FILI in ${FILS[@]}; do

MOD=`echo $FILI | rev | cut -d'_' -f3 | rev`


echo -e "Generating daily max. wetbulb temperature for [Model: $MOD] [Scenario: $RCP] \n"


cat > smme_wetbulb.tmp.ncl << EOF
; rasmussen; last updated: Fri Mar 21 11:02:41 PDT 2014

; Estimates maximum daily wetbulb temperature at the county level for 
; June, July, August only

external wet "f90/calc_wetbulb.so"

begin

          pattname = ${PATTNAME}
          pattmods = ${MODELS}

         ;tavgf = ${ROOT}+"/${MOD}/tas/smme_county_daily_${SET}_tas_${MOD}_${RCP}_19810101-21001231.nc"
          ;tavgf = "/home/dmr/RUCORE/SMME/${RCP^^}/TAS_DAILY/smme_county_daily_${SET}_tas_${MOD}_${RCP}_19810101-22001231.nc"
          tavgf = "$FILI"

          ; get temperature pathway associated with this model
          delim = "	"
          filName = "${WGT_DIR}/${RCP}_2090_SMME.tsv"
          print("Opening global MAGICC global mean temperature pathways: "+filName)
          
          lines = asciiread(filName,-1,"string")
          nlines = dimsizes(lines) - 1   ; First line is a header
          
          ; determine what MAGICC quants we are going to need
           models  = str_get_field(lines(1:),2,delim)
           do i=0,dimsizes(models)-1
            if str_get_cols(models(i),strlen(models(i))-1,strlen(models(i))-1) .eq."*" then
             models(i) = str_get_cols(models(i),0,strlen(models(i))-2)
             if isStrSubset(models(i), "_" ) then
               parts = str_split(models(i),"_")
               models(i) = parts(0)
             end if
            end if
           end do

           print(models)
           indx = ind(models.eq."${MOD}")
           if ismissing(indx) then
             print("I cannot find model: ${MOD} in ${WGT_DIR}/${RCP}_2090_SMME.tsv")
             exit
           end if
          
          
           ; Get MAGICC temperature pathway for this model
           nyear = (2200-${YR1})+1
           tas_path = new((/nyear/),"float",-999.999) ; (year)
           ii = 0
           do iy=$YR1, 2200
             tas_path(ii) = tofloat(str_get_field(lines(1+indx),4+(iy-1950),delim))
             ii = ii + 1
           end do
           delete(lines)


           ; Pattern scaling for forced change
           if (str_get_cols("${MOD}",0,3).eq."patt") then
              idx = ind(pattname.eq."${MOD}")
              pattern = pattmods(idx)
           else
              pattern = "${MOD}"
           end if

           diri = "${PATT_DIR}/tas/"

           ; BCSD grid to county centroid mapping
           lines = asciiread("bcsd_grid2cnty.csv",-1,"string")
           y = toint(str_get_field(lines(1:),3,","))
           x = toint(str_get_field(lines(1:),4,","))
           delete(lines)


           pattf = systemfunc("ls "+diri+"/seasonal_pattern_tas_bcsd_${RCP}_"+pattern+".nc")
           
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
           yr2 = toint(floor(xtime(dimsizes(xtime)-1)/1000))

           nyr = ((yr2-${YR1})+1)
           xtas = new((/nyr*92,ncnty/),"float",-999.99)
           xtas!0 = "time"
           xtas!1 = "county"

           outtime = new((/nyr*92/),"float",-999.99)
           outtime!0 = "time"

           it = 0
           ; We only want JJA days from Tavg record
           ;June 1 is 152, August 31 is 243
           do i=0,dimsizes(xtime)-1
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
           ncFileName = "smme_county_daily_${SET}_wetbulb_${MOD}_${RCP}_${YR1}0601-"+yr2+"0831.nc"

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



$ncl smme_wetbulb.tmp.ncl >> ncl.log
rm smme_wetbulb.tmp.ncl



done # each model
done # each RCP
