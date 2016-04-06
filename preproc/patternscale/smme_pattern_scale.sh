#!/bin/bash

ncl=/usr/local/ncl-6.3.0/bin/ncl

root="/home/dmr/jnk"

WGT_DIR="/home/dmr/acp_code/preproc/magicc/modelweights/"

# Output directory for seasonal means
seasonDir=$root"/merged/seasonal/"

# Start and End years for pattern scaling (must be greater than or equal to 1950)
YR1=1950
YR2=2100

rm ncl.log

for VAR in "tas" "pr" "tasmax" "tasmin" ; do
for RCP in "rcp26" "rcp45" "rcp60" "rcp85"; do

case $RCP in
  rcp26)
  # Models used for SMME pattern scaling
  models=("gfdl-esm2g" "giss-e2-r" "fgoals-g2" "mpi-esm-lr" "miroc-esm-chem" \
             "hadgem2-es" "miroc-esm-chem" )
  # Names for pattern scaled climate model output
  outname=("pattern2" "pattern4" "pattern5" "pattern7" \
              "pattern24" "pattern28" "pattern29")
  # Quantiles here must match pattern targets file
  quants=(4.005 10.000 16.000 30.000 90.000 98.995 98.995) 
  ;;
  rcp45)
  models=("gfdl-esm2g" "gfdl-esm2g" "miroc-esm-chem" "hadgem2-ao" \
              "miroc-esm-chem" "hadgem2-ao" "miroc-esm-chem" "hadgem2-ao")
  outname=("pattern1" "pattern4" "pattern38" "pattern39" \
              "pattern40" "pattern41" "pattern42" "pattern43")
  quants=(4.005 10.00 90.00 90.00 95.00 95.00 98.995 98.995) 
  ;;
  rcp60)
  models=( "gfdl-esm2m" "fio-esm" "gfdl-esm2m" "fio-esm" "miroc-esm-chem"  \
              "hadgem2-es" "miroc-esm-chem" "hadgem2-es" "miroc-esm-chem"  \
              "hadgem2-es")
  outname=( "pattern1" "pattern2" "pattern3" "pattern4"  \
               "pattern23" "pattern24"  \
              "pattern25" "pattern26" "pattern27" "pattern28" )
  quants=(4.005 4.005 10.00 16.00 90.00 90.00 95.00 95.00 98.995 98.995)
  ;;
  rcp85)
  models=( "giss-e2-r" "inmcm4" "giss-e2-r" "inmcm4" "miroc-esm-chem"  \
              "gfdl-cm3" "miroc-esm-chem" "gfdl-cm3" "miroc-esm-chem"  \
              "gfdl-cm3" "miroc-esm-chem" )
  outname=( "pattern1" "pattern2" "pattern3" "pattern4" 
            "pattern38" "pattern39" "pattern40"  \
            "pattern41" "pattern42" "pattern43"  \
            "pattern44" )
  quants=(4.005 4.005 10.000 10.000 84.000 90.000 90.000  \
          84.000 90.000 90.000 95.000 95.000 98.995 98.995)
  ;;
      *) echo "Not a valid RCP: $RCP" ;exit
esac


i=0
for MOD in ${models[@]}; do

# Seasonal pattern file
PATT_FIL=`ls ${root}/patterns/${VAR}/seasonal_pattern_${VAR}_bcsd_${RCP}_${MOD}.nc 2>/dev/null`
if [ -z $PATT_FIL ]; then continue; fi

# Monthly residual file
RESID_FIL=`ls ${root}/residuals/${VAR}/${VAR}_mon_residual_bcsd_${RCP}_${MOD}_*.nc 2>/dev/null`

if [ -z $RESID_FIL ]; then continue; fi

echo -e "Pattern scaling for [Variable: $VAR] [Scenario: $RCP] [Model: ${outname[$i]} (${models[$i]})] \n"

cat >smme_pattern_scale.tmp.ncl << EOF
; smme_pattern_scale.ncl
;
; Written by DJ Rasmussen, Fri Dec 27 12:32:29 PST 2013
; email: d-dot-m-dot-rasmussen-dot-jr-AT-gmail-dot-com
; Rhodium Group, LLC
; 
; This scripts pattern scales seasonal means to monthly using pre-calculated 
; patterns and then writes the model surrogate monthly mean values to a 
; netCDF file.

; References:

; Mitchell, T., Pattern Scaling: An Examination of the Accuracy of the 
; Technique for Describing Future Climates. (2003) Climatic Change.

; Rasmussen and co-authors, Probability-weighted ensembles of U.S. county-level 
; climate projections for climate risk analysis. (submitted) J. Appl. Met. Clim.

; Algorithm is as follows:

; For temperature:

; surrogate_month_Tanom(t2,i,j) = globT_anom_MAGICC(quantile,t1)*pattern(t3,i,j) 
;                                                             + residual(t2,i,j)

; Where t1 = year; t2 = month; t3 = season; i = latitude; j = longitude
; and surrogate_month_Tanom is with respect to each model's 1981-2010 
; period and the pattern and residual are from the same model (the models 
; used are defined below). The regression intercept is not included as it 
; is assumed that there is no local change under no global change.
;
; Note that pattern has units of:
;   local seasonal temperature (degree C/ K) per change in 
;                                      global mean annual temp (degree C/ K)

; For precipitation:
;
; surrogate_month_precip_anom(t2,i,j) = 
;        globT_anom_MAGICC(quantile,t1)*pattern(t3,i,j) + residual(t2,i,j)

; Where t1 = year; t2 = month; t3 = season; i = latitude; j = longitude
; and local_precip_anom is with respect to each model's 1981-2010 period 
; and the pattern and residual are from the same model (the models used 
; are defined below). The regression intercept is not included as it is 
; assumed that there is no local change under no global change.

; Note that pattern has units of:
;   local seasonal temperature ( mm per day) per change in 
;                              global mean annual temperature (degree C/ K)

; CMIP5 model paring that will be scaled to fill in the MAGICC6 
; global mean temperature distribution. 

;  RCP 8.5
;  Model Surrogate       CMIP5 model
; ------------------+---------------------------
;  pattern1_4.005   |  giss-e2-r ; cold, wet
;  pattern2_4.005   |  inmcm4 ; cold, dry
;  pattern3_10.000  |  giss-e2-r 
;  pattern4_10.000  |  inmcm4 
;  pattern38_84.000 |  miroc-esm-chem; hot, dry
;  pattern39_90.000 |  gfdl-cm3
;  pattern40_90.000 |  miroc-esm-chem
;  pattern41_95.000 |  gfdl-cm3
;  pattern42_95.000 |  miroc-esm-chem
;  pattern43_99.995 |  gfdl-cm3
;  pattern44_99.995 |  miroc-esm-chem

;  RCP 6.0
;  Model Surrogate       CMIP5 model
; ------------------+---------------------------
;  pattern1_4.005   |  gfdl-esm2m ; wet
;  pattern2_4.005   |  fio-esm  ; dry
;  pattern3_10.000  |  gfdl-esm2m
;  pattern4_10.000  |  fio-esm
;  pattern23_90.000 |  miroc-esm-chem
;  pattern24_90.000 |  hadgem2-es
;  pattern25_95.000 |  miroc-esm-chem
;  pattern26_95.000 |  hadgem2-es
;  pattern27_98.995 |  miroc-esm-chem
;  pattern28_98.995 |  hadgem2-es

;  RCP 4.5
;  Model Surrogate       CMIP5 model
; -------------------+--------------------------
;  pattern_1_4.005   |  gfdl-esm2g
;  pattern_4_10.000  |  gfdl-esm2m
;  pattern_38_90.000 |  miroc-esm-chem ; wet
;  pattern_39_90.000 |  hadgem2-ao   ; dry
;  pattern_40_95.000 |  miroc-esm-chem
;  pattern_41_95.000 |  hadgem2-ao
;  pattern_42_98.995 |  miroc-esm-chem
;  pattern_43_98.995 |  hadgem2-ao

;  RCP 2.6
;  Model Surrogate       CMIP5 model
; -------------------+--------------------------
;  pattern2_4.005    |  gfdl-esm2g
;  pattern4_10.000   |  giss-e2-r
;  pattern5_16.000   |  fgoals-g2
;  pattern7_30.000   |  mpi-esm-lr
;  pattern24_90.000  |  miroc-esm-chem
;  pattern28_98.995  |  hadgem2-es ; wet
;  pattern29_98.995  |  miroc-esm-chem ; dry

 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
 load "util.ncl"

begin

  month_abbr = (/"","Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep", \
                    "Oct","Nov","Dec"/)

 nyear = $YR2-$YR1+1

; open the global annual mean temperature from MAGICC
 filName = "${WGT_DIR}/${RCP}_2090_SMME.tsv"
 print("Opening global MAGICC global mean temperature pathways: "+filName)

 lines = asciiread(filName,-1,"string")
 nlines = dimsizes(lines) - 1   ; First line is a header

; determine what MAGICC quants we are going to need
  models  = str_get_field(lines(1:),2,"	")
  do i=0,dimsizes(models)-1
   if str_get_cols(models(i),strlen(models(i))-1,strlen(models(i))-1) .eq."*" then
    models(i) = str_get_cols(models(i),0,strlen(models(i))-2)
    if isStrSubset(models(i), "_" ) then
      parts = str_split(models(i),"_")
      models(i) = parts(0)
    end if
   end if
  end do
 indx = ind(models.eq."${outname[$i]}")
 if ismissing(indx) then
   print("I cannot find model: ${MOD} in ${WGT_DIR}/${RCP}_2090_SMME.tsv")
   exit
 end if

 ; Get MAGICC temperature pathway for this model
 magicc_tas = new((/nyear/),"float",-999.999) ; (year)

 ii = 0
 do iy=$YR1, $YR2
   magicc_tas(ii) = tofloat(str_get_field(lines(1+indx),4+(iy-1950),"	"))
   ii = ii + 1
 end do

; get MAGICC start and end years
 magicc_yrs = ispan($YR1,$YR2,1)
 print("MAGICC start year: "+magicc_yrs(0)+" MAGICC end year: "+magicc_yrs(dimsizes(magicc_yrs)-1))

; create arrays

 pattFil = "$PATT_FIL"
 print("Opening CONUS seasonal pattern file for "+str_upper("${MOD}")+" ..." +pattFil)
 fil1 = addfile(pattFil,"r")
 slope_conus = fil1->slope ; seasonal slopes
 global_temp =fil1->global_temp; 30-yr run-ave
 conus_int = fil1->intercept; seasonal y-intercept
 option = 1
 option@calendar = "standard"
 
 resFil = "$RESID_FIL"
 print("Opening CONUS BCSD monthly residual file for "+str_upper("${MOD}")+" ..." +resFil)
 fil2 = addfile(resFil,"r")
 residual = fil2->residual ; monthly residuals
 lat = fil2->lat
 lon = fil2->lon
 
 nmonth = 12
 nseason = 4
 nlon = dimsizes(lon)
 nlat = dimsizes(lat)
 nyear = dimsizes(fil2->time)/12

 surr_mon = new((/nyear*nmonth,nlat,nlon/),"float",1.E+20)
 surr_mon!0 = "time"
 surr_mon!1 = "lat"
 surr_mon!2 = "lon"
   
 surr_time = new((/nmonth*nyear/),"double",1.E+20)
 surr_time!0 = "time"
 surr_time@units = "days since 1950-01-01 00:00:00"
 surr_time@calendar = "standard"

 ii = 0
 im = 0 ; count months
 do iyear = 0, nyear-1
   do im = 1, 12
    if im .eq. 12 .or. im .le. 2 then
     season = "DJF"
     iseason = 0
    else if im .ge. 3 .and. im .le. 5 then
     season = "MAM"
     iseason = 1
    else if im .ge. 6 .and. im .le. 8 then
     season = "JJA"
     iseason = 2
    else if im .ge. 9 .and. im .le. 11 then
     season = "SON"
     iseason = 3
    else
     print("Error. im is "+im)
     exit
    end if
    end if
    end if
    end if
  
    surr_mon(ii,:,:) = magicc_tas(iyear)*slope_conus(iseason,:,:) + residual(ii,:,:)


    surr_time(ii) = cd_inv_calendar($YR1+iyear,im,15,12,0,0,"days since 1950-01-01 00:00:00",option)
    date_str = month_abbr(im) + " "  + sprinti("%0.4i", $YR1+iyear)
    print("Created pattern scaled data for: "+date_str)
    ii = ii + 1 
   end do; month loop
   print(magicc_yrs(iyear)+" "+magicc_tas(iyear)+" ${quants[$i]}")
 end do; year loop


; write synthetic data to disk as a netCDF file
 ncDir = "${root}/merged/${outname[$i]}/${RCP}/${VAR}"
 system("if ! test -d " + ncDir +" ; then mkdir -p " + ncDir + " ; fi")
 ncFil = "${VAR}_mon_bcsd_${RCP}_${outname[$i]}_"+$YR1+"01"+"-"+($YR1+iyear-1)+"12.nc"
 NCFILE = ncDir + "/" +ncFil
 
 print("Writing:" +NCFILE)
 setfileoption("nc","Format","NetCDF4") ; explicitly start file definition mode.
 ncdf = addfile(NCFILE,"c")
 
 ;fileattdef ( ncdf, fmod ) ; global attributes
 globAtt = True
 globAtt@Conventions = "None"
 globAtt@frequency = "daily"
 globAtt@creation_date = systemfunc ( "date" )
 fileattdef( ncdf, globAtt ) ; update attributes

; predefine the coordinate variables
 dimNames = (/"time","lat","lon"/)
 dimSizes = (/-1,nlat,nlon/)
 dimUnlim = (/True,False,False/)
 filedimdef( ncdf,dimNames,dimSizes,dimUnlim ) ; make time UNLIMITED dimension
 
 filevardef(ncdf,"time",typeof(surr_time),getvardims(surr_time))
 filevardef(ncdf,"lon",typeof(lon),getvardims(lon))
 filevardef(ncdf,"lat",typeof(lat),getvardims(lat))
 filevardef(ncdf,"${VAR}",typeof(surr_mon),getvardims(surr_mon))
;
 filevarattdef(ncdf,"time",surr_time)
 filevarattdef(ncdf,"lon",lon)
 filevarattdef(ncdf,"lat",lat)
 filevarattdef(ncdf,"${VAR}",surr_mon)
;
 varAtt = True
 varAtt@contents = "Gridded monthly pattern-scaled projections"
 varAtt@model = "Derived from model: ${MOD}"
 varAtt@experiment = "${RCP}"
 varAtt@pattern_name = "${outname[$i]}"
 varAtt@info = "q${quants[$i]}"
 varAtt@reference_period = "anomalized to 1981-2010"
 varAtt@history = "Monthly data produced by DJ Rasmussen (Rhodium Group, LLC); email: d.m.rasmussen.jr@gmail.com"
 varAtt@reference = "Method from Mitchell et al. (2003) from Climatic Change"
 varAtt@actual_range = (/ min(surr_mon), max(surr_mon) /)
 varAtt@time = surr_time(0)
;
 filevarattdef(ncdf,"${VAR}",varAtt)
;
 ncdf->time = (/surr_time/)
 ncdf->lat = (/lat/)
 ncdf->lon = (/lon/)
 ncdf->${VAR}  =  (/surr_mon/)
;   
end
EOF

$ncl smme_pattern_scale.tmp.ncl | tee ncl.log
rm smme_pattern_scale.tmp.ncl

i=`expr $i + 1`

done # model
done # RCP
done # variable


