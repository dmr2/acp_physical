
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"   

undef("clmDayTS")
function clmDayTS (x[*][*]:numeric, yyyyddd:integer)   
;
; calculate the mean Annual Cycle from daily data. 
; The return array will gave the raw climatology at each NCDC site
;
;              x(time,site)  <==== input dimension order
;              x!0 = "time"     <==== time is in days
;              x!1 = "site"
;
;    non-Leap  yyyyddd
;              1905001 =>  Jan  1, 1905
;              1905032 =>  Feb  1, 1905
;              1905059 =>  Feb 28, 1905
;              1905060 =>  Mar  1, 1905
;              1905365 =>  Dec 31, 1905
;
; Usage:  xClmDay = clmDAY_TLL (x, yyyyddd)
; -------

local dimx, ntim, nsite, ndys, days, clmDay, ndy, indx, year_day, nFill
begin

   dimx  = dimsizes (x)

   ntim  = dimx(0)
   nsite  = dimx(1)
   ndys  = 365                                    ; allow for leap year

   days  = yyyyddd - (yyyyddd/1000)*1000           ; strip year info [yyyy]

   clmDay= new((/ndys,nsite/),typeof(x), getFillValue(x) ) ; daily climatology
;
; Compute averages for each sequential day of the year. 
;
   do ndy=0,ndys-1;2                                 ; ndy=0->364 ==> day of year 1->365       
      indx = ind( days.eq.(ndy+1) )                ; indx:indx Mar 2012, handle 1 element
      if (.not.ismissing(indx(0))) then
          nindx = dimsizes(indx)
        if (nindx.eq.1) then                       ; force 3rd dimension
          clmDay(ndy,:) = dim_avg_n(x(indx:indx,:), 0) ; 5 Nov 2009 => dim_avg_n
        else
          clmDay(ndy,:) = dim_avg_n(x(indx,:), 0)  ; 5 Nov 2009 => dim_avg_n
        end if
      end if
      delete(indx)
   end do

  ;;For leap years only
  ;                                                 ; nominal day 366 
  ;                                                 ; ave(31 Dec + 1 Jan)=leap
  ; clmDay(ndys-1,:,:) = (clmDay(0,:,:) + clmDay(ndys-2,:,:))*0.5

   nFill = num(ismissing(clmDay))
   if (nFill.eq.0) then
       delete(clmDay@_FillValue)
   end if

   clmDay@long_name   = "Daily Climatology"
   if (isatt(x,"long_name")) then
       clmDay@long_name = clmDay@long_name +": "+x@long_name
   end if
   if (isatt(x,"units")) then
       clmDay@units     = x@units
   end if
   clmDay@information = "Raw daily averages across all years"
   clmDay@smoothing   = "None"

   year_day           = ispan(1,ndys,1)
   year_day@long_name = "day of year"
   year_day@units     = "ddd"

   clmDay!0           = "year_day"
   clmDay&year_day    =  year_day
   
   copy_VarCoords(x(0,:), clmDay(0,:))   ; trick
   delete(clmDay@year_day)                   ; clean up
 
   return (clmDay)
end
undef("calcDayAnomTS")
function calcDayAnomTS (x[*][*]:numeric, yyyyddd:integer, clmDay[*][*]:numeric)   

; Remove the Annual Cycle from daily data. 
; On return x will consist of deviations from each day's long term mean.
;
;              x(time,site)  <==== input dimension order
;              x!0 = "time"     <==== time is in days
;              x!1 = "site"
;
;    non-Leap  yyyyddd
;              1905001 =>  Jan  1, 1905
;              1905032 =>  Feb  1, 1905
;              1905059 =>  Feb 28, 1905
;              1905060 =>  Mar  1, 1905
;              1905365 =>  Dec 31, 1905

; Usage:  xAnom = calcDayAnomTS (x, yyyyddd, clmDay)

local dimx, ntim, nsite, ndys, days, xAnom, nt
begin

   dimx  = dimsizes (x)

   ntim  = dimx(0)
   nsite  = dimx(1)
   ndys  = 365                               

   days  = yyyyddd - (yyyyddd/1000)*1000      ; strip year info [yyyy]

; quick check [debug]
  ;if (.not.all(days(0:ndys-2).eq.clmDay&year_day(0:ndys-2))) then
  ;    print("calcDayAnomTLL: day mismatch")
  ;end if

; loop on every day and compute difference. 

   xAnom = (/ x  /)                                 ; create xAnom
   do nt=0,ntim-1
      xAnom(nt,:) = x(nt,:) - clmDay(days(nt)-1,:)   ; -1 for 0-based subscript
   end do

   if (isatt(x,"long_name")) then
       xAnom@long_name = "Anomalies: "+x@long_name
   else
       xAnom@long_name = "Anomalies from Daily Climatology"
   end if
   if (isatt(x,"units")) then
       xAnom@units = x@units
   end if

   copy_VarCoords(x, xAnom)

   return(xAnom)
end
undef("dlyClimAnomAdd")
function dlyClimAnomAdd (x[*][*]:numeric, yyyyddd:integer, norm[*][*]:numeric)   

; Add NCDC daily normals to the model 1981-2010 anomalies
; On return x will consist of additions from each day's 1981-2010 normal.
;
;              x(time,site)  <==== input dimension order
;              x!0 = "time"     <==== time is in days
;              x!1 = "site"
;
;    non-Leap  yyyyddd
;              1905001 =>  Jan  1, 1905
;              1905032 =>  Feb  1, 1905
;              1905059 =>  Feb 28, 1905
;              1905060 =>  Mar  1, 1905
;              1905365 =>  Dec 31, 1905

; Usage:  xAnom + NCDC Normals = dlyClimAnomAdd (x, yyyyddd, clmDay)

local dimx, ntim, nsite, ndys, days, xOut, nt
begin

   dimx  = dimsizes (x)

   ntim  = dimx(0)
   nsite  = dimx(1)
   ndys  = 365                               

   days  = yyyyddd - (yyyyddd/1000)*1000      ; strip year info [yyyy]

; quick check [debug]
  ;if (.not.all(days(0:ndys-2).eq.clmDay&year_day(0:ndys-2))) then
  ;    print("calcDayAnomTLL: day mismatch")
  ;end if

; loop on every day and compute difference. 

   xOut = (/ x  /)                                 ; create xAnom
   do nt=0,ntim-1
      xOut(nt,:) = x(nt,:) + norm(days(nt)-1,:)   ; -1 for 0-based subscript
   end do

   if (isatt(x,"long_name")) then
       xOut@long_name = "Anomalies + NCDC normals: "+x@long_name
   else
       xOut@long_name = "Anomalies + NCDC normals from Daily Climatology"
   end if
   if (isatt(x,"units")) then
       xOut@units = x@units
   end if

   copy_VarCoords(x, xOut)

   return(xOut)
end


;**********************************************************************
; Calculate long term monthly means (monthly climatology)
;
; The time dimension must be a multiple of 12
;
;    x(time,site)  <==== INPUT DIMENSION ORDER
;    x!0 = "time"
;    x!1 = "site"
;
undef("clmMonTS")
function clmMonTS (x[*][*]:numeric)
local dimx, ntim, nsite, i, nmo, nmos, monAveLLT
begin
  dimx = dimsizes(x)
  ntim = dimx(0)

  nmos  = 12
  modCheck ("clmMonTS", ntim, nmos)   ; error check

  nsite = dimx(1)
;
; Compute all 12 monthly averages first. 
;
   aveMonth = new((/nmos,nsite/),typeof(x) \
                                    ,getFillValue(x))
   do nmo=0,nmos-1
      aveMonth(nmo,:) = dim_avg_n(x(nmo:ntim-1:nmos,:), 0)
   end do

; copy attributes and add a couple of informational attributes
; only need to do this if plotting or writing to netCDF file

   copy_VarAtts (x,aveMonth)
   aveMonth@time_op_ncl = "Climatology: "+ (ntim/nmos) +" years"
   aveMonth@info    = "function clmMonLLT: contributed.ncl"

   aveMonth!0     = "month"           ; create a "month" named dim
   aveMonth&month = ispan(0,nmos-1,1) ; create a month coord var

   return (aveMonth)
end

;**********************************************************************
; Calculate anomalies from climatology [Remove annual cycle]
;
; Subtract the the long term means from each "month".
; On return x will consist of deviations from  each "month's" long term mean.
;
; The time dimension must be a multiple of 12
;
;     x(time,site)  <==== INPUT DIMENSION ORDER
;     x!0 = "time"
;     x!1 = "site"
;     xAve(12,site) <==== THE 12 is Jan, Feb, .., Dec
;
; Usage:  x     = calcMonAnomTS (x,xAve)  ; overwrites "x"
;         xAnom = calcMonAnomTS (x,xAve)  ; creates xAnom as new variable
; where   xAve  = clmMonTS (x)           

undef("calcMonAnomTS")
function calcMonAnomTS (x[*][*]:numeric, xAve[12][*]:numeric)
local dimx, ntim, yr, nmos, xAnom
begin
   dimx  = dimsizes (x)
   ntim  = dimx(0)

   nmos  = 12
   modCheck ("calcMonAnomTS", ntim, nmos)   ; error check

; Now loop on every year and compute difference. 
; The [yr:yr+nmos-1] strips out 12 months for each year. [array notation]

   xAnom = x          ; variable to variable copy [meta data]
   do yr=0,ntim-1,nmos
      xAnom(yr:yr+nmos-1,:) = (/ x(yr:yr+nmos-1,:)- xAve /)
   end do

; Create an informational attribute:  

   xAnom@anomaly_op_ncl  = "Anomalies from Annual Cycle: calcMonAnomTS: contributed.ncl"

   return (xAnom)
end

;***********************************************************************;
; Given a color map and the number of desired colors, this function 
; returns an array of color indexes that nicely span the full colormap.
;
; For a named colormap, the first two color values are not used,
; because these are the foregfloor/backgfloor colors.
;
; This function is very similar to the span_color_rgba function,
; which returns RGBA values. 
;
; The colormap can be a named colormap, like "rainbow", or an array
; of RGB (n,3) or RGBA (n,4).
;***********************************************************************
undef("span_color_indexes")
function span_color_indexes(cnlvls[*]:numeric,cmapt)
local ncols, lcount, fmin, fmax, fcols, icols, cmap
begin
  if(isstring(cmapt)) then
     cmap = read_colormap_file(cmapt)
  else if(isnumeric(cmapt)) then
    dims = dimsizes(cmapt)
    if(dims(0).lt.3.or.dims(0).gt.256.or..not.any(dims(1).ne.(/3,4/))) then
      print ("Error: span_color_indexex: cmap must be an n x 3 or n x 4 array of RGB or RGBA values, or a valid color map name")
      return(new(1,integer))   ; return missing
    end if
    cmap = cmapt
  else
    print ("Error: span_color_indexex: cmap must be an n x 3 or n x 4 array of RGB or RGBA values, or a valid color map name")
  end if
  end if

  ncols  = dimsizes(cmap(:,0))
  lcount = dimsizes(cnlvls)
;
; Start at index 0 and end at ncols-1 (the full range of the
; color map.
;
  minix = 0
  maxix = ncols-1

  fmin = new(1,float)    ; to make sure we get a missing value (?)
  fmax = new(1,float)
  fmin = minix
  fmax = maxix
  fcols = fspan(fmin,fmax,lcount+1)
  icols = tointeger(fcols + 0.5)
  return(icols)
end

undef("span_color_indexes2")
function span_color_indexes2(cmapt,ncolors)
local ncols, fmin, fmax, fcols, icols, cmap
begin
  if(isstring(cmapt)) then
     cmap = read_colormap_file(cmapt)
  else if(isnumeric(cmapt)) then
    dims = dimsizes(cmapt)
   if(dimsizes(dims).ne.2.or.dims(0).lt.3.or.dims(0).gt.256.or.\
       .not.any(dims(1).ne.(/3,4/))) then
      print ("Error: span_color_indexes: cmap must be an n x 3 or n x 4 array of RGB or RGBA values, or a valid color map name")
      return(new(1,integer))   ; return missing
    end if
    cmap = cmapt
  else
    print ("Error: span_color_indexes: cmap must be an n x 3 or n x 4 array of RGB or RGBA values, or a valid color map name")
  end if
  end if

  ncols  = dimsizes(cmap(:,0))
;
; Start at index 0 and end at ncols-1 (the full range of the
; color map.
;
  minix = 0
  maxix = ncols-1

  fmin = new(1,float)    ; to make sure we get a missing value (?)
  fmax = new(1,float)
  fmin = minix
  fmax = maxix
  fcols = fspan(fmin,fmax,ncolors)
  icols = tointeger(fcols + 0.5)
  if(isstring(cmapt)) then
    return(icols+2)
  else
    return(icols)
  end if
end


undef("QA_ndarray"); quickly QA's an array for missing or negative values
function QA_ndarray(xarr:numeric,ltz:logical)
local ndim, xflat
begin

  xflat=ndtooned(xarr) ;flatten array
  ndim=dimsizes(dimsizes(xarr))
  lqa=False

  if any(ismissing(xarr)) then
    print("Found missing values in array!")
    print("indicies are:" + ind(ismissing(xarr)) )
    print("Exiting...")
    exit
  end if

  if ltz then; check for values less than zero...
    if any(xarr.lt.0) then
     print("Found values less than zero in array!")
     print("indicies are:" + ind_resolve(ind(xflat.lt.0),ndim))
     print("Exiting...")
     exit
    end if
  end  if


  delete(xflat)
  return(lqa)
end
function closest_val_AnyOrder(xVal[1]:numeric, x:numeric)
local xAbsDif, xMinVal, iClose
begin
  xAbsDif = abs(xVal-x)
  iClose  = minind(xAbsDif)
  return(iClose)   ; original doc says "first occurence"
end

undef("getind_closest_latlon2d")
function getind_closest_latlon2d(lat2d[*][*]:numeric,lon2d[*][*]:numeric \
                        ,LAT:numeric, LON:numeric, xvar[*][*]:numeric \
                        ,sname:string, sstate:string)
local ij, lat1d, lon1d
begin

 ij = new ( (/2/) , "integer")

 lat1d  = ndtooned( lat2d )
 lon1d  = ndtooned( lon2d )
 n2d    = dimsizes( lat2d )

 dist  = gc_latlon(LAT,LON,lat1d,lon1d, 2,2)
 ind1d = dim_pqsort(dist,1)
 if (.not.ismissing(ind1d(0))) then
     i = 0
     ij = ind_resolve( ind1d(i), n2d )
     do while ismissing(xvar(ij(0),ij(1)))
       ij(:) = ind_resolve( ind1d(i), n2d )
       print("iteration: "+i+", x:"+ij(0)+" y:"+ij(1)+" lat:"+lat2d(ij(0),0)+" lon:"+lon2d(0,ij(1))+" val:"+xvar(ij(0),ij(1)))
       if i .gt. 300 then
         print("Over 300 iterations. Could not find any data nearby "+sname+","+sstate) 
         exit
       end if
       i=i+1
     end do
 else
     print("getind_latlon2d: lat="+ LAT(n)+"  lon="+ LON(n)+" problem")
 end if

 delete(dist)
 delete(ind1d)
 return( ij )
end


; ------------------------

; D Shea
; compute a user specified seasonal mean [all are three-month means]
; DJF,JFM,FMA,MAM,AMJ,MJJ,JJA,JAS,ASO,SON,OND,NDJ

; first (DJF=JF) /last (NDJ=ND) seasons are 2-month averages
;
; x(time,lat,lon),  x(time,lev,lat,lon)
;   ^^^^^^^^^^^^      ^^^^^^^^^^^^^^^^   
; must have named dim BUT can be ANY names
;
; The input "x" are assumed to contain monthly mean data
; The size of "time" MUST be divisible by 12.
; Also, it is assumed the "Jan" is the 1st month.
;
; xMon(time)  or  xMon(time,lat,lon)   or  xMon(time,lev,lat,lon)
; USAGE     xJJA = month2season (xMon, "JJA")
;
; RESULT    xJJA(time/12,lev,lat,lon)    xJJA(time/12,lat,lon)
;
; Note: this returns (7/2003) NMO as an attribute

undef ("month2season")
function month2season (xMon:numeric, SEASON:string)

local season,NMO,dimx,rank,ntim,nlat,mlon,nmos,nyrs,con \
    , nyrStrt,nyrLast,nyr,n,xSea, klev, dName,cv,xSea
begin
  season  =  (/"DJF","JFM","FMA","MAM","AMJ","MJJ" \
              ,"JJA","JAS","ASO","SON","OND","NDJ" /)

  NMO     = ind(season.eq.SEASON)  ; index corresponding to season
  if (ismissing(NMO)) then
      print ("contributed: month2season: bad season: SEASON="+SEASON)
      exit
  end if

  dimx    = dimsizes(xMon)
  rank    = dimsizes(dimx)
  if (rank.ne.2) then
      print ("contributed: month2season: rank="+rank)
      print ("----- rank currently not handled -----")
  end if

  nmos    = 12
  ntim    = dimx(0)
  modCheck ("month2season", ntim, nmos)

  if (rank.eq.2) then
      nsite    = dimx(rank-1)
  end if
  nyrs    = ntim/nmos
  con     = 1./3.

  nyrStrt = 0
  nyrLast = nyrs-1
  if (NMO.eq.0) then
      nyrStrt = 1
  end if
  if (NMO.eq.nmos-1) then
      nyrLast = nyrs-2
  end if

  if (rank.eq.2) then
      xSea = new ( (/nyrs,nsite/), typeof(xMon), getFillValue(xMon))
    do nyr=nyrStrt,nyrLast
         n = nyr*nmos + NMO
         xSea(nyr,:) = (xMon(n-1,:) + xMon(n,:) + xMon(n+1,:))*con
      end do
                                        ; special for beginning/end points
     if (NMO.eq.0) then
         n = 0
         xSea(0,:) = (xMon(n,:) + xMon(n+1,:))*0.5
     end if
     if (NMO.eq.nmos-1) then
         n = (nyrs-1)*nmos + NMO
         xSea(nyrs-1,:) = (xMon(n,:) + xMon(n-1,:))*0.5
     end if
  end if

  copy_VarAtts (xMon, xSea)
  if (isatt(xMon,"long_name") .or. isatt(xMon,"description") .or. \
      isatt(xMon,"standard_name") ) then
      xSea@long_name = SEASON+": "+getLongName(xMon)
  end if

  do n=1,rank-1                  ; copy spatial coordinates
     if (.not.ismissing(xMon!n)) then
         xSea!n = xMon!n
        if(iscoord(xMon,xMon!n))
           xSea&$xSea!n$ = xMon&$xMon!n$
        end if
     end if
  end  do

  dName        = xMon!0
  xSea!0       = dName

  if(iscoord(xMon,dName)) then
      cv = xMon&$dName$(NMO:ntim-1:nmos)
                                     ; possibly override
      xSea&$dName$ = cv
  end if

  xSea@NMO = NMO   ; for possible use in subscripting 
                   ; eg: nStrt= xSea@NMO ; time(nStrt:ntim-1,12)
  return (xSea)
  end


; ------------------------

; D Shea
; compute a user specified seasonal mean [all are three-month means]
; DJF,JFM,FMA,MAM,AMJ,MJJ,JJA,JAS,ASO,SON,OND,NDJ

; first (DJF=JF) /last (NDJ=ND) seasons are 2-month averages
;
; x(time,lat,lon),  x(time,lev,lat,lon)
;   ^^^^^^^^^^^^      ^^^^^^^^^^^^^^^^   
; must have named dim BUT can be ANY names
;
; The input "x" are assumed to contain monthly mean data
; The size of "time" MUST be divisible by 12.
; Also, it is assumed the "Jan" is the 1st month.
;
; xMon(time)  or  xMon(time,lat,lon)   or  xMon(time,lev,lat,lon)
; USAGE     xJJA = month2season_pr_total(xMon, "JJA")
;
; RESULT    xJJA(time/12,lev,lat,lon)    xJJA(time/12,lat,lon)
;
; Note: this returns (7/2003) NMO as an attribute

undef ("month2season_pr_total")
function month2season_pr_total (xMon:numeric, SEASON:string)

local season,NMO,dimx,rank,ntim,nlat,mlon,nmos,nyrs,con \
    , nyrStrt,nyrLast,nyr,n,xSea, klev, dName,cv,xSea
begin
  season  =  (/"DJF","JFM","FMA","MAM","AMJ","MJJ" \
              ,"JJA","JAS","ASO","SON","OND","NDJ" /)

  NMO     = ind(season.eq.SEASON)  ; index corresponding to season
  if (ismissing(NMO)) then
      print ("contributed: month2season_pr_total: bad season: SEASON="+SEASON)
      exit
  end if

  dimx    = dimsizes(xMon)
  rank    = dimsizes(dimx)
  if (rank.ne.2) then
      print ("contributed: month2season_pr_total: rank="+rank)
      print ("----- rank currently not handled -----")
  end if

  nmos    = 12
  ntim    = dimx(0)
  modCheck ("month2season_pr_total", ntim, nmos)

  if (rank.eq.2) then
      nsite    = dimx(rank-1)
  end if
  nyrs    = ntim/nmos

  nyrStrt = 0
  nyrLast = nyrs-1
  if (NMO.eq.0) then
      nyrStrt = 1
  end if
  if (NMO.eq.nmos-1) then
      nyrLast = nyrs-2
  end if

  if (rank.eq.2) then
      xSea = new ( (/nyrs,nsite/), typeof(xMon), getFillValue(xMon))
    do nyr=nyrStrt,nyrLast
         n = nyr*nmos + NMO
         xSea(nyr,:) = xMon(n-1,:) + xMon(n,:) + xMon(n+1,:)
      end do
                                        ; special for beginning/end points
     if (NMO.eq.0) then
         n = 0
         xSea(0,:) = xMon(n,:) + xMon(n+1,:)
     end if
     if (NMO.eq.nmos-1) then
         n = (nyrs-1)*nmos + NMO
         xSea(nyrs-1,:) = xMon(n,:) + xMon(n-1,:)
     end if
  end if

  copy_VarAtts (xMon, xSea)
  if (isatt(xMon,"long_name") .or. isatt(xMon,"description") .or. \
      isatt(xMon,"standard_name") ) then
      xSea@long_name = SEASON+": "+getLongName(xMon)
  end if

  do n=1,rank-1                  ; copy spatial coordinates
     if (.not.ismissing(xMon!n)) then
         xSea!n = xMon!n
        if(iscoord(xMon,xMon!n))
           xSea&$xSea!n$ = xMon&$xMon!n$
        end if
     end if
  end  do

  dName        = xMon!0
  xSea!0       = dName

  if(iscoord(xMon,dName)) then
      cv = xMon&$dName$(NMO:ntim-1:nmos)
                                     ; possibly override
      xSea&$dName$ = cv
  end if

  xSea@NMO = NMO   ; for possible use in subscripting 
                   ; eg: nStrt= xSea@NMO ; time(nStrt:ntim-1,12)
  return (xSea)
  end
