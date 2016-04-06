
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

undef("closest_val_AnyOrder")
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
     i=0
     ij = ind_resolve( ind1d(i), n2d )
     do while ismissing(xvar(ij(0),ij(1)))
       ij(:) = ind_resolve( ind1d(i), n2d )
       print("iteration: "+i+", "+ij(0)+" "+ij(1)+" "+lat2d(ij(0),0)+" "+lon2d(0,ij(1))+" "+xvar(ij(0),ij(1)))
       if i .gt. 300 then
         print("Could not find any data nearby "+sname+","+sstate) 
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
