      subroutine count(ncnty,nday,xvar,thres,icnt,l)
!*************************************************************
!  Written by DJ Rasmussen
!  Rhodium Group, LLC
!  Last Updated: Sun Aug 25 21:25:09 EDT 2013

! writes an array of residuals to a text file
!*************************************************************

      implicit none
      integer ncnty,nday
      real icnt(ncnty)
      real xvar(ncnty,nday),thres
      integer ic, id
      logical l

      if (l) then 
       do ic=1,ncnty
        do id=1,nday
          if(xvar(ic,id) .gt. thres)then
           icnt(ic)=icnt(ic)+1
          end if
        end do
       end do
      else
       do ic=1,ncnty
        do id=1,nday
          if(xvar(ic,id) .lt. thres)then
           icnt(ic)=icnt(ic)+1
          end if
        end do
       end do
      end if
      return
      end subroutine count
