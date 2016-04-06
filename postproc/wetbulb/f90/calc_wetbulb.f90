      subroutine estimwet(nparm,params,nday,n305,nyr,slope, &
                          gtas,xwet,xtas,hmin)
!  Written by DJ Rasmussen
!  Rhodium Group, LLC
!  Last Updated: Sat Mar 22 18:00:24 PDT 2014

       implicit none

! ..... in only
       integer, intent(in)  :: nparm, nday, nyr
       real,    intent(in)  :: params(nparm)

!      params(1)  .... slope of simple linear regression
!      params(2)  .... y-intercept of simple linear regression
!      params(3)  .... stddev of simple linear regression
!      params(4)  .... BIC of simple linear regression
!      params(5)  .... y-intercept of the piecewise regression
!      params(6)  .... parameter b1 of the piecewise regression
!      params(7)  .... parameter b2 of the piecewise regression
!      params(8)  .... breakpoint
!      params(9)  .... stddev of piecewise regression
!      params(10) .... BIC of the piecewise regression
 
       real,    intent(in)  :: gtas(nyr)
       real,    intent(in)  :: slope, hmin

! ..... in/out
       integer, intent(inout)  :: n305 ! count number of days > 305K
       real,    intent(inout)  :: xwet(nday)
       real,    intent(inout)  :: xtas(nday)

! ..... local variables 
       integer :: i, iy 
       real :: tdlt

       n305=0

! ..... 2-component model

          if (params(10) .eq. -999.99) then ! default to simple linear regression

             do i = 1, nday
               xwet(i) = params(2) + params(1)*xtas(i)
               if (xwet(i) .ge. 305.) then
                  n305 = n305 + 1
               end if
             end do

! ...... Use what ever model has smallest BIC
          else 
             
            if (params(10) .lt. params(4)) then

! ..... 4-component model
                iy=1
                do i=1,nday
   
                  tdlt = xtas(i) - slope*gtas(iy) ! forced local temp. change
                  xwet(i) = params(5) + params(6)*min(tdlt, params(8)) &
                            + params(7)*max(0., tdlt-params(8))  &
                            + params(1)*(slope*gtas(iy))
   
                  ! Where there is a negative slope, let's prevent wet bulb from 
                  ! being lower than the historical minimum
   
                  if (xwet(i) .lt. hmin) then 
                    xwet(i) = hmin
                  end if
   
                  if (xwet(i) .ge. 305.) then
                     n305 = n305 + 1
                  end if
   
                  if (mod(i, 92) .eq. 0) then ! 92 days in June + July + August
                    iy = iy + 1 ! use next record of global mean temperature
                  end if
   
                end do


            else 
! ..... 2-component model

                do i = 1, nday
                  xwet(i) = params(2) + params(1)*xtas(i)
                  if (xwet(i) .ge. 305.) then
                     n305 = n305 + 1
                  end if
                end do
    
            end if

          end if

       return
      end subroutine estimwet
