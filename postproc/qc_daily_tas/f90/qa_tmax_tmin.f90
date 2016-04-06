      subroutine qa_minmax(nsite,nday,icnt,xtmin,xtmax,xtas)
!  Written by DJ Rasmussen
!  Rhodium Group, LLC
!  Last Updated: Sat Mar 15 10:43:58 PDT 2014

       implicit none
       integer, intent(in)  :: nsite, nday
       real,    intent(inout)  :: xtmin(nday,nsite)
       real,    intent(inout)  :: xtmax(nday,nsite)
       real,    intent(in)  :: xtas(nday,nsite)
       integer,    intent(inout)  :: icnt

! ..... local variables 
       integer :: isite,i

       do isite=1,nsite
        do i=1,nday
          if (xtmin(i,isite) > xtas(i,isite) .and. &
                    xtmax(i,isite) > xtas(i,isite)) then

            xtmin(i,isite) = xtas(i,isite) - 2.5
            xtmax(i,isite) = xtas(i,isite) + 2.5
            icnt=icnt+1

          else if (xtmin(i,isite) > xtmax(i,isite)) then

            xtmin(i,isite) = xtas(i,isite) - 2.5
            xtmax(i,isite) = xtas(i,isite) + 2.5
            icnt=icnt+1

          else if (xtmin(i,isite) > xtas(i,isite)) then

            xtmin(i,isite) = xtas(i,isite) - 2.5
            xtmax(i,isite) = xtas(i,isite) + 2.5
            icnt=icnt+1
       
          else if (xtas(i,isite) > xtmax(i,isite)) then

            xtmin(i,isite) = xtas(i,isite) - 2.5
            xtmax(i,isite) = xtas(i,isite) + 2.5
            icnt=icnt+1

          end if
        end do ! each day
       end do ! each site

       return
      end subroutine qa_minmax
