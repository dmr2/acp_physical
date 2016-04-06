       subroutine ds(nsite,nday,obsday,f,modmm,xwork,icnt,iscnt,precip)
!  Written by DJ Rasmussen
!  Rhodium Group, LLC
!  Last Updated: Thu 03 Mar 2016 05:02:29 PM PST

! Does the precipitation downscaling including the spillover calculation

      implicit none
      logical, intent(in)  :: precip ! downscale precipitation (True, False)
      integer, intent(in)  :: nsite ! number of GHCND sites or counties
      integer, intent(in)  :: nday ! number of days in the month
      real,    intent(in)  :: obsday(nday,nsite) ! observed daily values
      real,    intent(in)  :: modmm(nsite) ! model monthly mean
      real,    intent(in)  :: f(nsite) ! ratio of model mon. mean and obs. mon. mean 
      real,    intent(inout)  :: xwork(nday,nsite) ! daily rainfall totals within the month
      integer, intent(out)  :: icnt ! total number of times spillover calc invoked
      integer, intent(out)  :: iscnt ! total number of sites where spillover calc invoked

! .... local variables
      integer :: isite,i,icheck
      real :: tenmax
      real ::twomax

      if ( precip ) then
          icnt=0
          do isite=1,nsite
           icheck=0
           ! threshold criteria
           tenmax=10*modmm(isite) ! 10x model monthly mean
           twomax=2*maxval(obsday(:,isite)) ! 2x max observed daily total for this month

           do i=1,nday
             xwork(i,isite)=obsday(i,isite)*f(isite) ! scale obs. dailies to match model monthly mean
           end do 

           do i=1,nday ! check if precip violates thresholds
             if ( xwork(i,isite) .gt. tenmax .and. &
                       xwork(i,isite) .gt. twomax ) then
               icheck=1
               icnt=icnt+1
               ! distribute precip amongst adjacent days
               if (i .eq. 1) then
                 xwork(i+1,isite)=xwork(i+1,isite)+xwork(i,isite)/3.
                 xwork(i+2,isite)=xwork(i+2,isite)+xwork(i,isite)/3.
               else if (i .eq. nday) then
                 xwork(i-1,isite)=xwork(i-1,isite)+xwork(i,isite)/3.
                 xwork(i-2,isite)=xwork(i-2,isite)+xwork(i,isite)/3.
               else 
                 xwork(i-1,isite)=xwork(i-1,isite)+xwork(i,isite)/3.
                 xwork(i+1,isite)=xwork(i+1,isite)+xwork(i,isite)/3.
               end if
               xwork(i,isite)=xwork(i,isite)/3.
             end if
           end do ! each day
           
           if ( icheck .eq. 1) then
             iscnt=iscnt+1
           end if
          end do ! each site
      else  ! for temperature
          do isite=1,nsite
           do i=1,nday
             xwork(i,isite) = modmm(isite) + (obsday(i,isite) - f(isite))
           end do ! each day in month
          end do ! each site
      end if


      return
      end subroutine ds
