!
! !MODULE: etrac_output.f90  --- NetCDF output
!
! !INTERFACE:
   MODULE mod_etrac_output
#ifndef TBNTonly_bulk_bud
!
! !DESCRIPTION:
!
! !USES:
   
   use netcdf
   
   use mod_etrac_common
   use mod_etrac_init, only: get_nc_field
   
   implicit none
!
!  default: all is private.
   private

   public init_tbnt_output
   public do_tbnt_output
!
! !LOCAL VARIABLES:
   ! relative fraction output
   integer                :: relFracRec, relFracTimeID
   ! absolute fraction output
   integer                :: absFracRec
   integer, dimension(:), allocatable :: absFracTimeID
   ! linked fluxes fraction output
   integer                :: lnkFracRec, lnkFracTimeID
!
!=======================================================================

   contains

!=======================================================================
!
! !INTERFACE:
   subroutine init_tbnt_output(ierr)
#define SUBROUTINE_NAME 'init_tbnt_output'
!
! !DESCRIPTION:
!  initialise output files
!
! !USES:
   implicit none
!
! !OUTPUT PARAMETERS:
   integer, intent(inout) :: ierr
!
! !LOCAL VARIABLES:
!
   integer                :: iTarget, ios, system
   logical                :: lExist, writeRel, writeLnk
   character(len=nameLen) :: csvFile
!
!-----------------------------------------------------------------------
#include "call-trace.inc"
   
   ! create output directory
   ios = 0
   lExist = .false.
   inquire(file=trim(output_dir),exist=lExist)
   if (.not.lExist) ios = system('mkdir -p '//trim(output_dir))
   if (ios/=0) then
      ierr = ierr + ios
      write(error_msg,'(a)')'Error creating output directory: '//trim(output_dir)
      call stop_tbnt(error_msg,ierr)
   end if
   
   ! initialize write switches (latet set to ".true." if new files are created)
   writeRel = .false.
   writeLnk = .false.
   
   ! create files and allocate fields for relative fractions
   if (TBNTrelFracOutStep>=0) then
      ! init nc file
      allocate ( TBNTrelVarID3D(TBNTnVars3D,TBNTnFractions) )
#ifndef TBNTnoVar2D
      allocate ( TBNTrelVarID2D(TBNTnVars2D,TBNTnFractions) )
#endif
      relFracFile = trim(output_dir)//trim(TBNTrun)//'_relative_fractions.nc'
      lExist = .false.
      inquire(file=trim(relFracFile),exist=lExist)
      if (.not.lExist.or..not.TBNTcontinue) then
         relFracRec = 0
         call nc_create(relFracFile, 0, 1, relFracTimeID, relFracFileID, ierr)
         writeRel = .true.
      else
         ! open existing file
         relFracRec = TBNTiStart-TBNTiOffset
         call nc_open(relFracFile, 0, 1, relFracTimeID, relFracFileID, ierr)
      end if
   end if

   ! allocate fields for absolute fractions
   if (TBNTabsFracOutStep>=0) then
      allocate ( TBNTabsFlxID2D(TBNTnFluxes2D,TBNTnFractions), TBNTabsFlxID3D(TBNTnFluxes3D,TBNTnFractions) )
      allocate ( TBNTabsVarID3D(TBNTnVars3D,TBNTnFractions) )
#ifndef TBNTnoVar2D
      allocate ( TBNTabsVarID2D(TBNTnVars2D,TBNTnFractions) )
#endif
      allocate ( absFracFile(TBNTnFractions), absFracFileID(TBNTnFractions), absFracTimeID(TBNTnFractions) )
   end if

   ! linked fraction fluxes
   if (TBNTlinkedFluxes) then
      ! init nc file and write initial state
      allocate ( TBNTlnkFlxID2D(TBNTnLinkedFluxes2D,TBNTnFractions), TBNTlnkFlxID3D(TBNTnLinkedFluxes3D,TBNTnFractions) )
      lnkFracFile = trim(output_dir)//trim(TBNTrun)//'_linked_fluxes.nc'
      lExist = .false.
      inquire(file=trim(lnkFracFile),exist=lExist)
      if (.not.lExist.or..not.TBNTcontinue) then
         ! create new file
         lnkFracRec = 0
         call nc_create(lnkFracFile, 0, 3, lnkFracTimeID, lnkFracFileID, ierr)
         writeLnk = .true.
      else
         ! open existing file
         lnkFracRec = TBNTiStart-TBNTiOffset
         call nc_open(lnkFracFile, 0, 3, lnkFracTimeID, lnkFracFileID, ierr)
      end if
   end if

   do iFrac = 1,TBNTnFractions
      ! write initial state for relative fractions (in case of new file)
      if (TBNTrelFracOutStep>=0.and.writeRel) call nc_write(relFracFile, iFrac, relFracFileID, relFractimeID, outTimeStart, 1, ierr)
      ! init nc files and write initial state for absolute fractions (in case of new file)
      if (TBNTabsFracOutStep>=0) then
         absFracFile(iFrac) = trim(output_dir)//trim(TBNTrun)//'_absolute_fractions_'//trim(TBNTfracName(iFrac))//'.nc'
         lExist = .false.
         inquire(file=trim(absFracFile(iFrac)),exist=lExist)
         if (.not.lExist.or..not.TBNTcontinue) then
            ! create new file and write initial state
            if (iFrac==1) absFracRec = 0
            call nc_create(absFracFile(iFrac), iFrac, 2, absFracTimeID(iFrac), absFracFileID(iFrac), ierr)
            call nc_write(absFracFile(iFrac), iFrac, absFracFileID(iFrac), absFractimeID(iFrac), outTimeStart, 2, ierr)
         else
            ! open existing file
            if (iFrac==1) absFracRec = TBNTiStart-TBNTiOffset
            call nc_open(absFracFile(iFrac), iFrac, 2, absFracTimeID(iFrac), absFracFileID(iFrac), ierr)
         end if
      end if
      ! write initial state for linked fraction fluxes (in case of new file)
      if (TBNTlinkedFluxes.and.writeLnk) call nc_write(lnkFracFile, iFrac, lnkFracFileID, lnkFractimeID, outTimeStart, 3, ierr)
   end do
   
   ! initialize target output
   if (TBNTtargetOut) then
      do iTarget = 1,TBNTnTargetVars
         csvFile = trim(output_dir)//trim(TBNTrun)//'_'//trim(TBNTtargetVars(iTarget)%code)//'_target_results.csv'
         lExist = .false.
         inquire(file=trim(csvFile),exist=lExist)
         if (.not.lExist.or..not.TBNTcontinue) then
            ! create new file
            open(target_unit+iTarget, file=trim(csvFile), status = 'replace')
            call asc_write_header(target_unit+iTarget,iTarget,ierr)
         else
            ! continue writing in existing file
            if (.not.lExist) then
               ierr = 1
               write(error_msg,'(a)')'Error opening existing file: '//trim(csvFile)
               call stop_tbnt(error_msg,ierr)
            end if
            open(target_unit+iTarget, file=trim(csvFile), status = 'old', position = 'append')
         end if
      end do
   end if

   end subroutine init_tbnt_output
!=======================================================================
!
! !INTERFACE:
   subroutine do_tbnt_output(ierr)
#define SUBROUTINE_NAME 'do_tbnt_output'
!
! !DESCRIPTION:
!  caller routine for all output routines
!
! !USES:
   implicit none
!
! !OUTPUT PARAMETERS:
   integer, intent(inout) :: ierr
!
! !LOCAL VARIABLES:
!
   integer  :: iTarget, ncStatus, ios
   real(dp) :: nc_time
!
!-----------------------------------------------------------------------
#include "call-trace.inc"

   ios = 0
   
   ! set time
   nc_time = outTimeStart + real(iStep,dp) * outTimeFac * real(TBNTstep,dp)/rMinPerDay
   
   ! write relative fractions
   if ((TBNTrelFracOutStep>0.and.mod(iStep,TBNTrelFracOutStep)==0).or. &
       (TBNTrelFracOutStep>=0.and.iStep==TBNTnSteps)) then
      do iFrac = 1,TBNTnFractions
         call nc_write(relFracFile, iFrac, relFracFileID, relFracTimeID, nc_time, 1, ierr)
      end do
      ncStatus = NF90_SYNC(relFracFileID)
      call nc_check(ncStatus, ios)
      if (ios/=0) then
         ierr = ierr + ios
         write(error_msg,'(a)')'NetCDF-error synchronizing data to disk: '//trim(relFracFile)//' - '//trim(NF90_STRERROR(ncStatus))
         call stop_tbnt(error_msg,ierr)
      end if
   end if
   ! write absolute fractions
   if ((TBNTabsFracOutStep>0.and.mod(iStep,TBNTabsFracOutStep)==0).or. &
       (TBNTabsFracOutStep>=0.and.iStep==TBNTnSteps)) then
      do iFrac = 1,TBNTnFractions
         call nc_write(absFracFile(iFrac), iFrac, absFracFileID(iFrac), absFracTimeID(iFrac), nc_time, 2, ierr)
         ncStatus = NF90_SYNC(absFracFileID(iFrac))
         call nc_check(ncStatus, ios)
         if (ios/=0) then
            ierr = ierr + ios
            write(error_msg,'(a)')'NetCDF-error synchronizing data to disk: '//trim(absFracFile(iFrac))//' - '//trim(NF90_STRERROR(ncStatus))
            call stop_tbnt(error_msg,ierr)
         end if
      end do
      TBNTflux3D = 0.e0
      TBNTflux2D = 0.e0
   end if
   ! write linked fraction fluxes
   if ((TBNTlnkFracOutStep>0.and.mod(iStep,TBNTlnkFracOutStep)==0).or. &
       (TBNTlnkFracOutStep>=0.and.iStep==TBNTnSteps)) then
      do iFrac = 1,TBNTnFractions
         call nc_write(lnkFracFile, iFrac, lnkFracFileID, lnkFracTimeID, nc_time, 3, ierr)
      end do
      ncStatus = NF90_SYNC(lnkFracFileID)
      call nc_check(ncStatus, ios)
      if (ios/=0) then
         ierr = ierr + ios
         write(error_msg,'(a)')'NetCDF-error synchronizing data to disk: '//trim(lnkFracFile)//' - '//trim(NF90_STRERROR(ncStatus))
         call stop_tbnt(error_msg,ierr)
      end if
      TBNTlinkedFlux3D = 0.e0
      TBNTlinkedFlux2D = 0.e0
   end if
   
   ! write daily target output
   if (TBNTtargetOut.and.mod(iStep*TBNTstep,1440)==0) then
      do iTarget = 1,TBNTnTargetVars
         call asc_write_data(target_unit+iTarget,iTarget,ierr)
      end do
   end if

   end subroutine do_tbnt_output
!=======================================================================
! !INTERFACE:
   subroutine nc_create(filename, iFrac_, iSwitch, timeID, fileID, ierr)
#define SUBROUTINE_NAME 'nc_create'
!
! !DESCRIPTION:
!  initialise a new NetCDF file
!
   implicit none
!
! !INPUT/OUTPUT PARAMETERS:
   character(len=fileLen), intent(in   ) :: filename     ! file name
   integer               , intent(in   ) :: iFrac_       ! current fraction
   integer               , intent(in   ) :: iSwitch      ! relative fractions (1), absolute fractions (2) or linked fluxes (3)
   integer               , intent(  out) :: timeID       ! time ID for current file
   integer               , intent(  out) :: fileID       ! file ID for current file
   integer               , intent(inout) :: ierr         ! error tracker
!
! LOCAL VARIABLES
   integer                            :: ios, ncStatus, j, k, i, iFrac1, iFrac2, iiFrac
   integer, dimension(:), allocatable :: dimid, dummyID
   character(len=nameLen)             :: varStr, varUnit, fluxStr, fluxUnit, varName, flxName
   logical                            :: found
!   
!-----------------------------------------------------------------------
#include "call-trace.inc"
!
   ios = 0
   
   ! initialise new file
   ncStatus = NF90_CREATE(trim(filename), ior(ior(NF90_CLOBBER, NF90_64BIT_OFFSET), NF90_SHARE), fileID)
   call nc_check(ncStatus, ios)
   if (ios/=0) then
      ierr = ierr + ios
      write(error_msg,'(a)')'NetCDF-error opening new file: '//trim(fileName)//' - '//trim(NF90_STRERROR(ncStatus))
      call stop_tbnt(error_msg,ierr)
   end if
   
   NCdims(4)%size = NF90_UNLIMITED
   NCgrid(4)%size = NF90_UNLIMITED
   
   ! define dimensions
   do i = 1,4
      ncStatus = NF90_DEF_DIM(fileID, NCdims(i)%name, NCdims(i)%size, NCdims(i)%ncid)
      call nc_check(ncStatus, ios)
      if (ios/=0) then
         ierr = ierr + ios
         write(error_msg,'(a)')'NetCDF-error defining dimension: '//trim(fileName)//' - '//trim(NCdims(i)%name) &
                               //' - '//trim(NF90_STRERROR(ncStatus))
         call stop_tbnt(error_msg,ierr)
      end if
   end do
   ! define grid variables
   do i = 1,4
      if (i<4) then
         found = .false.
         do j = 1,size(NCgrid(i)%size)
            do k = 1,3
               if (NCgrid(i)%size(j) == NCdims(k)%size) then
                  NCgrid(i)%dimid(j) = NCdims(k)%ncid
                  found = .true.
                  exit
               end if
            end do
            if (.not.found) then
               ierr = 1
               write(error_msg,'(a)')'NetCDF-error - output grid variable''s dimension not found: '//trim(NCgrid(i)%name)
               call stop_tbnt(error_msg,ierr)
            end if
         end do
      else
         NCgrid(i)%dimid(1) = NCdims(4)%ncid
      end if
      ncStatus = NF90_DEF_VAR(fileID, NCgrid(i)%name, NF90_DOUBLE, NCgrid(i)%dimid, NCgrid(i)%ncid)
      call nc_check(ncStatus, ios)
      if (ios/=0) then
         ierr = ierr + ios
         write(error_msg,*)'NetCDF-error defining grid variable: '//trim(fileName)//' - '//trim(NCgrid(i)%name) &
                               //' - '//trim(NF90_STRERROR(ncStatus)),NCgrid(i)%dimid
         call stop_tbnt(error_msg,ierr)
      end if
      ncStatus = NF90_PUT_ATT(fileID, NCgrid(i)%ncid, "units", NCgrid(i)%unit)
      call nc_check(ncStatus, ios)
      if (ios/=0) then
         ierr = ierr + ios
         write(error_msg,'(a)')'NetCDF-error adding grid variable unit: '//trim(fileName)//' - '//trim(NCgrid(i)%name) &
                               //' - '//trim(NF90_STRERROR(ncStatus))
         call stop_tbnt(error_msg,ierr)
      end if
   end do
   ! write time variable ID to output variable
   timeID = NCgrid(4)%ncid
   
   ! define units and name add-on
   select case (iSwitch)
      case (1)
         varStr   = 'relative fraction of'
         varUnit  = '%'
         fluxStr  = '' ! no flux output
         fluxUnit = ''
      case (2)
         varStr   = 'absolute fraction of'
         varUnit  = 'mmol'
         fluxStr  = 'fraction flux of'
         fluxUnit = 'mmol'
     case (3)
         varStr   = '' ! no variable output
         varUnit  = ''
         fluxStr  = 'fraction flux of'
         fluxUnit = 'mmol'
      case default
         ierr = ierr + 1
         write(error_msg,'(a)')'NetCDF-error in nc_create: Invalid ''iSwitch'' value.'
         call stop_tbnt(error_msg,ierr)
   end select
   
   if (iFrac_==0) then
      iFrac1 = 1
      iFrac2 = TBNTnFractions
   else
      iFrac1 = iFrac_
      iFrac2 = iFrac_
   end if
   
   if (iSwitch<=2) then  
      do iiFrac = iFrac1,iFrac2
#ifndef TBNTnoVar2D
         ! 2D variables
         allocate ( dimid(3), dummyID(TBNTnVars2D) )
         dimid = (/ NCdims(1)%ncID, NCdims(2)%ncid, NCdims(4)%ncid /)
         do iVar = 1,TBNTnVars2D
            varName = trim(TBNTvar(TBNTnVars3D+iVar)%name)//'-'//trim(TBNTfracName(iiFrac))
            ncStatus = NF90_DEF_VAR(fileID, trim(varName), NF90_DOUBLE, dimid, dummyID(iVar))
            call nc_check(ncStatus, ios)
            if (ios/=0) then
               ierr = ierr + ios
               write(error_msg,'(a)')'NetCDF-error defining variable: '//trim(fileName)//' - '//trim(varName)// &
                                     ' - '//trim(NF90_STRERROR(ncStatus))
               call stop_tbnt(error_msg,ierr)
            end if
            ! add attributes
            ncStatus = NF90_PUT_ATT(fileID, dummyID(iVar), 'long_name', trim(varStr)//' '//trim(varName))
            call nc_check(ncStatus, ios)
            if (ios/=0) then
               ierr = ierr + ios
               write(error_msg,'(a)')'NetCDF-error adding variable "long_name": '//trim(fileName)//' - '// &
                                     trim(varName)//' - '//trim(NF90_STRERROR(ncStatus))
               call stop_tbnt(error_msg,ierr)
            end if
            ncStatus = NF90_PUT_ATT(fileID, dummyID(iVar), 'units', trim(varUnit))
            call nc_check(ncStatus, ios)
            if (ios/=0) then
               ierr = ierr + ios
               write(error_msg,'(a)')'NetCDF-error adding variable "unit": '//trim(fileName)//' - '// &
                                     trim(varName)//' - '//trim(NF90_STRERROR(ncStatus))
               call stop_tbnt(error_msg,ierr)
            end if
            ncStatus = NF90_PUT_ATT(fileID, dummyID(iVar), '_FillValue', fail)
            call nc_check(ncStatus, ios)
            if (ios/=0) then
               ierr = ierr + ios
               write(error_msg,'(a)')'NetCDF-error adding variable "FillValue": '//trim(fileName)//' - '// &
                                     trim(varName)//' - '//trim(NF90_STRERROR(ncStatus))
               call stop_tbnt(error_msg,ierr)
            end if
         end do
         if (iSwitch==1) then
            TBNTrelVarID2D(:,iiFrac) = dummyID
         else
            TBNTabsVarID2D(:,iiFrac) = dummyID
         end if
         deallocate ( dummyID, dimid )
#endif
         ! 3D variables
         allocate ( dimid(4), dummyID(TBNTnVars3D) )
         dimid = (/ NCdims(:)%ncid /)
         ! variables
         do iVar = 1,TBNTnVars3D
            varName = trim(TBNTvar(iVar)%name)//'-'//trim(TBNTfracName(iiFrac))
            ncStatus = NF90_DEF_VAR(fileID, trim(varName), NF90_DOUBLE, dimid, dummyID(iVar))
            call nc_check(ncStatus, ios)
            if (ios/=0) then
               ierr = ierr + ios
               write(error_msg,'(a)')'NetCDF-error defining variable: '//trim(fileName)//' - '//trim(varName)// &
                                     ' - '//trim(NF90_STRERROR(ncStatus))
               call stop_tbnt(error_msg,ierr)
            end if
            ! add attributes
            ncStatus = NF90_PUT_ATT(fileID, dummyID(iVar), 'long_name', trim(varStr)//' '//trim(varName))
            call nc_check(ncStatus, ios)
            if (ios/=0) then
               ierr = ierr + ios
               write(error_msg,'(a)')'NetCDF-error adding variable "long_name": '//trim(fileName)//' - '// &
                                  trim(varName)//' - '//trim(NF90_STRERROR(ncStatus))
               call stop_tbnt(error_msg,ierr)
            end if
            ncStatus = NF90_PUT_ATT(fileID, dummyID(iVar), 'units', trim(varUnit))
            call nc_check(ncStatus, ios)
            if (ios/=0) then
               ierr = ierr + ios
               write(error_msg,'(a)')'NetCDF-error adding variable "unit": '//trim(fileName)//' - '// &
                                     trim(varName)//' - '//trim(NF90_STRERROR(ncStatus))
               call stop_tbnt(error_msg,ierr)
            end if
            ncStatus = NF90_PUT_ATT(fileID, dummyID(iVar), '_FillValue', fail)
            call nc_check(ncStatus, ios)
            if (ios/=0) then
               ierr = ierr + ios
               write(error_msg,'(a)')'NetCDF-error adding variable "FillValue": '//trim(fileName)//' - '// &
                                     trim(varName)//' - '//trim(NF90_STRERROR(ncStatus))
               call stop_tbnt(error_msg,ierr)
            end if
         end do
         if (iSwitch==1) then
            TBNTrelVarID3D(:,iiFrac) = dummyID
         else
            TBNTabsVarID3D(:,iiFrac) = dummyID
         end if
         deallocate ( dummyID, dimid )
      end do
   end if
   
   if (iSwitch==2) then ! create fraction fluxes only in case of linked fluxes or absolute fraction output
      ! 2D fluxes
      allocate ( dimid(3) )
      dimid = (/ NCdims(1)%ncid, NCdims(2)%ncid, NCdims(4)%ncid /)
      iOff = TBNTnFluxes3D+TBNTnLinkedFluxes3D
      do iFlux = 1,TBNTnFluxes2D
         flxName = trim(TBNTflux(iOff+iFlux)%name)//'-'//trim(TBNTfracName(iFrac_))
         ncStatus = NF90_DEF_VAR(fileID, trim(flxName), NF90_DOUBLE, dimid, TBNTabsFlxID2D(iFlux,iFrac_))
         call nc_check(ncStatus, ios)
         if (ios/=0) then
            ierr = ierr + ios
            write(error_msg,'(a)')'NetCDF-error defining variable: '//trim(fileName)//' - '//trim(flxName)// &
                                  ' - '//trim(NF90_STRERROR(ncStatus))
            call stop_tbnt(error_msg,ierr)
         end if
         ! add attributes
         ncStatus = NF90_PUT_ATT(fileID, TBNTabsFlxID2D(iFlux,iFrac_), 'long_name', trim(fluxStr)//' '//trim(flxName))
         call nc_check(ncStatus, ios)
         if (ios/=0) then
            ierr = ierr + ios
            write(error_msg,'(a)')'NetCDF-error adding variable "long_name": '//trim(fileName)//' - '// &
                                  trim(flxName)//' - '//trim(NF90_STRERROR(ncStatus))
            call stop_tbnt(error_msg,ierr)
         end if
         ncStatus = NF90_PUT_ATT(fileID, TBNTabsFlxID2D(iFlux,iFrac_), 'units', trim(fluxUnit))
         call nc_check(ncStatus, ios)
         if (ios/=0) then
            ierr = ierr + ios
            write(error_msg,'(a)')'NetCDF-error adding variable "unit": '//trim(fileName)//' - '// &
                                  trim(flxName)//' - '//trim(NF90_STRERROR(ncStatus))
            call stop_tbnt(error_msg,ierr)
         end if
         ncStatus = NF90_PUT_ATT(fileID, TBNTabsFlxID2D(iFlux,iFrac_), '_FillValue', fail)
         call nc_check(ncStatus, ios)
         if (ios/=0) then
            ierr = ierr + ios
            write(error_msg,'(a)')'NetCDF-error adding variable "FillValue": '//trim(fileName)//' - '// &
                                  trim(flxName)//' - '//trim(NF90_STRERROR(ncStatus))
            call stop_tbnt(error_msg,ierr)
         end if
      end do
      deallocate ( dimid )
      ! 3D fluxes
      allocate ( dimid(4) )
      dimid = (/ NCdims(:)%ncid /)
      iOff = 0
      do iFlux = 1,TBNTnFluxes3D
         flxName = trim(TBNTflux(iOff+iFlux)%name)//'-'//trim(TBNTfracName(iFrac_))
         ncStatus = NF90_DEF_VAR(fileID, trim(flxName), NF90_DOUBLE, dimid, TBNTabsFlxID3D(iFlux,iFrac_))
         call nc_check(ncStatus, ios)
         if (ios/=0) then
            ierr = ierr + ios
            write(error_msg,'(a)')'NetCDF-error defining flux: '//trim(fileName)//' - '//trim(flxName)// &
                                  ' - '//trim(NF90_STRERROR(ncStatus))
            call stop_tbnt(error_msg,ierr)
         end if
         ! add attributes
         ncStatus = NF90_PUT_ATT(fileID, TBNTabsFlxID3D(iFlux,iFrac_), 'long_name', trim(fluxStr)//' '//trim(flxName))
         call nc_check(ncStatus, ios)
         if (ios/=0) then
            ierr = ierr + ios
            write(error_msg,'(a)')'NetCDF-error adding flux "long_name": '//trim(fileName)//' - '// &
                                  trim(flxName)//' - '//trim(NF90_STRERROR(ncStatus))
            call stop_tbnt(error_msg,ierr)
         end if
         ncStatus = NF90_PUT_ATT(fileID, TBNTabsFlxID3D(iFlux,iFrac_), 'units', trim(fluxUnit))
         call nc_check(ncStatus, ios)
         if (ios/=0) then
            ierr = ierr + ios
            write(error_msg,'(a)')'NetCDF-error adding flux "unit": '//trim(fileName)//' - '// &
                                  trim(flxName)//' - '//trim(NF90_STRERROR(ncStatus))
            call stop_tbnt(error_msg,ierr)
         end if
         ncStatus = NF90_PUT_ATT(fileID, TBNTabsFlxID3D(iFlux,iFrac_), '_FillValue', fail)
         call nc_check(ncStatus, ios)
         if (ios/=0) then
            ierr = ierr + ios
            write(error_msg,'(a)')'NetCDF-error adding flux "FillValue": '//trim(fileName)//' - '// &
                                  trim(flxName)//' - '//trim(NF90_STRERROR(ncStatus))
            call stop_tbnt(error_msg,ierr)
         end if
      end do
      deallocate ( dimid )
   else if (iSwitch==3) then
      ! 2D fluxes
      allocate ( dimid(3) )
      dimid = (/ NCdims(1)%ncid, NCdims(2)%ncid, NCdims(4)%ncid /)
      iOff = TBNTnFluxes3D+TBNTnLinkedFluxes3D+TBNTnFluxes2D
      allocate ( dummyID(TBNTnLinkedFluxes2D) )
      do iiFrac = iFrac1,iFrac2
         do iFlux = 1,TBNTnLinkedFluxes2D
            flxName = trim(TBNTflux(iOff+iFlux)%name)//'-'//trim(TBNTfracName(iiFrac))
            ncStatus = NF90_DEF_VAR(fileID, trim(flxName), NF90_DOUBLE, dimid, dummyID(iFlux))
            call nc_check(ncStatus, ios)
            if (ios/=0) then
               ierr = ierr + ios
               write(error_msg,'(a)')'NetCDF-error defining variable: '//trim(fileName)//' - '//trim(flxName)// &
                                     ' - '//trim(NF90_STRERROR(ncStatus))
               call stop_tbnt(error_msg,ierr)
            end if
            ! add attributes
            ncStatus = NF90_PUT_ATT(fileID, dummyID(iFlux), 'long_name', trim(fluxStr)//' '//trim(flxName))
            call nc_check(ncStatus, ios)
            if (ios/=0) then
               ierr = ierr + ios
               write(error_msg,'(a)')'NetCDF-error adding variable "long_name": '//trim(fileName)//' - '// &
                                     trim(flxName)//' - '//trim(NF90_STRERROR(ncStatus))
               call stop_tbnt(error_msg,ierr)
            end if
            ncStatus = NF90_PUT_ATT(fileID, dummyID(iFlux), 'units', trim(fluxUnit))
            call nc_check(ncStatus, ios)
            if (ios/=0) then
               ierr = ierr + ios
               write(error_msg,'(a)')'NetCDF-error adding variable "unit": '//trim(fileName)//' - '// &
                                     trim(flxName)//' - '//trim(NF90_STRERROR(ncStatus))
               call stop_tbnt(error_msg,ierr)
            end if
            ncStatus = NF90_PUT_ATT(fileID, dummyID(iFlux), '_FillValue', fail)
            call nc_check(ncStatus, ios)
            if (ios/=0) then
               ierr = ierr + ios
               write(error_msg,'(a)')'NetCDF-error adding variable "FillValue": '//trim(fileName)//' - '// &
                                     trim(flxName)//' - '//trim(NF90_STRERROR(ncStatus))
               call stop_tbnt(error_msg,ierr)
            end if
         end do
         TBNTlnkFlxID2D(:,iiFrac) = dummyID
      end do  
      deallocate ( dimid, dummyID )
      ! 3D fluxes
      allocate ( dimid(4) )
      dimid = (/ NCdims(:)%ncid /)
      iOff = TBNTnFluxes3D
      allocate ( dummyID(TBNTnLinkedFluxes3D) )
      do iiFrac = iFrac1,iFrac2
         do iFlux = 1,TBNTnLinkedFluxes3D
            flxName = trim(TBNTflux(iOff+iFlux)%name)//'-'//trim(TBNTfracName(iiFrac))
            ncStatus = NF90_DEF_VAR(fileID, trim(flxName), NF90_DOUBLE, dimid, dummyID(iFlux))
            call nc_check(ncStatus, ios)
            if (ios/=0) then
               ierr = ierr + ios
               write(error_msg,'(a)')'NetCDF-error defining flux: '//trim(fileName)//' - '//trim(flxName)// &
                                     ' - '//trim(NF90_STRERROR(ncStatus))
               call stop_tbnt(error_msg,ierr)
            end if
            ! add attributes
            ncStatus = NF90_PUT_ATT(fileID, dummyID(iFlux), 'long_name', trim(fluxStr)//' '//trim(flxName))
            call nc_check(ncStatus, ios)
            if (ios/=0) then
               ierr = ierr + ios
               write(error_msg,'(a)')'NetCDF-error adding flux "long_name": '//trim(fileName)//' - '// &
                                     trim(flxName)//' - '//trim(NF90_STRERROR(ncStatus))
               call stop_tbnt(error_msg,ierr)
            end if
            ncStatus = NF90_PUT_ATT(fileID, dummyID(iFlux), 'units', trim(fluxUnit))
            call nc_check(ncStatus, ios)
            if (ios/=0) then
               ierr = ierr + ios
               write(error_msg,'(a)')'NetCDF-error adding flux "unit": '//trim(fileName)//' - '// &
                                     trim(flxName)//' - '//trim(NF90_STRERROR(ncStatus))
               call stop_tbnt(error_msg,ierr)
            end if
            ncStatus = NF90_PUT_ATT(fileID, dummyID(iFlux), '_FillValue', fail)
            call nc_check(ncStatus, ios)
            if (ios/=0) then
               ierr = ierr + ios
               write(error_msg,'(a)')'NetCDF-error adding flux "FillValue": '//trim(fileName)//' - '// &
                                     trim(flxName)//' - '//trim(NF90_STRERROR(ncStatus))
               call stop_tbnt(error_msg,ierr)
            end if
         end do
         TBNTlnkFlxID3D(:,iiFrac) = dummyID
      end do
      deallocate ( dimid, dummyID )
   end if
   
   ! end define mode   
   ncStatus = NF90_ENDDEF(fileID)
   call nc_check(ncStatus, ios)
   if (ios/=0) then
      ierr = ierr + ios
      write(error_msg,'(a)')'NetCDF-error ending define mode: '//trim(fileName)
      call stop_tbnt(error_msg,ierr)
   endif
   
   ! write grid coordinates data
   do i = 1,3
      if (allocated(NCgrid(i)%data1D)) then
         ncStatus = NF90_PUT_VAR(fileID, NCgrid(i)%ncid, NCgrid(i)%data1D)
      else
         ncStatus = NF90_PUT_VAR(fileID, NCgrid(i)%ncid, NCgrid(i)%data2D)
      end if
      call nc_check(ncStatus, ios)
      if (ios/=0) then
         ierr = ierr + ios
         write(error_msg,'(a)')'NetCDF-error writing grid information: '//trim(fileName)//' - '// &
                               trim(NCgrid(i)%name)//' - '//trim(NF90_STRERROR(ncStatus))
         call stop_tbnt(error_msg,ierr)
      end if
   end do
   
   end subroutine nc_create
!=======================================================================
! !INTERFACE:
   subroutine nc_open(filename, iFrac_, iSwitch, timeID, fileID, ierr)
#define SUBROUTINE_NAME 'nc_open'
!
! !DESCRIPTION:
!  open existing NetCDF file and get variable IDs
!
   implicit none
!
! !INPUT/OUTPUT PARAMETERS:
   character(len=fileLen), intent(in   ) :: filename     ! file name
   integer               , intent(in   ) :: iFrac_       ! current fraction
   integer               , intent(in   ) :: iSwitch      ! relative fractions (1), absolute fractions (2) or linked fluxes (3)
   integer               , intent(  out) :: timeID       ! time ID for current file
   integer               , intent(  out) :: fileID       ! file ID for current file
   integer               , intent(inout) :: ierr         ! error tracker
!
! LOCAL VARIABLES
   integer                            :: ios, ncStatus, iFrac1, iFrac2, iiFrac
   integer, dimension(:), allocatable :: dummyID
   character(len=nameLen)             :: varName, flxName
   logical                            :: lExist
!   
!-----------------------------------------------------------------------
#include "call-trace.inc"
!
   ios = 0
   
   ! open existing NetCDF file
   lexist = .false.
   inquire(file=trim(filename),exist=lExist)
   if (.not.lExist) then
      ierr = 1
      write(error_msg,'(a)')'Error - file for continued writing does not exist: '//trim(filename)
      call stop_tbnt(error_msg,ierr)
   end if
   ncStatus = NF90_OPEN(trim(filename), ior(NF90_WRITE, NF90_SHARE), fileID)
   call nc_check(ncStatus, ios)
   if (ios/=0) then
      ierr = ierr + ios
      write(error_msg,'(a)')'NetCDF-error opening existing file: '//trim(fileName)//' - '//trim(NF90_STRERROR(ncStatus))
      call stop_tbnt(error_msg,ierr)
   end if
   
   NCdims(4)%size = NF90_UNLIMITED
   NCgrid(4)%size = NF90_UNLIMITED
   
   ! get time dimension's ID (= unlimited dimension's ID)
   ncStatus = NF90_INQUIRE(fileID, unlimitedDimId = timeID)
   call nc_check(ncStatus, ios)
   if (ios/=0) then
      ierr = ierr + ios
      write(error_msg,'(a)')'NetCDF-error inquiring time ID from file: '//trim(filename)
      call stop_tbnt(error_msg,ierr)
   end if
   
   ! get output variable IDs
   if (iFrac_==0) then
      iFrac1 = 1
      iFrac2 = TBNTnFractions
   else
      iFrac1 = iFrac_
      iFrac2 = iFrac_
   end if
   
   if (iSwitch<=2) then  
      do iiFrac = iFrac1,iFrac2
#ifndef TBNTnoVar2D
         ! 2D variables
         allocate ( dummyID(TBNTnVars2D) )
         do iVar = 1,TBNTnVars2D
            varName = trim(TBNTvar(TBNTnVars3D+iVar)%name)//'-'//trim(TBNTfracName(iiFrac))
            ncStatus = NF90_INQ_VARID(fileID, trim(varName), dummyID(iVar))
            call nc_check(ncStatus, ios)
            if (ios/=0) then
               ierr = ierr + ios
               write(error_msg,'(a)')'NetCDF-error inquiring variable ID: '//trim(fileName)//' - '//trim(varName)// &
                                     ' - '//trim(NF90_STRERROR(ncStatus))
               call stop_tbnt(error_msg,ierr)
            end if
         end do
         if (iSwitch==1) then
            TBNTrelVarID2D(:,iiFrac) = dummyID
         else
            TBNTabsVarID2D(:,iiFrac) = dummyID
         end if
         deallocate ( dummyID )
#endif
         ! 3D variables
         allocate ( dummyID(TBNTnVars3D) )
         ! variables
         do iVar = 1,TBNTnVars3D
            varName = trim(TBNTvar(iVar)%name)//'-'//trim(TBNTfracName(iiFrac))
            ncStatus = NF90_INQ_VARID(fileID, trim(varName), dummyID(iVar))
            call nc_check(ncStatus, ios)
            if (ios/=0) then
               ierr = ierr + ios
               write(error_msg,'(a)')'NetCDF-error inquiring variable ID: '//trim(fileName)//' - '//trim(varName)// &
                                     ' - '//trim(NF90_STRERROR(ncStatus))
               call stop_tbnt(error_msg,ierr)
            end if
         end do
         if (iSwitch==1) then
            TBNTrelVarID3D(:,iiFrac) = dummyID
         else
            TBNTabsVarID3D(:,iiFrac) = dummyID
         end if
         deallocate ( dummyID )
      end do
   end if
   
   ! get fraction flux IDs only in case of linked fluxes or absolute fraction output
   if (iSwitch==2) then ! fraction fluxes
      ! 2D fluxes
      iOff = TBNTnFluxes3D+TBNTnLinkedFluxes3D
      do iFlux = 1,TBNTnFluxes2D
         flxName = trim(TBNTflux(iOff+iFlux)%name)//'-'//trim(TBNTfracName(iFrac_))
         ncStatus = NF90_INQ_VARID(fileID, trim(flxName), TBNTabsFlxID2D(iFlux,iFrac_))
         call nc_check(ncStatus, ios)
         if (ios/=0) then
            ierr = ierr + ios
            write(error_msg,'(a)')'NetCDF-error inquiring flux ID: '//trim(fileName)//' - '//trim(flxName)// &
                                  ' - '//trim(NF90_STRERROR(ncStatus))
            call stop_tbnt(error_msg,ierr)
         end if
      end do
      ! 3D fluxes
      do iFlux = 1,TBNTnFluxes3D
         flxName = trim(TBNTflux(iFlux)%name)//'-'//trim(TBNTfracName(iFrac_))
         ncStatus = NF90_INQ_VARID(fileID, trim(flxName), TBNTabsFlxID3D(iFlux,iFrac_))
         call nc_check(ncStatus, ios)
         if (ios/=0) then
            ierr = ierr + ios
            write(error_msg,'(a)')'NetCDF-error inquiring flux ID: '//trim(fileName)//' - '//trim(flxName)// &
                                  ' - '//trim(NF90_STRERROR(ncStatus))
            call stop_tbnt(error_msg,ierr)
         end if
      end do
   else if (iSwitch==3) then ! linked fluxes
      ! 2D fluxes
      iOff = TBNTnFluxes3D+TBNTnLinkedFluxes3D+TBNTnFluxes2D
      allocate ( dummyID(TBNTnLinkedFluxes2D) )
      do iiFrac = iFrac1,iFrac2
         do iFlux = 1,TBNTnLinkedFluxes2D
            flxName = trim(TBNTflux(iOff+iFlux)%name)//'-'//trim(TBNTfracName(iiFrac))
            ncStatus = NF90_INQ_VARID(fileID, trim(flxName), dummyID(iFlux))
            call nc_check(ncStatus, ios)
            if (ios/=0) then
               ierr = ierr + ios
               write(error_msg,'(a)')'NetCDF-error inquiring flux ID: '//trim(fileName)//' - '//trim(flxName)// &
                                     ' - '//trim(NF90_STRERROR(ncStatus))
               call stop_tbnt(error_msg,ierr)
            end if
         end do
         TBNTlnkFlxID2D(:,iiFrac) = dummyID
      end do  
      deallocate ( dummyID )
      ! 3D fluxes
      iOff = TBNTnFluxes3D
      allocate ( dummyID(TBNTnLinkedFluxes3D) )
      do iiFrac = iFrac1,iFrac2
         do iFlux = 1,TBNTnLinkedFluxes3D
            flxName = trim(TBNTflux(iOff+iFlux)%name)//'-'//trim(TBNTfracName(iiFrac))
            ncStatus = NF90_INQ_VARID(fileID, trim(flxName), dummyID(iFlux))
            call nc_check(ncStatus, ios)
            if (ios/=0) then
               ierr = ierr + ios
               write(error_msg,'(a)')'NetCDF-error inquiring flux ID: '//trim(fileName)//' - '//trim(flxName)// &
                                     ' - '//trim(NF90_STRERROR(ncStatus))
               call stop_tbnt(error_msg,ierr)
            end if
         end do
         TBNTlnkFlxID3D(:,iiFrac) = dummyID
      end do
      deallocate ( dummyID )
   end if
   
   end subroutine nc_open
!=======================================================================
!
! !INTERFACE:
   subroutine nc_write(filename, iFrac_, fileID, timeID, tbnt_time, iSwitch, ierr)
#define SUBROUTINE_NAME 'nc_write'
!
! !DESCRIPTION:
!  write netcdf file
!
! !USES:
   implicit none
!
! !OUTPUT PARAMETERS:
   character(len=fileLen), intent(in   ) :: filename
   integer               , intent(in   ) :: iFrac_, fileID, timeID, iSwitch
   integer               , intent(inout) :: ierr
   real(dp)              , intent(in   ) :: tbnt_time
!
! !LOCAL VARIABLES:
   integer                            :: ios, recOut, j, k, i, jSwitch, ncStatus
   integer                            :: ncStart2D(3), ncCount2D(3), ncStart3D(4), ncCount3D(4)
   real(dp)                           :: buffer2D(n_x,n_y), buffer3D(n_x,n_y,n_z)
   character(len=nameLen)             :: varName, flxName
   integer, dimension(:), allocatable :: dummyID
!
!-----------------------------------------------------------------------
#include "call-trace.inc"
   
   ios = 0
   recOut = 0
   
   select case (iSwitch)
      case (1)
         if (iFrac_==1) relFracRec = relFracRec + 1
         recOut = relFracRec
      case (2)
         if (iFrac_==1) absFracRec = absFracRec + 1
         recOut = absFracRec
      case (3)
         if (iFrac_==1) lnkFracRec = lnkFracRec + 1
         recOut = lnkFracRec
      case default
         ierr = ierr + 1
         write(error_msg,'(a)')'NetCDF-error: Invalid ''iSwitch'' value.'
         call stop_tbnt(error_msg,ierr)
   end select
   
   ! write time
   ncStatus = NF90_PUT_VAR(fileID, timeID, tbnt_time, (/recOut/))
   call nc_check(ncStatus, ios)
   if (ios/=0) then
      ierr = ierr + ios
      write(error_msg,'(a)')'NetCDF-error writing time: '//trim(filename)//' - '//trim(NF90_STRERROR(ncStatus))
      call stop_tbnt(error_msg,ierr)
   end if
   
   ! initialise output fields variables
   ncStart2D = (/   1,   1, recOut /)
   ncCount2D = (/ n_x, n_y,      1 /)
   buffer2D = fail
   ncStart3D = (/   1,   1,   1, recOut /)
   ncCount3D = (/ n_x, n_y, n_z,      1 /)
   buffer3D = fail
      
   if (iSwitch<=2) then
#ifndef TBNTnoVar2D
      ! 2D variables
      allocate ( dummyID(TBNTnVars2D) )
      if (iSwitch==1) then
         dummyID = TBNTrelVarID2D(:,iFrac_)
      else
         dummyID = TBNTabsVarID2D(:,iFrac_)
      end if
      do iVar = 1,TBNTnVars2D
         iVarPnt = TBNTvar2Dpnt(iVar,iFrac_)
         do j = TBNTyStart,TBNTyEnd
#ifdef TBNTswitchNS
            jSwitch = n_y-j+1
#else
            jSwitch = j
#endif
            do i = TBNTxStart,TBNTxEnd
               if (k_index(i,j)<=0) cycle
               iPut = index2Dto1D(i,j)
               if (excludeCell2D(iPut)) cycle

               if (iSwitch==1) then
                  buffer2D(i,jSwitch) = TBNTrelVar2D(iPut,iVarPnt)*100.e0 ! percent
               else
                  buffer2D(i,jSwitch) = TBNTvar2D(iPut,iVarPnt)
               end if
            end do
         end do
         ncStatus = NF90_PUT_VAR(fileID, dummyID(iVar), buffer2D, ncStart2D, ncCount2D)
         call nc_check(ncStatus, ios)
         if (ios/=0) then
            ierr = ierr + ios
            varName = trim(TBNTvar(TBNTnVars3D+iVar)%name)//'-'//trim(TBNTfracName(iFrac_))
            write(error_msg,'(a)')'NetCDF-error writing variable: '//trim(filename)//' - '// &
                                  trim(varName)//' - '//trim(NF90_STRERROR(ncStatus))
            call stop_tbnt(error_msg,ierr)
         end if
      end do
      deallocate ( dummyID )
#endif
      ! write 3D fields
      allocate ( dummyID(TBNTnVars3D) )
      if (iSwitch==1) then
         dummyID = TBNTrelVarID3D(:,iFrac_)
      else
         dummyID = TBNTabsVarID3D(:,iFrac_)
      end if
      ! variables
      do iVar = 1,TBNTnVars3D
         iVarPnt = TBNTvar3Dpnt(iVar,iFrac_)
         do j = TBNTyStart,TBNTyEnd
#ifdef TBNTswitchNS
            jSwitch = n_y-j+1
#else
            jSwitch = j
#endif
            do i = TBNTxStart,TBNTxEnd
               if (k_index(i,j)<=0) cycle
               do k = 1,k_index(i,j)
                  iPut = index3Dto1D(i,j,k)
                  if (excludeCell3D(iPut)) cycle
                  if (iSwitch==1) then
                     buffer3D(i,jSwitch,k) = TBNTrelVar3D(iPut,iVarPnt)*100.e0 ! percent
                  else
                     buffer3D(i,jSwitch,k) = TBNTvar3D(iPut,iVarPnt)
                  end if
               end do
            end do
         end do
         ncStatus = NF90_PUT_VAR(fileID, dummyID(iVar), buffer3D, ncStart3D, ncCount3D)
         call nc_check(ncStatus, ios)
         if (ios/=0) then
            ierr = ierr + ios
            varName = trim(TBNTvar(iVar)%name)//'-'//trim(TBNTfracName(iFrac_))
            write(error_msg,'(a)')'NetCDF-error writing variable: '//trim(filename)//' - '// &
                                  trim(varName)//' - '//trim(NF90_STRERROR(ncStatus))
            call stop_tbnt(error_msg,ierr)
         end if
      end do
      deallocate ( dummyID )
   end if
   
   if (iSwitch==1) return ! write fluxes only in case of absolute fraction output or linked fluxes output
   
   ! 2D fluxes
   if (iSwitch==2) then ! absolute fractions
      do iFlux = 1,TBNTnFluxes2D
         iFlxPnt = TBNTflux2Dpnt(iFlux,iFrac_)
         do j = TBNTyStart,TBNTyEnd
#ifdef TBNTswitchNS
            jSwitch = n_y-j+1
#else
            jSwitch = j
#endif
            do i = TBNTxStart,TBNTxEnd
               if (k_index(i,j)<=0) cycle
               iPut = index2Dto1D(i,j)
               if (excludeCell2D(iPut)) cycle
               buffer2D(i,jSwitch) = TBNTflux2D(iPut,iFlxPnt)
            end do
         end do
         ncStatus = NF90_PUT_VAR(fileID, TBNTabsFlxID2D(iFlux,iFrac_), buffer2D, ncStart2D, ncCount2D)
         call nc_check(ncStatus, ios)
         if (ios/=0) then
            ierr = ierr + ios
            flxName = trim(TBNTflux(TBNTnFluxes3D+TBNTnLinkedFluxes3D+iFlux)%name)//'-'//trim(TBNTfracName(iFrac_))
            write(error_msg,'(a)')'NetCDF-error writing variable: '//trim(filename)//' - '// &
                                  trim(flxName)//' - '//trim(NF90_STRERROR(ncStatus))
            call stop_tbnt(error_msg,ierr)
         end if
      end do
      ! 3D fluxes
      do iFlux = 1,TBNTnFluxes3D
         iFlxPnt = TBNTflux3Dpnt(iFlux,iFrac_)
         do j = TBNTyStart,TBNTyEnd
#ifdef TBNTswitchNS
            jSwitch = n_y-j+1
#else
            jSwitch = j
#endif
            do i = TBNTxStart,TBNTxEnd
               if (k_index(i,j)<=0) cycle
               do k = 1,k_index(i,j)
                  iPut = index3Dto1D(i,j,k)
                  if (excludeCell3D(iPut)) cycle
                  buffer3D(i,jSwitch,k) = TBNTflux3D(iPut,iFlxPnt)
               end do
            end do
         end do
         ncStatus = NF90_PUT_VAR(fileID, TBNTabsFlxID3D(iFlux,iFrac_), buffer3D, ncStart3D, ncCount3D)
         call nc_check(ncStatus, ios)
         if (ios/=0) then
            ierr = ierr + ios
            flxName = trim(TBNTflux(iFlux)%name)//'-'//trim(TBNTfracName(iFrac_))
            write(error_msg,'(a)')'NetCDF-error writing variable: '//trim(filename)//' - '// &
                                  trim(flxName)//' - '//trim(NF90_STRERROR(ncStatus))
            call stop_tbnt(error_msg,ierr)
         end if
      end do
   else ! iSwitch==3: linked fluxes
      ! 2D fluxes
      do iFlux = 1,TBNTnLinkedFluxes2D
         iFlxPnt = TBNTlinkedFlux2Dpnt(iFlux,iFrac_)
         do j = TBNTyStart,TBNTyEnd
#ifdef TBNTswitchNS
            jSwitch = n_y-j+1
#else
            jSwitch = j
#endif
            do i = TBNTxStart,TBNTxEnd
               if (k_index(i,j)<=0) cycle
               iPut = index2Dto1D(i,j)
               if (excludeCell2D(iPut)) cycle
               buffer2D(i,jSwitch) = TBNTlinkedFlux2D(iPut,iFlxPnt)
            end do
         end do
         ncStatus = NF90_PUT_VAR(fileID, TBNTlnkFlxID2D(iFlux,iFrac_), buffer2D, ncStart2D, ncCount2D)
         call nc_check(ncStatus, ios)
         if (ios/=0) then
            ierr = ierr + ios
            flxName = trim(TBNTflux(TBNTnFluxes3D+TBNTnLinkedFluxes3D+TBNTnFluxes2D+iFlux)%name)//'-'//trim(TBNTfracName(iFrac_))
            write(error_msg,'(a)')'NetCDF-error writing variable: '//trim(filename)//' - '// &
                                  trim(flxName)//' - '//trim(NF90_STRERROR(ncStatus))
            call stop_tbnt(error_msg,ierr)
         end if
      end do
      ! 3D fluxes
      do iFlux = 1,TBNTnLinkedFluxes3D
         iFlxPnt = TBNTlinkedFlux3Dpnt(iFlux,iFrac_)
         do j = TBNTyStart,TBNTyEnd
#ifdef TBNTswitchNS
            jSwitch = n_y-j+1
#else
            jSwitch = j
#endif
            do i = TBNTxStart,TBNTxEnd
               if (k_index(i,j)<=0) cycle
               do k = 1,k_index(i,j)
                  iPut = index3Dto1D(i,j,k)
                  if (excludeCell3D(iPut)) cycle
                  buffer3D(i,jSwitch,k) = TBNTlinkedFlux3D(iPut,iFlxPnt)
               end do
            end do
         end do
         ncStatus = NF90_PUT_VAR(fileID, TBNTlnkFlxID3D(iFlux,iFrac_), buffer3D, ncStart3D, ncCount3D)
         call nc_check(ncStatus, ios)
         if (ios/=0) then
            ierr = ierr + ios
            flxName = trim(TBNTflux(TBNTnFluxes3D+iFlux)%name)//'-'//trim(TBNTfracName(iFrac_))
            write(error_msg,'(a)')'NetCDF-error writing variable: '//trim(filename)//' - '// &
                                  trim(flxName)//' - '//trim(NF90_STRERROR(ncStatus))
            call stop_tbnt(error_msg,ierr)
         end if
      end do
   end if
   
   end subroutine nc_write
!=======================================================================
!
! !INTERFACE:
   subroutine asc_write_header(fUnit,iTargetVar,ierr)
#define SUBROUTINE_NAME 'asc_write_header'
!
! !DESCRIPTION:
!  write header to target output file
!
! !USES:
   implicit none
!
! !OUTPUT PARAMETERS:
!
   integer               , intent(in   ) :: fUnit, iTargetVar
   integer               , intent(inout) :: ierr
!
! !LOCAL VARIABLES:
!
   integer         , parameter :: maxLen = 10*lineLen
   character(len=1), parameter :: separator = ';'
   
   integer                     :: iArea
   character(len=maxLen)       :: headerStr
!
!-----------------------------------------------------------------------
#include "call-trace.inc"

   ! create string header string with area names
   headerStr = trim(TBNTtargetVars(iTargetVar)%name)//separator//trim(TBNTtargetAreas(1)%name)
   do iArea = 2,TBNTnTargetAreas
      headerStr = trim(headerStr)//separator//trim(TBNTtargetAreas(iArea)%name)
   end do
   if (len_trim(headerStr)>maxLen) then
      ierr = ierr + 1
      write(error_msg,'(2(a,i5),a)')'Header line length exceeds maximum allowed length: ',len_trim(headerStr), &
                                    '> ',maxLen, new_line('a')//'Check length of target area names.'
      call stop_tbnt(error_msg,ierr)
   end if
   write(fUnit,'(a)')trim(headerStr)

   end subroutine asc_write_header
!=======================================================================
!
! !INTERFACE:
   subroutine asc_write_data(fUnit,iTargetVar,ierr)
#define SUBROUTINE_NAME 'asc_write'
!
! !DESCRIPTION:
!  write target output file
!
! !USES:
   implicit none
!
! !OUTPUT PARAMETERS:
!
   integer     , intent(in   ) :: fUnit, iTargetVar
   integer     , intent(inout) :: ierr
!
! !LOCAL VARIABLES:
!
   integer         , parameter :: maxLen = 10*lineLen
   character(len=1), parameter :: separator = ';'
   
   integer                     :: ios
   character(len=formatLen)    :: formatStr
!
!-----------------------------------------------------------------------
#include "call-trace.inc"

   formatStr = '(a,";",    (E20.10,";"),E20.10)'
   write(formatStr(8:11),'(i4)') TBNTnTargetAreas-1
   
   ! write data to file
   write(fUnit,trim(formatStr)) 'totMass', TBNTtargetWeights(:,iTargetVar)
   do iFrac = 1,TBNTnFractions-1
      write(fUnit,trim(formatStr),iostat=ios) trim(TBNTfracName(iFrac)), TBNTtargetVals(:,iTargetVar,iFrac)
      if (ios/=0) then
         ierr = ierr + ios
         write(error_msg,'(a)')'Error writing target area data.'
         call stop_tbnt(error_msg,ierr)
      end if
   end do

   end subroutine asc_write_data
!=======================================================================
#endif
end module mod_etrac_output
