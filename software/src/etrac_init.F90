!
! !MODULE: etrac_init.f90  --- initializing ETRAC 
!
! !INTERFACE:
   MODULE mod_etrac_init
!
! !DESCRIPTION:
!
! !USES:
   use netcdf
   use mod_etrac_common
   
   implicit none
!
!  default: all is private.
   private
   
   public init_tbnt, get_nc_field
!
! !LOCAL VARIABLES:
!
   character(len=fileLen), parameter :: work_dir = './'
   
   integer                :: n_total_flx, n_total_var, n_progn_var, n_diagn_var, n_dummy_var
   character(len=fileLen) :: fluxes_vars_file, dummy_vars_file, rivers_file
   character(len=fileLen) :: bulk_nc_dir, bulk_nc_file, init_nc_dir, init_nc_file
#ifndef TBNTonly_bulk_bud
   character(len=fileLen) :: target_area_file, target_vars_file, linked_fluxes_file
#endif
   logical                :: includeAtmos, useDummySrc
!=======================================================================

   contains

!=======================================================================
!
! !INTERFACE:
   subroutine init_tbnt(ierr)
#define SUBROUTINE_NAME 'init_tbnt'
!
! !DESCRIPTION:
! initialize TBNT
!
! !USES:
   implicit none
!
! !OUTPUT PARAMETERS:
   integer, intent(inout) :: ierr
!
! !LOCAL VARIABLES:
   integer                :: ios, system
#ifdef TBNTonly_bulk_bud
   integer                :: idum
#endif
   character(len=yearLen) :: yearStr
   logical                :: lExist
!
! !FOR NAMELISTS
!
   integer                :: isWarmStart, continueWrite, year, timeStep, startStep, endStep, offsetStep, subTimeStep
   integer                :: nTracedElements, nRiverFractions, nOpenBFractions, nAtmosFractions, initialFraction
   character(len=nameLen) :: runID
#ifndef TBNTonly_bulk_bud
   integer                :: relFracOutStep, absFracOutStep, targetOutput, linkedFlxOutStep
   character(len=fileLen) :: linkedFlx_file
#endif

   namelist /tbnt_run_nml/ runID, isWarmstart, continueWrite, tbnt_set_dir, bulk_nc_dir, bulk_nc_file, init_nc_dir, init_nc_file

   namelist /time_setup_nml/ year, timeStep, startStep, endStep, offsetStep, subTimeStep
   
   namelist /area_vol_nml/ areaVar, volVar
   
   namelist /fraction_numbers_nml/ nTracedElements, nRiverFractions, nOpenBFractions, nAtmosFractions, initialFraction
   
#ifndef TBNTonly_bulk_bud
   namelist /linked_fluxes_nml/ linkedFlxOutStep, linkedFlx_file

   namelist /output_nml/ output_dir, relFracOutStep, absFracOutStep, targetOutput, target_area_file, target_vars_file
#endif

!-----------------------------------------------------------------------
#include "call-trace.inc"

   ! open TBNT log file
   open(tbnt_log_unit,file=trim(work_dir)//trim(tbnt_log_file),action='write',status='replace')
#ifndef TBNTonly_bulk_bud
   write(tbnt_log_unit, '(a)') 'initialize TBNT'
   call breakLine(100)
#endif

   ! open TBNT settings file: tbnt_set.nml
   init_nc_file = ' ' ! to get rid of 'unused'  warning
   open(tbnt_settings_unit,file=trim(work_dir)//trim(tbnt_settings_file),action='read',status='old')
   ! run ID & warm start
   read(tbnt_settings_unit,nml=tbnt_run_nml)
   TBNTrun = runID
   if (isWarmstart==1) then
      TBNTwarm = .true.
   else
      TBNTwarm = .false.
   end if
   ! continue writing results to existing files
   if (continueWrite==1) then
      TBNTcontinue = .true.
   else
      TBNTcontinue = .false.
   end if
   ! time setup
   read(tbnt_settings_unit,nml=time_setup_nml)
   TBNTyear = year
   ! check for leap year
   if ((mod(year,4)==0.and.mod(year,100)/=0).or.mod(year,400)==0) then
      TBNTyearDays = 366
   else
      TBNTyearDays = 365
   end if
   TBNTstep = timeStep
   ! select whole year if startStep or endStep = 0
   if (startStep==0.or.endStep==0) then
      startStep = 1
      endStep = TBNTyearDays
   else
      startStep = max(1,startStep)
      endStep = min(TBNTyearDays,endStep)
   end if
   TBNTstart = timeStep*startStep
   TBNTend = timeStep*endStep
   TBNTiOffset = offsetStep
   TBNToffset = timeStep*offsetStep
#ifndef TBNTonly_bulk_bud
   if (mod(timeStep,subTimeStep)/=0) then
      ierr = ierr + 1
      write(error_msg,'(a)')'Invalid TBNT time setup: subTimeStep must be a proper divider of timeStep. Check time_setup_nml.'
      call stop_tbnt(error_msg,ierr)
   end if     
   TBNTnSubSteps = timeStep/subTimeStep
   TBNTtimeFac = 1.e0/real(TBNTnSubSteps,dp)
#else
   idum = subTimeStep
#endif
   if (mod(iMinPerDay,TBNTstep)/=0) then
      ierr = ierr + 1
      write(error_msg,'(a)')'Invalid TBNT time setup: timeStep must be a proper divider '// &
                            'of 1440 (minutes per day). Check time_setup_nml.'
      call stop_tbnt(error_msg,ierr)
   end if
   rMinPerDay = real(iMinPerDay,dp)
   TBNTiStart = int(real(TBNTstart,dp)/real(TBNTstep,dp))
   TBNTiEnd   = int(real(TBNTend,dp)/real(TBNTstep,dp))
   TBNTnSteps = TBNTiEnd - TBNTiStart + 1
   TBNTnDays  = ceiling(real(TBNTnSteps,dp)*real(TBNTstep,dp)/rMinPerDay)
   
   ! create run ID consisting of run name and year
   write(yearStr,'(i4)') TBNTyear
   TBNTrun = trim(TBNTrun)//'_'//yearStr
   
   ! get names of area and volume state variables to be read from NetCDF files
   areaVar = ' '
   volVar  = ' '
   read(tbnt_settings_unit,nml=area_vol_nml)

   ! initialise the model input files and grid
   call init_model(ierr)
   
   ! read basic TBNT setup (number of different tracers)
   read(tbnt_settings_unit,nml=fraction_numbers_nml)
   TBNTiniFrac = initialFraction
   TBNTnTracedElements = nTracedElements
   TBNTnRiverFractions = nRiverFractions
   TBNTnOpenBFractions = nOpenBFractions
   TBNTnAtmosFractions = nAtmosFractions
   TBNTnFractions = nRiverFractions + nOpenBFractions + nAtmosFractions 
   if (TBNTnFractions == 0) then
      ierr = ierr + 1
      write(error_msg,'(a)')'No tracer sources defined.'
      call stop_tbnt(error_msg,ierr)
   end if

   ! initialise the list of traced fluxes and variables
   call init_tbnt_fluxes_vars(ierr)

#ifndef TBNTonly_bulk_bud
   ! initialise linked fluxes
   TBNTlinkedFluxes = .false.
   TBNTnLinkedFluxFractions3D = 0
   TBNTnLinkedFluxFractions2D = 0
   TBNTlnkFracOutStep = -1
   read(tbnt_settings_unit,nml=linked_fluxes_nml)
   if (linkedFlxOutStep>=0) then
      TBNTlinkedFluxes = .true.
      linked_fluxes_file = linkedFlx_file
      TBNTlnkFracOutStep = linkedFlxOutStep
      call init_tbnt_linked_fluxes(ierr)
   end if
   
   ! initialize tracer sources
   call init_tbnt_sources(ierr)
   
   ! initialise output
   output_dir = ' '
   TBNTrelFracOutStep = -1
   TBNTabsFracOutStep = -1
   
   read(tbnt_settings_unit,nml=output_nml)   
   
   TBNTrelFracOutStep = relFracOutStep
   TBNTabsFracOutStep = absFracOutStep
   
   ! create output directory
   if (TBNTrelFracOutStep>=0.or.TBNTabsFracOutStep>=0.or.linkedFlxOutStep>=0.or.targetOutput==1) then
      ios = 0
      lExist = .false.
      inquire(file=trim(output_dir),exist=lExist)
      if (.not.lExist) ios = system('mkdir -p '//trim(output_dir))
      if (ios/=0) then
         ierr = ierr + ios
         write(error_msg,'(a)')'Error creating output directory: '//trim(output_dir)
         call stop_tbnt(error_msg,ierr)
      end if
   end if
   
   ! initialise target areas & variables
   if (targetOutput==1) then
      TBNTtargetOut = .true.
      call init_tbnt_target_areas(ierr)
      call init_tbnt_target_vars(ierr)
   else if (targetOutput==0) then
      TBNTtargetOut = .false.
   else
      ierr = ierr + 1
      write(error_msg,'(a)')'Invalid value for switch: targetOutput. Check output_nml.'
      call stop_tbnt(error_msg,ierr)
   end if

   ! initialise fraction names if warm start file is used or output will be created
   if (TBNTwarm.or.TBNTlinkedFluxes.or.TBNTrelFracOutStep>=0.or.TBNTabsfracOutStep>=0) call init_tbnt_fraction_names

   ! create empty fields and pointer matrix for ALL variable fractions
   call init_tbnt_fraction_fields(ierr)
#else
   ! initialize bulk variable and fluxes to be balanced
   call init_bulk_budget(ierr)
#endif
   
   ! open hydrodynamic files
   call init_cell_size(ierr)

   ! open netcdf file with bulk fluxes
   call init_bulk_fields
   
   ! close TBNT namelist
   close(tbnt_settings_unit)
   
   end subroutine init_tbnt

!=======================================================================
!
! !INTERFACE:
   subroutine init_model(ierr)
#define SUBROUTINE_NAME 'init_model'
!
! !DESCRIPTION:
! initialise the underlying model
!
! !USES:

   implicit none
!
! !OUTPUT PARAMETERS:
   integer, intent(inout) :: ierr
   integer                :: ios, ncStatus
   character(len=fileLen) :: filename
   logical                :: lExist
! !FOR NAME LIST
   character(len=fileLen) :: model_grid_file, model_rivers_file, model_fluxes_file, model_dummy_vars_file
   
   namelist /model_nml/ model_set_dir, model_grid_file, model_rivers_file, model_fluxes_file, model_dummy_vars_file
!
!-----------------------------------------------------------------------
#include "call-trace.inc"

   read(tbnt_settings_unit, nml=model_nml)
   ! check existence of input files
   ! grid file
   model_grid_file = trim(model_set_dir)//trim(model_grid_file)
   inquire(file = trim(model_grid_file), exist = lExist)
   if (.not.lExist) then
      ierr = 1
      write(error_msg,'(2a)')'File does not exist: ',trim(model_grid_file)
      call stop_tbnt(error_msg,ierr)
   end if
   
   ! river list
   model_rivers_file = trim(model_set_dir)//trim(model_rivers_file)
   inquire(file = trim(model_rivers_file), exist = lExist)
   if (.not.lExist) then
      ierr = 1
      write(error_msg,'(2a)')'File does not exist: ',trim(model_rivers_file)
      call stop_tbnt(error_msg,ierr)
   end if
   ! flux list
   model_fluxes_file = trim(model_set_dir)//trim(model_fluxes_file)
   inquire(file = trim(model_fluxes_file), exist = lExist)
   if (.not.lExist) then
      ierr = 1
      write(error_msg,'(2a)')'File does not exist: ',trim(model_fluxes_file)
      call stop_tbnt(error_msg,ierr)
   end if
   
   fluxes_vars_file = model_fluxes_file ! used by 'init_tbnt_fluxes_vars'
   rivers_file      = model_rivers_file ! used by 'init_tbnt_rivers'
   
   dummy_vars_file = model_dummy_vars_file ! used by 'init_tbnt_fluxes_vars' in case of dummy variables existing in current set-up
   
   ! initialise grid
   call init_model_grid(model_grid_file,ierr)

   filename = trim(bulk_nc_dir)//trim(bulk_nc_file)
   ncStatus = NF90_OPEN(trim(filename), NF90_NOWRITE, bulk_nc_unit)
   ios = 0
   call nc_check(ncStatus, ios)
   if (ios/=0) then
      ierr = ierr + ios
      write(error_msg,'(2a)')'NetCDF-error opening bulk input file: ',trim(filename)
      call stop_tbnt(error_msg,ierr)
   end if
   
   ! get grid dimensions from bulk NetCDF file
   ! dimensions are later on copied to output files
   call get_nc_dimensions(filename, bulk_nc_unit, ierr)

   end subroutine init_model
!=======================================================================
!
! !INTERFACE:
   subroutine init_model_grid(filename, ierr)
#define SUBROUTINE_NAME 'init_model_grid'
!
! !DESCRIPTION:
! initialise the underlying model grid
!
! !USES:

   implicit none
!
! !OUTPUT PARAMETERS:
   integer               , intent(inout) :: ierr
   character(len=fileLen), intent(in   ) :: filename
!
! !LOCAL VARIABLES:
!
   integer :: nSkip, i, ios
   integer, dimension(:), allocatable :: iSurface, iBottom

#ifdef TBNTconvert_3Dto1D
   integer :: iiPut, j, k, switch_iDepNS
   logical :: lExist

   character(len=formatLen) :: frmt

  !FOR NAMELIST
   integer                :: offset_xStart, offset_xEnd, offset_yStart, offset_yEnd
   character(len=nameLen) :: dimName_x, dimName_y, dimName_z, gridName_x, gridName_y, gridName_z
   character(len=fileLen) :: model_idep_file

   namelist /indexing_nml/ dimName_x, dimName_y, dimName_z, gridName_x, gridName_y, gridName_z,   &
                           n_x, n_y, n_z, offset_xStart, offset_xEnd, offset_yStart, offset_yEnd, &
                           model_idep_file
#endif
!-----------------------------------------------------------------------
#include "call-trace.inc"
   
   ios = 0
   
   open(input_unit, file = trim(filename), action = 'read')
   ! check header length
   call get_headLen(input_unit, filename, nSkip, ierr)
   ! skip header
   call skipLines(input_unit,nSkip)
   ! read grid neighbours for each directional component
   read(input_unit,*, iostat = ierr) nWetCells, nComponents
   if (ios/=0) then
      ierr = ierr + ios
      write(error_msg,'(2a)')'Error reading: ',trim(filename)
      call stop_tbnt(error_msg,ierr)
   end if
   allocate ( neighbours(nWetCells,nComponents), compSign(nComponents) )
   neighbours = 0
   compSign = 0
   do i = 1,nComponents
      read(input_unit, *, iostat = ierr) compSign(i), neighbours(:,i)
      if (ios/=0) exit
   end do
   if (ios/=0) then
      ierr = ierr + ios
      write(error_msg,'(2a)')'Error reading: ',trim(filename)
      call stop_tbnt(error_msg,ierr)
   end if
   ! read surface cell indices
   read(input_unit,*, iostat = ierr) nSurfCells
   if (ios/=0) then
      ierr = ierr + ios
      write(error_msg,'(2a)')'Error reading: ',trim(filename)
      call stop_tbnt(error_msg,ierr)
   end if
   allocate ( iSurface(nSurfCells) )
   iSurface = 0
   read(input_unit, *, iostat = ierr) iSurface
   if (ios/=0) then
      ierr = ierr + ios
      write(error_msg,'(2a)')'Error reading: ',trim(filename)
      call stop_tbnt(error_msg,ierr)
   end if
   ! create surface mask and pairs of pelagic/surface cells
   allocate ( surfaceCell(nWetCells), surf2pelag(nSurfCells) )
   surfaceCell = .false.
   do i = 1,nSurfCells
      surfaceCell(iSurface(i)) = .true.
      surf2pelag(i) = iSurface(i)
   end do
   deallocate ( iSurface )
   ! read bottom cell indices
   read(input_unit,*, iostat = ierr) nBotCells
   if (ios/=0) then
      ierr = ierr + ios
      write(error_msg,'(2a)')'Error reading: ',trim(filename)
      call stop_tbnt(error_msg,ierr)
   end if
   allocate ( iBottom(nBotCells) )
   iBottom = 0
   read(input_unit, *, iostat = ierr) iBottom
   if (ios/=0) then
      ierr = ierr + ios
      write(error_msg,'(2a)')'Error reading: ',trim(filename)
      call stop_tbnt(error_msg,ierr)
   end if
   ! create bottom mask and pairs of pelagic/bottom cells
   allocate ( bottomCell(nWetCells), bottom2pelag(nBotCells) )
   bottomCell = .false.
   do i = 1,nBotCells
      bottomCell(iBottom(i)) = .true.
      bottom2pelag(i) = iBottom(i)
   end do
   deallocate ( iBottom )   
   
   close(input_unit)
   
   ! initialise masks to exclude cells from calculation
   allocate ( excludeCell3D(nWetCells), excludeCell2D(nBotCells) )
   excludeCell3D = .false.
   excludeCell2D = .false.
   
#ifdef TBNTconvert_3Dto1D
   dimName_x = ' '
   dimName_y = ' '
   dimName_z = ' '
   gridName_x = ' '
   gridName_y = ' '
   gridName_z = ' '
   n_x = 0
   n_y = 0
   n_z = 0
   offset_xStart = 0
   offset_xEnd   = 0
   offset_yStart = 0
   offset_yEnd   = 0
   model_idep_file = ' '

   read(tbnt_settings_unit, nml=indexing_nml)

   ! set index range in x- and y-direction
   TBNTxStart = 1 + offset_xStart
   TBNTxEnd = n_x - offset_xEnd
   TBNTyStart = 1 + offset_yStart
   TBNTyEnd = n_y - offset_yEnd
   
   ! set dimension names for NetCDF I/O
   NCdims(1)%name = trim(dimName_x)
   NCdims(2)%name = trim(dimName_y)
   NCdims(3)%name = trim(dimName_z)
   
   NCgrid(1)%name = trim(gridName_x)
   NCgrid(2)%name = trim(gridName_y)
   NCgrid(3)%name = trim(gridName_z)

   ! read map with maximum depth indices
   lExist = .false.
   inquire(file = trim(model_set_dir)//trim(model_idep_file), exist = lExist)
   if (.not.lExist) then
      ierr = 1
      write(error_msg,'(2a)')'File does not exist: ',trim(model_set_dir)//trim(model_idep_file)
      call stop_tbnt(error_msg,ierr)
   end if

   open(input_unit, file=trim(model_set_dir)//trim(model_idep_file), action='read', status='old')
   ! check header length
   call get_headLen(input_unit, trim(model_set_dir)//trim(model_idep_file), nSkip, ierr)
   ! skip header
   call skipLines(input_unit,nSkip)

   allocate ( k_index(n_x,n_y) )
   k_index = 0
   switch_iDepNS = 0

   frmt = '(    i4)'
   write(frmt(2:5), '(i4)') n_x
   read(input_unit,*) switch_iDepNS
   do j = 1,n_y
      if (switch_iDepNS==0) then
         read(input_unit, trim(frmt), iostat=ios) k_index(:,j)
      else
         read(input_unit, trim(frmt), iostat=ios) k_index(:,n_y-j+1)
      end if
      if (ios/=0) then
         ierr = 1
         write(error_msg,'(2a)')'Error reading from file: ',trim(model_idep_file)
         call stop_tbnt(error_msg,ierr)
      end if
   end do
   close(input_unit)
   
   ! create index fields pointing from 1D to 3D/2D and vice versa
   allocate ( index3Dto1D(n_x, n_y, n_z), index2Dto1D(n_x,n_y), &
              index1Dto3D(nWetCells,3), index1Dto2D(nBotCells,2) )
   index3Dto1D = 0
   index2Dto1D = 0
   index1Dto3D = 0
   index1Dto2D = 0
   
   iPut = 0
   iiPut = 0
   do j = 1,n_y
      do i = 1,n_x
         if (k_index(i,j)<=0) cycle
         do k = 1,k_index(i,j)
            iPut = iPut + 1
            index3Dto1D(i,j,k) = iPut
            index1Dto3D(iPut,1) = i
            index1Dto3D(iPut,2) = j
            index1Dto3D(iPut,3) = k
         end do
         iiPut = iiPut + 1
         index2Dto1D(i,j) = iiPut
         index1Dto2D(iiPut,1) = i
         index1Dto2D(iiPut,2) = j
      end do
   end do
#endif
   
   end subroutine init_model_grid
!=======================================================================
#ifdef TBNTconvert_3Dto1D
!
! !INTERFACE:
subroutine get_nc_dimensions(fileName, fileUnit, ierr)
#define SUBROUTINE_NAME 'get_nc_dimensions'
!
! !DESCRIPTION:
! get dimension names and information from bulk NetCDF file for use in output files
!
   implicit none
!
! !OUTPUT PARAMETERS:
   character(len=nameLen), intent(in   ) :: fileName
   integer               , intent(in   ) :: fileUnit
   integer               , intent(inout) :: ierr
!
! !LOCAL VARIABLES:
!
   integer                             :: ncStatus, ios, i, j, k,timeID, nDims, nVarDims
   integer , dimension(:), allocatable :: dimIDs
   real(dp), dimension(:), allocatable :: outTime
   character(len=nameLen)              :: unitStr
   logical                             :: found
!
!-----------------------------------------------------------------------
#include "call-trace.inc"
   
   ios = 0
   nDims = 0
   timeID = 0
   
   ! get number of dimensions and time dimension ID
   ncStatus = NF90_INQUIRE(fileUnit, nDimensions = nDims, unlimitedDimId = timeID)
   call nc_check(ncStatus, ios)
   if (ios/=0) then
      ierr = ierr + ios
      write(error_msg,'(a)')'NetCDF-error inquiring number of dimensions and time ID: '//trim(fileName)
      call stop_tbnt(error_msg,ierr)
   end if
   if (nDims<4) then
      ierr = 1
      write(error_msg,'(a,i2,a)')'NetCDF-error regarding number of dimensions: '//trim(fileName)//  &
                                 new_line('a')//'Minimum is 4, only ', nDims, ' found.'
      call stop_tbnt(error_msg,ierr)
   end if
   if (timeID<=0) then
      ierr = 1
      write(error_msg,'(a)')'NetCDF-error finding time dimension index: '//trim(fileName)
      call stop_tbnt(error_msg,ierr)
   end if
   
   ! get dimensions' IDs and sizes
   NCdims%ncid = 0
   NCdims%size = 0
   ! spatial dimensions
   do i = 1,3
      ncStatus = NF90_INQ_DIMID(fileUnit, trim(NCdims(i)%name), NCdims(i)%ncid)
      call nc_check(ncStatus, ios)
      if (ios/=0) then
         ierr = ierr + ios
         write(error_msg,'(a)')'NetCDF-error inquiring dimension ID for '//trim(NCdims(i)%name)//': '//trim(fileName)
         call stop_tbnt(error_msg,ierr)
      end if
      if (NCdims(i)%ncid<=0) then
         ierr = 1
         write(error_msg,'(a)')'NetCDF-error - invalid dimension ID: '//trim(NCdims(i)%name)
         call stop_tbnt(error_msg,ierr)
      end if
      ncStatus = NF90_INQUIRE_DIMENSION(fileUnit, NCdims(i)%ncid, len = NCdims(i)%size)
      call nc_check(ncStatus, ios)
      if (ios/=0) then
         ierr = ierr + ios
         write(error_msg,'(a)')'NetCDF-error inquiring dimension length: '//trim(NCdims(i)%name)
         call stop_tbnt(error_msg,ierr)
      end if
      if (NCdims(i)%size<=0) then
         ierr = 1
         write(error_msg,'(a)')'NetCDF-error - invalid dimension length: '//trim(NCdims(i)%name)
         call stop_tbnt(error_msg,ierr)
      end if
   end do
   ! time dimension
   NCdims(4)%ncid = timeID
   NCdims(4)%name = ' '
   ncStatus = NF90_INQUIRE_DIMENSION(fileUnit, timeID, NCdims(4)%name, NCdims(4)%size)
   call nc_check(ncStatus, ios)
   if (ios/=0) then
      ierr = ierr + ios
      write(error_msg,'(a)')'NetCDF-error inquiring unlimited dimension name and length'
      call stop_tbnt(error_msg,ierr)
   end if
   if (len_trim(NCdims(i)%name)==0) then
      ierr = 1
      write(error_msg,'(a)')'NetCDF-error - invalid name for unlimited dimension'
      call stop_tbnt(error_msg,ierr)
   end if
   if (NCdims(i)%size<=0) then
      ierr = 1
      write(error_msg,'(a)')'NetCDF-error - invalid length for unlimited dimension'
      call stop_tbnt(error_msg,ierr)
   end if
   
   ! get grid variables' IDs and sizes
   NCgrid%unit = ' '
   NCgrid%ncid = 0
   ! spatial information
   do i = 1,3
      ncStatus = NF90_INQ_VARID(fileUnit, trim(NCgrid(i)%name), NCgrid(i)%ncid)
      call nc_check(ncStatus, ios)
      if (ios/=0) then
         ierr = ierr + ios
         write(error_msg,'(a)')'NetCDF-error inquiring grid variable ID: '//trim(NCgrid(i)%name)
         call stop_tbnt(error_msg,ierr)
      end if
      if (NCgrid(i)%ncid<=0) then
         ierr = 1
         write(error_msg,'(a)')'NetCDF-error - invalid grid variable ID: '//trim(NCgrid(i)%name)
         call stop_tbnt(error_msg,ierr)
      end if
      ncStatus = NF90_INQUIRE_VARIABLE(fileUnit, NCgrid(i)%ncid, ndims = nVarDims)
      call nc_check(ncStatus, ios)
      if (ios/=0) then
         ierr = ierr + ios
         write(error_msg,'(a)')'NetCDF-error inquiring grid variable''s number of dimension: '//trim(NCgrid(i)%name)
         call stop_tbnt(error_msg,ierr)
      end if
      if (nVarDims<1.or.nVarDims>2) then
         ierr = 1
         write(error_msg,'(a)')'NetCDF-error - invalid number of dimensions for grid variable: '//trim(NCgrid(i)%name)
         call stop_tbnt(error_msg,ierr)
      end if
      ! check if 'units' attribute exists
      ncStatus = NF90_INQUIRE_ATTRIBUTE(fileUnit, NCgrid(i)%ncid, 'units')
      call nc_check(ncStatus, ios)
      if (ios/=0.and.ncStatus==NF90_ENOTATT) then
         NCgrid(i)%unit = 'dimensionless'
      elseif (ios/=0.and.ncStatus/=NF90_ENOTATT) then
         ierr = ierr + ios
         write(error_msg,'(a)')'NetCDF-error inquiring variable unit: '//trim(NCgrid(i)%name)
         call stop_tbnt(error_msg,ierr)
      else
         ncStatus = NF90_GET_ATT(fileUnit, NCgrid(i)%ncid, 'units', unitStr)
         call nc_check(ncStatus, ios)
         if (ios/=0) then
            ierr = ierr + ios
            write(error_msg,'(a)')'NetCDF-error reading grid variable unit: '//trim(NCgrid(i)%name)
            call stop_tbnt(error_msg,ierr)
         end if
         NCgrid(i)%unit = trim(unitStr)
      end if
      
      ! note that information gathered here refer to input file. dimension IDs stored to NCgrid do not match
      ! those of the output files and will be overwritten when output files are created
      allocate ( NCgrid(i)%size(nVarDims), NCgrid(i)%dimid(nVarDims), dimIDs(nVarDims) )
      NCgrid(i)%size = 0
      dimIDs = 0
      ncStatus = NF90_INQUIRE_VARIABLE(fileUnit, NCgrid(i)%ncid, dimids = dimIDs)
      call nc_check(ncStatus, ios)
      if (ios/=0) then
         ierr = ierr + ios
         write(error_msg,'(a,i3)')'NetCDF-error inquiring grid variable''s dimension IDs: '//trim(NCgrid(i)%name)
         call stop_tbnt(error_msg,ierr)
      end if
      do j = 1,nVarDims
         if (dimIDs(j)<=0) then
            ierr = 1
            write(error_msg,'(a)')'NetCDF-error - invalid grid variable''s dimension ID: '//trim(NCgrid(i)%name)
            call stop_tbnt(error_msg,ierr)
         end if
         found = .false.
         do k = 1,3
            if (dimIDs(j)==NCdims(k)%ncid) then
               NCgrid(i)%size(j) = NCdims(k)%size
               NCgrid(i)%dimid(j) = NCdims(k)%ncid
               found = .true.
               exit
            end if
         end do
         if (.not.found) then
            ierr = 1
            write(error_msg,'(a)')'NetCDF-error - grid variable''s dimension not found: '//trim(NCgrid(i)%name)
            call stop_tbnt(error_msg,ierr)
         end if
      end do
      deallocate ( dimIDs )
      if (nVarDims==1) then
         allocate ( NCgrid(i)%data1D(NCgrid(i)%size(1)) )
         ncStatus = NF90_GET_VAR(fileUnit, NCgrid(i)%ncid, NCgrid(i)%data1D)
      else
         allocate ( NCgrid(i)%data2D(NCgrid(i)%size(1), NCgrid(i)%size(2)) )
         ncStatus = NF90_GET_VAR(fileUnit, NCgrid(i)%ncid, NCgrid(i)%data2D)
      end if
      call nc_check(ncStatus, ios)
      if (ios/=0) then
         ierr = ierr + ios
         write(error_msg,'(a)')'NetCDF-error reading grid variable data: '//trim(NCgrid(i)%name)
         call stop_tbnt(error_msg,ierr)
      end if
   end do
   
   ! time information
   allocate ( NCgrid(4)%size(1), NCgrid(4)%dimid(1) )
   NCgrid(4)%name = trim(NCdims(4)%name)
   NCgrid(4)%size(1) = NCdims(4)%size
   NCgrid(4)%dimid(1) = timeID
   ncStatus = NF90_INQ_VARID(fileUnit, trim(NCgrid(4)%name), NCgrid(4)%ncid)
   call nc_check(ncStatus, ios)
   if (ios/=0) then
      ierr = ierr + ios
      write(error_msg,'(a)')'NetCDF-error inquiring time variable ID: '//trim(NCgrid(4)%name)
      call stop_tbnt(error_msg,ierr)
   end if
   if (NCgrid(4)%ncid<=0) then
      ierr = 1
      write(error_msg,'(a)')'NetCDF-error - invalid grid variable ID: '//trim(NCgrid(4)%name)
      call stop_tbnt(error_msg,ierr)
   end if
   ncStatus = NF90_INQUIRE_ATTRIBUTE(fileUnit, NCgrid(4)%ncid, 'units')
   call nc_check(ncStatus, ios)
   if (ios/=0.and.ncStatus/=NF90_ENOTATT) then
      ierr = ierr + ios
      write(error_msg,'(a)')'NetCDF-error inquiring time unit: '//trim(NCgrid(4)%name)
      call stop_tbnt(error_msg,ierr)
   end if
   ncStatus = NF90_GET_ATT(fileUnit, NCgrid(i)%ncid, 'units', unitStr)
   call nc_check(ncStatus, ios)
   if (ios/=0) then
      ierr = ierr + ios
      write(error_msg,'(a)')'NetCDF-error reading grid variable unit: '//trim(NCgrid(i)%name)
      call stop_tbnt(error_msg,ierr)
   end if
   NCgrid(i)%unit = trim(unitStr)
   i = index(unitStr,'since')
   if (i==0) then
      ierr = 1
      write(error_msg,'(a)')'Error reading reference year from time unit.'
      call stop_tbnt(error_msg,ierr)
   end if
   ! determine time unit (i.e., seconds or minutes etc.)
   if (index(unitStr(1:i),'second')>0) then
      outTimeFac = 86400.e0
   elseif (index(unitStr(1:i),'minute')>0) then
      outTimeFac = 1440.e0
   elseif (index(unitStr(1:i),'hour')>0) then
      outTimeFac = 24.e0
   elseif (index(unitStr(1:i),'day')>0) then
      outTimeFac = 1.e0
   else
      ierr = 1
      write(error_msg,'(a)')'Error reading time unit from time unit attribute. No valid unit.'
      call stop_tbnt(error_msg,ierr)
   end if
   ! read initial time value
   allocate ( outTime(NCgrid(4)%size(1)) )
   ncStatus = NF90_GET_VAR(fileUnit, NCgrid(4)%ncid, outTime)
   call nc_check(ncStatus, ios)
   if (ios/=0) then
      ierr = ierr + ios
      write(error_msg,'(a)')'NetCDF-error reading starting time: '//trim(NCgrid(i)%name)
      call stop_tbnt(error_msg,ierr)
   end if
   outTimeStart = outTime(TBNTiStart)
   deallocate ( outTime )
   
   end subroutine get_nc_dimensions
#endif
!=======================================================================
!
! !INTERFACE:
subroutine init_tbnt_fluxes_vars(ierr)
#define SUBROUTINE_NAME 'init_tbnt_fluxes_vars'
!
! !DESCRIPTION:
! get fluxes and variables related to traced elements
!
   implicit none
!
! !OUTPUT PARAMETERS:
   integer, intent(inout) :: ierr
!
! !LOCAL VARIABLES:
!
   integer                :: ios, iElem, ind, ind1, ind2, nEntry, iiVar
   integer                :: iFlux, iRivVar
   real(dp)               :: facVal
   character(len=nameLen) :: elem, facStr, readStr, dummyVar, auxVar
   character(len=lineLen) :: str
   logical                :: useFlux, useVar, lExist

   type(TBNTfluxProps)                              :: flx
   type(TBNTfluxProps)  , dimension(:), allocatable :: flxList2D, flxList3D
   type(TBNTmodVarProps), dimension(:), allocatable :: varList2D, varList3D
!
! !FOR NAMELISTS
!
   character(len=nameLen) :: tracedElement

   namelist /element_nml/ tracedElement
!
!-----------------------------------------------------------------------
#include "call-trace.inc"

   ios = 0
   open(input_unit, file=trim(fluxes_vars_file), &
                    form='formatted', action='read', iostat=ios)
   if (ios/=0) then
      ierr = ierr + ios
      write(error_msg,'(a)')'Error opening '//trim(fluxes_vars_file)//'.'
      call stop_tbnt(error_msg,ierr)
   end if

   ! initialize temporary flux list
   allocate ( flxList2D(TBNTnMaxFluxes), flxList3D(TBNTnMaxFluxes) )
   ! main flux
   flxList2D%name = '======'
   flxList2D%type = 0
   flxList2D%i_comp = 0
   flxList2D%iTBNT = 0
   flxList2D%varIn%name = '======'
   flxList2D%varIn%type = 0
   flxList2D%varIn%iTBNT = 0
   flxList2D%varIn%auxVar%name = '======'
   flxList2D%varIn%auxVar%type = 0
   flxList2D%varIn%auxVar%iTBNT = 0
   flxList2D%varIn%auxVar%fac = 1.e0
   flxList2D%varOut%name = '======'
   flxList2D%varOut%type = 0
   flxList2D%varOut%iTBNT = 0
   flxList2D%varOut%auxVar%name = '======'
   flxList2D%varOut%auxVar%type = 0
   flxList2D%varOut%auxVar%iTBNT = 0
   flxList2D%varOut%auxVar%fac = 1.e0
   ! auxiliary flux
   flxList2D%auxFlux%name = '======='
   flxList2D%auxFlux%type = 0
   flxList2D%auxFlux%fac  = 1.e0
   
   flxList3D = flxList2D
   
   n_total_flx = 0
   TBNTnFluxes2D = 0
   TBNTnFluxes3D = 0
   
   includeAtmos = .false.
   
   do iElem = 1,TBNTnTracedElements
      read(tbnt_settings_unit,nml=element_nml)
      ios = 0
      ! read flux entry
      do while (ios==0)
         read(input_unit,'(a)',iostat=ios)str ! read line
         if (index(str(1:10),'!+')>0) cycle
         ind = index(str, ';')
         elem = trim(adjustl(str(1:ind-1)))
         useFlux = .false.
         call string_check(elem, tracedElement, useFlux)
         if (.not.useFlux) cycle
         str = str(ind+1:len_trim(str))
         ind = index(str, ';')
         flx%name = trim(adjustl(str(1:ind-1)))
         ! check if current flux is already in the list
         useFlux = .false.
         do iFlux = 1,TBNTnFluxes3D
            call string_check(trim(flx%name),trim(flxList3D(iFlux)%name),useFlux)
            if (useFlux) exit
         end do
         if (useFlux) continue
         do iFlux = 1,TBNTnFluxes2D
            call string_check(trim(flx%name),trim(flxList2D(iFlux)%name),useFlux)
            if (useFlux) exit
         end do
         if (useFlux) continue
         
         ! initialize temporary flux
         ! main flux
         flx%type = 0
         flx%i_comp = 0
         flx%iTBNT = 0
         flx%varIn%name = '======'
         flx%varIn%type = 0
         flx%varIn%iTBNT = 0
         flx%varIn%auxVar%name = '======'
         flx%varIn%auxVar%type = 0
         flx%varIn%auxVar%iTBNT = 0
         flx%varIn%auxVar%fac = 0.e0
         flx%varOut%name = '======'
         flx%varOut%type = 0
         flx%varOut%iTBNT = 0
         flx%varOut%auxVar%name = '======'
         flx%varOut%auxVar%type = 0
         flx%varOut%auxVar%iTBNT = 0
         flx%varOut%auxVar%fac = 0.e0
         ! auxiliary flux
         flx%auxFlux%name = '======='
         flx%auxFlux%type = 0
         flx%auxFlux%fac  = 1.e0
         
         n_total_flx = n_total_flx + 1
         ! read flux information
         nEntry = 0
         str = str(ind+1:len_trim(str))
         ind = index(str, ';')
         do while (ind > 0)
            nEntry = nEntry + 1
            readStr = adjustl(str(1:ind-1))
            select case (nEntry)
               case (1)
                  read(readStr,*) flx%type
               case (2)
                  read(readStr,*) flx%i_comp
               case (3)
                  flx%varIn%name = trim(readStr)
               case (4)
                  read(readStr,*) flx%varIn%type
               case (5)
                  flx%varOut%name = trim(readStr)
               case (6)
                  read(readStr,*) flx%varOut%type
               case (7)
                  flx%auxFlux%name = trim(readStr)
               case (8)
                  read(readStr,*) flx%auxFlux%type
               case default
                  ierr = 1
                  write(error_msg,'(a)')'Error reading model_fluxes.txt. Check entry for flux: '//trim(flx%name)
                  call stop_tbnt(error_msg,ierr)
            end select
            str = str(ind+1:len_trim(str))
            ind = index(str, ';')
         end do
         if (((flx%type==TBNTder2DFlux.or.flx%type==TBNTder3DFlux).and.                           &
              (flx%auxFlux%type==TBNTder2DFlux.or.flx%auxFlux%type==TBNTder3DFlux)).or.           &
             ((flx%type==TBNTadvecFlux.or.flx%auxFlux%type==TBNTadvecFlux.or.                     &
               flx%type==TBNTdiffuFlux.or.flx%auxFlux%type==TBNTdiffuFlux).and.flx%i_comp<=0).or. &
             ((flx%type/=TBNTadvecFlux.and.flx%auxFlux%type/=TBNTadvecFlux.and.                   &
               flx%type/=TBNTdiffuFlux.and.flx%auxFlux%type/=TBNTdiffuFlux).and.flx%i_comp>0)) then
            ierr = 1
            write(error_msg,'(a)')'Invalid flux entry for flux: '//trim(flx%name)
            call stop_tbnt(error_msg,ierr)
         end if
         readStr = adjustl(str(1:len_trim(str)))
         if (nEntry==5) then
            read(readStr,*) flx%varOut%type
         elseif (nEntry==8) then
            facStr = trim(readStr)
            call interpret_mathString(facStr,flx%auxFlux%fac)
         else
            ierr = 1
            write(error_msg,'(a)')'Error reading model_fluxes.txt.'
            call stop_tbnt(error_msg,ierr)
         end if
         if (flx%type==TBNTder2DFlux.or.flx%type==TBNTsed2DFlux) then
            ! add current flux to 2D flux list
            TBNTnFluxes2D = TBNTnFluxes2D + 1
            flx%iTBNT = TBNTnFluxes2D
            flxList2D(TBNTnFluxes2D) = flx
         else
            ! add current flux to final flux list
            TBNTnFluxes3D = TBNTnFluxes3D + 1
            flx%iTBNT = TBNTnFluxes3D
            flxList3D(TBNTnFluxes3D) = flx
         end if
         if (.not.includeAtmos) then ! check if fluxes at atmosphere-ocean interface are involved
            if (flx%type==TBNTa2s3DFlux.or.flx%auxFlux%type==TBNTa2s3DFlux.or. &
                flx%type==TBNTa2s2DFlux.or.flx%auxFlux%type==TBNTa2s2DFlux.or. &
                flx%type==TBNTatm3DFlux.or.flx%auxFlux%type==TBNTatm3DFlux.or. &
                flx%type==TBNTatm2DFlux.or.flx%auxFlux%type==TBNTatm2DFlux) includeAtmos = .true.
         end if
      end do
      ! check if end of file is reached or a read error occured
      if (ios<0) then 
         if (iElem<TBNTnTracedElements) rewind(input_unit)
      else
         ierr = ierr + ios
         write(error_msg,'(a)')'Error reading model_fluxes.txt.'
         call stop_tbnt(error_msg,ierr)
      end if
   end do
   
   close(input_unit)

   ! copy flux list to final flux collection: first 3D fluxes, second 2D fluxes
   allocate( TBNTflux(n_total_flx) )
   TBNTflux(1:TBNTnFluxes3D) = flxList3D(1:TBNTnFluxes3D)
   TBNTflux(TBNTnFluxes3D+1:n_total_flx) = flxList2D(1:TBNTnFluxes2D)
   deallocate( flxList2D, flxList3D )
   
#ifndef TBNTonly_bulk_bud

   ! check if number of fractions > 0, if atmospheric deposition is removed
   if (.not.includeAtmos) then
      TBNTnFractions = TBNTnFractions - TBNTnAtmosFractions
      TBNTnAtmosFractions = 0
      if (TBNTnFractions == 0) then
         ierr = ierr + 1
         write(error_msg,'(a)')'Atmospheric deposition source turned off. No sources left.'
         call stop_tbnt(error_msg,ierr)
      end if      
   end if

   ! write flux list to log file
   write(tbnt_log_unit,'(a)'     ) 'LIST OF FLUXES'
   write(tbnt_log_unit,'(a,i3)'  ) 'Total number of involved fluxes: ',n_total_flx
   write(tbnt_log_unit,'(a,i3)'  ) 'Number of involved 3D fluxes   : ',TBNTnFluxes3D
   write(tbnt_log_unit,'(a,i3)'  ) 'Number of involved 2D fluxes   : ',TBNTnFluxes2D
   write(tbnt_log_unit,'(a)'     ) ''
   write(tbnt_log_unit,'(4(a25))') 'Flux name', 'Flux type', 'Input variable', 'Output variable'
   do iFlux = 1,n_total_flx
      write(tbnt_log_unit,'(a25,i25,2a25)') trim(TBNTflux(iFlux)%name), TBNTflux(iFlux)%type, &
                                            trim(TBNTflux(iFlux)%varIn%name), trim(TBNTflux(iFlux)%varOut%name)
   end do
   call breakLine(100)
#endif
      
   ! create list of variables
   allocate ( varList2D(2*n_total_flx), varList3D(2*n_total_flx) )
   varList2D%name = '======'
   varList2D%type = 0
   varList2D%iTBNT = 0
   varList2D%auxVar%name = '======'
   varList2D%auxVar%type = 0
   varList2D%auxVar%iTBNT = 0
   varList2D%auxVar%fac = 0.e0
   
   varList3D = varList2D
   
   n_total_var = 0
   n_progn_var = 0
   n_diagn_var = 0
   n_dummy_var = 0
   TBNTnVars3D = 0
   TBNTnVars2D = 0
   
   do iFlux = 1,n_total_flx
      ! check if input variable is already in the list
      useVar = .false.
      do iVar = 1,TBNTnVars2D
         call string_check(TBNTflux(iFlux)%varIn%name, varList2D(iVar)%name, useVar)
         if (useVar) exit
      end do
      if (.not.useVar) then
         do iVar = 1,TBNTnVars3D
            call string_check(TBNTflux(iFlux)%varIn%name, varList3D(iVar)%name, useVar)
            if (useVar) exit
         end do
      end if
      ! add input variable to list
      if (.not.useVar) then
         if (TBNTflux(iFlux)%varIn%type==TBNTpro2dVar.or. &
             TBNTflux(iFlux)%varIn%type==TBNTdia2dVar.or. &
             TBNTflux(iFlux)%varIn%type==TBNTdum2dVar) then
            TBNTnVars2D = TBNTnVars2D + 1
            varList2D(TBNTnVars2D) = TBNTflux(iFlux)%varIn
            varList2D(TBNTnVars2D)%iTBNT = TBNTnVars2D
            if (TBNTflux(iFlux)%varIn%type==TBNTpro2dVar) then
               n_progn_var = n_progn_var + 1
            elseif (TBNTflux(iFlux)%varIn%type==TBNTdia2dVar) then
               n_diagn_var = n_diagn_var + 1
            elseif (TBNTflux(iFlux)%varIn%type==TBNTdum2dVar) then
               n_dummy_var = n_dummy_var + 1
            end if
         else
            TBNTnVars3D = TBNTnVars3D + 1
            varList3D(TBNTnVars3D) = TBNTflux(iFlux)%varIn
            varList3D(TBNTnVars3D)%iTBNT = TBNTnVars3D
            if (TBNTflux(iFlux)%varIn%type==TBNTpro3dVar) then
               n_progn_var = n_progn_var + 1
            elseif (TBNTflux(iFlux)%varIn%type==TBNTdia3dVar) then
               n_diagn_var = n_diagn_var + 1
            elseif (TBNTflux(iFlux)%varIn%type==TBNTdum3dVar) then
               n_dummy_var = n_dummy_var + 1
            end if
         end if
      end if
      ! check if output variable == input variable
      call string_check(TBNTflux(iFlux)%varOut%name, TBNTflux(iFlux)%varIn%name, useVar)
      if (useVar) then
         TBNTflux(iFlux)%varOut = TBNTflux(iFlux)%varIn
         cycle
      end if
      ! check if output variable is already in the list
      do iVar = 1,TBNTnVars2D
         call string_check(TBNTflux(iFlux)%varOut%name, varList2D(iVar)%name, useVar)
         if (useVar) exit
      end do
      if (useVar) cycle
      do iVar = 1,TBNTnVars3D
         call string_check(TBNTflux(iFlux)%varOut%name, varList3D(iVar)%name, useVar)
         if (useVar) exit
      end do
      if (useVar) cycle
      ! add output variable to list
      if (TBNTflux(iFlux)%varOut%type==TBNTpro2dVar.or. &
          TBNTflux(iFlux)%varOut%type==TBNTdia2dVar.or. &
          TBNTflux(iFlux)%varOut%type==TBNTdum2dVar) then
         TBNTnVars2D = TBNTnVars2D + 1
         varList2D(TBNTnVars2D) = TBNTflux(iFlux)%varOut
         varList2D(TBNTnVars2D)%iTBNT = TBNTnVars2D
         if (TBNTflux(iFlux)%varOut%type==TBNTpro2dVar) then
            n_progn_var = n_progn_var + 1
         elseif (TBNTflux(iFlux)%varOut%type==TBNTdia2dVar) then
            n_diagn_var = n_diagn_var + 1
         elseif (TBNTflux(iFlux)%varOut%type==TBNTdum2dVar) then
            n_dummy_var = n_dummy_var + 1
         end if
      else
         TBNTnVars3D = TBNTnVars3D + 1
         varList3D(TBNTnVars3D) = TBNTflux(iFlux)%varOut
         varList3D(TBNTnVars3D)%iTBNT = TBNTnVars3D
         if (TBNTflux(iFlux)%varOut%type==TBNTpro3dVar) then
            n_progn_var = n_progn_var + 1
         elseif (TBNTflux(iFlux)%varOut%type==TBNTdia3dVar) then
            n_diagn_var = n_diagn_var + 1
         elseif (TBNTflux(iFlux)%varOut%type==TBNTdum3dVar) then
            n_dummy_var = n_dummy_var + 1
         end if
      end if
   end do

#ifdef TBNTnoVar2D
   if (TBNTnVars2D>0) then
      ! variable is 2D despite CPP switch TBNTnoVar2D
      ierr = ierr + 1
      write(error_msg,'(a)')'Number of 2D variables greater zero despite CPP switch TBNTnoVar2D'// &
                            new_line('a')//'  >  Check your setup file: '//trim(tbnt_settings_file)
      call stop_tbnt(error_msg,ierr)
   end if
#endif
   
   n_total_var = TBNTnVars2D + TBNTnVars3D
   
   ! get auxiliary information about dummy variables
   if (n_dummy_var>0) then
      ! check existence of dummy variable file
      lExist = .false.
      inquire(file = trim(model_set_dir)//trim(dummy_vars_file), exist = lExist)
      if (.not.lExist) then
         ierr = 1
         write(error_msg,'(2a)')'File does not exist: ',trim(model_set_dir)//trim(dummy_vars_file)
         call stop_tbnt(error_msg,ierr)
      end if
      dummy_vars_file = trim(model_set_dir)//trim(dummy_vars_file)
      ! open dummy variable file
      ios = 0
      open(input_unit, file=trim(dummy_vars_file), &
                       form='formatted', action='read', iostat=ios)
      if (ios/=0) then
         ierr = ierr + ios
         write(error_msg,'(a)')'Error opening '//trim(dummy_vars_file)//'.'
         call stop_tbnt(error_msg,ierr)
      end if
      ! read dummy variable information
      nEntry = 0
      do while (ios==0)
         read(input_unit,'(a)',iostat=ios) str ! read line
         if (index(str(1:10),'!+')>0) cycle
         ind1 = index(str, ';')
         ind2 = index(str, ';', .true.)
         dummyVar = trim(adjustl(str(1:ind1-1)))
         auxVar = trim(adjustl(str(ind1+1:ind2-1)))
         facStr = trim(adjustl(str(ind2+1:len_trim(str))))
         call interpret_mathString(facStr,facVal)
         useVar = .false.
         ! check if dummy and auxiliary variable are identical
         call string_check(dummyVar, auxVar, useVar)
         if (useVar.or.facVal==0.) cycle ! if identical or factor = 0: dummy variable does not need to be created
         ! loop over 3D variable
         varLoop3D: do iVar = 1,TBNTnVars3D
            if (varList3D(iVar)%type/=TBNTdum3dVar) cycle
            call string_check(varList3D(iVar)%name, dummyVar, useVar)
            if (useVar) then
               varList3D(iVar)%auxVar%name = trim(auxVar)
               varList3D(iVar)%auxVar%fac = facVal
               do iiVar = 1,TBNTnVars3D
                  if (varList3D(iiVar)%type==TBNTdum3dVar) cycle
                  call string_check(varList3D(iiVar)%name, auxVar, useVar)
                  if (useVar) then
                     varList3D(iVar)%auxVar%iTBNT = varList3D(iiVar)%iTBNT
                     varList3D(iVar)%auxVar%type = varList3D(iiVar)%type
                     exit varLoop3D
                  end if
               end do
            end if
         end do varLoop3D
#ifndef TBNTnoVar2D
         if (.not.useVar) then
            ! loop over 2D variable
            varLoop2D: do iVar = 1,TBNTnVars2D
               if (varList2D(iVar)%type/=TBNTdum2dVar) cycle
               call string_check(varList2D(iVar)%name, dummyVar, useVar)
               if (useVar) then
                  varList2D(iVar)%auxVar%name = trim(auxVar)
                  varList2D(iVar)%auxVar%fac = facVal
                  do iiVar = 1,TBNTnVars2D
                     if (varList2D(iiVar)%type==TBNTdum2dVar) cycle
                     call string_check(varList2D(iiVar)%name, auxVar, useVar)
                     if (useVar) then
                        varList2D(iVar)%auxVar%iTBNT = varList2D(iiVar)%iTBNT
                        varList2D(iVar)%auxVar%type = varList2D(iiVar)%type
                        exit varLoop2D
                     end if
                  end do
               end if
            end do varLoop2D
         end if
#endif
      end do
   end if 
   
   ! copy variables to final variable collection: first 3D variables, second 2D variables
   allocate ( TBNTvar(n_total_var) )
   TBNTvar(1:TBNTnVars3D) = varList3D(1:TBNTnVars3D)
#ifndef TBNTnoVar2D
   TBNTvar(TBNTnVars3D+1:n_total_var) = varList2D(1:TBNTnVars2D)
#endif
   deallocate ( varList2D, varList3D )
   
   ! complete variable entries for all fluxes
   do iFlux = 1,n_total_flx
      if (TBNTflux(iFlux)%varIn%iTBNT==0) then
         do iVar = 1,n_total_var
            call string_check(TBNTflux(iFlux)%varIn%name, TBNTvar(iVar)%name, useVar)
            if (useVar) TBNTflux(iFlux)%varIn = TBNTvar(iVar)
         end do
      end if
      if (TBNTflux(iFlux)%varOut%iTBNT==0) then
         do iVar = 1,n_total_var
            call string_check(TBNTflux(iFlux)%varOut%name, TBNTvar(iVar)%name, useVar)
            if (useVar) TBNTflux(iFlux)%varOut = TBNTvar(iVar)
         end do
      end if
   end do
   
   ! get number of different river fluxes
   nRivFluxes = 0
   iRivVar = 0
   do iFlux = 1,n_total_flx
      if (TBNTflux(iFlux)%type/=TBNTrvdisFlux.and.TBNTflux(iFlux)%auxFlux%type/=TBNTrvdisFlux) cycle
      if (iRivVar==0) then
         iRivVar = TBNTflux(iFlux)%varOut%iTBNT
         nRivFluxes = 1
      elseif (TBNTflux(iFlux)%varOut%iTBNT==iRivVar) then
         nRivFluxes = nRivFluxes + 1
      end if
   end do
   
#ifndef TBNTonly_bulk_bud
   ! write variable list to log file
   write(tbnt_log_unit,'(a)'   )  'LIST OF VARIABLES'
   write(tbnt_log_unit,'(a,i3)')  'Total number of traced variables     : ',n_total_var
   write(tbnt_log_unit,'(a,i3)')  'Number of traced 3D variables        : ',TBNTnVars3D
   write(tbnt_log_unit,'(a,i3)')  'Number of traced 2D variables        : ',TBNTnVars2D
   write(tbnt_log_unit,'(a,i3)')  'Number of traced prognostic variables: ',n_progn_var
   write(tbnt_log_unit,'(a,i3)')  'Number of traced derived variables   : ',n_diagn_var
   write(tbnt_log_unit,'(a,i3)')  'Number of traced dummy variables     : ',n_dummy_var
   write(tbnt_log_unit,'(a)'   )  ''
   write(tbnt_log_unit,'(2(a25))')'Variable name', 'Variable type'
   do iVar = 1,n_total_var
      write(tbnt_log_unit,'(a25,i25)') trim(TBNTvar(iVar)%name), TBNTvar(iVar)%type
   end do
   call breakLine(100)
#endif

   end subroutine init_tbnt_fluxes_vars
!=======================================================================
#ifndef TBNTonly_bulk_bud
!
! !INTERFACE:
   subroutine init_tbnt_sources(ierr)
#define SUBROUTINE_NAME 'init_tbnt_sources'
!
! !DESCRIPTION:
! initialize sources
!
! !USES:

   implicit none
!
! !OUTPUT PARAMETERS:
   integer, intent(inout) :: ierr
!
! !LOCAL VARIABLES
   integer :: i, nInputs
!
!-----------------------------------------------------------------------
#include "call-trace.inc"

   ! write number of fractions and fraction definition to log file
   write(tbnt_log_unit,'(a)'     )'LIST OF SOURCES'
   write(tbnt_log_unit,'(a,i3)'  )'Total number of sources        : ', TBNTnFractions
   write(tbnt_log_unit,'(a,i3)'  )'Number of river sources        : ', TBNTnRiverFractions
   write(tbnt_log_unit,'(a,i3)'  )'Number of open boundary sources: ', TBNTnOpenBFractions
   write(tbnt_log_unit,'(a,i3)'  )'Number of atmospheric sources  : ', TBNTnAtmosFractions
   write(tbnt_log_unit,'(a)'     )''
   write(tbnt_log_unit,'(3(a25))')'Source group name', 'Source group index', '# of grid points'

   ! adding fraction for untraced substance
   useDummySrc = .true.
   TBNTnFractions = TBNTnFractions + 1

   ! initialise source list   
   allocate ( TBNTsource(TBNTnFractions) )
   TBNTsource%name = ''
   TBNTsource%type = 0
   TBNTsource%nSubsources = 0
   TBNTnDummyRivers = 0
   TBNTnDummyAtmos  = 0
   ! read river sources and add them to source list (always called because untraced rivers are added to dummy tracer)
   call init_tbnt_rivers(ierr)
   ! read open boundary sources and add them to source list
   if (TBNTnOpenBFractions > 0) call init_tbnt_openb(ierr)
   ! read atmospheric deposition sources and add them to source list
   if (includeAtmos) call init_tbnt_atmos(ierr)
   
   ! check if any dummy sources exist
   if (TBNTnDummyRivers+TBNTnDummyAtmos==0) then
      useDummySrc = .false.
      TBNTnFractions = TBNTnFractions - 1
   end if
   
   do iFrac = 1,TBNTnFractions
      nInputs = 0
      do i = 1,TBNTsource(iFrac)%nSubsources
         nInputs = nInputs + TBNTsource(iFrac)%subsource(i)%nInputs
      end do
      write(tbnt_log_unit,'(a25,2(i25))')trim(TBNTsource(iFrac)%name), iFrac, nInputs
   end do
   
   call breakLine(100)

   end subroutine init_tbnt_sources
#endif
!=======================================================================
#ifndef TBNTonly_bulk_bud
!
! !INTERFACE:
   subroutine init_tbnt_rivers(ierr)
#define SUBROUTINE_NAME 'init_tbnt_rivers'
!
! !DESCRIPTION:
! initialize river sources
!
! !USES:

   implicit none
!
! !OUTPUT PARAMETERS:
   integer, intent(inout) :: ierr
!
! !LOCAL VARIABLES:
!   
   integer                            :: i, ir, nModelRivers, nInputs, nMaxIn, flxType, auxFlxType, iRivFlux
   integer, dimension(:), allocatable :: nRivs
   logical                            :: isIncluded
   logical, dimension(:), allocatable :: isTracedRiv
   character(len=nameLen), dimension(:), allocatable :: riverName
   
   integer                              :: nSkip, ind
   integer, dimension(:,:), allocatable :: riverInputs
   character(len=lineLen)               :: inLine
   logical, dimension(:)  , allocatable :: isUsed
!
! !FOR NAMELISTS
!
   integer                            :: nTracedRivers  , rivGroupID
   character(len=nameLen)             :: tracedGroupName, tracedRivName

   namelist /river_groups_nml/ tracedGroupName, nTracedRivers
   namelist /river_fraction_nml/ tracedRivName , rivGroupID

!-----------------------------------------------------------------------
#include "call-trace.inc"

   ! read river list
   open(input_unit, file = trim(rivers_file))
   ! check header length
   call get_headLen(input_unit, rivers_file, nSkip, ierr)
   ! skip header
   call skipLines(input_unit,nSkip)
   ! read river information
   read(input_unit,*, iostat = ierr) nModelRivers, nMaxIn
   allocate ( riverInputs(nModelRivers,nMaxIn), riverName(nModelRivers) )
   do i = 1,nModelRivers
      read(input_unit, '(a)', iostat = ierr) inLine
      ind = index(inLine,'-')
      if (ind==0) ierr = 1      
      if (ierr/=0) exit
      read(inLine(1:ind-1), *) ir, riverInputs(ir,:)
      riverName(ir) = trim(adjustl(inLine(ind+1:len_trim(inLine))))
      call capitalise_string(riverName(ir))
   end do
   close(input_unit)
   if (ierr/=0) then
      write(error_msg,'(2a)')'Error reading: ',trim(rivers_file)
      call stop_tbnt(error_msg,ierr)
   end if
   
   allocate ( isUsed(nWetCells) )

   isUsed = .false.   
   allocate ( isTracedRiv(nModelRivers) )
   
   ! get river sources
   isTracedRiv = .false.
   if (TBNTnRiverFractions>0) then
      do iFrac = 1,TBNTnRiverFractions
         read(tbnt_settings_unit, nml=river_groups_nml)
         TBNTsource(iFrac)%name = trim(tracedGroupName)
         TBNTsource(iFrac)%type = TBNTriverSource
         TBNTsource(iFrac)%nSubsources = nTracedRivers
         allocate( TBNTsource(iFrac)%subsource(nTracedRivers) )
         TBNTsource(iFrac)%subsource%name = '==='
         TBNTsource(iFrac)%subsource%type = TBNTriverSource
         TBNTsource(iFrac)%subsource%index = 0
         TBNTsource(iFrac)%subsource%nInputs = 0
      end do
      
      allocate ( nRivs(TBNTnRiverFractions) )
      nRivs = 0
      do iSource = 1,sum(TBNTsource(1:TBNTnRiverFractions)%nSubsources)
         read(tbnt_settings_unit, nml=river_fraction_nml)
         if (rivGroupID<0 .or. rivGroupID>TBNTnRiverFractions) then
            ierr = ierr + 1
            write(error_msg,'(2(a))')'River source index out of bounds, check tbnt_set.nml: ',trim(tracedRivName)
            call stop_tbnt(error_msg,ierr)
         end if
         call capitalise_string(tracedRivName)
         isIncluded = .false.
         do ir = 1,nModelRivers
            if (index(riverName(ir), trim(tracedRivName)) > 0) isIncluded = .true.
            if (isIncluded) then
               isTracedRiv(ir) = .true.
               if (isUsed(riverInputs(ir,1))) then
                  TBNTsource(rivGroupID)%nSubsources = TBNTsource(rivGroupID)%nSubsources - 1
                  exit
               end if
               nRivs(rivGroupID) = nRivs(rivGroupID) + 1
               TBNTsource(rivGroupID)%subsource(nRivs(rivGroupID))%name = trim(tracedRivName)
               TBNTsource(rivGroupID)%subsource(nRivs(rivGroupID))%index = ir
               isUsed(riverInputs(ir,1)) = .true.
               nInputs = count(riverInputs(ir,:)>0)
               TBNTsource(rivGroupID)%subsource(nRivs(rivGroupID))%nInputs = nInputs
               allocate ( TBNTsource(rivGroupID)%subsource(nRivs(rivGroupID))%i_from(nInputs), &
                          TBNTsource(rivGroupID)%subsource(nRivs(rivGroupID))%i_to(nInputs)  , &
                          TBNTsource(rivGroupID)%subsource(nRivs(rivGroupID))%i_comp(nInputs), &
                          TBNTsource(rivGroupID)%subsource(nRivs(rivGroupID))%sign(nInputs)   )
               TBNTsource(rivGroupID)%subsource(nRivs(rivGroupID))%i_from = riverInputs(ir,1:nInputs)
               TBNTsource(rivGroupID)%subsource(nRivs(rivGroupID))%i_to   = riverInputs(ir,1:nInputs)
               TBNTsource(rivGroupID)%subsource(nRivs(rivGroupID))%i_comp = 0
               TBNTsource(rivGroupID)%subsource(nRivs(rivGroupID))%sign   = 0.e0
               exit
            end if
         end do
         if (.not.isIncluded) then
            ierr = ierr + 1
            write(error_msg,'(2(a))')'River is not included: ',trim(tracedRivName)
            call stop_tbnt(error_msg,ierr)
         end if
      end do
      do iFrac = 1,TBNTnRiverFractions
         if (TBNTsource(iFrac)%nSubsources - nRivs(iFrac)/=0) then
            ierr = ierr + 1
            write(error_msg,'(a,i3)')'Number of river subsources does not fit - river group index: ', iFrac
            call stop_tbnt(error_msg,ierr)
         end if
      end do
      deallocate ( nRivs )
   end if
   
   ! add untraced rivers to dummy source
   TBNTnDummyRivers = nModelRivers - count(isTracedRiv)
   
   iFrac = TBNTnFractions
   TBNTsource(iFrac)%name = 'UNTRACED'
   TBNTsource(iFrac)%type = TBNTdummySource
   TBNTsource(iFrac)%nSubsources = TBNTnDummyRivers
   
   if (TBNTnDummyRivers>0) then
      allocate ( TBNTsource(iFrac)%subsource(TBNTnDummyRivers) )
      TBNTsource(iFrac)%subsource%type = TBNTriverSource
      TBNTsource(iFrac)%subsource%nInputs = 0     
      iPut = 0
      do ir = 1,nModelRivers
         if (isTracedRiv(ir)) cycle
         if (isUsed(riverInputs(ir,1))) then
            TBNTsource(iFrac)%nSubsources = TBNTsource(iFrac)%nSubsources - 1
            TBNTnDummyRivers = TBNTnDummyRivers - 1
            cycle
         end if
         iPut = iPut + 1
         TBNTsource(iFrac)%subsource(iPut)%name = trim(riverName(ir))
         TBNTsource(iFrac)%subsource(iPut)%index = ir
         isUsed(riverInputs(ir,1)) = .true.
         nInputs = count(riverInputs(ir,:)>0)
         TBNTsource(iFrac)%subsource(iPut)%nInputs = nInputs
         allocate ( TBNTsource(iFrac)%subsource(iPut)%i_from(nInputs), TBNTsource(iFrac)%subsource(iPut)%i_to(nInputs), &
                    TBNTsource(iFrac)%subsource(iPut)%i_comp(nInputs), TBNTsource(iFrac)%subsource(iPut)%sign(nInputs) )
         TBNTsource(iFrac)%subsource(iPut)%i_from = riverInputs(ir,1:nInputs)
         TBNTsource(iFrac)%subsource(iPut)%i_to = riverInputs(ir,1:nInputs)
         TBNTsource(iFrac)%subsource(iPut)%i_comp = 0
         TBNTsource(iFrac)%subsource(iPut)%sign = 0.e0
      end do
   end if
   deallocate ( isTracedRiv, isUsed, riverInputs )
         
   ! get pointer for river discharge flux of different variables
   allocate ( iRvDisFlux(TBNTnVars3D,nRivFluxes) )
   iRvDisFlux = 0
   do iVar = 1,TBNTnVars3D
      iRivFlux = 0
      if (TBNTvar(iVar)%type==TBNTdum3dVar.and.TBNTvar(iVar)%auxVar%fac==0.e0) cycle
      do iFlux = 1,TBNTnFluxes3D
         flxType = TBNTflux(iFlux)%type
         auxFlxType = TBNTflux(iFlux)%auxFlux%type
         if ((flxType==TBNTrvdisFlux.or.auxFlxType==TBNTrvdisFlux).and.TBNTflux(iFlux)%varIn%iTBNT==iVar) then
            iRivFlux = iRivFlux + 1
            iRvDisFlux(iVar,iRivFlux) = iFlux
            if (iRivFlux==nRivFluxes) exit
         end if
      end do
      if (sum(iRvDisFlux(iVar,:))==0) then
         ierr = ierr + 1
         write(error_msg,'(2a)')'No river flux found: ',trim(TBNTvar(iVar)%name)
         call stop_tbnt(error_msg,ierr)
      end if
   end do
   
   deallocate ( riverName )

   end subroutine init_tbnt_rivers
#endif
!=======================================================================
#ifndef TBNTonly_bulk_bud
!
! !INTERFACE:
   subroutine init_tbnt_openb(ierr)
#define SUBROUTINE_NAME 'init_tbnt_openb'
!
! !DESCRIPTION:
! initialize open boundary sources
!
! !USES:

   implicit none
!
! !OUTPUT PARAMETERS:
   integer, intent(inout) :: ierr
!
! !LOCAL VARIABLES:
!
   integer                              :: i, nInputs
   integer                              :: nSkip, ios, nExclude, iCheck
   integer, dimension(:  ), allocatable :: iExcluded
   integer, dimension(:,:), allocatable :: sourceInfo
   logical                              :: lExist
   character(len=formatLen)             :: numStr
   character(len=fileLen)               :: filename
!
! !FOR NAMELISTS
!
   character(len=fileLen) :: openb_file

   namelist /openb_fraction_nml/ openb_file

!-----------------------------------------------------------------------
#include "call-trace.inc"

   ! read file with open boundary information
   read(tbnt_settings_unit, nml=openb_fraction_nml)
   filename = trim(tbnt_set_dir)//trim(openb_file)
   inquire(file=trim(filename),exist=lExist)
   if (.not.lExist) then
      ierr = ierr + 1
      write(error_msg,'(a)')'Open boundary input file does not exist: '//trim(filename)
      call stop_tbnt(error_msg,ierr)
   end if
   open(input_unit, file = trim(filename), action = 'read')
   ios = 0
   ! check header length
   call get_headLen(input_unit, filename, nSkip, ierr)
   ! skip header
   call skipLines(input_unit, nSkip)
   ! read open boundary information
   read(input_unit,*,iostat=ios) nInputs
   if (nInputs/=TBNTnOpenBFractions) then
      ierr = ierr + 1
      write(error_msg,'(a)')'Mismatch in number of open boundary sources: check tbnt_set.nml and '//trim(openb_file)//'.'
      call stop_tbnt(error_msg,ierr)
   end if
   iOff = TBNTnRiverFractions
   TBNTsource(iOff+1:iOff+TBNTnOpenBFractions)%nSubsources = 0
   do iFrac = iOff+1,iOff+TBNTnOpenBFractions
      read(input_unit,*,iostat=ios) iCheck, nInputs
      if (ios/=0) then
         ierr = ierr + ios
         write(error_msg,'(a)')'Error reading open boundary source information.'
         call stop_tbnt(error_msg,ierr)
      end if
      if (iCheck>TBNTnOpenBFractions.or.iCheck<0.or.iCheck/=iFrac-iOff) then
         ierr = ierr + 1
         write(error_msg,'(a)')'Invalid open boundary source identifier.'
         call stop_tbnt(error_msg,ierr)
      end if
      write(numStr,'(i1)') iCheck
      TBNTsource(iFrac)%name = 'OpenBoundary-'//trim(numStr)
      TBNTsource(iFrac)%type = TBNTopenbSource
      TBNTsource(iFrac)%nSubsources = 1
      allocate ( sourceInfo(nInputs,4), TBNTsource(iFrac)%subsource(1) )
      do i = 1,4
         read(input_unit,*,iostat=ios) sourceInfo(:,i)
         if (ios/=0) then
            ierr = ierr + ios
            write(error_msg,'(a)')'Error reading open boundary source information 2.'
            call stop_tbnt(error_msg,ierr)
         end if
      end do
      TBNTsource(iFrac)%subsource(1)%name    = trim(TBNTsource(iFrac)%name)
      TBNTsource(iFrac)%subsource(1)%type    = TBNTopenbSource
      TBNTsource(iFrac)%subsource(1)%index   = 0
      TBNTsource(iFrac)%subsource(1)%nInputs = nInputs
      allocate ( TBNTsource(iFrac)%subsource(1)%i_from(nInputs), TBNTsource(iFrac)%subsource(1)%i_to(nInputs), &
                 TBNTsource(iFrac)%subsource(1)%i_comp(nInputs), TBNTsource(iFrac)%subsource(1)%sign(nInputs) )
      TBNTsource(iFrac)%subsource(1)%i_from  = sourceInfo(:,1)
      TBNTsource(iFrac)%subsource(1)%i_to    = sourceInfo(:,2)
      TBNTsource(iFrac)%subsource(1)%i_comp  = sourceInfo(:,3)
      TBNTsource(iFrac)%subsource(1)%sign    = real(sourceInfo(:,4),dp)
      do i = 1,nInputs
         if (sourceInfo(i,1)==sourceInfo(i,2).or.sourceInfo(i,3)<=0.or. &
             sourceInfo(i,3)>nComponents.or.abs(sourceInfo(i,4))/=1.e0) then
            ierr = ierr + 1
            write(error_msg,'(a,5i6)')trim(TBNTsource(iFrac)%name)//' - invalid set for input index: ', i, sourceInfo(i,:)
            call stop_tbnt(error_msg,ierr)
         end if
      end do
      deallocate ( sourceInfo )
   end do
   ! read exclude list (cells to be excluded from TBNT calculation)
   read(input_unit,*,iostat=ios) nExclude
   if (ios/=0) then
      ierr = ierr + ios
      write(error_msg,'(a)')'Error reading exclude list information.'
      call stop_tbnt(error_msg,ierr)
   end if
   allocate ( iExcluded(nExclude) )
   read(input_unit,*,iostat=ios) iExcluded
   if (ios/=0) then
      ierr = ierr + ios
      write(error_msg,'(a)')'Error reading exclude list information.'
      call stop_tbnt(error_msg,ierr)
   end if
   ! 3D exlude list
   do i = 1,nExclude
      excludeCell3D(iExcluded(i)) = .true.
   end do
   deallocate ( iExcluded )
   ! 2D exclude list
   iPut = 0
   do i = 1,nWetCells
      if (.not.bottomCell(i)) cycle
      iPut = iPut + 1
      if (excludeCell3D(i)) excludeCell2D(iPut) = .true. 
   end do
   
   end subroutine init_tbnt_openb
#endif
!=======================================================================
#ifndef TBNTonly_bulk_bud
!
! !INTERFACE:
   subroutine init_tbnt_atmos(ierr)
#define SUBROUTINE_NAME 'init_tbnt_atmos'
!
! !DESCRIPTION:
! initialize atmospheric deposition sources
!
! !USES:

   implicit none
!
! !OUTPUT PARAMETERS:
   integer, intent(inout) :: ierr
!
! !LOCAL VARIABLES:
!
   integer                               :: i, nInputs
   integer                               :: ios, nSkip, iCheck
   integer, dimension(:), allocatable    :: sourceInfo, i_out
   logical                               :: lExist
   logical, dimension(:), allocatable    :: isDummy
   character(len=formatLen)              :: numStr
   character(len=fileLen)                :: filename
   type(TBNTsourceProps)                 :: dummySrc
!
! !FOR NAMELISTS
!
   character(len=fileLen) :: atmos_file

   namelist /atmos_fraction_nml/ atmos_file

!-----------------------------------------------------------------------
#include "call-trace.inc"

   ! mask for atmospheric dummy sources
   allocate ( isDummy(nWetCells) )
   where (.not.surfaceCell.or.excludeCell3D)
      isDummy = .false.
   else where
      isDummy = .true.
   end where
   ! get atmospheric sources
   if (TBNTnAtmosFractions>0) then
      ! read file with atmospheric sources information
      read(tbnt_settings_unit, nml=atmos_fraction_nml)
      filename = trim(tbnt_set_dir)//trim(atmos_file)
      inquire(file=trim(filename),exist=lExist)
      if (.not.lExist) then
         ierr = ierr + 1
         write(error_msg,'(a)')'Atmospheric sources file does not exist: '//trim(filename)
         call stop_tbnt(error_msg,ierr)
      end if
      open(input_unit, file = trim(filename), action = 'read')
      ! check header length
      call get_headLen(input_unit, filename, nSkip, ierr)
      ! skip header
      call skipLines(input_unit, nSkip)
      ! read atmospheric sources information
      read(input_unit,*,iostat=ios) nInputs
      if (nInputs/=TBNTnAtmosFractions) then
         ierr = ierr + 1
         write(error_msg,'(a)')'Mismatch in number of atmospheric sources: check tbnt_set.nml and '//trim(atmos_file)//'.'
         call stop_tbnt(error_msg,ierr)
      end if
      iOff = TBNTnRiverFractions + TBNTnOpenBFractions
      do iFrac = iOff+1,iOff+TBNTnAtmosFractions
         read(input_unit, *, iostat = ios) iCheck, nInputs
         if (ios/=0) then
            ierr = ierr + ios
            write(error_msg,'(a)')'Error reading atmospheric source information.'
            call stop_tbnt(error_msg,ierr)
         end if
         if (iCheck>TBNTnAtmosFractions.or.iCheck<0.or.iCheck/=iFrac-iOff) then
            ierr = ierr + 1
            write(error_msg,'(a)')'Invalid atmospheric source identifier.'
            call stop_tbnt(error_msg,ierr)
         end if
         write(numStr,'(i1)') iCheck
         TBNTsource(iFrac)%name = 'Atmosphere-'//trim(numStr)
         TBNTsource(iFrac)%type = TBNTatmosSource
         TBNTsource(iFrac)%nSubsources = 1
         allocate ( TBNTsource(iFrac)%subsource(1) )
         TBNTsource(iFrac)%subsource(1)%name    = trim(TBNTsource(iFrac)%name)
         TBNTsource(iFrac)%subsource(1)%type    = TBNTatmosSource
         TBNTsource(iFrac)%subsource(1)%index   = 0
         TBNTsource(iFrac)%subsource(1)%nInputs = nInputs
         allocate ( TBNTsource(iFrac)%subsource(1)%i_from(nInputs), TBNTsource(iFrac)%subsource(1)%i_to(nInputs), &
                    TBNTsource(iFrac)%subsource(1)%i_comp(nInputs), TBNTsource(iFrac)%subsource(1)%sign(nInputs), &
                    sourceInfo(nInputs) )
         TBNTsource(iFrac)%subsource(1)%i_comp = 0
         TBNTsource(iFrac)%subsource(1)%sign   = 0.e0
         read(input_unit, *, iostat = ios) sourceInfo
         if (ios/=0) then
            ierr = ierr + ios
            write(error_msg,'(a)')'Error reading atmospheric source information.'
            call stop_tbnt(error_msg,ierr)
         end if
         TBNTsource(iFrac)%subsource(1)%i_from = sourceInfo
         TBNTsource(iFrac)%subsource(1)%i_to   = sourceInfo
         write(numStr,'(i6.6)')i
         TBNTsource(iFrac)%subsource(1)%name = trim(TBNTsource(iFrac)%name)
         isDummy(sourceInfo) = .false.
         deallocate ( sourceInfo )
      end do
      close(input_unit)
   end if
   
   ! add untraced atmospheric deposition to dummy source
   if (count(isDummy)>0) then
      iOff = TBNTnDummyRivers
      TBNTnDummyAtmos = 1
      iFrac = TBNTnFractions
      dummySrc%name = 'UNTRACED'
      dummySrc%type = TBNTdummySource
      dummySrc%nSubsources = iOff + 1
      allocate ( dummySrc%subsource(iOff+1) )
      ! copy existing untraced river sources to dummy source
      if (iOff>0) dummySrc%subsource(1:iOff) = TBNTsource(iFrac)%subsource(1:iOff)
      write(numStr,'(i1)') iOff+1
      dummySrc%subsource(iOff+1)%name  = trim(dummySrc%name)//'_'//trim(numStr)
      dummySrc%subsource(iOff+1)%type  = TBNTatmosSource
      dummySrc%subsource(iOff+1)%index = 0
      allocate ( i_out(nWetCells) )
      nInputs = 0
      do i = 1,nWetCells
         if (.not.isDummy(i).or.excludeCell3D(i).or..not.surfaceCell(i)) cycle
         nInputs = nInputs + 1
         i_out(nInputs) = i
      end do
      dummySrc%subsource(iOff+1)%nInputs = nInputs
      allocate ( dummySrc%subsource(iOff+1)%i_from(nInputs), dummySrc%subsource(iOff+1)%i_to(nInputs), &
                 dummySrc%subsource(iOff+1)%i_comp(nInputs), dummySrc%subsource(iOff+1)%sign(nInputs) )
      dummySrc%subsource(iOff+1)%i_comp = 0
      dummySrc%subsource(iOff+1)%sign   = 0.e0
      dummySrc%subsource(iOff+1)%i_from = i_out(1:nInputs)
      dummySrc%subsource(iOff+1)%i_to   = i_out(1:nInputs)
      TBNTsource(iFrac) = dummySrc ! copy updated dummy source
      deallocate ( dummySrc%subsource, i_out )
   end if

   end subroutine init_tbnt_atmos
#endif
!=======================================================================
#ifndef TBNTonly_bulk_bud
!
! !INTERFACE:
   subroutine init_tbnt_linked_fluxes(ierr)
#define SUBROUTINE_NAME 'init_tbnt_linked_fluxes'
!
! !DESCRIPTION:
! read linked fluxes from file add fluxes to flux list
!
! !USES:

   implicit none
!
! !OUTPUT PARAMETERS:
   integer, intent(inout) :: ierr
!
! !LOCAL VARIABLES
   integer                :: ios, ind, flxType, i, nEntries
   character(len=5)       :: frmtStr
   character(len=lineLen) :: str
   character(len=nameLen) :: flxName, varName
   logical                :: found
   type(TBNTfluxProps), dimension(:), allocatable :: flxList2D, flxList3D
!
!-----------------------------------------------------------------------
#include "call-trace.inc"

   ios = 0

   open(input_unit, file = trim(tbnt_set_dir)//trim(linked_fluxes_file), action = 'read', status = 'old',iostat=ios)
   if (ios/=0) then
      ierr = ierr + ios
      ierr = ierr + 1
      write(error_msg,'(a)') 'Error opening linked fluxes file.'
      call stop_tbnt(error_msg,ierr)
   end if
   
   ! copy existing flux lists to temporary lists
   allocate ( flxList2D(TBNTnMaxFluxes), flxList3D(TBNTnMaxFluxes) )
   flxList3D(1:TBNTnFluxes3D) = TBNTflux(1:TBNTnFluxes3D)
   flxList2D(1:TBNTnFluxes2D) = TBNTflux(TBNTnFluxes3D+1:TBNTnFluxes3D+TBNTnFluxes2D)
   deallocate ( TBNTflux )
   
   nEntries = 3 ! entries per flux
   
   found = .false.
   
   TBNTnLinkedFluxes2D = 0
   TBNTnLinkedFluxes3D = 0
   read(input_unit,'(a)',iostat=ios)str ! read first line
   do while (ios==0)
      if (index(str(1:10),'!+')>0) then
         read(input_unit,'(a)',iostat=ios)str ! read next line
         cycle
      end if
      do i = 1,nEntries
         if (i<nEntries) ind = index(str, ';')
         select case (i)
            case (1)
               flxName = trim(adjustl(str(1:ind-1)))
            case (2)
               frmtStr ='(iXX)'
               write(frmtStr(3:4),'(i2)') ind-1
               read(str(1:ind-1),frmtStr) flxType
            case (3)
               varName = trim(adjustl(str))
         end select
         if (i<nEntries) str = str(ind+1:len_trim(str))
      end do
      if (flxType==TBNTsed2DFlux.or.flxType==TBNTder2DFlux) then
         TBNTnLinkedFluxes2D = TBNTnLinkedFluxes2D + 1
         flxList2D(TBNTnFluxes2D+TBNTnLinkedFluxes2D)%name = flxName
         flxList2D(TBNTnFluxes2D+TBNTnLinkedFluxes2D)%type = flxType
         flxList2D(TBNTnFluxes2D+TBNTnLinkedFluxes2D)%iTBNT = TBNTnFluxes2D+TBNTnLinkedFluxes2D
         do iVar = 1,n_total_var
            call string_check(varName,TBNTvar(iVar)%name,found)
            if (found) then
               flxList2D(TBNTnFluxes2D+TBNTnLinkedFluxes2D)%varIn  = TBNTvar(iVar)
               flxList2D(TBNTnFluxes2D+TBNTnLinkedFluxes2D)%varOut = TBNTvar(iVar)
               exit
            end if
         end do
      else
         TBNTnLinkedFluxes3D = TBNTnLinkedFluxes3D + 1
         flxList3D(TBNTnFluxes3D+TBNTnLinkedFluxes3D)%name = flxName
         flxList3D(TBNTnFluxes3D+TBNTnLinkedFluxes3D)%type = flxType
         flxList3D(TBNTnFluxes3D+TBNTnLinkedFluxes3D)%iTBNT = TBNTnFluxes3D+TBNTnLinkedFluxes3D
         do iVar = 1,n_total_var
            call string_check(varName,TBNTvar(iVar)%name,found)
            if (found) then
               flxList3D(TBNTnFluxes3D+TBNTnLinkedFluxes3D)%varIn  = TBNTvar(iVar)
               flxList3D(TBNTnFluxes3D+TBNTnLinkedFluxes3D)%varOut = TBNTvar(iVar)
               exit
            end if
         end do
      end if
      read(input_unit,'(a)',iostat=ios)str ! read next line
   end do
   
   if (ios>0) then
      ierr = ierr + ios
      write(error_msg,'(a)')'Error reading linked fluxes file.'
      call stop_tbnt(error_msg,ierr)
   end if
   
   allocate ( TBNTflux(TBNTnFluxes3D+TBNTnLinkedFluxes3D+TBNTnFluxes2D+TBNTnLinkedFluxes2D) )
   iOff = TBNTnFluxes3D+TBNTnLinkedFluxes3D
   TBNTflux(1:iOff) = flxList3D(1:iOff)
   TBNTflux(iOff+1:iOff+TBNTnFluxes2D+TBNTnLinkedFluxes2D) = flxList2D(1:TBNTnFluxes2D+TBNTnLinkedFluxes2D)
   deallocate (flxList3D, flxList2D )
   
   ! write variable list to log file
   write(tbnt_log_unit,'(a)'   )  'LIST OF LINKED FLUXES'
   write(tbnt_log_unit,'(a,i3)')  'Total number of linked fluxes: ',TBNTnLinkedFluxes3D+TBNTnLinkedFluxes2D
   write(tbnt_log_unit,'(a,i3)')  'Number of linked 3D fluxes   : ',TBNTnLinkedFluxes3D
   write(tbnt_log_unit,'(a,i3)')  'Number of linked 2D fluxes   : ',TBNTnLinkedFluxes2D
   write(tbnt_log_unit,'(a)'   )  ''
   write(tbnt_log_unit,'(4(a25))')'Flux name', 'Flux type', 'Linking var name', 'Linking var type'
   do iFlux = TBNTnFluxes3D+1,TBNTnFluxes3D+TBNTnLinkedFluxes3D
      write(tbnt_log_unit,'(2(a25,i25))') trim(TBNTflux(iFlux)%name)      , TBNTflux(iFlux)%type,      &
                                          trim(TBNTflux(iFlux)%varIn%name), TBNTflux(iFlux)%varIn%type
   end do
   iOff = TBNTnFluxes3D+TBNTnLinkedFluxes3D+TBNTnFluxes2D
   do iFlux = 1,TBNTnLinkedFluxes2D
      write(tbnt_log_unit,'(2(a25,i25))') trim(TBNTflux(iOff+iFlux)%name)      , TBNTflux(iOff+iFlux)%type,      &
                                          trim(TBNTflux(iOff+iFlux)%varIn%name), TBNTflux(iOff+iFlux)%varIn%type
   end do
   call breakLine(100)

   end subroutine init_tbnt_linked_fluxes
#endif
!=======================================================================
#ifdef TBNTonly_bulk_bud
!
! !INTERFACE:
   subroutine init_bulk_budget(ierr)
#define SUBROUTINE_NAME 'init_bulk_budget'
!
! !DESCRIPTION:
! reads nice map
!
! !USES:

   implicit none

!
! !OUTPUT PARAMETERS:
   integer, intent(inout) :: ierr
!
! !LOCAL VARIABLES:
   integer                                        :: ind, iPut, iiFlux
#ifdef TBNTconvert_3Dto1D
   integer                                        :: ii
#endif

   logical                                        :: isFound
   character(len=nameLen)                         :: vari, varo
   type(TBNTfluxProps), dimension(:), allocatable :: flux
!
! !FOR NAMELIST
   character(len=nameLen) :: balVarName
   integer                :: iBal, jBal, kBal

   namelist /bulk_bud_nml/ balVarName, iBal, jBal, kBal
   
!-----------------------------------------------------------------------
#include "call-trace.inc"

   write(tbnt_log_unit,'(a)'     )'initialize BULK BUDGET'
   call breakLine(100)
   read(tbnt_settings_unit, nml=bulk_bud_nml)
   isFound = .false.
   ! find variable to be balanced in TBNT variable selection
   do iVar = 1,n_total_var
      call string_check(TBNTvar(iVar)%name, balVarName, isFound)
      if (isFound) then
         BULKbudVar = TBNTvar(iVar)
         exit
      end if
   end do
   if (.not.isFound) then
      ierr = ierr + 1
      write(error_msg,'(a)') 'Budget variable not found: '//trim(balVarName)
      call stop_tbnt(error_msg,ierr)
   end if
   if ((BULKbudVar%type==TBNTdum3dVar.or.BULKbudVar%type==TBNTdum2dVar).and.BULKbudVar%auxVar%fac==0.e0) then
      ierr = ierr + 1
      write(error_msg,'(a)') 'Budget variable is dummy variable w/o auxiliary variable: '//trim(balVarName)// &
                             '. Budget calculation not possible.'
      call stop_tbnt(error_msg,ierr)
   end if
#ifdef TBNTnoVar2D
   if (BULKbudVar%type==TBNTpro2dVar.or.BULKbudVar%type==TBNTdia2dVar.or.BULKbudVar%type==TBNTdum2dVar) then
      ! variable is 2D despite CPP switch TBNTnoVar2D
      ierr = ierr + 1
      write(error_msg,'(a)')'Variable is 2D variable despite CPP switch TBNTnoVar2D: '//trim(BULKbudVar%name)// &
                            new_line('a')//'  >  Check your setup file: '//trim(tbnt_settings_file)
      call stop_tbnt(error_msg,ierr)
   end if
#endif
#ifdef TBNTconvert_3Dto1D
   if (k_index(iBal,jBal)<=0) then
      ierr = ierr + 1
      write(error_msg,'(a,2i4)')'Selected budget position (i,j) is on land: ',iBal, jBal
      call stop_tbnt(error_msg,ierr)
   end if
   if (BULKbudVar%type/=TBNTpro2dVar.and.BULKbudVar%type/=TBNTdia2dVar.and.BULKbudVar%type/=TBNTdum2dVar) then
      iBud = index3Dto1D(iBal,jBal,kBal)
      if (iBud==0) then
         ii = 1
         iBud = index3Dto1D(iBal,jBal,ii)
         do while (index3Dto1D(iBal,jBal,ii+1)>0)
            ii = ii + 1
            iBud = index3Dto1D(iBal,jBal,ii)
         end do
      end if
   else
      iBud = index2Dto1D(iBal,jBal)
   end if
#else
   ierr = ierr + 1
   write(error_msg,'(a)') 'Error using compiler flags:'//new_line('a')//                                       &
                          ' ==> Currently conversion from 3D to 1D fields is still required.'//new_line('a')// &
                          ' ==> Set TBNTconvert_3Dto1D.'
   call stop_tbnt(error_msg,ierr)
#endif
   
   allocate ( flux(n_total_flx) )
   iPut = 0
   ! check for 3D fluxes
   BULKnBudFluxes3D = 0
   do iFlux = 1,TBNTnFluxes3D
      vari = TBNTflux(iFlux)%varIn%name
      varo = TBNTflux(iFlux)%varOut%name
      ind = max(index(trim(vari), trim(balVarName)),index(trim(varo), trim(balVarName)))
      if (ind > 0) then
         iPut = iPut + 1
         flux(iPut) = TBNTflux(iFlux)
         BULKnBudFluxes3D = BULKnBudFluxes3D + 1
      end if
   end do
   ! check for 2D fluxes
   BULKnBudFluxes2D = 0
   do iFlux = 1,TBNTnFluxes2D
      iiFlux = TBNTnFluxes3D + iFlux
      vari = TBNTflux(iiFlux)%varIn%name
      varo = TBNTflux(iiFlux)%varOut%name
      ind = max(index(trim(vari), trim(balVarName)),index(trim(varo), trim(balVarName)))
      if (ind > 0) then
         iPut = iPut + 1
         flux(iPut) = TBNTflux(iiFlux)
         BULKnBudFluxes2D = BULKnBudFluxes2D + 1
      end if
   end do
   
   BULKnBudFluxes = BULKnBudFluxes3D + BULKnBudFluxes2D
   if (BULKnBudFluxes==0) then
      ierr = ierr + 1
      write(error_msg,'(a)') 'No budget fluxes found for: '//trim(balVarName)
      call stop_tbnt(error_msg,ierr)
   end if
   
   allocate ( BULKbudFlux(BULKnBudFluxes) )
   BULKbudFlux = flux(1:BULKnBudFluxes)
   
   deallocate ( flux )
   
   write(tbnt_log_unit,'(a)'     )'BALANCED VARIABLE           : '//trim(balVarName)
   write(tbnt_log_unit,'(a)'     )''
   write(tbnt_log_unit,'(a)'     )'LIST OF BALANCED FLUXES'
   write(tbnt_log_unit,'(a,i3)'  )'Number of involved fluxes   : ', BULKnBudFluxes
   write(tbnt_log_unit,'(a,i3)'  )'Number of involved 3D fluxes: ', BULKnBudFluxes3D
   write(tbnt_log_unit,'(a,i3)'  )'Number of involved 2D fluxes: ', BULKnBudFluxes2D
   write(tbnt_log_unit,'(a)'     )''
   write(tbnt_log_unit,'(4(a25))') 'Flux name', 'Flux type', 'Input variable', 'Output variable'
   do iFlux = 1,BULKnBudFluxes
      write(tbnt_log_unit,'(a25,i25,2a25)') trim(BULKbudFlux(iFlux)%name), BULKbudFlux(iFlux)%type, &
                                            trim(BULKbudFlux(iFlux)%varIn%name), trim(BULKbudFlux(iFlux)%varOut%name)
   end do
   call breakLine(100)
   
   end subroutine init_bulk_budget
#endif
!=======================================================================
#ifndef TBNTonly_bulk_bud
!
! !INTERFACE:
   subroutine init_tbnt_fraction_names
!
! !DESCRIPTION:
!  create fraction suffixes for output
!
! !USES:
   implicit none
!
! !OUTPUT PARAMETERS:
!
! !LOCAL VARIABLES:
!
  character(len=nameLen) :: riverFormat, openBFormat, atmosFormat, outName
!
!-----------------------------------------------------------------------
   
   if (TBNTnRiverFractions>0) then
      riverFormat = '("River-",iX.X)'
      if (TBNTnRiverFractions>99) then
         write(riverFormat(12:14),'(a)') '3.3'
      else
         if (TBNTnRiverFractions>9) then
            write(riverFormat(12:14),'(a)') '2.2'
         else
            write(riverFormat(12:14),'(a)') '1.1'
         end if
      end if
   end if
   if (TBNTnOpenBFractions>0) then
      openBFormat = '("OpenB-",iX.X)'
      if (TBNTnOpenBFractions>99) then
         write(openBFormat(12:14),'(a)') '3.3'
      else
         if (TBNTnOpenBFractions>9) then
            write(openBFormat(12:14),'(a)') '2.2'
         else
            write(openBFormat(12:14),'(a)') '1.1'
         end if
      end if
   end if
   if (TBNTnAtmosFractions>0) then
      atmosFormat = '("Atmos-",iX.X)'
      if (TBNTnAtmosFractions>99) then
         write(atmosFormat(12:14),'(a)') '3.3'
      else
         if (TBNTnatmosFractions>9) then
            write(atmosFormat(12:14),'(a)') '2.2'
         else
            write(atmosFormat(12:14),'(a)') '1.1'
         end if
      end if
   end if
      
   ! variables
   allocate ( TBNTfracName(TBNTnFractions) )
   do iFrac = 1,TBNTnFractions
      outName = ''
      if (TBNTsource(iFrac)%type==TBNTriverSource) then
         write(outName,riverFormat) iFrac
      elseif (TBNTsource(iFrac)%type==TBNTopenBSource) then
         write(outName,openBFormat) iFrac - TBNTnRiverFractions
      elseif (TBNTsource(iFrac)%type==TBNTatmosSource) then
         write(outName,atmosFormat) iFrac - TBNTnRiverFractions - TBNTnOpenBFractions
      elseif (TBNTsource(iFrac)%type==TBNTdummySource) then
         outName = 'Untraced'
      end if
      TBNTfracName(iFrac) = trim(outName)
   end do

   end subroutine init_tbnt_fraction_names
#endif
!=======================================================================
#ifndef TBNTonly_bulk_bud
!
! !INTERFACE:
   subroutine init_tbnt_target_areas(ierr)
#define SUBROUTINE_NAME 'init_tbnt_target_areas'
!
! !DESCRIPTION:
!  initialise target areas
!
! !USES:
   implicit none
!
! !OUTPUT PARAMETERS:
!
   integer, intent(inout) :: ierr
!
! !LOCAL VARIABLES:
!
   integer                            :: ios, nSkip, i1, i2, iArea, nCells, i, iCell
   integer, dimension(:), allocatable :: iCells
   character(len=lineLen)             :: inLine
   character(len=fileLen)             :: filename
   logical                            :: lExist, found
!
!-----------------------------------------------------------------------
#include "call-trace.inc"
   
   TBNTnTargetAreas = 0
   
   filename = trim(tbnt_set_dir)//trim(target_area_file)
   inquire(file = trim(filename), exist = lExist)
   if (.not.lExist) then
      ierr = 1
      write(error_msg,'(2a)')'File does not exist: ',trim(filename)
      call stop_tbnt(error_msg,ierr)
   end if
   
   ios = 0
   
   open(input_unit, file = trim(filename), action = 'read')
   ! check header length
   call get_headLen(input_unit, filename, nSkip, ierr)
   ! skip header
   call skipLines(input_unit,nSkip)
   ! create target area list
   read(input_unit,*,iostat=ios) TBNTnTargetAreas
   if (ios/=0) then
      ierr = ierr + ios
      write(error_msg,'(2a)')'Error reading number of target areas from file: ',trim(filename)
      call stop_tbnt(error_msg,ierr)
   end if
   if (TBNTnTargetAreas<=0) then
      ierr = 1
      write(error_msg,'(2a)')'Number of target areas less or equal zero. Check file: ',trim(filename)
      call stop_tbnt(error_msg,ierr)
   end if
   allocate ( TBNTtargetAreas(TBNTnTargetAreas) )
   do iArea = 1,TBNTnTargetAreas
      read(input_unit,'(a)',iostat=ios) inLine
      if (ios/=0) then
         ierr = ierr + ios
         write(error_msg,'(a,i5)')'Error reading header for target area #: ',iArea
         call stop_tbnt(error_msg,ierr)
      end if
      i2 = scan(inLine,'0123456789',.true.)
      i1 = index(inLine(1:i2),' ',.true.)
      if (i1==0.or.i2==0) then
         ierr = 1
         write(error_msg,'(a,i5)')'Error reading header for target area #: ',iArea
         call stop_tbnt(error_msg,ierr)
      end if
      TBNTtargetAreas(iArea)%name = trim(adjustl(inLine(1:i1)))
      read(inLine(i1+1:i2),*) nCells
      allocate ( TBNTtargetAreas(iArea)%isPart3D(nWetCells), TBNTtargetAreas(iArea)%isPart2D(nBotCells) )
      TBNTtargetAreas(iArea)%isPart3D = .false.
      TBNTtargetAreas(iArea)%isPart2D = .false.
      allocate ( iCells(nCells) )
      read(input_unit,*,iostat=ios) iCells
      if (ios/=0) then
         ierr = ierr + ios
         write(error_msg,'(a,i5)')'Error reading cells for target area #: ',iArea
         call stop_tbnt(error_msg,ierr)
      end if
      do iCell = 1,nCells
         if (excludeCell3D(iCells(iCell))) cycle
         TBNTtargetAreas(iArea)%isPart3D(iCells(iCell)) = .true.
         if (.not.bottomCell(iCells(iCell))) cycle
         found = .false.
         do i = 1,nBotCells
            if (bottom2pelag(i)==iCells(iCell)) then
               found = .true.
               exit
            end if
         end do
         if (.not.found) then
            ierr = 1
            write(error_msg,'(a,i5,a,i10)')'Bottom cell not found for target area #: ',iArea,', cell #: ',iCell
            call stop_tbnt(error_msg,ierr)
         end if
         if (.not.excludeCell2D(i)) TBNTtargetAreas(iArea)%isPart2D(i) = .true.
      end do
      deallocate ( iCells )
      
   end do
   close(input_unit)
   
   ! write target areas to log file
   write(tbnt_log_unit,'(a)'     )'LIST OF TARGET AREAS'
   write(tbnt_log_unit,'(a,i3)'  )'Total number of areas: ', TBNTnTargetAreas
   write(tbnt_log_unit,'(a)'     )''
   write(tbnt_log_unit,'(3(a25))')'Area name', '# of 3D grid points', '# of 2D grid points'
   do iArea = 1,TBNTnTargetAreas
      write(tbnt_log_unit,'(a25,(2(i25)))') trim(TBNTtargetAreas(iArea)%name), &
                                            count(TBNTtargetAreas(iArea)%isPart3D), &
                                            count(TBNTtargetAreas(iArea)%isPart2D)
   end do
   call breakLine(100)

   end subroutine init_tbnt_target_areas
#endif
!=======================================================================
#ifndef TBNTonly_bulk_bud
!
! !INTERFACE:
   subroutine init_tbnt_target_vars(ierr)
#define SUBROUTINE_NAME 'init_tbnt_target_vars'
!
! !DESCRIPTION:
!  initialise target variables
!
! !USES:
   implicit none
!
! !OUTPUT PARAMETERS:
!
   integer, intent(inout) :: ierr
!
! !LOCAL VARIABLES:
!
   integer                            :: ios, nSkip, nVars, iTargetVar, iiVar, iDay1, iDay2
   real(dp)                           :: nTargetDays
   character(len=lineLen)             :: inLine
   character(len=fileLen)             :: filename
   logical                            :: lExist, found
!
!-----------------------------------------------------------------------
#include "call-trace.inc"
   
   TBNTnTargetVars = 0
   
   filename = trim(tbnt_set_dir)//trim(target_vars_file)
   inquire(file = trim(filename), exist = lExist)
   if (.not.lExist) then
      ierr = 1
      write(error_msg,'(2a)')'File does not exist: ',trim(filename)
      call stop_tbnt(error_msg,ierr)
   end if
   
   ios = 0

   open(input_unit, file = trim(filename), action = 'read')
   ! check header length
   call get_headLen(input_unit, filename, nSkip, ierr)
   ! skip header
   call skipLines(input_unit,nSkip)
   ! create target area list
   read(input_unit,*,iostat=ios) TBNTnTargetVars
   if (ios/=0) then
      ierr = ierr + ios
      write(error_msg,'(2a)')'Error reading number of target variables from file: ',trim(filename)
      call stop_tbnt(error_msg,ierr)
   end if
   if (TBNTnTargetVars<=0) then
      ierr = 1
      write(error_msg,'(a)')'Number of target variables less or equal zero.'//new_line('a')// &
                            ' ==>  check file: '//trim(filename)
      call stop_tbnt(error_msg,ierr)
   end if
   allocate ( TBNTtargetVars(TBNTnTargetVars) )
   do iTargetVar = 1,TBNTnTargetVars
      read(input_unit,'(a)',iostat=ios) inLine
      if (ios/=0) then
         ierr = ierr + ios
         write(error_msg,'(a,i5)')'Error reading name for target variable #: ',iTargetVar
         call stop_tbnt(error_msg,ierr)
      end if
      TBNTtargetVars(iTargetVar)%name = trim(adjustl(inLine))
      read(input_unit,'(a)',iostat=ios) inLine
      if (ios/=0) then
         ierr = ierr + ios
         write(error_msg,'(a,i5)')'Error reading code for target variable #: ',iTargetVar
         call stop_tbnt(error_msg,ierr)
      end if
      TBNTtargetVars(iTargetVar)%code = trim(adjustl(inLine))
      read(input_unit,'(a)',iostat=ios) inLine
      if (ios/=0) then
         ierr = ierr + ios
         write(error_msg,'(a,i5)')'Error reading code for target variable #: ',iTargetVar
         call stop_tbnt(error_msg,ierr)
      end if
      read(inLine,*) nVars
      if (nVars<=0) then
         ierr = ierr + 1
         write(error_msg,'(a)')'Invalid number of contributing variables.'//new_line('a')//  &
                               ' ==> check file    : '//trim(filename)//new_line('a')//      &
                               ' ==> check variable: '//trim(TBNTtargetVars(iTargetVar)%name)
         call stop_tbnt(error_msg,ierr)
      end if
      TBNTtargetVars(iTargetVar)%nVars = nVars
      allocate ( TBNTtargetVars(iTargetVar)%vars(nVars) )
      do iiVar = 1,nVars
         read(input_unit,*,iostat=ios) inLine
         if (ios/=0) then
            ierr = ierr + ios
            write(error_msg,'(a,i5,a)')'Error reading contributing variable #: ',iiVar, new_line('a')// &
                                       ' ==> check file    : '//trim(filename)//new_line('a')//         &
                                       ' ==> check variable: '//trim(TBNTtargetVars(iTargetVar)%name)
                                        
            call stop_tbnt(error_msg,ierr)
         end if
         inLine = adjustl(inLine)
         found = .false.
         do iVar = 1,n_total_var
            call string_check(inLine,TBNTvar(iVar)%name,found)
            if (found) exit
         end do
         if (.not.found) then
            ierr = ierr + 1
            write(error_msg,'(a)')'Contributing variable not found: '//trim(inLine)//new_line('a')// &
                                  ' ==> check file    : '//trim(filename)//new_line('a')//           &
                                  ' ==> check variable: '//trim(TBNTtargetVars(iTargetVar)%name)
            call stop_tbnt(error_msg,ierr)
         end if
         TBNTtargetVars(iTargetVar)%vars(iiVar) = TBNTvar(iVar)
      end do
      read(input_unit,'(a)',iostat=ios) inLine
      if (ios/=0) then
         ierr = ierr + ios
         write(error_msg,'(a,i5)')'Error reading starting day for target variable #: ',iTargetVar
         call stop_tbnt(error_msg,ierr)
      end if
      read(inLine,*) iDay1
      read(input_unit,'(a)',iostat=ios) inLine
      if (ios/=0) then
         ierr = ierr + ios
         write(error_msg,'(a,i5)')'Error reading ending day for target variable #: ',iTargetVar
         call stop_tbnt(error_msg,ierr)
      end if
      read(inLine,*) iDay2
      if (iDay1<0.or.iDay2<0) then
         ierr = ierr + 1
         write(error_msg,'(a,i5,a)')'Invalid period definition for target variable #:',iTargetVar, &
                                    new_line('a')//' ==> check file: '//trim(filename)
         call stop_tbnt(error_msg,ierr)
      end if
      if (iDay1==0.or.iDay2==0) then
         ! select whole year
         iDay1 = 1
         iDay2 = TBNTyearDays
      end if
      ! create mask for target period
      allocate ( TBNTtargetVars(iTargetVar)%targetPeriod(TBNTnSteps) )
      TBNTtargetVars(iTargetVar)%targetPeriod = .false.
      do iStep = 1,TBNTnSteps
         if (iDay1<=iDay2) then
            if (TBNTstart+iStep*TBNTstep>iMinPerDay*iDay1.and.TBNTstart+(iStep-1)*TBNTstep<=iMinPerDay*iDay2) then
               TBNTtargetVars(iTargetVar)%targetPeriod(iStep) = .true.
            end if
         else
            if (TBNTstart+iStep*TBNTstep>iMinPerDay*iDay1.or.TBNTstart+(iStep-1)*TBNTstep<=iMinPerDay*iDay2) then
               TBNTtargetVars(iTargetVar)%targetPeriod(iStep) = .true.
            end if
         end if
      end do
   end do
   close(input_unit)
   
   ! allocate empty fields
   allocate ( TBNTtargetVals(TBNTnTargetAreas, TBNTnTargetVars, TBNTnFractions-1) )
   allocate ( TBNTtargetWeights(TBNTnTargetAreas, TBNTnTargetVars) )
   
   ! write target variables to log file
   write(tbnt_log_unit,'(a)') 'LIST OF TARGET VARIABLES'
   write(tbnt_log_unit,'(a,i3)') 'Total number of variables: ', TBNTnTargetVars
   write(tbnt_log_unit,'(a)')''
   write(tbnt_log_unit,'(a30,a20,2a25)')'Variable name', 'Variable code', '# of model variables', 'Target period [days]'
   do iVar = 1,TBNTnTargetVars
      nTargetDays = real(count(TBNTtargetVars(iVar)%targetPeriod),dp)*real(TBNTstep,dp)/rMinPerDay
      write(tbnt_log_unit,'(a30,a20,i25,f25.5)') trim(TBNTtargetVars(iVar)%name), trim(TBNTtargetVars(iVar)%code), &
                                                 TBNTtargetVars(iVar)%nVars, nTargetDays
   end do
   call breakLine(100)

   end subroutine init_tbnt_target_vars
#endif
!=======================================================================
#ifndef TBNTonly_bulk_bud
!
! !INTERFACE:
   subroutine init_tbnt_fraction_fields(ierr)
#define SUBROUTINE_NAME 'init_tbnt_fraction_fields'
!
! !DESCRIPTION:
!  initialise arrays for fractions fluxes, relative variable fractions, fraction pointers
!
   implicit none
!
! !OUTPUT PARAMETERS:
!
  integer, intent(inout) :: ierr
!
! !LOCAL VARIABLES:
!
!-----------------------------------------------------------------------
#include "call-trace.inc"

   ! create pointer for fraction fluxes
   ! pointer - sequence: all fluxes on variables of source 1, all fluxes on variables of source 2, ...
   allocate ( TBNTflux3Dpnt(TBNTnFluxes3D, TBNTnFractions) )
   allocate ( TBNTflux2Dpnt(TBNTnFluxes2D, TBNTnFractions) )
   do iFrac = 1,TBNTnFractions
      do iFlux = 1,TBNTnFluxes3D
         TBNTflux3Dpnt(iFlux,iFrac) = iFlux + (iFrac-1)*TBNTnFluxes3D
      end do
      do iFlux = 1,TBNTnFluxes2D
         TBNTflux2Dpnt(iFlux,iFrac) = iFlux + (iFrac-1)*TBNTnFluxes2D
      end do
   end do
   TBNTnFluxFractions3D = TBNTnFractions*TBNTnFluxes3D
   TBNTnFluxFractions2D = TBNTnFractions*TBNTnFluxes2D
   
   ! create pointer for linked fraction fluxes
   ! pointer - sequence: all fluxes on variables of source 1, all fluxes on variables of source 2, ...
   allocate ( TBNTlinkedFlux3Dpnt(TBNTnLinkedFluxes3D, TBNTnFractions) )
   allocate ( TBNTlinkedFlux2Dpnt(TBNTnLinkedFluxes2D, TBNTnFractions) )
   do iFrac = 1,TBNTnFractions
      do iFlux = 1,TBNTnLinkedFluxes3D
         TBNTlinkedFlux3Dpnt(iFlux,iFrac) = iFlux + (iFrac-1)*TBNTnLinkedFluxes3D
      end do
      do iFlux = 1,TBNTnLinkedFluxes2D
         TBNTlinkedFlux2Dpnt(iFlux,iFrac) = iFlux + (iFrac-1)*TBNTnLinkedFluxes2D
      end do
   end do
   TBNTnLinkedFluxFractions3D = TBNTnFractions*TBNTnLinkedFluxes3D
   TBNTnLinkedFluxFractions2D = TBNTnFractions*TBNTnLinkedFluxes2D
   
   ! initialise pointer for variables
   ! pointer - sequence: all variables of source 1, all variables of source 2, ...
   allocate ( TBNTvar3Dpnt(TBNTnVars3D,TBNTnFractions) )
#ifndef TBNTnoVar2D
   allocate ( TBNTvar2Dpnt(TBNTnVars2D,TBNTnFractions) )
#endif
   do iFrac = 1,TBNTnFractions
      do iVar = 1,TBNTnVars3D
         TBNTvar3Dpnt(iVar,iFrac) = iVar + (iFrac-1)*TBNTnVars3D
      end do
#ifndef TBNTnoVar2D
      do iVar = 1,TBNTnVars2D
         TBNTvar2Dpnt(iVar,iFrac) = iVar + (iFrac-1)*TBNTnVars2D
      end do
#endif
   end do
   
   ! create empty fields used during calculation
   allocate ( tempFlxFrac2D(nBotCells,TBNTnFluxFractions2D) )
   allocate ( tempFlxFrac3D(nWetCells,TBNTnFluxFractions3D) )
   tempFlxFrac2D = 0.e0
   tempFlxFrac3D = 0.e0
   
   if (TBNTlinkedFluxes) then
      allocate ( lnkFlxFrac2D(nBotCells,TBNTnLinkedFluxFractions2D) )
      allocate ( lnkFlxFrac3D(nWetCells,TBNTnLinkedFluxFractions3D) )
      lnkFlxFrac2D = 0.e0
      lnkFlxFrac3D = 0.e0
   end if
   
   ! create empty variable fractions
   TBNTnVarFractions3D = TBNTnVars3D*TBNTnFractions
   allocate ( TBNTvar3D(nWetCells,TBNTnVarFractions3D), TBNTrelVar3D(nWetCells,TBNTnVarFractions3D) )
   allocate ( varChange3D(nWetCells,TBNTnVarFractions3D) )
   TBNTvar3D = 0.e0
   TBNTrelVar3D = 0.e0
   varChange3D = 0.e0
#ifndef TBNTnoVar2D
   TBNTnVarFractions2D = TBNTnVars2D*TBNTnFractions
   allocate ( TBNTvar2D(nBotCells,TBNTnVarFractions2D), TBNTrelVar2D(nBotCells,TBNTnVarFractions2D) )
   allocate ( varChange2D(nBotCells,TBNTnVarFractions2D) )
   TBNTvar2D = 0.e0
   TBNTrelVar2D = 0.e0
   varChange2D = 0.e0
# endif
   
   ! create empty fields for cumulated fluxes and initialise with ZERO
   if (TBNTabsFracOutStep>=0) then
      allocate ( TBNTflux3D(nWetCells, TBNTnFluxFractions3D) )
      allocate ( TBNTflux2D(nBotCells, TBNTnFluxFractions2D) )
      TBNTflux3D = 0.e0
      TBNTflux2D = 0.e0
   end if
   
   ! create empty fields for cumulated fluxes (incl. linked fluxes) and initialise with ZERO
   if (TBNTlinkedFluxes) then
      allocate ( TBNTlinkedFlux3D(nWetCells, TBNTnLinkedFluxFractions3D) )
      allocate ( TBNTlinkedFlux2D(nBotCells, TBNTnLinkedFluxFractions2D) )
      TBNTlinkedFlux3D = 0.e0
      TBNTlinkedFlux2D = 0.e0
   end if

   ! initialise variable fractions
   if (TBNTwarm) then ! use warmstart file
      call get_tbnt_warmstart(ierr)
   else
      if (TBNTiniFrac<=0) then
         ! no specific source selected as initial source
         if (useDummySrc) then
            ! UNTRACED source exists
            ! => set traced fractions to zero, untraced to 1
            iOff = TBNTnVars3D*(TBNTnFractions-1)
            do iFrac = iOff+1,iOff+TBNTnVars3D
               where (.not.excludeCell3D) TBNTrelVar3D(:,iFrac) = 1.e0
            end do
#ifndef TBNTnoVar2D
            iOff = TBNTnVars2D*(TBNTnFractions-1)
            do iFrac = iOff+1,iOff+TBNTnVars2D
               where (.not.excludeCell2D) TBNTrelVar2D(:,iFrac) = 1.e0
            end do
#endif
         else
            ! UNTRACED source does not exist
            ! => equally distribute initial mass among all sources
            do iFrac = 1,TBNTnVarFractions3D
               where (.not.excludeCell3D) TBNTrelVar3D(:,iFrac) = 1.e0/real(TBNTnFractions,dp)
            end do
#ifndef TBNTnoVar2D
            do iFrac = 1,TBNTnVarFractions2D
               where (.not.excludeCell2D) TBNTrelVar2D(:,iFrac) = 1.e0/real(TBNTnFractions,dp)
            end do
#endif
         end if
      else
         ! specific source selected as initial source
         if (TBNTiniFrac>TBNTnFractions) then
            ierr = ierr + 1
            write(error_msg,'(a)')'Index of selected initial fraction greater than number of sources.'
            call stop_tbnt(error_msg,ierr)
         end if 
         ! attribute initial mass to selected source
         iOff = TBNTnVars3D*(TBNTiniFrac-1)
         do iFrac = iOff+1,iOff+TBNTnVars3D
            where (.not.excludeCell3D) TBNTrelVar3D(:,iFrac) = 1.e0
         end do
#ifndef TBNTnoVar2D
         iOff = TBNTnVars2D*(TBNTiniFrac-1)
         do iFrac = iOff+1,iOff+TBNTnVars2D
            where (.not.excludeCell2D) TBNTrelVar2D(:,iFrac) = 1.e0
         end do
#endif
      end if
   end if

   end subroutine init_tbnt_fraction_fields
#endif
!=======================================================================
#ifndef TBNTonly_bulk_bud
!
! !INTERFACE:
   subroutine get_tbnt_warmstart(ierr)
#define SUBROUTINE_NAME 'get_tbnt_warmstart'
!
! !DESCRIPTION:
! get initial relative variable fractions from warmstart file
!
! !USES:

   implicit none
!
! !OUTPUT PARAMETERS:
!
  integer, intent(inout) :: ierr
!
! !LOCAL VARIABLES:
!
   integer                :: ios, ncStatus, warm_unit, timeID, iRecord
   character(len=nameLen) :: varName
   character(len=fileLen) :: filename
!
!-----------------------------------------------------------------------
#include "call-trace.inc"
    
   ios = 0

   ! open warmstart file
   filename = trim(init_nc_dir)//trim(init_nc_file)
   ncStatus = NF90_OPEN(trim(filename), NF90_NOWRITE, warm_unit)
   call nc_check(ncStatus, ios)
   if (ios/=0) then
      ierr = ierr + ios
      write(error_msg,'(a)')'NetCDF-error opening warmstart file: '//trim(filename)
      call stop_tbnt(error_msg,ierr)
   end if
   
   ! get correct time record for initialisation
   if (TBNTcontinue) then
      ! continue calculation for selected year
      iRecord = (TBNTstart-TBNToffset)/TBNTstep
   else
      ! TBNT calculation for new year => get last record of warmstart file
      ncStatus = NF90_INQUIRE(warm_unit, unlimitedDimId = timeID)
      call nc_check(ncStatus, ios)
      if (ios/=0) then
         ierr = ierr + ios
         write(error_msg,'(a)')'NetCDF-error inquiring time ID from file: '//trim(filename)
         call stop_tbnt(error_msg,ierr)
      end if
      ncStatus = NF90_INQUIRE_DIMENSION(warm_unit, timeID, len=iRecord)
      call nc_check(ncStatus, ios)
      if (ios/=0) then
         ierr = ierr + ios
         write(error_msg,'(a)')'NetCDF-error inquiring last record number from file: '//trim(filename)
         call stop_tbnt(error_msg,ierr)
      end if
   endif
   
   ! read initial values for relative variable fractions
   do iFrac = 1,TBNTnFractions
      do iVar = 1,TBNTnVars3D
         iVarPnt = TBNTvar3Dpnt(iVar,iFrac)
         varName = trim(TBNTvar(iVar)%name)//'-'//trim(TBNTfracName(iFrac))
         call get_nc_field(warm_unit, trim(varName), iRecord, ierr, fieldData3D=TBNTrelVar3D(:,iVarPnt))
         where (.not.excludeCell3D)
            TBNTrelVar3D(:,iVarPnt) = real(0.01,dp)*TBNTrelVar3D(:,iVarPnt)
         else where
            TBNTrelVar3D(:,iVarPnt) = real(0.0,dp)
         end where
      end do
#ifndef TBNTnoVar2D
      do iVar = 1,TBNTnVars2D
         iVarPnt = TBNTvar2Dpnt(iVar,iFrac)
         varName = trim(TBNTvar(TBNTnVars3D+iVar)%name)//'-'//trim(TBNTfracName(iFrac))
         call get_nc_field(warm_unit, trim(varName), iRecord, ierr, fieldData2D=TBNTrelVar2D(:,iVarPnt))
         where (.not.excludeCell2D)
            TBNTrelVar2D(:,iVarPnt) = real(0.01,dp)*TBNTrelVar2D(:,iVarPnt)
         else where
            TBNTrelVar2D(:,iVarPnt) = real(0.0,dp)
         end where
      end do
#endif
   end do
   
   ! close warmstart file
   call nc_close(filename, warm_unit, ierr, .true.)

   end subroutine get_tbnt_warmstart
#endif
!=======================================================================
!
! !INTERFACE:
   subroutine init_cell_size(ierr)
!
! !DESCRIPTION:
!  open hydrodynamic files and allocate empty arrays
!
! !USES:

   implicit none
!
! !OUTPUT PARAMETERS:
   integer, intent(inout) :: ierr
!
! !LOCAL VARIABLES:
   integer                :: i, j, k
!
!-----------------------------------------------------------------------

   ! to supress compiler warning
   i = 0
   j = 0
   k = 0

   ! allocate 3D volume buffers
   allocate ( volOldIn(nWetCells), volNewIn(nWetCells) )
   allocate ( volOld(nWetCells), volNew(nWetCells) )
   volOldIn = 0.e0
   volNewIn = 0.e0
   volOld   = 0.e0
   volNew   = 0.e0
   
#if !defined TBNTmass_fluxes || !defined TBNTnoVar2D
   ! read area fields
   allocate ( area2D(nBotCells) )
   call get_nc_field(bulk_nc_unit, trim(areaVar), 1, ierr, fieldData2D=area2D, isTimeResolved=.false.)
#ifndef TBNTmass_fluxes
#ifdef TBNTconvert_3Dto1D
   ! create 3D area array for easier multiplication of 3D fluxes
   allocate ( area3D(nWetCells) )
   do j = 1,n_y
      do i = 1,n_x
         if (k_index(i,j)<=0) cycle
         do k = 1,k_index(i,j)
            area3D(index3Dto1D(i,j,k)) = area2D(index2Dto1D(i,j))
         end do
      end do
   end do
#else
   ierr = ierr + 1
   write(error_msg,'(a)') 'Error using compiler flags:'//new_line('a')//                                       &
                          ' ==> Currently conversion from 3D to 1D fields is still required.'//new_line('a')// &
                          ' ==> Set TBNTconvert_3Dto1D.'
   call stop_tbnt(error_msg,ierr)
#endif
#endif
#endif            

   end subroutine init_cell_size
!=======================================================================
!
! !INTERFACE:
   subroutine init_bulk_fields
!
! !DESCRIPTION:
! open netcdf file and initialize bulk quantities (variables and fluxes)
!
! !USES:
!
! !OUTPUT PARAMETERS:
!
! !LOCAL VARIABLES:
!
!----------------------------------------------------------------------
#ifndef TBNTonly_bulk_bud
   ! TBNT
   ! BULK variables
   allocate ( BULKflux3D(nWetCells, TBNTnFluxes3D+TBNTnLinkedFluxes3D) )
   allocate ( BULKflux2D(nBotCells, TBNTnFluxes2D+TBNTnLinkedFluxes2D) )
   BULKflux3D = 0.e0
   BULKflux2D = 0.e0

   allocate ( BULKvar3D(nWetCells, TBNTnVars3D) )
   BULKvar3D  = 0.e0
#ifndef TBNTnoVar2D
   allocate ( BULKvar2D(nBotCells, TBNTnVars2D) )
   BULKvar2D  = 0.e0
#endif
#else
   ! BUDGET ONLY
   ! BULK variable
#ifndef TBNTnoVar2D
   if (BULKbudVar%type==TBNTpro2dVar.or.BULKbudVar%type==TBNTdia2dVar) then
      allocate ( BULKvar2D(nBotCells,1) )
      BULKvar2D = 0.e0
   else
#endif
      allocate ( BULKvar3D(nWetCells,1) )
      BULKVar3D = 0.e0
#ifndef TBNTnoVar2D
   end if
#endif
   ! BULK fluxes
   allocate ( BULKflux3D(nWetCells,BULKnBudFluxes3D), BULKflux2D(nBotCells,BULKnBudFluxes2D) )
   BULKflux3D = 0.e0
   BULKflux2D = 0.e0
#endif

   end subroutine init_bulk_fields
!=======================================================================
!
! !INTERFACE:
   subroutine get_nc_field(fileUnit, fieldName, iRecord, ierr, isTimeResolved, fieldData2D, fieldData3D)
!
! !DESCRIPTION:
! get data from NetCDF file
!
! !USES:

   implicit none
!
! !OUTPUT PARAMETERS:
   integer         , intent(in   )           :: fileUnit, iRecord
   character(len=*), intent(in   )           :: fieldName
   integer         , intent(inout)           :: ierr
   logical         , intent(in   ), optional :: isTimeResolved
   real(dp)        , intent(  out), optional :: fieldData2D(nBotCells), fieldData3D(nWetCells)
!
! !LOCAL VARIABLES:
!
   integer                                   :: fieldUnit, fieldType, fieldDimensionality, nFieldDims, dimLength
   integer                                   :: ios, i, iLon, iLat, iDep, n, nNCfieldDims, ncStatus
   integer,    dimension(:)    , allocatable :: fieldDimIDs, ncStart, ncCount
#ifdef TBNTconvert_3Dto1D
   integer                                   :: jSwitch, j, k
#endif
   real   ,    dimension(:,:  ), allocatable :: fieldData2Din 
   real   ,    dimension(:,:,:), allocatable :: fieldData3Din
   character(len=nameLen)                    :: dimName
!-----------------------------------------------------------------------
!
   
   fieldDimensionality = 0
   if (present(fieldData2D)) then
      fieldDimensionality = 2
   elseif (present(fieldData3D)) then
      fieldDimensionality = 3
   elseif (present(fieldData2D).and.present(fieldData3D)) then
      ierr = ierr + 1
      write(error_msg,'(a)')'Two output parameters ''fieldData2D'' and ''fieldData3D'' set. Only one allowed.'
      call stop_tbnt(error_msg,ierr)
   else
      ierr = ierr + 1
      write(error_msg,'(a)')'No output parameter set. Set ''fieldData2D'' or ''fieldData3D''.'
      call stop_tbnt(error_msg,ierr)
   end if
   
   if (.not.present(isTimeResolved).or.(present(isTimeResolved).and.isTimeResolved)) then
      nFieldDims = fieldDimensionality + 1
   else
      nFieldDims = fieldDimensionality
   end if

   allocate ( fieldDimIDs(nFieldDims), ncStart(nFieldDims), ncCount(nFieldDims) )
   fieldDimIDs = -1
   ncStart = 1
   ncCount = 1
   
   ! inquire field ID
   ios = 0
   ncStatus = NF90_INQ_VARID(fileUnit, trim(fieldName), fieldUnit)
   call nc_check(ncStatus, ios)
   if (ios/=0) then
      ierr = ierr + ios
      deallocate ( fieldDimIDs, ncStart, ncCount )
      write(error_msg,'(a)')'NetCDF-error reading field ID: '//trim(fieldName)//' - '//trim(nf90_strerror(ncStatus))
      call stop_tbnt(error_msg,ierr)
   end if
   ! inquire field metadata
   ncStatus = NF90_INQUIRE_VARIABLE(fileUnit, fieldUnit, xtype=fieldType, nDims=nNCfieldDims, dimids=fieldDimIDs)
   call nc_check(ncStatus, ios)
   if (ios/=0) then
      ierr = ierr + ios
      deallocate ( fieldDimIDs, ncStart, ncCount )
      write(error_msg,'(a)')'NetCDF-error reading field meta data: '//trim(fieldName)//' - '//trim(nf90_strerror(ncStatus))
      call stop_tbnt(error_msg,ierr)
   end if
   if (nNCfieldDims /= nFieldDims) then
      ierr = ierr + 1
      write(error_msg,'(a)')'NetCDF dimensions do not fit field dimensions for: '//trim(fieldName)
      call stop_tbnt(error_msg,ierr)
   end if
   ! inquire field dimension lengths and set start and count for reading
   iLon = 0
   iLat = 0
   iDep = 0
   do n = 1,nFieldDims
      ncStatus = NF90_INQUIRE_DIMENSION(fileUnit,fieldDimIDs(n),dimName,dimLength)
      call nc_check(ncStatus, ios)
      if (ios/=0) then
         ierr = ierr + ios
         deallocate ( fieldDimIDs, ncStart, ncCount )
         write(error_msg,'(a,i5,a)')'NetCDF-error reading dimension meta data: ',fieldDimIDs(n),' - '//trim(nf90_strerror(ncStatus))
         call stop_tbnt(error_msg,ierr)
      end if
      if     (trim(dimName) == trim(NCdims(1)%name)) then
         iLon = n
         ncCount(n) = dimLength
      elseif (trim(dimName) == trim(NCdims(2)%name)) then
         iLat = n
         ncCount(n) = dimLength
      elseif (trim(dimName) == trim(NCdims(3)%name)) then
         iDep = n
         ncCount(n) = dimLength
      elseif (n==nFieldDims) then
         ncStart(n) = iRecord
      end if
   end do
   if (fieldDimensionality==2) then
      if (iLon==0.or.iLat==0) ierr = ierr + 1
   else
      if (iLon==0.or.iLat==0.or.iDep==0) ierr = ierr + 1
   end if
   if (ierr/=0) then
      deallocate ( fieldDimIDs, ncStart, ncCount )
      write(error_msg,'(a)')'NetCDF-error: At least one required dimension was not found.'
      call stop_tbnt(error_msg,ierr)
   end if
   
   if (fieldDimensionality==2) then ! get 2D field
      allocate ( fieldData2Din(n_x,n_y) )
      ncStatus = NF90_GET_VAR(fileUnit, fieldUnit, fieldData2Din, start=ncStart, count=ncCount)
#ifdef TBNTconvert_3Dto1D
      do j = 1,n_y
# ifdef TBNTswitchNS
         jSwitch = n_y-j+1
# else
         jSwitch = j
# endif
         do i = 1,n_x
            iPut = index2Dto1D(i,jSwitch)
            if (iPut<=0) cycle
            fieldData2D(iPut) = fieldData2Din(i,j) 
         end do
      end do
#endif
      deallocate ( fieldData2Din )
   else if (fieldDimensionality==3) then ! get 3D field
      allocate ( fieldData3Din(n_x,n_y,n_z) )
      ncStatus = NF90_GET_VAR(fileUnit, fieldUnit, fieldData3Din, start=ncStart, count=ncCount)
#ifdef TBNTconvert_3Dto1D
      do j = 1,n_y
# ifdef TBNTswitchNS
         jSwitch = n_y-j+1
# else
         jSwitch = j
# endif
         do i = 1,n_x
            if (k_index(i,j)<=0) cycle
            do k = 1,k_index(i,j)
               iPut = index3Dto1D(i,jSwitch,k)
               if (iPut<=0) cycle
               fieldData3D(iPut) = fieldData3Din(i,j,k)
            end do
         end do
      end do
#endif
      deallocate ( fieldData3Din )
   end if
   call nc_check(ncStatus, ios)
   if (ios/=0) then
      ierr = ierr + ios
      deallocate ( fieldDimIDs, ncStart, ncCount )
      write(error_msg,'(3a,i5)')'NetCDF-error reading ',trim(fieldName),' at iStep = ',iStep
      call stop_tbnt(error_msg,ierr)
   end if
   
   deallocate ( fieldDimIDs, ncStart, ncCount )

   end subroutine get_nc_field
!=======================================================================
end module mod_etrac_init
