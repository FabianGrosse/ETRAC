!
! !MODULE: etrac_common.f90 --- common definitions for ETRAC

! !INTERFACE:
   MODULE mod_etrac_common
!
! !DESCRIPTION:
! This module is the container for TBNT variables intensively used
!
! !USES:
  use netcdf
!
   implicit none
!
!  default: all is public.
   public
!
! !PARAMETERS
!
! character/string lengths
   integer               , parameter                      :: yearLen       = 4
   integer               , parameter                      :: formatLen     = 50
   integer               , parameter                      :: nameLen       = 50
   integer               , parameter                      :: fileLen       = 200
   integer               , parameter                      :: lineLen       = 300
   integer               , parameter                      :: errorLen      = 999
! maximum number of BULK fluxes
   integer               , parameter                      :: TBNTnMaxFluxes = 800
! real kind (double precision)
   integer               , parameter                      :: dp            = selected_real_kind(15,307)
! time related
   integer               , parameter                      :: iMinPerDay    = 1440 ! integer "minutes per day"
! values for variable types
   integer               , parameter                      :: TBNTpro3dVar  = 1    ! 3D prognostic variable
   integer               , parameter                      :: TBNTpro2dVar  = 2    ! 2D prognostic variable
   integer               , parameter                      :: TBNTdia3dVar  = 3    ! 3D diagnostic variable (i.e. derived)
   integer               , parameter                      :: TBNTdia2dVar  = 4    ! 2D diagnostic variable (i.e. derived)
   integer               , parameter                      :: TBNTdum3dVar  = 5    ! 3D dummy variable (no real state variable in underlying model; created internally during TBNT calculation)
   integer               , parameter                      :: TBNTdum2dVar  = 6    ! 2D dummy variable (no real state variable in underlying model; created internally during TBNT calculation)
! values for flux types
   integer               , parameter                      :: TBNTpelagFlux = 1    ! pelagic flux
   integer               , parameter                      :: TBNTatm2DFlux = 2    ! 2D atmospheric deposition
   integer               , parameter                      :: TBNTatm3DFlux = 3    ! 3D atmospheric deposition
   integer               , parameter                      :: TBNTsed2DFlux = 4    ! 2D flux between sediment and ocean
   integer               , parameter                      :: TBNTsed3DFlux = 5    ! 3D flux between sediment and ocean
   integer               , parameter                      :: TBNTadvecFlux = 6    ! advective transport flux
   integer               , parameter                      :: TBNTdiffuFlux = 7    ! diffusive transport flux
   integer               , parameter                      :: TBNTrvdisFlux = 8    ! river discharge
   integer               , parameter                      :: TBNTprevaFlux = 9    ! precipitation minus evaporation
   integer               , parameter                      :: TBNTa2s2DFlux =10    ! 2D air-sea flux (o2o/o2c)
   integer               , parameter                      :: TBNTa2s3DFlux =11    ! 3D air-sea flux (o2o/o2c)
   integer               , parameter                      :: TBNTder2DFlux =12    ! 2D derived flux (i.e. sediment fluxes on derived variables, e.g. bap_sed)
   integer               , parameter                      :: TBNTder3DFlux =13    ! 3D derived flux (i.e. physical and sediment fluxes on derived variables, e.g. tu_bap)
! values for source types
   integer               , parameter                      :: TBNTriverSource = 1  ! river source
   integer               , parameter                      :: TBNTopenbSource = 2  ! open boundary source
   integer               , parameter                      :: TBNTatmosSource = 3  ! atmospheric source
   integer               , parameter                      :: TBNTdummySource = 4  ! dummy for all sources being not traced being added to the untraced substance
! units and file names for namelist and log file
   integer               , parameter                      :: tbnt_log_unit      = 1
   integer               , parameter                      :: tbnt_settings_unit = 2
   character(len=fileLen), parameter                      :: tbnt_log_file      = 'etrac_logfile.dat'
   character(len=fileLen)                                 :: tbnt_settings_file = 'etrac_set.nml'
   character(len=fileLen)                                 :: model_set_dir        ! directory with model setup files
   character(len=fileLen)                                 :: tbnt_set_dir         ! directory with TBNT setup files
! input unit for model input files
   integer               , parameter                      :: input_unit         = 10
! output unit for target information
   integer               , parameter                      :: target_unit        = 11
! fill value for NetCDF files
   real(dp)          , parameter                          :: fail               = -9.999e10
   
! !TYPES needed for NetCDF file dimensions and grid
   type TBNTdimProps
      character(len=nameLen)                              :: name                 ! dimension name
      integer                                             :: ncid                 ! dimension NetCDF ID
      integer                                             :: size                 ! dimension size
   end type TBNTdimProps
   
   type TBNTgridVarProps
      character(len=nameLen)                              :: name                 ! variable name
      character(len=nameLen)                              :: unit                 ! variable unit
      integer                                             :: ncid                 ! variable NetCDF ID
      integer , dimension(:)  , allocatable               :: dimid                ! dimension NetCDF IDs related to variable
      integer , dimension(:)  , allocatable               :: size                 ! variable size
      real(dp), dimension(:)  , allocatable               :: data1D               ! variable values (if 1D variable)
      real(dp), dimension(:,:), allocatable               :: data2D               ! variable values (if 2D variable)
   end type TBNTgridVarProps
   
! !TYPES needed for variables and flux collection
   type TBNTauxVarProps
      character(len=nameLen)                              :: name                 ! auxiliary variable name (e.g. "n4n_n3n")
      integer                                             :: type                 ! auxiliary variable type
      integer                                             :: iTBNT                ! BULK variable index
      real(dp)                                            :: fac                  ! variable multiplier to calculate dummy variable from auxiliary state variable (default = 1., i.e. no auxiliary flux needed)
   end type TBNTauxVarProps
   
   type TBNTmodVarProps
      character(len=nameLen)                              :: name                 ! variable name (e.g. "n3n")
      integer                                             :: type                 ! variable type: 1 = prognostic, 2 = derived, 3 = sediment, 4 = dummy
      integer                                             :: iTBNT                ! TBNT variable index
      type(TBNTauxVarProps)                               :: auxVar               ! auxiliary variable, if variable is dummy variable
   end type TBNTmodVarProps

   type TBNTauxFluxProps
      character(len=nameLen)                              :: name                 ! flux name (e.g. "n4n_n3n")
      integer                                             :: type                 ! flux type
      real(dp)                                            :: fac                  ! flux multiplier to get right flux for derived variable (default = 1., i.e. no auxiliary flux needed)
   end type TBNTauxFluxProps

   type TBNTfluxProps
      character(len=nameLen)                              :: name                 ! flux name (e.g. "n4n_n3n")
      integer                                             :: type                 ! flux type
      integer                                             :: i_comp               ! directional component (ZERO = default, NON-ZERO = transport flux)
      integer                                             :: iTBNT                ! BULK flux index
      type(TBNTmodVarProps)                               :: varIn                ! input variable
      type(TBNTmodVarProps)                               :: varOut               ! output variable
      type(TBNTauxFluxProps)                              :: auxFlux              ! auxiliary flux if flux is derived flux
   end type TBNTfluxProps

! !TYPES needed for source collection
   type TBNTsubsourceProps
      character(len=nameLen)                              :: name                 ! name of subsource (e.g. "Elbe")
      integer                                             :: type                 ! type of sub source: 1 = river, 2 = open boundary, 3 = atmospheric
      integer                                             :: index                ! if river source: index = ir, else: index = 0
      integer                                             :: nInputs              ! number of source points
      integer , dimension(:), allocatable                 :: i_from               ! 1D indices of source points
      integer , dimension(:), allocatable                 :: i_to                 ! 1D indices of target points
      integer , dimension(:), allocatable                 :: i_comp               ! index of the directional component (only relevant for open boundary sources)
      real(dp), dimension(:), allocatable                 :: sign                 ! in + or - direction of the directional component? (only relevant for open boundary sources)
   end type TBNTsubsourceProps

   type TBNTsourceProps
      character(len=nameLen)                              :: name                 ! source name (e.g. "Elbe and Weser")
      integer                                             :: type                 ! type of source: 1 = river, 2 = open boundary, 3 = atmospheric, 4 = dummy
      integer                                             :: nSubsources          ! number of different subsources
      type(TBNTsubsourceProps), dimension(:), allocatable :: subsource            ! subsources
   end type TBNTsourceProps

#ifndef TBNTonly_bulk_bud
! !TYPES needed for target areas
   type TBNTtarAreaProps
      character(len=nameLen)                              :: name                 ! target area name (e.g. "NL-C2")
      logical, dimension(:    ), allocatable              :: isPart2D             ! logical index array for target area (2D)
      logical, dimension(:    ), allocatable              :: isPart3D             ! logical index array for target area (3D)
   end type TBNTtarAreaProps
! !TYPES needed for target variables
   type TBNTtarVarProps
      character(len=nameLen)                              :: name                 ! target variable name (e.g. "total nitrogen")
      character(len=nameLen)                              :: code                 ! target variable code (e.g. "TN")
      logical, dimension(:), allocatable                  :: targetPeriod         ! mask for target period
      integer                                             :: nVars                ! number of model variables associated with the target variable
      type(TBNTmodVarProps), dimension(:), allocatable    :: vars                 ! list of model variables associated with the target variable
   end type TBNTtarVarProps
#endif

! !VARIABLES
   ! MODEL RELATED: GRID
   integer                                                :: nWetCells            ! number of WET grid cells
   integer                                                :: nComponents          ! number of directional components
   integer                                                :: nSurfCells           ! number of WET surface grid cells
   integer                                                :: nBotCells            ! number of WET bottom grid cells
   integer , dimension(:  ), allocatable                  :: surf2pelag           ! translator of surface cell indices to adjacent pelagic cells
   integer , dimension(:  ), allocatable                  :: bottom2pelag         ! translator of bottom cell indices to adjacent pelagic cells
   integer , dimension(:,:), allocatable                  :: neighbours           ! list of grid cell neighbours
   real(dp), dimension(:  ), allocatable                  :: compSign             ! sign of transport flux relative to its location of definition (i.e. source or sink)
   logical , dimension(:  ), allocatable                  :: surfaceCell          ! mask of surface cells
   logical , dimension(:  ), allocatable                  :: bottomCell           ! mask of bottom cells
   logical , dimension(:  ), allocatable                  :: excludeCell2D        ! mask of bottom cells excluded from calculations
   logical , dimension(:  ), allocatable                  :: excludeCell3D        ! mask of pelagic cells excluded from calculations
   character(len=nameLen)                                 :: areaVar              ! name of area variable to be read from NetCDF file
   character(len=nameLen)                                 :: volVar               ! name of volume variable to be read from NetCDF file
#ifdef TBNTconvert_3Dto1D
   integer                                                :: n_x                  ! number of indices in horizontal x-direction
   integer                                                :: n_y                  ! number of indices in horizontal y-direction
   integer                                                :: n_z                  ! number of indices in vertical direction
   integer                                                :: TBNTxStart           ! start index in x-direction
   integer                                                :: TBNTxEnd             ! end index in x-direction
   integer                                                :: TBNTyStart           ! start index in y-direction
   integer                                                :: TBNTyEnd             ! end index in y-direction
   integer               , dimension(:,:,:), allocatable  :: index3Dto1D          ! 1D index for 3D indices
   integer               , dimension(:,:)  , allocatable  :: index2Dto1D          ! 1D index for 2D indices
   integer               , dimension(:,:)  , allocatable  :: index1Dto3D          ! 3D indices for 1D index vector
   integer               , dimension(:,:)  , allocatable  :: index1Dto2D          ! 2D indices for 1D index vector
   integer               , dimension(:,:)  , allocatable  :: k_index              ! 2D field with maximum depth index
   type(TBNTdimProps)    , dimension(4)                   :: NCdims               ! list of NetCDF file dimensions needed for output
   type(TBNTgridVarProps), dimension(4)                   :: NCgrid               ! NetCDF grid data needed for output
#endif
   ! run related
   character(len=nameLen)                                 :: TBNTrun              ! run ID
   logical                                                :: TBNTwarm             ! use warm start/initialisation file?
   logical                                                :: TBNTcontinue         ! continue writing into existing files (in case of restart during year)
   ! time related
   integer                                                :: TBNTyear             ! TBNT year
   integer                                                :: TBNTstep             ! TBNT time step in minutes
   integer                                                :: TBNTstart            ! TBNT start time in minutes
   integer                                                :: TBNTend              ! TBNT end time in minutes
   integer                                                :: TBNToffset           ! TBNT start offset in minutes
   integer                                                :: TBNTnSubSteps        ! # of substeps per TBNT step
   integer                                                :: TBNTyearDays         ! days of year
   integer                                                :: TBNTiStart           ! start index of time loop
   integer                                                :: TBNTiEnd             ! end index of time loop
   integer                                                :: TBNTiOffset          ! start offset as number of time steps
   integer                                                :: TBNTnSteps           ! steps of TBNT run
   integer                                                :: TBNTnDays            ! days of TBNT run
   real(dp)                                               :: TBNTtimeFac          ! time scaling factor per substep (=1/TBNTnSubSteps)
   real(dp)                                               :: rMinPerDay           ! floating point "minutes per day"
   real(dp)                                               :: outTimeStart         ! initial time of input/output NetCDF files
   real(dp)                                               :: outTimeFac           ! conversion factor for NetCDF output time
   ! source related
   integer                                                :: TBNTnTracedElements  ! number of traced elements
   integer                                                :: TBNTnRiverFractions  ! number of river source groups
   integer                                                :: TBNTnOpenBFractions  ! number of open boundary source groups
   integer                                                :: TBNTnAtmosFractions  ! number of atmospheric deposition source groups
   integer                                                :: TBNTnFractions       ! total number of source groups/fractions
   integer                                                :: TBNTiniFrac          ! index of source fraction to which initial mass is atrributed
   integer                                                :: TBNTnDummyRivers     ! number of rivers going into dummy fraction
   integer                                                :: TBNTnDummyAtmos      ! number of atmospheric deposition sources going into dummy fraction
   type(TBNTsourceProps), dimension(:)      , allocatable :: TBNTsource           ! collection of all fraction sources
   ! variable + flux related
   integer                                                :: TBNTnFluxes3D        ! total number of 3D bulk fluxes
   integer                                                :: TBNTnFluxFractions3D ! total number of 3D fractions fluxes
   integer              , dimension(:,:)    , allocatable :: TBNTflux3Dpnt        ! pointer for 3D fraction fluxes
   real(dp)             , dimension(:,:)    , allocatable :: TBNTflux3D           ! array of 3D fraction fluxes for output
   real(dp)             , dimension(:,:)    , allocatable :: tempFlxFrac3D        ! array of 3D fraction fluxes used for calculation
   integer                                                :: TBNTnFluxes2D        ! total number of 2D bulk fluxes
   integer                                                :: TBNTnFluxFractions2D ! total number of 2D fraction fluxes
   integer              , dimension(:,:)    , allocatable :: TBNTflux2Dpnt        ! pointer for 2D fraction fluxes
   real(dp)             , dimension(:,:)    , allocatable :: TBNTflux2D           ! array of 2D fraction fluxes used for output
   real(dp)             , dimension(:,:)    , allocatable :: tempFlxFrac2D        ! array of 2D fraction fluxes used for calculation
   integer                                                :: TBNTnVars3D          ! total number of 3D bulk variables
   integer                                                :: TBNTnVarFractions3D  ! total number of 3D variable fractions
   integer              , dimension(:,:)    , allocatable :: TBNTvar3Dpnt         ! pointer for 3D variable fractions
   real(dp)             , dimension(:,:)    , allocatable :: TBNTvar3D            ! array of 3D absolute variable fractions
   real(dp)             , dimension(:,:)    , allocatable :: TBNTrelVar3D         ! array of 3D absolute variable fractions
   real(dp)             , dimension(:,:)    , allocatable :: varChange3D          ! array of integrated change in 3D fraction variables per time step
   integer                                                :: TBNTnVars2D          ! total number of 2D bulk variables
   integer                                                :: TBNTnVarFractions2D  ! total number of 2D variable fractions
   integer              , dimension(:,:)    , allocatable :: TBNTvar2Dpnt         ! pointer for 2D variable fractions
   real(dp)             , dimension(:,:)    , allocatable :: TBNTvar2D            ! array of 2D absolute variable fractions
   real(dp)             , dimension(:,:)    , allocatable :: TBNTrelVar2D         ! array of 2D absolute variable fractions
   real(dp)             , dimension(:,:)    , allocatable :: varChange2D          ! array of integrated change in 2D fraction variables per time step
   type(TBNTfluxProps)  , dimension(:)      , allocatable :: TBNTflux             ! collection of all fluxes on fractionated variables
   type(TBNTmodVarProps), dimension(:)      , allocatable :: TBNTvar              ! collection of all fractionated variables
   ! linked fluxes related
   logical                                                :: TBNTlinkedFluxes           ! switch for linked fluxes (1=on, 0=off)
   integer                                                :: TBNTnLinkedFluxes3D        ! total number of linked 3D bulk fluxes
   integer                                                :: TBNTnLinkedFluxFractions3D ! total number of linked 3D fractions fluxes
   integer                                                :: TBNTnLinkedFluxes2D        ! total number of linked 2D bulk fluxes
   integer                                                :: TBNTnLinkedFluxFractions2D ! total number of linked 2D fraction fluxes
   integer              , dimension(:,:)    , allocatable :: TBNTlinkedFlux3Dpnt        ! pointer for 3D fraction fluxes
   integer              , dimension(:,:)    , allocatable :: TBNTlinkedFlux2Dpnt        ! pointer for 2D fraction fluxes
   real(dp)             , dimension(:,:)    , allocatable :: TBNTlinkedFlux3D           ! array of 3D linked fraction fluxes used for output
   real(dp)             , dimension(:,:)    , allocatable :: TBNTlinkedFlux2D           ! array of 2D linked fraction fluxes used for output
   real(dp)             , dimension(:,:)    , allocatable :: lnkFlxFrac3D               ! array of 3D linked fraction fluxes used for calculation
   real(dp)             , dimension(:,:)    , allocatable :: lnkFlxFrac2D               ! array of 32 linked fraction fluxes used for calculation
#ifndef TBNTonly_bulk_bud
   ! output related
   character(len=fileLen)                                 :: output_dir           ! output directory
   integer                                                :: relFracFileID        ! netCDF file ID for relative fractions
   character(len=fileLen)                                 :: relFracFile          ! netCDF filename for relative fractions
   integer               , dimension(:)     , allocatable :: absFracFileID        ! netCDF file IDs for absolute fractions
   character(len=fileLen), dimension(:)     , allocatable :: absFracFile          ! netCDF filenames for absolute fractions
   integer                                                :: lnkFracFileID        ! netCDF file ID for linked fraction fluxes
   character(len=fileLen)                                 :: lnkFracFile          ! netCDF filename for linked fraction fluxes
   integer                                                :: TBNTrelFracOutStep   ! output interval for relative fractions
   integer                                                :: TBNTabsFracOutStep   ! output interval for fraction fluxes and variables
   integer                                                :: TBNTlnkFracOutStep   ! output interval for linked fractions fluxes
   logical                                                :: TBNTtargetOut        ! output switch for target areas (1=on, 0=off)
   integer                                                :: TBNTnTargetAreas     ! number of target areas
   type(TBNTtarAreaProps), dimension(:)     , allocatable :: TBNTtargetAreas      ! list with target areas
   integer                                                :: TBNTnTargetVars      ! number of target variables
   type(TBNTtarVarProps) , dimension(:)     , allocatable :: TBNTtargetVars       ! list with target variables
   real(dp)              , dimension(:,:,:) , allocatable :: TBNTtargetVals       ! array with output values for target areas and variable fractions
   real(dp)              , dimension(:,:)   , allocatable :: TBNTtargetWeights    ! array with weights (i.e., total mass) of target variables for target areas
   character(len=nameLen), dimension(:)     , allocatable :: TBNTfracName         ! output suffix for fractions
   integer               , dimension(:,:)   , allocatable :: TBNTabsFlxID3D       ! output IDs for 3D absolute fraction fluxes
   integer               , dimension(:,:)   , allocatable :: TBNTlnkFlxID3D       ! output IDs for 3D linked fraction fluxes
   integer               , dimension(:,:)   , allocatable :: TBNTabsFlxID2D       ! output IDs for 2D absolute fraction fluxes
   integer               , dimension(:,:)   , allocatable :: TBNTlnkFlxID2D       ! output IDs for 2D linked fraction fluxes
   integer               , dimension(:,:)   , allocatable :: TBNTrelVarID3D       ! output IDs for 3D relative fraction variables
   integer               , dimension(:,:)   , allocatable :: TBNTabsVarID3D       ! output IDs for 3D absolute fraction variables
#ifndef TBNTnoVar2D
   integer               , dimension(:,:)   , allocatable :: TBNTrelVarID2D       ! output IDs for 2D relative fraction variables
   integer               , dimension(:,:)   , allocatable :: TBNTabsVarID2D       ! output IDs for 2D absolute fraction variables
#endif
   character(len=formatLen)                               :: TBNTrelFracFormat    ! output format for warmstart file
#endif
   ! auxiliary variables
   integer                                                :: bulk_nc_unit         ! unit for NetCDF file with bulk quantities
   integer                                                :: ierr                 ! error counter
   character(len=errorLen)                                :: error_msg            ! error message
   integer                                                :: iStep, iSubStep      ! auxiliary index variables
   integer                                                :: iFrac, iVar, iFlux   ! auxiliary index variables
   integer                                                :: iSource, iPut, iOff  ! auxiliary index variables
   integer                                                :: iVarPnt, iFlxPnt     ! auxiliary pointer variables
   ! auxiliary variables to compute amounts of substances and fluxes on amounts of substances
   real(dp)             , dimension(:)      , allocatable :: volOldIn, volNewIn   ! 1D array for daily volumes from vol.bin
   real(dp)             , dimension(:)      , allocatable :: volOld, volNew       ! 1D arrays for volumes interpolated to start & end of time step
#if !defined TBNTnoVar2D || !defined TBNTmass_fluxes
   real(dp)             , dimension(:)      , allocatable :: area2D               ! area field for sediment cells
#endif
#ifndef TBNTmass_fluxes
   real(dp)             , dimension(:)      , allocatable :: area3D               ! area field for pelagic cells
#endif
   ! auxiliary variables to handle rivers
   integer                                                :: nRivFluxes           ! number of river fluxes
   integer              , dimension(:,:)    , allocatable :: iRvDisFlux           ! pointer to discharge flux for different variables
   ! BULK variables and fluxes
   real(dp)             , dimension(:,:)    , allocatable :: BULKvar3D            ! 3D bulk variable
   real(dp)             , dimension(:,:)    , allocatable :: BULKvar2D            ! 2D bulk variable
   real(dp)             , dimension(:,:)    , allocatable :: BULKflux3D           ! 3D bulk flux
   real(dp)             , dimension(:,:)    , allocatable :: BULKflux2D           ! 2D bulk flux

#ifdef TBNTonly_bulk_bud
   ! variables for ONLY BUDGET MODE
   integer                                                :: iBud                 ! budget location
   type(TBNTmodVarProps)                                  :: BULKbudVar           ! bulk variable to be balanced
   real(dp)                                               :: BULKbudVarNew        ! time-integrated balance variable
   integer                                                :: BULKnBudFluxes       ! total number of balance fluxes
   integer                                                :: BULKnBudFluxes2D     ! number of 2D balance fluxes
   integer                                                :: BULKnBudFluxes3D     ! number of 3D balance fluxes
   type(TBNTfluxProps)  , dimension(:)      , allocatable :: BULKbudFLux          ! list of bulk fluxes to be balanced
#ifdef TBNTbulk_bud_out
   integer, parameter                                     :: BULKbudUnit = 500    ! unit for bulk flux output file
#endif
#endif
   
!==============================================================================================================================================

   contains
   
!==============================================================================================================================================
!
! !INTERFACE:
   subroutine breakLine(length)
!
! !DESCRIPTION:
!  write separation line to log file
!
   implicit none
!
! !OUTPUT PARAMETERS:
   integer, intent(in) :: length
!
!-----------------------------------------------------------------------
   
   write(tbnt_log_unit,'(a)') new_line('a')//repeat('=', length)//new_line('a')

   end subroutine breakLine
!=======================================================================
!
! !INTERFACE:
   subroutine get_headLen(fUnit, fName, nSkip, ierr)
!
! !DESCRIPTION:
! get length of header
!
   implicit none
!
! !OUTPUT PARAMETERS:
   integer               , intent(in   ) :: fUnit
   character(len=fileLen), intent(in   ) :: fName
   integer               , intent(  out) :: nSkip
   integer               , intent(inout) :: ierr
!
! !LOCAL VARIABLES:
   integer                :: ios
   character(len=lineLen) :: str
!
!-----------------------------------------------------------------------

   read(fUnit, '(a)') str
   nSkip = 0
   do while (index(str,'!+') > 0)
      nSkip = nSkip + 1
      read(fUnit, '(a)', iostat = ios) str
      if (ios/=0) exit
   end do
   if (ios/=0) then
      ierr = ierr + ios
      write(error_msg,'(2a)')'Error reading: ',trim(fName)
      call stop_tbnt(error_msg,ierr)
   end if
   rewind(fUnit)

   end subroutine get_headLen
!=======================================================================
!
! !INTERFACE:
   subroutine skiplines(fUnit, nLines)
!
! !DESCRIPTION:
! skips nLines of file fUnit
!
   implicit none
!
! !LOCAL VARIABLES:
   integer            :: fUnit, nLines, n
   character(len=lineLen) :: str
!
!-----------------------------------------------------------------------

   do n =1,nLines
      read(fUnit,'(a)')str
   end do

   end subroutine skiplines
!=======================================================================
!
! !INTERFACE:
   subroutine string_check(str1, str2, identical)
!
! !DESCRIPTION:
! checks, if two strings match
!
! !USES:

   implicit none
!
! !OUTPUT PARAMETERS:
   character(len=*)     :: str1, str2
   logical, intent(out) :: identical

!-----------------------------------------------------------------------

   identical = .false.
   if (trim(str1) == trim(str2)) identical = .true.

   end subroutine string_check
! ======================================================================
!
! !INTERFACE:
   subroutine capitalise_string(str)
!
! !DESCRIPTION:
! capitalize a string
!
   implicit none
!
! !OUTPUT PARAMETERS:
   character(len=*), intent(inout) :: str
!
! !LOCAL VARIABLES:
   character(len=*), parameter :: ABCc = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
   character(len=*), parameter :: ABCs = 'abcdefghijklmnopqrstuvwxyz'
   
   integer :: n, nFirst, nLast, is
!
!-----------------------------------------------------------------------

   nFirst = scan(str,ABCs)
   if (nFirst<=0) return
   nLast = scan(str,ABCs,.true.)
   do n = nFirst,nLast
      is = scan(ABCs, str(n:n))
      if (is>0) str(n:n) = ABCc(is:is)
   end do

   end subroutine capitalise_string
!=======================================================================
!
! !INTERFACE:
   subroutine interpret_mathString(str, num)
!
! !DESCRIPTION:
! interpret string containing mathematical expressions: '*' and '/'
!
   implicit none
!
! !OUTPUT PARAMETERS:
  character(len=nameLen), intent(inout) :: str
  real(dp)              , intent(  out) :: num
!
! !LOCAL VARIABLES:
   character(len=2), parameter :: ops = '*/'
   
   integer          :: ind
   real(dp)         :: fac
   character(len=1) :: op
!
!-----------------------------------------------------------------------
   
   ind = scan(str,ops)
   if (ind <= 0) then
      read(str(1:len_trim(str)),*) num
      return
   end if
   
   read(str(1:ind-1),*)num
   op = str(ind:ind)
   str = str(ind+1:len_trim(str))
   ind = scan(str,ops)
   do while (ind > 0)
      read(str(1:ind-1),*) fac
      if (op==ops(1:1)) then
         num = num*fac
      else
         num = num/fac
      end if
      op = str(ind:ind)
      str = str(ind+1:len_trim(str))
      ind = scan(str,ops)
   end do
   
   read(str(1:len_trim(str)),*) fac
   if (op==ops(1:1)) then
      num = num*fac
   else
      num = num/fac
   end if

   end subroutine interpret_mathString
!=======================================================================
!
! !INTERFACE:
   subroutine nc_check(ncStatusIn,iret)
!
! !DESCRIPTION:
! checks NetCDF-action status
!
   implicit none
!
! !OUTPUT PARAMETERS:   
   integer, intent(in ) :: ncStatusIn
   integer, intent(out) :: iret
!   
!-----------------------------------------------------------------------
!
   iret = 0
   if (ncStatusIn /= NF90_NOERR) iret = 1

   end subroutine nc_check
!=======================================================================
! !INTERFACE:
   subroutine nc_close(filename, fileID, ierr, stop_flag)
#define SUBROUTINE_NAME 'nc_close'
!
! !DESCRIPTION:
!  initialise a new NetCDF file
!
   implicit none
!
! !OUTPUT PARAMETERS:
   character(len=fileLen), intent(in   ) :: filename
   integer               , intent(in   ) :: fileID
   integer               , intent(inout) :: ierr
   logical               , intent(in   ) :: stop_flag
!
! LOCAL VARIABLES
   integer                               :: ios, ncStatus
!   
!-----------------------------------------------------------------------
#include "call-trace.inc"

   ios = 0
   
   ! close netcdf file
   ncStatus = nf90_close(fileID)
   call nc_check(ncStatus, ios)
   if (ios/=0) then
      if (stop_flag) then
         ierr = ierr + ios
         write(error_msg,'(a)')'NetCDF-error closing file: '//trim(filename)//' - '//trim(nf90_strerror(ncStatus))
         call stop_tbnt(error_msg,ierr)
      else
         write(tbnt_log_unit,'(a)')'NetCDF-error closing file (during program abortion): '//trim(filename)//' - '//trim(nf90_strerror(ncStatus))
      end if
   end if
   
   end subroutine nc_close
!=======================================================================
#ifndef TBNTonly_bulk_bud
!
! !INTERFACE:
   subroutine finalize_tbnt_output(ierr, stop_flag)
#define SUBROUTINE_NAME 'finalize_tbnt_output'
!
! !DESCRIPTION:
!  close TBNT output
!
! !USES:
   implicit none
!
! !OUTPUT PARAMETERS:
   integer, intent(inout) :: ierr
   logical, intent(in   ) :: stop_flag
!
! !LOCAL VARIABLES:
!
   integer :: iTarget
!
!-----------------------------------------------------------------------
#include "call-trace.inc"
   
   if (TBNTrelFracOutStep>=0) call nc_close(relFracFile, relFracFileID, ierr, stop_flag)

   if (TBNTabsFracOutStep>=0) then
      do iFrac = 1,TBNTnFractions
         call nc_close(absFracFile(iFrac), absFracFileID(iFrac), ierr, stop_flag)
      end do
      deallocate ( TBNTflux2D, TBNTflux3D )
   end if
   
   if (TBNTlinkedFluxes) then
      call nc_close(lnkFracFile, lnkFracFileID, ierr, stop_flag)
      deallocate ( TBNTlinkedFlux2D, TBNTlinkedFlux3D )
   end if
   
   if (TBNTtargetOut) then
      do iTarget = 1,TBNTnTargetVars
         close(target_unit+iTarget)
      end do
      deallocate ( TBNTtargetVals, TBNTtargetWeights )
   end if   

   end subroutine finalize_tbnt_output
#endif
!=======================================================================
!
! !INTERFACE:
   subroutine stop_tbnt(msg,ierr)
!
! !DESCRIPTION:
! stops tbnt run
!
! !USES:

   implicit none
!
! !LOCAL VARIABLES:
   character(len=errorLen) :: msg
   integer, intent(inout)  :: ierr
!
!-----------------------------------------------------------------------

   write(tbnt_log_unit, '(i3,2x,a)') ierr, trim(msg)
#ifndef TBNTonly_bulk_bud
   write(tbnt_log_unit, '(a)')'==========================================='//new_line('a')// &
                              '================ STOP TBNT ================'//new_line('a')// &
                              '==========================================='
   ! close all open output files
   call finalize_tbnt_output(ierr, .false.)
#else
   write(tbnt_log_unit, '(a)')'==========================================='//new_line('a')// &
                              '============ STOP BULK BUDGET ============='//new_line('a')// &
                              '==========================================='
#endif
   close(tbnt_log_unit)
   stop

   end subroutine stop_tbnt
!=======================================================================
   end module mod_etrac_common
!=======================================================================
