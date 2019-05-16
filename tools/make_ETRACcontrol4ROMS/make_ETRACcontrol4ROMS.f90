! this program reads the ROMS grid file, river forcing file and ASCII nice maps for TBNT sources
! and creates standardised input files for the ETRAC application
program make_ETRACcontrol

   use netcdf

   implicit none
   
   integer, parameter :: inUnit = 1
   integer, parameter :: outUnit = 2
   integer, parameter :: nComponents = 3         ! number of directional components (usually 3)

   real(8), parameter :: eps = 1.0e-10
      
   character(len=2), parameter :: com = '!+'
   
   type boxProps
      character(len=20)                  :: name
      character(len= 2)                  :: code
      integer                            :: nCells
      integer, dimension(:), allocatable :: iCells
   end type
   
   integer                                :: n_xi, n_eta, n_sigma, i, j, k, iSep, iRiv, iOut
   integer                                :: iBox, nBoxes, nBoxesFinal, iOff, ncUnit, ncStatus
   integer                                :: iiwet, iiwet2, ios, counter, nExclude, nEdge, nDigits
   integer                                :: nSrc, iMinSrc, iMaxSrc, iSrc, ic
   integer                                :: n_riv, n_times, nCells
   integer, dimension(2)                  :: dims2D
   integer, dimension(:    ), allocatable :: indSrc, openb_src_fr1D, openb_src_to1D, openb_src_comp, openb_src_sign, iRiver
   integer, dimension(:    ), allocatable :: exclude, rivIndex, riv_xi, riv_eta, riv_kMax, compSign, nRiver
   integer, dimension(:,:  ), allocatable :: i_fr, atmos_src, openb_src_fr, openb_src_to, exclude2D, i_surf_bot, kDiv, k_index
   integer, dimension(:,:,:), allocatable :: index1D
   
   real(8)                                :: hinv, cff_r, cff1_r
   real(8), dimension(2)                  :: boxLon, boxDepth
   real(8), dimension(:)    , allocatable :: sc_r, Cs_r, riv_xi_r, riv_eta_r, riv_dir
   real(8), dimension(:,:)  , allocatable :: h, rivQ, riv_shape
   real(8), dimension(:,:,:), allocatable :: z_rho, rivMaskTmp
   real(8), dimension(:,:,:,:), allocatable :: rivMask
   
   real(8), dimension(:,:)  , allocatable :: lat_rho, lon_rho, mask_rho
   
   character(len=10)  :: boxName, boxCode
   character(len=20)  :: formatStr
   character(len=100) :: breakLine
   character(len=100), dimension(:), allocatable :: riverVar, riverNames
   character(len=500) :: line, nice_mapFile
   character(len=500) :: rivName, boxFile(2), boxLabelFile
   character(len=500) :: inDir, outDir, mapDir, boxDir, etracMapDir
   character(len=2), dimension(:,:), allocatable :: boxes2D
   
   logical :: lExist, isSrc, doBoxes, doRivMask
   logical, dimension(:    ), allocatable :: useRiv
   logical, dimension(:,:  ), allocatable :: isAtmos, isOpenb, boxMask
   logical, dimension(:,:,:), allocatable :: openb_dir, isRiv
   
   type(boxProps), dimension(:), allocatable :: boxList

   ! variables for name lists
   integer            :: n_boundary, boxMaxXi, boxMaxEta, riverIndex, river_xi_off, river_eta_off, nRivers
   real(8)            :: boxSepDepth, hcrit
   character(len=50)  :: domain, gridID, riverID, boxID, boxInfoType, riverMask, riverName, n_xi_var, n_eta_var
   character(len=500) :: gridFile, hisFile, riverFile
   logical            :: closedNB, closedSB, closedWB, closedEB, includeAtmos

   namelist / etrac_setup_nml / domain, gridFile, gridID, hisFile, n_xi_var, n_eta_var,        &
                                hcrit, n_boundary, closedNB, closedSB, closedWB, closedEB,     &
                                includeAtmos, riverFile, riverID, river_xi_off, river_eta_off, &
                                nRivers, boxID, boxInfoType, boxSepDepth, boxMaxXi, boxMaxEta
                               
   namelist / river_mask_nml / riverName, riverMask

   namelist / river_list_nml / riverName, riverIndex
   
   !===========================================================================================================
   ! get user setup from setup.nml
   
   open(inUnit, file='setup.nml', action='read', status='old')
   read(inUnit, nml=etrac_setup_nml)
   
   ! check domain identifier
   if (trim(domain)=='') then
      domain = 'NGoMex'
      write(*,'("No domain name provided, use default: ",a)') domain
   end if
   
   ! check name of TBNT file
   if (trim(hisFile)=='') then
      hisFile = '/misc/2/output/grosse/roms854_TBNT-N_NGoMex/tbnt_mch_2000_0001.nc'
      write(*,'("No TBNT file location provided, use default: ",a)') hisFile
   end if
   
   ! check name of grid file
   if (trim(gridFile)=='') then
      gridFile = '/misc/1/input/grosse/ROMS_forcing/roms854/mch_grd.nc'
      write(*,'("No grid file location provided, use default: ",a)') gridFile
   endif
   
   ! check river information
   if (trim(riverFile)=='') then
      riverFile = '/misc/1/input/grosse/ROMS_forcing/roms854/mch_river_1982_2016_Rdet_prov.nc'
      write(*,'("No river file location provided, use default: ",a)') riverFile
   end if
   if (trim(riverID)=='') then
      riverID = 'Mississippi_Atchafalaya'
      write(*,'("No river ID provided, use default: ",a)') riverID
   end if
   if (nRivers>0) then
      doRivMask = .true.
      allocate ( riverNames(nRivers), riverVar(nRivers) )
      do i = 1,nRivers
         read(inUnit, nml=river_mask_nml)
         riverNames(i) = trim(riverName)
         riverVar(i) = trim(riverMask)
      end do
   elseif (nRivers<0) then
      doRivMask = .false.
      nRivers = abs(nRivers)
      allocate ( riverNames(nRivers), iRiver(nRivers) )
      do i = 1,nRivers
         read(inUnit, nml=river_list_nml)
         riverNames(i) = trim(riverName)
         iRiver(i) = riverIndex
      end do
   else
      write(*,'("No valid number of rivers provided (nRivers must be unequal ZERO). Abort program execution")')
      stop
   end if
   
   ! check box information
   if (trim(boxID)=='none'.or.trim(boxID)=='off') then
      write(*,'("Box aggregation disabled.")')
      doBoxes = .false.
   elseif (trim(boxID)=='') then
      boxID = 'NGoMex_4+8boxes'
      write(*,'("No boxID provided, use default: ",a)') boxID
      doBoxes = .true.
   else
      doBoxes = .true.
   end if
   if (trim(boxInfoType)==''.and.doBoxes) then
      boxInfoType = 'depth_lon'
      write(*,'("No boxInfoType provided, use default: ",a)') boxInfoType
   end if
   close(inUnit)
   
   !==============================================================================
   
   ios = 0
 
   inDir  = './Input/'
   outDir = './Output/'
   
   mapDir = trim(inDir)//'nice-maps/'//trim(domain)//'/'
   if (trim(boxInfoType)=='depth_lon') then
     boxDir = trim(inDir)//'box-defs/'
   elseif (trim(boxInfoType)=='nicemap') then
     boxDir = trim(mapDir)//'box-defs/'
   else
      write(*,'("No valid box file type provided. Abort program.")')
      stop
   end if
   etracMapDir = trim(mapDir)//'etrac-maps/' 
   
   breakLine = com
   do i = 3,100
      breakLine(i:i) = '='
   end do
   
   ! read grid information
   call open_nc_file(trim(gridFile), ncUnit)
   call get_nc_dim(ncUnit, trim(n_xi_var), n_xi)
   call get_nc_dim(ncUnit, trim(n_eta_var), n_eta)

   nEdge   = floor(log10(real(n_eta,8))) + 1;
   nDigits = floor(log10(real(n_xi,8)));

   ! allocate fields for grid variables
   !allocate ( lat_rho(n_xi,n_eta), lon_rho(n_xi,n_eta), mask_rho(n_xi,n_eta) )
   dims2D(1) = n_xi
   dims2D(2) = n_eta
   call get_nc_var(ncUnit, 'lon_rho' , 2, (/ n_xi, n_eta /), var2D = lon_rho)
   call get_nc_var(ncUnit, 'lat_rho' , 2, (/ n_xi, n_eta /), var2D = lat_rho)
   call get_nc_var(ncUnit, 'mask_rho', 2, (/ n_xi, n_eta /), var2D = mask_rho)

   ncStatus = NF90_CLOSE(ncUnit)
   call nc_check(ncStatus, ios)
   if (ios/=0) then
      write(*,'("NetCDF-error closing ",a,".")') trim(gridFile)
      stop
   end if

   ! close boundaries if applicable
   if (closedNB) mask_rho(:,n_eta-n_boundary+1:n_eta) = 0.e0
   if (closedSB) mask_rho(:,1:1+n_boundary-1) = 0.e0
   if (closedEB) mask_rho(n_xi-n_boundary+1:n_xi,:) = 0.e0
   if (closedWB) mask_rho(1:1+n_boundary-1,:) = 0.e0

   iiwet2 = nint(sum(mask_rho)) ! number of water columns

   ! get depth-related parameters from history file and calculate depth levels at rho points (z_rho)
   call open_nc_file(trim(hisFile),ncUnit)
   call get_nc_dim(ncUnit, 's_rho', n_sigma)
   
   allocate( sc_r(n_sigma), Cs_r(n_sigma), h(n_xi,n_eta), z_rho(n_xi,n_eta,n_sigma) )

   call get_nc_var(ncUnit, 's_rho', 1, (/n_sigma/)      , var1D = sc_r)
   call get_nc_var(ncUnit, 'Cs_r' , 1, (/n_sigma/)      , var1D = Cs_r)
   call get_nc_var(ncUnit, 'h'    , 2, (/ n_xi, n_eta /), var2D = h)
   
   z_rho = 0.e0
   do j = 1,n_eta
      do i = 1,n_xi
         if (mask_rho(i,j)<0.5) cycle
         hinv = 1.e0/(hcrit+h(i,j))
         do k = 1,n_sigma
            cff_r = hcrit * sc_r(k)
            cff1_r = (cff_r + Cs_r(k)*h(i,j)) * hinv
            z_rho(i,j,k) = h(i,j) * cff1_r
         end do
      end do
   end do

   deallocate ( sc_r, Cs_r )

   ncStatus = NF90_CLOSE(ncUnit)
   call nc_check(ncStatus, ios)
   if (ios/=0) then
      write(*,'("NetCDF-error closing ",a,".")') trim(hisFile)
      stop
   end if

   iiwet = iiwet2 * n_sigma ! number of wet points

   call system("mkdir -p "//trim(outDir))

   ! CREATE MAP FILE WITH MAXIMUM k INDICES
   allocate ( k_index(n_xi,n_eta) )
   where (mask_rho >= 0.5)
      k_index = n_sigma
   else where (mask_rho < 0.5)
      k_index = 0
   end where
   formatStr = '(    i4)'
   write(formatStr(2:5), '(i4)') n_xi
   open(outUnit, file = trim(outDir)//'model_iDep_'//trim(gridID)//'.txt')
   write(outUnit,'(a)')'!+=================================================================================================='
   write(outUnit,'(a)')'!+================= 2D MAP OF MAXIMUM DEPTH INDEX FOR EACH HORIZONTAL GRID POINT  =================='
   write(outUnit,'(a)')'!+=================================================================================================='
   write(outUnit,'(a)')'!+'
   write(outUnit,'(a)')'!+ - THE FIRST LINE CONTAINS AN INTEGER INDICATING WHETHER THE NORTH-SOUTH INDEXING OF THE MAP IS'
   write(outUnit,'(a)')'!+   IN DIRECTION OF THE NetCDF INDEXING (value == 0) OR NOT (else)'
   write(outUnit,'(a)')'!+'
   write(outUnit,'(a)')'!+ - THE FOLLOWING LINES CONTAIN THE DEPTH INDEX MAP'
   write(outUnit,'(a)')'!+'
   write(outUnit,'(a)')'!+ - THE NUMBER OF ROWS MUST EQUAL THE NUMBER OF LATITUDINAL INDICES OF THE MODEL DOMAIN'
   write(outUnit,'(a)')'!+'
   write(outUnit,'(a)')'!+ - THE NUMBER OF ENTRIES PER ROW MUST EQUAL THE NUMBER OF LATITUDINAL INDICES OF THE MODEL DOMAIN'
   write(outUnit,'(a)')'!+'
   write(outUnit,'(a)')'!+=================================================================================================='
   write(outUnit,'(a)')'!+========================================== MAP START ============================================='
   write(outUnit,'(a)')'!+=================================================================================================='
   write(outUnit,'(i4)') 0
   ! write k_index
   do i = 1,n_eta
      write(outUnit,trim(formatStr)) k_index(:,i)
   end do
   deallocate ( k_index )
   close(outUnit)
   
   ! create 3D vector with 1D grid cell index and 1D vectors with surface and bottom cell indices
   allocate ( index1D(n_xi,n_eta,n_sigma) )
   index1D = 0
   counter = 0
   do j = 1,n_eta
      do i = 1,n_xi
         if (mask_rho(i,j)<0.5) cycle
         do k = 1,n_sigma
            counter = counter + 1
            index1D(i,j,k) = counter
         end do
      end do
   end do
   if (counter/=iiwet) then
      write(*,'(a)') 'Grid dimension mismatch. Check consistency of file: '// &
                     new_line('a')//' ==> '//trim(gridFile)//'.'
      stop
   end if
   
   ! calculate neighbour affected by flux at (i,j,k) for each directional component
   ! (physical fluxes are defined on interface in negative index direction)
   allocate ( i_fr(iiwet,nComponents), i_surf_bot(iiwet2,2) )
   counter = 0
   i_fr = 0
   do j = 1,n_eta
      do i = 1,n_xi
         if (mask_rho(i,j)<0.5) cycle
         counter = counter + 1
         i_surf_bot(counter,1) = index1D(i,j,n_sigma)          ! surface cell
         i_surf_bot(counter,2) = index1D(i,j,1)                ! bottom cell
         do k = 1,n_sigma
            if (i>1) i_fr(index1D(i,j,k),1) = index1D(i-1,j,k) ! western neighbour
            if (j>1) i_fr(index1D(i,j,k),2) = index1D(i,j-1,k) ! southern neighbour
            if (k>1) i_fr(index1D(i,j,k),3) = index1D(i,j,k-1) ! lower neighbour
         end do
      end do
   end do
 
   ! CREATE NEIGHBOUR FILE
   open(outUnit, file = trim(outDir)//'model_grid_'//trim(gridID)//'.txt')
   ! write header
   write(outUnit,'(a)')'!+=================================================================================================='
   write(outUnit,'(a)')'!+============== LIST OF NEIGHBOURS FOR EACH GRID CELL AND DIRECTIONAL COMPONENT ==================='
   write(outUnit,'(a)')'!+=================================================================================================='
   write(outUnit,'(a)')'!+'
   write(outUnit,'(a)')'!+ - THE LIST BELOW MUST CONTAIN *N* COMPLETE SETS OF NEIGHBOURS FOR EACH WET POINT OF THE MODEL AND'
   write(outUnit,'(a)')'!+   THE LIST OF SURFACE AND BOTTOM CELL INDICES'
   write(outUnit,'(a)')'!+'
   write(outUnit,'(a)')'!+ - THE FIRST LINE CONTAINS THE NUMBER OF WET POINTS AND THE NUMBER OF DIRECTIONAL COMPONENTS *N*,'
   write(outUnit,'(a)')'!+   e.g. a regular 3D grid has 3 directional components'
   write(outUnit,'(a)')'!+'
   write(outUnit,'(a)')'!+ - BELOW THE FIRST LINE, THE LIST MUST CONTAIN ONE COMPLETE SET OF INFORMATION FOR EACH'
   write(outUnit,'(a)')'!+   DIRECTIONAL COMPONENT'
   write(outUnit,'(a)')'!+'
   write(outUnit,'(a)')'!+ - A COMPLETE SET OF INFORMATION CONSISTS OF THE SIGN OF A FLUX RELATIVE TO ITS LOCATION IN THE'
   write(outUnit,'(a)')'!+   FIRST LINE FOLLOWED BY THE LIST OF NEIGHBOURS IN NEGATIVE DIRECTION OF THE COMPONENT'
   write(outUnit,'(a)')'!+'
   write(outUnit,'(a)')'!+ - EACH SET MUST CONTAIN ONLY *ONE SINGLE NEIGHBOUR IN NEGATIVE DIRECTION* PER EACH GRID POINT'
   write(outUnit,'(a)')'!+'
   write(outUnit,'(a)')'!+ - EACH SET MUST BE A CONTINUOUS SERIES OF INTEGER VALUES'
   write(outUnit,'(a)')'!+'
   write(outUnit,'(a)')'!+ - THERE MUST BE NO SEPARATOR (e.g. blank line, comment etc.) BETWEEN THE DIFFERENT SETS'
   write(outUnit,'(a)')'!+'
   write(outUnit,'(a)')'!+ - IN THE CASE OF AN UNSTRUCTURED GRID, THE NUMBER OF NEIGHBOURS PER CELL MAY VARY, THEN *N* MUST'
   write(outUnit,'(a)')'!+   BE THE MAXIMUM NUMBER OF NEIGHBOURS A CELL HAS WITHIN THE UNDERLYING MODEL'
   write(outUnit,'(a)')'!+'
   write(outUnit,'(a)')'!+ - THE 1st SET CONTAINS THE NEIGHBOURS IN NEGATIVE DIRECTION OF THE 1st DIRECTIONAL COMPONENT,'
   write(outUnit,'(a)')'!+   THE 2nd SET CONTAINS THE NEIGHBOURS IN NEGATIVE DIRECTION OF THE 2nd DIRECTIONAL COMPONENT, ...'
   write(outUnit,'(a)')'!+'
   write(outUnit,'(a)')'!+ - AFTER THE LAST SET OF NEIGHBOURS THE LIST OF SURFACE CELL INDICES AND BOTTOM CELL INDICES MUST'
   write(outUnit,'(a)')'!+   FOLLOW'
   write(outUnit,'(a)')'!+'
   write(outUnit,'(a)')'!+ - THE FIRST LINE OF THE SURFACE LIST CONTAINS THE NUMBER OF SURFACE CELLS, THEN THE CORRESPONDING'
   write(outUnit,'(a)')'!+   INDICES FOLLOW'
   write(outUnit,'(a)')'!+'
   write(outUnit,'(a)')'!+ - THE FIRST LINE OF THE BOTTOM LIST CONTAINS THE NUMBER OF BOTTOM CELLS, THEN THE CORRESPONDING'
   write(outUnit,'(a)')'!+   INDICES FOLLOW'
   write(outUnit,'(a)')'!+'
   write(outUnit,'(a)')'!+=================================================================================================='
   write(outUnit,'(a)')'!+========================================== LIST START ============================================'
   write(outUnit,'(a)')'!+=================================================================================================='

   ! define sign of flux relative to its own grid cell
   ! ("1" means source for mass in cell in case of positive flux, "-1" means sink)
   allocate ( compSign(nComponents) )
   compSign(1) = 1
   compSign(2) = 1
   compSign(3) = 1

   ! write index sets
   write(outUnit, '(5i10)') iiwet, nComponents
   do i = 1,nComponents
      write(outUnit, '(i10)') compSign(i)
      write(outUnit, '(10i10)') i_fr(:,i)
   end do
   do i = 1,2
      write(outUnit, '(i10)') iiwet2
      write(outUnit, '(10i10)') i_surf_bot(:,i)
   end do
   close(outUnit)
   deallocate ( i_fr, i_surf_bot, compSign )

   ! GET OPEN BOUNDARY SOURCES
   allocate ( openb_src_fr(n_xi,n_eta), openb_src_to(n_xi,n_eta), isOpenb(n_xi,n_eta), openb_dir(n_xi,n_eta,4) )
   nice_mapFile = trim(etracMapDir)//'etrac_openb_source_map_from_'//trim(domain)//'.txt'
   call read_integer_map(nice_mapFile, openb_src_fr, 1, n_eta, 1, n_xi, nEdge, nDigits)
   nice_mapFile = trim(etracMapDir)//'etrac_openb_source_map_to_'//trim(domain)//'.txt'
   call read_integer_map(nice_mapFile, openb_src_to, 1, n_eta, 1, n_xi, nEdge, nDigits)
   ! open output file
   open(outUnit, file = trim(outDir)//'etrac_openb_source_'//trim(gridID)//'_'//trim(domain)//'.txt')
   ! write header
   write(outUnit,'(a)')'!+=================================================================================================='
   write(outUnit,'(a)')'!+=============== LIST OF OPEN BOUNDARY SOURCES AND CELLS EXCLUDED FROM CALCULATION ================'
   write(outUnit,'(a)')'!+=================================================================================================='
   write(outUnit,'(a)')'!+'
   write(outUnit,'(a)')'!+ - THE VERY FIRST LINE OF THE LIST BELOW MUST CONTAIN THE NUMBER OF OPEN BOUNDARY SOURCES, THEN'
   write(outUnit,'(a)')'!+   THE DIFFERENT ENTRIES FOR EACH OF THESE SOURCES MUST FOLLOW'
   write(outUnit,'(a)')'!+'
   write(outUnit,'(a)')'!+ - THE LIST BELOW MUST CONTAIN ON COMPLETE SET OF INFORMATION PER OPEN BOUNDARY SOURCE AND THE'
   write(outUnit,'(a)')'!+   LIST OF CELLS EXCLUDED FROM ETRAC CALCULATION DUE TO THE LOCATION OF THE OPEN BOUNDARIES'
   write(outUnit,'(a)')'!+'
   write(outUnit,'(a)')'!+ - OPEN BOUNDARY SOURCES ARE USING THE EXCHANGE FLUXES FROM A SOURCE CELL INTO A TARGET CELL,'
   write(outUnit,'(a)')'!+   THEREFORE EACH POINT OF AN OPEN BOUNDARY SOURCE CONSISTS OF A SOURCE AND TARGET CELL, THE'
   write(outUnit,'(a)')'!+   CORRESPONDIONG DIRECTIONAL COMPONENT AND THE SIGN OF THE FLUX'
   write(outUnit,'(a)')'!+'
   write(outUnit,'(a)')'!+ - THE SIGN OF THE FLUX IS ''+1'' IF THE FLUX FROM THE SOURCE TO THE TARGET CELL IS IN POSITIVE'
   write(outUnit,'(a)')'!+   FLUX DIRECTION OF THE CORRESPONDING DIRECTIONAL COMPONENT, OTHERWISE THE SIGN IS ''-1'''
   write(outUnit,'(a)')'!+'
   write(outUnit,'(a)')'!+ - THE FIRST LINE OF EACH SET CONTAINS THE IDENTIFIER OF THE OPEN BOUNDARY SOURCE (INTEGER) AND'
   write(outUnit,'(a)')'!+   THE NUMBER OF GRID CELLS *N* APPERTAINING TO THIS SOURCE'
   write(outUnit,'(a)')'!+'
   write(outUnit,'(a)')'!+ - BELOW THIS FIRST LINE, THE LIST FOR THE *N* SOURCE CELLS, TARGET CELLS, DIRECTIONAL COMPONENTS'
   write(outUnit,'(a)')'!+   AND SIGNS MUST FOLLOW'
   write(outUnit,'(a)')'!+'
   write(outUnit,'(a)')'!+ - EACH OF THESE LISTS MUST BE A CONSECUTIVE SERIES OF INTEGER VALUES'
   write(outUnit,'(a)')'!+'
   write(outUnit,'(a)')'!+ - THERE MUST BE NO SEPARATOR (e.g. blank line, comment etc.) BETWEEN THE DIFFERENT LISTS AND SETS'
   write(outUnit,'(a)')'!+'
   write(outUnit,'(a)')'!+ - THE 1st SET CONTAINS THE DIFFERENT INFORMATION OF THE 1st OPEN BOUNDARY SOURCE, THE 2nd SET'
   write(outUnit,'(a)')'!+   CONTAINS THE DIFFERENT INFORMATION OF THE 2nd OPEN BOUNDARY SOURCE, ...'
   write(outUnit,'(a)')'!+'
   write(outUnit,'(a)')'!+ - THE LIST OF EXCLUDED CELLS CONSISTS OF A HEADER LINE CONTAINING THE NUMBER OF EXCLUDED CELLS'
   write(outUnit,'(a)')'!+   AND THE SUBSEQUENT CONSECUTIVE LIST OF CELL INDICES TO BE EXCLUDED'
   write(outUnit,'(a)')'!+'
   write(outUnit,'(a)')'!+ - EXCLUDED CELLS ARE DEFINED AS THOSE CELLS LYING OUTSIDE THE DEFINED OPEN BOUNDARIES'
   write(outUnit,'(a)')'!+'
   write(outUnit,'(a)')'!+=================================================================================================='
   write(outUnit,'(a)')'!+========================================== LIST START ============================================'
   write(outUnit,'(a)')'!+=================================================================================================='

   ! determine number of valid open boundary points
   iMinSrc = minval(openb_src_fr, openb_src_fr>0)
   iMaxSrc = maxval(openb_src_fr, openb_src_fr>0)
   
   write(outUnit,'(i10)') iMaxSrc-iMinSrc+1
   do iSrc = iMinSrc,iMaxSrc
      isOpenb = .false.
      openb_dir = .false.
      nSrc = 0
      do j= 1,n_eta
         do i = 1,n_xi
            if (openb_src_fr(i,j)==iSrc) then
               isSrc = .false.
               if (openb_src_to(i+1,j)==iSrc) then ! west => east
                  isOpenb(i,j) = .true.
                  openb_dir(i,j,1) = .true.
                  nSrc = nSrc + n_sigma
                  isSrc = .true.
               end if
               if (openb_src_to(i-1,j)==iSrc) then ! east => west
                  isOpenb(i,j) = .true.
                  openb_dir(i,j,2) = .true.
                  nSrc = nSrc + n_sigma
                  isSrc = .true.
               end if
               if (openb_src_to(i,j-1)==iSrc) then ! north => south
                  isOpenb(i,j) = .true.
                  openb_dir(i,j,3) = .true.
                  nSrc = nSrc + n_sigma
                  isSrc = .true.
               end if
               if (openb_src_to(i,j+1)==iSrc) then ! south => north
                  isOpenb(i,j) = .true.
                  openb_dir(i,j,4) = .true.
                  nSrc = nSrc + n_sigma
                  isSrc = .true.
               end if
               if (.not.isSrc) then
                  write(*,'("Invalid open boundary point: i = ",i3,", j = ",i3)') i, j
                  deallocate ( openb_src_fr, openb_src_to, index1D )
                  stop
               end if
            end if
         end do
      end do
      ! get pairs of open boundary points
      allocate ( openb_src_fr1D(nSrc), openb_src_to1D(nSrc), openb_src_comp(nSrc), openb_src_sign(nSrc) )
      openb_src_fr1D = 0
      openb_src_to1D = 0
      openb_src_comp = 0
      openb_src_sign = 0
      ic = 0
      do j = 1,n_eta
         do i = 1,n_xi
            if (.not.isOpenb(i,j)) cycle
            do k = 1,n_sigma
               if (openb_dir(i,j,1)) then ! west => east
                 ic = ic + 1
                 openb_src_fr1D(ic) = index1D(i  ,j,k)
                 openb_src_to1D(ic) = index1D(i+1,j,k)
                 openb_src_comp(ic) = 1
                 openb_src_sign(ic) = 1
               end if
               if (openb_dir(i,j,2)) then ! east => west
                  ic = ic + 1
                  openb_src_fr1D(ic) = index1D(i  ,j,k)
                  openb_src_to1D(ic) = index1D(i-1,j,k)
                  openb_src_comp(ic) = 1
                  openb_src_sign(ic) = -1
               end if
               if (openb_dir(i,j,3)) then ! north => south
                  ic = ic + 1
                  openb_src_fr1D(ic) = index1D(i,j  ,k)
                  openb_src_to1D(ic) = index1D(i,j-1,k)
                  openb_src_comp(ic) = 2
                  openb_src_sign(ic) = -1
               end if
               if (openb_dir(i,j,4)) then ! south => north
                  ic = ic + 1
                  openb_src_fr1D(ic) = index1D(i,j  ,k)
                  openb_src_to1D(ic) = index1D(i,j+1,k)
                  openb_src_comp(ic) = 2
                  openb_src_sign(ic) = 1
               end if
            end do
         end do
      end do
      if (ic/=nSrc) then
         write(*,'(a,3i10)') 'Mismatch in number of open boundary points: ', ic, nSrc, iSrc
         deallocate ( openb_src_fr1D, openb_src_to1D, openb_src_comp, openb_src_sign  )
         deallocate ( openb_src_fr, openb_src_to, index1D )
         stop
      end if
      
      ! write information to file
      write(outUnit,'(2i10)') iSrc, nSrc
      write(outUnit,'(10i10)') openb_src_fr1D
      write(outUnit,'(10i10)') openb_src_to1D
      write(outUnit,'(10i10)') openb_src_comp
      write(outUnit,'(10i10)') openb_src_sign
      
      deallocate ( openb_src_fr1D, openb_src_to1D, openb_src_comp, openb_src_sign )
   end do

   deallocate ( openb_src_fr, openb_src_to, openb_dir )
   
   ! GET EXCLUDED CELL INDICES
   allocate ( exclude2D(n_xi,n_eta) )
   nice_mapFile = trim(etracMapDir)//'etrac_exclude_map_'//trim(domain)//'.txt'
   call read_integer_map(nice_mapFile, exclude2D, 1, n_eta, 1, n_xi, nEdge, nDigits)
   nExclude = count(exclude2D==1) * n_sigma
   allocate ( exclude(nExclude) )
   exclude = 0
   ic = 0
   do j = 1,n_eta
      do i = 1,n_xi
         if (exclude2D(i,j)/=1) cycle
         do k = 1,n_sigma
            ic = ic + 1
            exclude(ic) = index1D(i,j,k)
         end do
      end do
   end do
   if (ic/=nExclude) then
      write(*,'(a,2i10)') 'Mismatch in number of excluded points: ',ic, nExclude
      deallocate ( index1D, exclude, exclude2D )
      stop
   end if
   write(outUnit,'(i10)') nExclude
   write(outUnit,'(10i10)') exclude
   close(outUnit)
   
   if (includeAtmos) then
      ! GET ATMOSPHERIC BOUNDARY SOURCES
      allocate ( atmos_src(n_xi,n_eta), isAtmos(n_xi,n_eta) )
      nice_mapFile = trim(etracMapDir)//'etrac_atmos_source_map.txt'
      call read_integer_map(nice_mapFile, atmos_src, 1, n_eta, 1, n_xi, nEdge, nDigits)
      ! write file with 1D indices for atmospheric sources
      open(outUnit, file = trim(outDir)//'etrac_atmos_source_'//trim(gridID)//'_'//trim(domain)//'.txt')
      ! write header
      write(outUnit,'(a)')'!+=================================================================================================='
      write(outUnit,'(a)')'!+===================== LIST OF INPUT POINTS FOR ATMOSPHERIC N DEPOSITION =========================='
      write(outUnit,'(a)')'!+=================================================================================================='
      write(outUnit,'(a)')'!+'
      write(outUnit,'(a)')'!+ - THE VERY FIRST LINE OF THE LIST BELOW MUST CONTAIN THE NUMBER OF DIFFERENT ATMOSPHERIC SOURCES'
      write(outUnit,'(a)')'!+'
      write(outUnit,'(a)')'!+ - THE LIST BELOW MUST CONTAIN ONE COMPLETE SET OF INPUT POINTS PER ATMOSPHERIC SOURCE'
      write(outUnit,'(a)')'!+'
      write(outUnit,'(a)')'!+ - THE FIRST LINE OF EACH SET CONTAINS THE IDENTIFIER OF THE ATMOSPHERIC SOURCE (INTEGER) AND THE'
      write(outUnit,'(a)')'!+   NUMBER OF GRID CELLS *N* ASSOCIATED WITH THIS SOURCE'
      write(outUnit,'(a)')'!+'
      write(outUnit,'(a)')'!+ - BELOW THIS FIRST LINE, THE *N* MEMBERS OF THE CORRESPONDING SOURCE MUST BE LISTED'
      write(outUnit,'(a)')'!+'
      write(outUnit,'(a)')'!+ - EACH SET MUST BE A CONSECUTIVE SERIES OF INTEGER VALUES'
      write(outUnit,'(a)')'!+'
      write(outUnit,'(a)')'!+ - THERE MUST BE NO SEPARATOR (e.g. blank line, comment etc.) BETWEEN THE DIFFERENT SETS'
      write(outUnit,'(a)')'!+'
      write(outUnit,'(a)')'!+=================================================================================================='
      write(outUnit,'(a)')'!+========================================== LIST START ============================================'
      write(outUnit,'(a)')'!+=================================================================================================='
      iMinSrc = minval(atmos_src, atmos_src>0)
      iMaxSrc = maxval(atmos_src, atmos_src>0)
      write(outUnit,'(i10)') iMaxSrc-iMinSrc+1
      do iSrc = iMinSrc,iMaxSrc
         isAtmos = .false.
         nSrc = count(atmos_src==iSrc)
         do j = 1,n_eta
            do i = 1,n_xi
               if (atmos_src(i,j)==iSrc) then
                  if (exclude2D(i,j)==0) then
                     isAtmos(i,j) = .true.
                  else
                     nSrc = nSrc - 1
                  end if
               end if
            end do
         end do
         allocate ( indSrc(nSrc) )
         indSrc = pack(index1D(:,:,n_sigma), isAtmos)
         write(outUnit, '(2i10)') iSrc, nSrc
         write(outUnit,'(10i10)') indSrc
         deallocate ( indSrc )
      end do
      close(outUnit)
      deallocate ( atmos_src, isAtmos )
   end if
   
   deallocate ( exclude, exclude2D )
   

   ! CREATE 1D RIVER FILE
   ! get information from NetCDF file
   call open_nc_file(trim(riverFile),ncUnit)
   call get_nc_dim(ncUnit, 'river'     , n_riv)
   call get_nc_dim(ncUnit, 'river_time', n_times)
   
   !allocate( riv_xi_r(n_riv), riv_eta_r(n_riv), riv_dir(n_riv), rivQ(n_riv,n_times), riv_shape(n_riv,n_sigma) )
   call get_nc_var(ncUnit, 'river_Xposition', 1, (/ n_riv /)                  , var1D = riv_xi_r)
   call get_nc_var(ncUnit, 'river_Eposition', 1, (/ n_riv /)                  , var1D = riv_eta_r)
   call get_nc_var(ncUnit, 'river_direction', 1, (/ n_riv /)                  , var1D = riv_dir)
   call get_nc_var(ncUnit, 'river_transport', 2, (/ n_riv, n_times /)         , var2D = rivQ)
   call get_nc_var(ncUnit, 'river_Vshape'   , 2, (/ n_riv, n_sigma /)         , var2D = riv_shape)
   if (doRivMask) then
      allocate ( rivMask(n_riv,n_sigma,n_times,nRivers) )
      do i = 1,nRivers
         call get_nc_var(ncUnit, trim(riverVar(i)), 3, (/ n_riv, n_sigma, n_times /), var3D = rivMaskTmp)
         rivMask(:,:,:,i) = rivMaskTmp
         deallocate ( rivMaskTmp )
      end do
   else
      
   end if

   ncStatus = NF90_CLOSE(ncUnit)
   call nc_check(ncStatus, ios)
   if (ios/=0) then
      write(*,'("NetCDF-error closing ",a,".")') trim(riverFile)
      stop
   end if

   allocate ( riv_xi(n_riv), riv_eta(n_riv), riv_kMax(n_riv), useRiv(n_riv) )
   riv_xi  = nint(riv_xi_r) - river_xi_off
   riv_eta = nint(riv_eta_r) - river_eta_off
   useRiv = .true.
   riv_kMax = 0
   
   ! get correct position of river input rho (i,j)-point and get maximum input depth of river
   ! What needs to be considered?
   ! 1. output data have an index offset of +1 in Xi and Eta direction compared to locations from river file (due to halo lines)
   ! 2. negative discharge (rivQ) implies reduction of actual grid index by 1
   ! 3. river direction (riv_dir) defines whether river goes in Xi (riv_dir=0) or Eta direction (riv_dir=1)
   ! NOTE: (1) and (2) balance out in case of rivQ<0
   do i = 1,n_riv
      if (riv_dir(i)<0.5e0) then ! Xi direction
         riv_eta(i) = riv_eta(i) + 1
         if (sum(rivQ(i,:)) > 0.0e0) riv_xi(i) = riv_xi(i) + 1
      else                       ! Eta direction
         riv_xi(i) = riv_xi(i) + 1
         if (sum(rivQ(i,:)) > 0.0e0) riv_eta(i) = riv_eta(i) + 1
      end if
      if (riv_xi(i)<1.or.riv_xi(i)>n_xi.or.riv_eta(i)<1.or.riv_eta(i)>n_eta) then
         useRiv(i) = .false.
         cycle
      end if
      riv_kMax(i) = count(riv_shape(i,:)>eps)
   end do
   
   if (doRivMask) then
      allocate ( isRiv(n_riv,n_sigma,nRivers) )
      isRiv = .false.
      do i = 1,nRivers
         if (.not.useRiv(i)) cycle
         where (rivMask(:,:,1,i) > 0.5) isRiv(:,:,i) = .true.
      end do
      do i = 1,n_sigma
         if (count(isRiv(:,i,:)) /= count(useRiv)) then
            write(*,'(a)') 'Number of rivers does not match river flags.'
            stop
         end if
      end do
      deallocate ( rivMask )
   end if

   deallocate ( riv_xi_r, riv_eta_r, riv_dir, rivQ, riv_shape )
   
   

   formatStr = '(i3,XXi10,a,a)'
   write(formatStr(5:6),'(i2)') n_sigma
   open(outUnit, file = trim(outDir)//'model_rivers_'//trim(gridID)//'_'//trim(riverID)//'.txt')
   ! write header to file
   write(outUnit,'(a)')'!+=================================================================================================='
   write(outUnit,'(a)')'!+================================ LIST OF INPUT POINTS FOR RIVERS ================================='
   write(outUnit,'(a)')'!+=================================================================================================='
   write(outUnit,'(a)')'!+'
   write(outUnit,'(a)')'!+ - THE LIST BELOW MUST CONTAIN THE INFORMATION ABOUT ALL RIVERS ENTERING THE MODEL DOMAIN'
   write(outUnit,'(a)')'!+'
   write(outUnit,'(a)')'!+ - THE FIRST LINE CONTAINS THE NUMBER OF RIVERS AND THE MAXIMUM NUMBER OF GRID CELLS A SINGLE'
   write(outUnit,'(a)')'!+   RIVER CAN ENTER'
   write(outUnit,'(a)')'!+'
   write(outUnit,'(a)')'!+ - BELOW THIS FIRST LINE, THE LIST OF RIVERS MUST FOLLOW'
   write(outUnit,'(a)')'!+'
   write(outUnit,'(a)')'!+ - EACH LINE CONTAINS THE INFORMATION FOR ONE RIVER CONSISTING OF, THE RIVER NUMBER, THE INDICES'
   write(outUnit,'(a)')'!+   OF THE GRID CELLS IN WHICH THE RIVER ENTERS, AND THE RIVER NAME'
   write(outUnit,'(a)')'!+'
   write(outUnit,'(a)')'!+ - A DASH (-) MUST BE USED TO SEPARATE THE RIVER NAME FROM THE INDEX LIST'
   write(outUnit,'(a)')'!+'
   write(outUnit,'(a)')'!+ - THE RIVER NAME MUST BE WRITTEN IN CAPITAL LETTERS'
   write(outUnit,'(a)')'!+'
   write(outUnit,'(a)')'!+=================================================================================================='
   write(outUnit,'(a)')'!+========================================== LIST START ============================================'
   write(outUnit,'(a)')'!+=================================================================================================='
   write(outUnit,'(2i10)') count(useRiv), n_sigma
   
   ! write river information
   allocate ( rivIndex(n_sigma) )
   
   do i = 1,nRivers
      call capitalise_string(riverNames(i))
      call replace_char(riverNames(i),' ','_')
   end do
   
   iOut = 0
   
   if (doRivMask) then
      allocate ( nRiver(nRivers) )
      nRiver = 0
      do i = 1,n_riv
         if (.not.useRiv(i)) cycle
         rivIndex = 0
         rivIndex(1:riv_kMax(i)) = index1D(riv_xi(i),riv_eta(i),n_sigma-riv_kMax(i)+1:n_sigma)
         if (count(rivIndex>0)==0) then
            write(*,'(a)') 'Error: river is on land.'
            stop
         end if
         do iRiv = 1,nRivers
            if (isRiv(i,1,iRiv)) then
               nRiver(iRiv) = nRiver(iRiv) + 1
               write(rivName,'(a,"-",i2.2)') trim(riverNames(iRiv)), nRiver(iRiv)
               exit
            end if
         end do
         iOut = iOut+1
         write(outUnit, trim(formatStr)) iOut, rivIndex, ' - ', trim(rivName)
      end do
      deallocate ( isRiv )
   else
      do i = 1,n_riv
         if (.not.useRiv(i)) cycle
         rivIndex = 0
         iRiv = iRiver(i)
         rivIndex(1:riv_kMax(iRiv)) = index1D(riv_xi(iRiv),riv_eta(iRiv),n_sigma-riv_kMax(iRiv)+1:n_sigma)
         if (count(rivIndex>0)==0) then
            write(*,'(a)') 'Error: river is on land.'
            stop
         end if
         iOut = iOut+1
         write(outUnit, trim(formatStr)) iOut, rivIndex, ' - ', trim(riverNames(i))
      end do
   end if
   
   close(outUnit)
   deallocate (rivIndex, riverNames, riv_xi, riv_eta, riv_kMax, useRiv )
   
   
   
   
   ! finish script execution if no box treatment
   if (.not.doBoxes) stop
   
   ! CREATE BOX SET-UP FILE
   call capitalise_string(boxInfoType)
   if (trim(boxInfoType)=='NICEMAP') THEN
      boxLabelFile = trim(boxDir)//'box-table_'//trim(gridID)//'_'//trim(boxID)//'.txt'
      boxFile(1) = trim(boxDir)//'nice-map_'//trim(gridID)//'_'//trim(boxID)//'_top.txt'
      boxFile(2) = trim(boxDir)//'nice-map_'//trim(gridID)//'_'//trim(boxID)//'_bot.txt'
      inquire(file=trim(boxFile(2)), exist=lExist)
      allocate ( kDiv(n_xi,n_eta) )
      kDiv = n_sigma
      if (lExist.and.boxSepDepth>0.0e0) then ! subdivision between surface and bottom boxes at "boxSepDepth"
         do j = 1,n_eta
            do i = 1,n_xi
               kDiv(i,j) = count(z_rho(i,j,:)<=boxSepDepth)
            end do
         end do
      end if
      ! get number of boxes + box names and codes
      call system('wc -l '//trim(boxLabelFile)//' > scratch.txt')
      open(inUnit,file='scratch.txt')
      read(inUnit,*) nBoxes
      close(inUnit)
      call system('rm scratch.txt')
      nBoxes = nBoxes - 1
      allocate ( boxList(nBoxes) )
      boxList%nCells = 0
      open(inUnit, file = trim(boxLabelFile), action = 'read', status = 'old')
      call skipLines(inUnit,1)
      do i = 1,nBoxes
         read(inUnit, '(2a10)') boxName, boxCode
         boxList(i)%name = trim(adjustl(boxName))
         boxList(i)%code = trim(adjustl(boxCode))
      end do
      close(inUnit)
      ! get box cells
      allocate ( boxes2D(n_xi,n_eta) )
      ! read upper boxes
      call read_character_map(boxFile(1), boxes2D, 1, n_eta, 1, n_xi, nEdge)
      nBoxesFinal = nBoxes
      do iBox = 1,nBoxes
         do j = 1,n_eta
            if (j>boxMaxEta) cycle
            do i = 1,n_xi
               if (boxes2D(i,j)=='##'.or.i>boxMaxXi) cycle
               if (boxes2D(i,j)==boxList(iBox)%code) then
                  boxList(iBox)%nCells = boxList(iBox)%nCells + min(n_sigma,kDiv(i,j))
               end if
            end do
         end do
         if (boxList(iBox)%nCells==0) then
            nBoxesFinal = nBoxesFinal - 1
            cycle
         end if
         allocate ( boxList(iBox)%iCells(boxList(iBox)%nCells) )
         iOff = 0
         do j = 1,n_eta
            if (j>boxMaxEta) cycle
            do i = 1,n_xi
               if (boxes2D(i,j)=='##'.or.i>boxMaxXi) cycle
               if (boxes2D(i,j)==boxList(iBox)%code) then
                  do k = 1,min(n_sigma,kDiv(i,j))
                     boxList(iBox)%iCells(iOff+k) = index1D(i,j,k)
                  end do
                  iOff = iOff + min(n_sigma,kDiv(i,j))
               end if
            end do
         end do
      end do
      ! read lower boxes
      if (lExist) then
         call read_character_map(boxFile(2), boxes2D, 1, n_eta, 1, n_xi, nEdge)
         do iBox = 1,nBoxes
            if (boxList(iBox)%nCells>0) cycle
            do j = 1,n_eta
            if (j>boxMaxEta) cycle
               do i = 1,n_xi
                  if (boxes2D(i,j)=='##'.or.i>boxMaxXi) cycle
                  if (boxes2D(i,j)==boxList(iBox)%code) then
                     boxList(iBox)%nCells = boxList(iBox)%nCells + n_sigma - kDiv(i,j)
                  end if
               end do
            end do
            if (boxList(iBox)%nCells==0) then
               nBoxesFinal = nBoxesFinal - 1
               cycle
            end if
            allocate ( boxList(iBox)%iCells(boxList(iBox)%nCells) )
            iOff = 0
            do j = 1,n_eta
            if (j>boxMaxEta) cycle
               do i = 1,n_xi
                  if (boxes2D(i,j)=='##'.or.i>boxMaxXi) cycle
                  if (boxes2D(i,j)==boxList(iBox)%code) then
                     do k = kDiv(i,j)+1,n_sigma
                        boxList(iBox)%iCells(iOff+k-kDiv(i,j)) = index1D(i,j,k)
                     end do
                     iOff = iOff + n_sigma - kDiv(i,j)
                  end if
               end do
            end do
         end do
      end if
   else if (boxInfoType=='DEPTH_LON') then
      boxFile(1) = trim(boxDir)//trim(boxID)//'.txt'
      call system('wc -l '//trim(boxFile(1))//' > scratch.txt')
      open(inUnit,file='scratch.txt')
      read(inUnit,*) nBoxes
      close(inUnit)
      call system('rm scratch.txt')
      allocate ( boxList(nBoxes) )
      boxLon = 1000.e0
      boxDepth = -1.e0
      open(inUnit, file=trim(boxFile(1)), action = 'read', status = 'old')
      nBoxesFinal = 0
      allocate ( boxMask(n_xi,n_eta) )
      do iBox = 1,nBoxes
         boxList(iBox)%code = '##'
         boxList(iBox)%nCells = 0
         read(inUnit, '(a)') line
         if (len_trim(line)==0) cycle
         iSep = index(line,':')
         if (iSep==0) then
            write(*,'("Error reading box file. Line separator (:) not found. Abort program.")')
            deallocate ( h, boxList )
            close(inUnit)
            stop
         end if
         if (iSep>1) boxName = trim(line(1:iSep-1))
         read(line(iSep+1:len_trim(line)),*) boxDepth, boxLon
         if (count(boxLon==1000.e0)>0.or.count(boxDepth<0.e0)>0) then
            write(*,'("Error reading box file. Invalid data entries. Abort program.")')
            deallocate ( h, boxList, boxMask )
            close(inUnit)
            stop
         end if
         boxMask = lon_rho>=boxLon(1).and.lon_rho<boxLon(2).and.h>boxDepth(1).and.h<=boxDepth(2).and.mask_rho>0.5
         if (boxMaxXi <=n_xi)  boxMask(min(n_xi,boxMaxXi+1):n_xi,:)    = .false.
         if (boxMaxEta<=n_eta) boxMask(:,min(n_eta,boxMaxEta+1):n_eta) = .false.
         nCells = count(boxMask)*n_sigma
         if (nCells==0) cycle
         nBoxesFinal = nBoxesFinal + 1
         boxList(iBox)%nCells = nCells
         if (iSep>1) then
            boxList(iBox)%name = trim(boxName)
         else
            write(boxList(iBox)%name, '("box"i3.3)') nBoxesFinal
         endif
         allocate ( boxList(iBox)%iCells(nCells) )
         ic = 0
         do j = 1,n_eta
            do i = 1,n_xi
               if (.not.boxMask(i,j)) cycle
               do k = 1,n_sigma
                  ic = ic + 1
                  boxList(iBox)%iCells(ic) = index1D(i,j,k)
               end do
            end do
         end do
         if (ic.ne.nCells) then
            write(*,'("Mismatch in grid cell indices for box: ",a)') boxList(iBox)%name
            deallocate ( h, boxMask )
            close(inUnit)
            stop
         end if
      end do
      deallocate ( boxMask )
      close(inUnit)
   end if
      
   open(outUnit, file=trim(outdir)//'target_areas_'//trim(gridID)//'_'//trim(boxID)//'.txt')
   ! write header to file
   write(outUnit,'(a)')'!+=================================================================================================='
   write(outUnit,'(a)')'!+====================================== LIST OF TARGET AREAS ======================================'
   write(outUnit,'(a)')'!+=================================================================================================='
   write(outUnit,'(a)')'!+'
   write(outUnit,'(a)')'!+ - THE LIST BELOW MUST CONTAIN ONE COMPLETE SET OF INFORMATION PER TARGET AREA'
   write(outUnit,'(a)')'!+'
   write(outUnit,'(a)')'!+ - A TARGET AREA IS DEFINED BY ITS NAME, NUMBER OF WET CELLS ASSOCIATED WITH THE AREA AND THE LIST'
   write(outUnit,'(a)')'!+   OF CELL INDICES'
   write(outUnit,'(a)')'!+'
   write(outUnit,'(a)')'!+ - THE FIRST LINE OF THE LIST CONTAINS THE NUMBER OF TARGET AREAS'
   write(outUnit,'(a)')'!+'
   write(outUnit,'(a)')'!+ - THE FIRST LINE OF EACH SET CONTAINS THE NAME OF THE TARGET AREA AND THE NUMBER OF GRID CELLS'
   write(outUnit,'(a)')'!+   *N* ASSOCIATED WITH THIS AREA'
   write(outUnit,'(a)')'!+'
   write(outUnit,'(a)')'!+ - BELOW THIS FIRST LINE, THE LIST OF THE *N* CELLS OF THE TARGET AREA MUST FOLLOW'
   write(outUnit,'(a)')'!+'
   write(outUnit,'(a)')'!+ - THIS LIST MUST BE A CONSECUTIVE SERIES OF INTEGER VALUES'
   write(outUnit,'(a)')'!+'
   write(outUnit,'(a)')'!+ - THERE MUST BE NO SEPARATOR (e.g. blank line, comment etc.) BETWEEN THE DIFFERENT LISTS AND SETS'
   write(outUnit,'(a)')'!+'
   write(outUnit,'(a)')'!+ - THE 1st SET CONTAINS THE DIFFERENT INFORMATION ON THE 1st TARGET AREA, THE 2nd SET CONTAINS THE'
   write(outUnit,'(a)')'!+   DIFFERENT INFORMATION ON THE 2nd TARGET AREA, ...'
   write(outUnit,'(a)')'!+'
   write(outUnit,'(a)')'!+=================================================================================================='
   write(outUnit,'(a)')'!+========================================== LIST START ============================================'
   write(outUnit,'(a)')'!+=================================================================================================='
   write(outUnit,'(i10)') nBoxesFinal
   ! write box list
   do iBox = 1,nBoxes
      if (boxList(iBox)%nCells==0) cycle
      write(outUnit,'(a10,i10)') boxList(iBox)%name, boxList(iBox)%nCells
      write(outUnit,'(10i10)') boxList(iBox)%iCells
   end do
   close(outUnit)
   
   deallocate ( index1D, boxList )
   
contains
   
! ================================ SUBROUTINES ===================================
! !INTERFACE:
   subroutine read_integer_map(nice_map, nice_array, jStart, jEnd, iStart, iEnd, nEdge, nDigits)
!
! !DESCRIPTION:
! reads nice map containing integers
!
   implicit none
!
! !OUTPUT PARAMETERS:
   integer, intent(inout) :: nice_array(n_xi,n_eta)
!
! !LOCAL VARIABLES:
   integer, parameter :: nice_map_unit = 101

   integer           , intent(in) :: iStart, iEnd, jStart, jEnd, nEdge, nDigits
   character(len=150), intent(in) :: nice_map
   integer                        :: i, j, i1, i2
   character(len=4)               :: formatStr
   character(len=nEdge)           :: str
   character(len=nDigits*n_xi)    :: nice_str
   logical                        :: lExist
!
!---------------------------------------------------------------------------------
   
   lExist = .false.
   inquire(file=trim(nice_map), exist=lExist)
   if (.not.lExist) then
      write(*,'("File does not exist: ",a,".")') trim(nice_map)
      stop
   end if
   
   formatStr = '(iX)'
   write(formatStr(3:3),'(i1)') nDigits
   nice_array = -1
   
   open(nice_map_unit, file = trim(nice_map), action = 'read', status = 'old')
   call skiplines(nice_map_unit,2)
   do j = jStart,jEnd
      read(nice_map_unit,'(3(a))')str,nice_str,str
      do i = iStart,iEnd
         i1 = (i-1)*nDigits + 1
         i2 = i*nDigits
         if ((nice_str(i1:i1)=='.'.and.nice_str(i2:i2)=='.').or. &
             (nice_str(i1:i1)=='#'.and.nice_str(i2:i2)=='#')) cycle
         read(nice_str(i1:i2),formatStr) nice_array(i,n_eta-j+1)
      end do
   end do
   close(nice_map_unit)

   end subroutine read_integer_map
! ================================================================================
! !INTERFACE:
   subroutine read_character_map(nice_map, nice_array, jStart, jEnd, iStart, iEnd, nEdge)
!
! !DESCRIPTION:
! reads nice map containing characters
!
   implicit none
!
! !OUTPUT PARAMETERS:
   character(len=2), intent(inout) :: nice_array(n_xi,n_eta)
!
! !LOCAL VARIABLES:
   integer, parameter :: nice_map_unit = 101

   integer           , intent(in) :: iStart, iEnd, jStart, jEnd, nEdge
   character(len=150), intent(in) :: nice_map
   integer                        :: i, j, i1, i2
   character(len=nEdge)           :: str
   character(len=2*n_xi)          :: nice_str
   logical                        :: lExist
!
!---------------------------------------------------------------------------------
   
   lExist = .false.
   inquire(file=trim(nice_map), exist=lExist)
   if (.not.lExist) then
      write(*,'("File does not exist: ",a,".")') trim(nice_map)
      stop
   end if
   
   nice_array = '##'
   
   open(nice_map_unit, file = trim(nice_map), action = 'read', status = 'old')
   call skiplines(nice_map_unit,2)
   do j = jStart,jEnd
      read(nice_map_unit,'(3(a))')str,nice_str,str
      do i = iStart+2*(iStart-1),iEnd
         i1 = 2*i-1
         i2 = 2*i
         if (nice_str(i1:i2)=='..'.or.nice_str(i1:i2)=='##') cycle
         nice_array(i,n_eta-j+1)(1:2) = nice_str(i1:i2)
      end do
   end do
   close(nice_map_unit)

   end subroutine read_character_map
! ================================================================================
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
   integer, intent(in) :: fUnit, nLines
   integer             :: n
   character(len=200)  :: str
!
!---------------------------------------------------------------------------------

   do n =1,nLines
      read(fUnit,'(a)')str
   end do

   end subroutine skiplines
! ================================================================================
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
!---------------------------------------------------------------------------------

   nFirst = scan(trim(str),ABCs)
   if (nFirst<=0) return
   nLast = scan(trim(str),ABCs,.true.)
   do n = nFirst,nLast
      is = scan(ABCs, str(n:n))
      if (is>0) str(n:n) = ABCc(is:is)
   end do

   end subroutine capitalise_string
! ================================================================================
!
! !INTERFACE:
   subroutine replace_char(str,charOld,charNew)
!
! !DESCRIPTION:
! replace all occurences of character 'charOld' in string 'str' by character 'charNew'
!
   implicit none
!
! !OUTPUT PARAMETERS:
   character(len=*), intent(inout) :: str
   character(len=1), intent(in   ) :: charOld, charNew
!
! !LOCAL VARIABLES:
   
   integer :: n, nFirst, nLast
!
!---------------------------------------------------------------------------------

   nFirst = scan(trim(str),charOld)
   if (nFirst<=0) return
   nLast = scan(trim(str),charOld,.true.)
   do n = nFirst,nLast
      if (str(n:n)==charOld) str(n:n) = charNew
   end do

   end subroutine replace_char
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
   integer, intent(in   ) :: ncStatusIn
   integer, intent(inout) :: iret
!   
!-----------------------------------------------------------------------
!
   if (ncStatusIn /= nf90_noerr) iret = 1

   end subroutine nc_check
!=======================================================================
!
! !INTERFACE:
   subroutine open_nc_file(ncFile, ncUnit)
!
! !DESCRIPTION:
! read dimension from NetCDF file
!
   implicit none
!
! !OUTPUT PARAMETERS:
   character(len=*), intent(in)  :: ncFile
   integer         , intent(out) :: ncUnit
!
! !LOCAL VARIABLES:
   integer :: ncStatus, ios
   logical :: lExist
!
!---------------------------------------------------------------------------------

   lExist = .false.
   ios = 0

   inquire(file=trim(ncFile), exist=lExist)
   if (.not.lExist) then
      write(*,'(a)') 'File does not exist: '//trim(ncFile)//'.'
      stop
   end if
   ncStatus = NF90_OPEN(trim(ncFile), NF90_NOWRITE, ncUnit)
   call nc_check(ncStatus, ios)
   if (ios/=0) then
      write(*,'("NetCDF-error opening ",a,".")') trim(ncFile)
      stop
   end if

   end subroutine open_nc_file
!=======================================================================
!
! !INTERFACE:
   subroutine get_nc_dim(fileUnit, dimName, dimVal)
!
! !DESCRIPTION:
! read dimension from NetCDF file
!
   implicit none
!
! !OUTPUT PARAMETERS:
   integer         , intent(in)  :: fileUnit
   character(len=*), intent(in)  :: dimName
   integer         , intent(out) :: dimVal
!
! !LOCAL VARIABLES:
   integer :: ios, ncStatus, dimID
!
!---------------------------------------------------------------------------------

   ios = 0

   ncStatus = NF90_INQ_DIMID(fileUnit, trim(dimName), dimID)
   call nc_check(ncStatus, ios)
   if (ios/=0) then
      write(*,'("NetCDF-error inquiring dimension ID for ",a,".")') trim(dimName)
      stop
   end if

   ncStatus = NF90_INQUIRE_DIMENSION(fileUnit, dimID, len=dimVal)
   call nc_check(ncStatus, ios)
   if (ios/=0) then
      write(*,'("NetCDF-error inquiring ",a,".")') trim(dimName)
      stop
   end if

   end subroutine get_nc_dim
!=======================================================================
!
! !INTERFACE:
   subroutine get_nc_var(fileUnit, varName, nDims, dimSize, var1D, var2D, var3D)
!
! !DESCRIPTION:
! read variable from NetCDF file
!
   implicit none
!
! !OUTPUT PARAMETERS:
   integer          , intent(in)    :: fileUnit, nDims
   integer          , intent(in)    :: dimSize(nDims)
   character(len=*) , intent(in)    :: varName
   real(8), optional, allocatable, intent(out) :: var1D(:)
   real(8), optional, allocatable, intent(out) :: var2D(:,:)
   real(8), optional, allocatable, intent(out) :: var3D(:,:,:)
!
! !LOCAL VARIABLES:
   integer :: ios, ncStatus, varID
!
!---------------------------------------------------------------------------------

   ios = 0
   
   ncStatus = NF90_INQ_VARID(fileUnit, trim(varName), varID)
   call nc_check(ncStatus, ios)
   if (ios/=0) then
      write(*,'("NetCDF-error inquiring variable ID for ",a,".")') trim(varName)
      stop
   end if
   if (nDims==1.and.present(var1D)) then
      allocate ( var1D(dimSize(1)) )
      ncStatus = NF90_GET_VAR(fileUnit, varID, var1D)
   else if (nDims==2.and.present(var2D)) then
      allocate ( var2D(dimSize(1), dimSize(2)) )
      ncStatus = NF90_GET_VAR(fileUnit, varID, var2D)
   else if (nDims==3.and.present(var3D)) then
      allocate ( var3D(dimSize(1), dimSize(2), dimSize(3)) )
      ncStatus = NF90_GET_VAR(fileUnit, varID, var3D)
   else
      write(*,'("NetCDF-error reading ",a,". Invalid dimensions.")') trim(varName)
      stop
   end if

   call nc_check(ncStatus, ios)
   if (ios/=0) then
      write(*,'("NetCDF-error reading ",a,". ",a)') trim(varName), NF90_STRERROR(ncStatus)
      stop
   end if

   end subroutine get_nc_var
!=======================================================================
   
end program make_ETRACcontrol
