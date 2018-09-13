! this program reads the ROMS grid file and creates an ASCII land-sea mask
! used as template for the creation of the ETRAC input control files
program make_nicemap4ROMS

   use netcdf

   implicit none
   
   integer, parameter :: inUnit = 1
   integer, parameter :: outUnit = 2
   
   integer :: i, j, ncUnit, ios, i1, i2, nEdge, nDigits
   integer :: ncStatus, n_xi, n_eta, n_sig, value, ioff, o_j
   
   real(8)              :: eps, divider
   real(8), allocatable :: mask_rho(:,:), mask_excl(:,:), var_excl(:,:,:)
   
   character(len=   4) :: frmts
   character(len=   6) :: frmti
   character(len=  20) :: exclStr
   character(len= 150) :: outDir, mapDir
   character(len=1000), allocatable :: headerlines(:), mapLines(:)
   
   logical :: lExist, doExclude

   ! variables for name lists
   character(len=50)  :: domain, excludeVar, n_xi_var, n_eta_var, n_sig_var
   character(len=500) :: gridFile, excludeFile
   integer            :: n_boundary
   logical            :: closedNB, closedSB, closedWB, closedEB

   namelist / grid_nml / domain, gridFile, n_xi_var, n_eta_var, n_sig_var, n_boundary,   &
                         closedNB, closedSB, closedWB, closedEB, excludeFile, excludeVar
   
   !==============================================================================
   ! get user setup from setup.nml
   
   open(inUnit, file='grid.nml', action='read', status='old')
   read(inUnit, nml=grid_nml)
   if (trim(domain)=='') then
      domain = 'NGoMex'
      write(*,'("No domain name provided, use default: ",a)') domain
   end if
   if (trim(gridFile)=='') then
      gridFile = '/misc/1/input/grosse/ROMS_forcing/roms854/mch_grd.nc'
      write(*,'("No grid file location provided, use default: ",a)') gridFile
   endif
   ! check information on exluded cells
   if (trim(excludeFile)=='none'.or.trim(excludeFile)=='off') then
      write(*,'("Cell exclusion disabled.")')
      doExclude = .false.
   elseif (trim(excludeFile)=='') then
      write(*,'("No cell exclusion file provided, disable cell exclusion.")')
      doExclude = .false.
   else
      doExclude = .true.
   end if
   if (trim(excludeVar)==''.and.doExclude) then
      excludeVar = 'NO3_NudgeCoef'
      write(*,'("No excludeVar provided, use default: ",a)') excludeVar
   end if
   
   close(inUnit)
   
   !==============================================================================
   
   ios = 0
   lExist = .true.
 
   outDir  = './Input/'
   mapDir = trim(outDir)//'nice-maps/'//trim(domain)//'/etrac-maps/'

   inquire(file=trim(mapDir)//'.', exist=lExist)
   if (.not.lExist) call system('mkdir -p '//trim(mapDir))

   ! read grid information
   inquire(file=trim(gridFile), exist=lExist)
   if (.not.lExist) then
      write(*,'(a)') 'File does not exist: '//trim(gridFile)//'.'
      stop
   end if
   ncStatus = NF90_OPEN(trim(gridFile), NF90_NOWRITE, ncUnit)
   call nc_check(ncStatus, ios)
   if (ios/=0) then
      write(*,'("NetCDF-error opening ",a,".")') trim(gridFile)
      stop
   end if

   call get_nc_dim(ncUnit, trim(n_xi_var), n_xi)
   call get_nc_dim(ncUnit, trim(n_eta_var), n_eta)
   
   nEdge   = floor(log10(real(n_eta,8))) + 1;
   nDigits = floor(log10(real(n_xi,8)));
   divider = 10**real(nDigits,8)
   
   ! allocate and read land-sea mask
   allocate ( mask_rho(n_xi,n_eta) )
   call get_nc_var(ncUnit, 'mask_rho', (/ n_xi, n_eta /), mask_rho)
   
   ! close grid file
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
   
   ! initialize nice map array
   allocate ( headerLines(2), mapLines(n_eta) )
   headerLines = repeat('_', n_xi*nDigits+2*nEdge)
   mapLines = repeat('_', nEdge)//repeat('#', n_xi*nDigits)//repeat('_', nEdge)

   ! fill array with actual values
   frmti = '(iX.X)'
   write(frmti(3:3), '(i1)') nDigits
   write(frmti(5:5), '(i1)') nDigits

   do i = 1,n_xi
      i1 = nEdge + (i-1)*nDigits + 1
      i2 = nEdge + i*nDigits
      value = nint(mod(real(i,8),divider))
      if (value==0) then
         write(headerLines(1)(i1:i1),'(i1)') nint(real(i,8)/divider)
      end if
      write(headerLines(2)(i1:i2), frmti) value
   end do
   
   frmti = '(iX)'
   frmts = '(aX)'
   write(frmts(3:3), '(i1)') nDigits

   ! create template nice map file
   open(outUnit, file=trim(mapDir)//'ROMS_'//trim(domain)//'_nice-map.txt', &
                 form='formatted', action='write', iostat=ios)
   if (ios/=0) then
      write(*,'("Error opening file: ",a)') trim(mapDir)//'ROMS_'//trim(domain)//'_nice-map.txt'
      stop
   end if

   do i = 1,2
      write(outUnit,'(a)') trim(headerLines(i))
   end do

   iOff = nEdge + nDigits*n_xi
   do j = 1,n_eta
      o_j = floor(log10(real(n_eta-j+1,8)))
      i1 = nEdge - o_j
      i2 = nEdge
      write(frmti(3:3), '(i1)') o_j+1
      write(mapLines(j)(i1:i2), trim(frmti)) n_eta-j+1
      write(mapLines(j)(iOff+i1:iOff+i2), trim(frmti)) n_eta-j+1
      do i = 1,n_xi
         i1 = nEdge + (i-1)*nDigits + 1
         i2 = nEdge + i*nDigits
         if (mask_rho(i,n_eta-j+1)>0.5) cycle
         write(mapLines(j)(i1:i2), frmts) repeat('.', i2-i1+1)
      end do
      write(outUnit,'(a)') trim(mapLines(j))
   end do

   do i = 2,1,-1
      write(outUnit,'(a)') trim(headerLines(i))
   end do
   
   close(outUnit)
   
   ! finish execution if no exclude file has to be read
   if (.not.doExclude) stop
   
   ! read netCDF exclude file and create nice map with excluded cells marked
   inquire(file=trim(excludeFile), exist=lExist)
   if (.not.lExist) then
      write(*,'(a)') 'File does not exist: '//trim(excludeFile)//'.'
      stop
   end if
   ncStatus = NF90_OPEN(trim(excludeFile), NF90_NOWRITE, ncUnit)
   call nc_check(ncStatus, ios)
   if (ios/=0) then
      write(*,'("NetCDF-error opening ",a,".")') trim(excludeFile)
      stop
   end if
   
   call get_nc_dim(ncUnit, trim(n_sig_var), n_sig)
   
   allocate ( var_excl(n_xi,n_eta,n_sig), mask_excl(n_xi,n_eta) )
   call get_nc_var(ncUnit, trim(excludeVar), (/ n_xi, n_eta, n_sig /), var_excl)
   mask_excl = mask_rho
   where (var_excl(:,:,1)>eps.and.mask_rho>0.5) mask_excl = -1.e0
   
   open(outUnit, file=trim(mapDir)//'ROMS_'//trim(domain)//'_exclude_cells_nice-map.txt', &
                 form='formatted', action='write', iostat=ios)
   if (ios/=0) then
      write(*,'("Error opening file: ",a)') trim(mapDir)//'ROMS_'//trim(domain)//'_exclude_cell_nice-map.txt'
      stop
   end if

   if (nDigits==1) then
      exclStr = '1'
   else if (nDigits==2) then
      exclStr = '01'
   end if
   
   do i = 1,2
      write(outUnit,'(a)') trim(headerLines(i))
   end do

   iOff = nEdge + nDigits*n_xi
   do j = 1,n_eta
      o_j = floor(log10(real(n_eta-j+1,8)))
      i1 = nEdge - o_j
      i2 = nEdge
      write(frmti(3:3), '(i1)') o_j+1
      write(mapLines(j)(i1:i2), trim(frmti)) n_eta-j+1
      write(mapLines(j)(iOff+i1:iOff+i2), trim(frmti)) n_eta-j+1
      do i = 1,n_xi
         i1 = nEdge + (i-1)*nDigits + 1
         i2 = nEdge + i*nDigits
         if (mask_excl(i,n_eta-j+1)>0.5) then
            cycle
         else if (mask_excl(i,n_eta-j+1)>-0.5.and.mask_excl(i,n_eta-j+1)<0.5) then
            write(mapLines(j)(i1:i2), frmts) repeat('.', i2-i1+1)
         else
            write(mapLines(j)(i1:i2), frmts) trim(exclStr)
         end if
      end do
      write(outUnit,'(a)') trim(mapLines(j))
   end do

   do i = 2,1,-1
      write(outUnit,'(a)') trim(headerLines(i))
   end do
   
   close(outUnit)

contains
!=======================================================================
! SUBROUTINES
!=======================================================================
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
   subroutine get_nc_var(fileUnit, varName, dimSize, varVal)
!
! !DESCRIPTION:
! read variable from NetCDF file
!
   implicit none
!
! !OUTPUT PARAMETERS:
   integer         , intent(in)    :: fileUnit
   integer         , intent(in)    :: dimSize(2)
   character(len=*), intent(in)    :: varName
   real(8)         , intent(inout) :: varVal(dimSize(1), dimSize(2))
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

   ncStatus = NF90_GET_VAR(fileUnit, varID, varVal)
   call nc_check(ncStatus, ios)
   if (ios/=0) then
      write(*,'("NetCDF-error reading ",a,".")') trim(varName)
      stop
   end if

   end subroutine get_nc_var
!=======================================================================
   
end program make_nicemap4ROMS
