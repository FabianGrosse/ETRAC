!
! !MODULE: etrac_main.f90  --- ETRAC main routine
!
! !INTERFACE:
   MODULE mod_etrac_main
!
! !DESCRIPTION:
!
! !USES:
   use mod_etrac_common
   use mod_etrac_init, only: get_nc_field
   
   implicit none
!
!  default: all is private.
   private
   
   public do_tbnt, get_bulk_vars, get_vols
#ifndef TBNTonly_bulk_bud
   public do_tbnt_targets, get_tbnt_vars, final_bulk_check
#endif
!
! !LOCAL VARIABLES:
!#ifndef TBNTonly_bulk_bud
!   real(dp), dimension(:,:)    , allocatable :: tempFlxFrac3D
!   real(dp), dimension(:,:)    , allocatable :: tempFlxFrac2D
!   real(dp), dimension(:,:)    , allocatable :: lnkFlxFrac3D
!   real(dp), dimension(:,:)    , allocatable :: lnkFlxFrac2D
!   real(dp), dimension(:,:)    , allocatable :: varChange3D
!#ifndef TBNTnoVar2D
!   real(dp), dimension(:,:)    , allocatable :: varChange2D
!#endif
!#endif
!
!=======================================================================

   contains

!=======================================================================
!
! !INTERFACE:
   subroutine do_tbnt(iReadRec, ierr)
#define SUBROUTINE_NAME 'do_tbnt'
!
! !DESCRIPTION:
! do TBNT calculations
!
! !USES:
   implicit none
!
! !OUTPUT PARAMETERS:
   integer, intent(inout) :: ierr
   integer, intent(in   ) :: iReadRec
!
! !LOCAL VARIABLES:
#ifdef TBNTonline_budget_check
   integer :: i
#endif
#ifndef TBNTonly_bulk_bud
   integer :: level
   real(dp):: flxFac, cumFac
#endif
!
!-----------------------------------------------------------------------
#include "call-trace.inc"

   write(tbnt_log_unit, '(a,i5)')'  MAIN STEP:', iStep

#ifdef TBNTonly_bulk_bud
   ! get volumes
   call get_vols(iReadRec, ierr)
#endif

   ! get bulk balance fluxes
   call get_bulk_fluxes(iReadRec+1, ierr)

#ifndef TBNTonly_bulk_bud
   ! ======================================================
   ! ================= DO TBNT SIMULATION =================
   ! ======================================================

!   allocate ( tempFlxFrac2D(nBotCells,TBNTnFluxFractions2D) )
!   allocate ( tempFlxFrac3D(nWetCells,TBNTnFluxFractions3D) )
!#ifndef TBNTnoVar2D
!   allocate ( varChange2D(nBotCells,TBNTnVarFractions2D) )
!#endif
!   allocate ( varChange3D(nWetCells,TBNTnVarFractions3D) )
!   if (TBNTlinkedFluxes) then
!      allocate ( lnkFlxFrac2D(nBotCells,TBNTnLinkedFluxFractions2D) )
!      allocate ( lnkFlxFrac3D(nWetCells,TBNTnLinkedFluxFractions3D) )
!   end if
   
   do iSubStep = 1,TBNTnSubSteps
      cumFac = 0.e0
      flxFac = 1.e0
      level = 1
      call tbnt_recursive_step(iReadRec, cumFac, flxFac, level, ierr)
      
      ! check bulk variable for positivity
      call check_bulk_vars(.true.,ierr)
   end do

#ifdef TBNTonline_budget_check
   ! check calculated bulk variable versus bulk variable from NetCDF file
   call get_bulk_vars(iReadRec+1, ierr)
   do iVar = 1,TBNTnVars3D
      BULKvar3D(:,iVar) = BULKvar3D(:,iVar)/volOld*volNew
   end do
#ifndef TBNTnoVar2D
   do iVar = 1,TBNTnVars2D
      do i = 1,nBotCells
         iPut = bottom2pelag(i)
         if (.not.surfaceCell(iPut)) cycle
         BULKvar2D(i,iVar) = BULKvar2D(i,iVar)/volOld(iPut)*volNew(iPut)
      end do
   end do
#endif
   
   call check_bulk_vars(.false.,ierr)
#endif
   
!   deallocate ( tempFlxFrac2D, tempFlxFrac3D, varChange3D )
!#ifndef TBNTnoVar2D
!   deallocate ( varChange2D )
!#endif
!
!   if (TBNTlinkedFluxes) deallocate ( lnkFlxFrac2D, lnkFlxFrac3D )

#else
   ! ======================================================
   ! =================== DO BULK BUDGET ===================
   ! ======================================================

   ! get bulk variables at beginning of each day
   call get_bulk_vars(iReadRec, ierr)
   
   ! calculate bulk budget
   call calc_bulk_bud(iReadRec, ierr)
   
#endif
   
   end subroutine do_tbnt
!=======================================================================
!
! !INTERFACE:
   subroutine get_vols(iReadRec, ierr)
!
! !DESCRIPTION:
! read volumes
!
! !USES:
   implicit none
!
! !OUTPUT PARAMETERS:
   integer, intent(inout) :: ierr
   integer, intent(in   ) :: iReadRec
!
! !LOCAL VARIABLES:
!
   integer                :: iVolRec
   real(dp)               :: alpha, beta
!-----------------------------------------------------------------------

   ! read volume
   !if (mod(TBNTstep*(iReadRec-1), iMinPerDay)==0.or.iReadRec==TBNTiStart) then
      iVolRec = floor(real(TBNTstep*iReadRec,dp)/rMinPerDay)! + 1 ! should be correct this way (w/o +1)
      if (TBNTstep/=iMinPerDay) iVolRec = iVolRec + 1
      iVolRec = min(iVolRec,TBNTyearDays)
      call get_nc_field(bulk_nc_unit, trim(volVar), iVolRec, ierr, fieldData3D=volOldIn)
      iVolRec = min(iVolRec+1,TBNTyearDays+1)
      call get_nc_field(bulk_nc_unit, trim(volVar), iVolRec, ierr, fieldData3D=volNewIn)
      !volOld = volOldIn
      !volNew = volNewIn
   !end if
   
   ! interpolate volume
   alpha = real(mod(TBNTstep*(iReadRec-1),iMinPerDay),dp)/rMinPerDay
   beta = 1.e0 - alpha
   volOld = beta*volOldIn + alpha*volNewIn
   alpha = real(mod(TBNTstep*iReadRec,iMinPerDay),dp)/rMinPerDay
   if (alpha==0.e0) alpha = 1.e0
   beta = 1.e0 - alpha
   volNew = beta*volOldIn + alpha*volNewIn
   
   end subroutine get_vols
!=======================================================================
!
! !INTERFACE:
   subroutine get_bulk_vars(iReadRec, ierr)
#define SUBROUTINE_NAME 'get_bulk_vars'
!
! !DESCRIPTION:
! read bulk variables
!
! !USES:
   implicit none
!
! !OUTPUT PARAMETERS:
   integer, intent(inout) :: ierr
   integer, intent(in   ) :: iReadRec
!
! !LOCAL VARIABLES:
#if !defined TBNTonly_bulk_bud & !defined TBNTnoVar2D
   integer                :: iiVar
#endif
!-----------------------------------------------------------------------
#include "call-trace.inc"

#ifndef TBNTonly_bulk_bud
! TBNT
   iVarPnt = 0
   do iVar = 1,TBNTnVars3D ! 3D variables
      if (TBNTvar(iVar)%type==TBNTdum3dVar.and.TBNTvar(iVar)%auxVar%fac==0.e0) cycle ! skip dummy variables w/o auxiliary variable
      iVarPnt = TBNTvar(iVar)%iTBNT
      if (TBNTvar(iVar)%type==TBNTdum3dVar) then
         call get_nc_field(bulk_nc_unit, trim(TBNTvar(iVar)%auxVar%name), iReadRec, ierr, fieldData3D=BULKvar3D(:,iVarPnt))
         BULKvar3D(:,iVarPnt) = BULKvar3D(:,iVarPnt) * TBNTvar(iVar)%auxVar%fac
      else
         call get_nc_field(bulk_nc_unit, trim(TBNTvar(iVar)%name), iReadRec, ierr, fieldData3D=BULKvar3D(:,iVarPnt))
      end if
      BULKvar3D(:,iVarPnt) = BULKvar3D(:,iVarPnt)*volOld ! calculate amount of substance
   end do
#ifndef TBNTnoVar2D
   do iVar = 1,TBNTnVars2D ! 2D variables
      iiVar = TBNTnVars3D+iVar
      if (TBNTvar(iiVar)%type==TBNTdum2dVar.and.TBNTvar(iiVar)%auxVar%fac==0.e0) cycle ! skip dummy variables w/o auxiliary variable
      iVarPnt = TBNTvar(iiVar)%iTBNT
      if (TBNTvar(iiVar)%type==TBNTdum2dVar) then
         call get_nc_field(bulk_nc_unit, trim(TBNTvar(iiVar)%auxVar%name), iReadRec, ierr, fieldData2D=BULKvar2D(:,iVarPnt))
         BULKvar2D(:,iVarPnt) = BULKvar2D(:,iVarPnt) * TBNTvar(iiVar)%auxVar%fac
      else
         call get_nc_field(bulk_nc_unit, trim(TBNTvar(iiVar)%name), iReadRec, ierr, fieldData2D=BULKvar2D(:,iVarPnt))
      end if
      BULKvar2D(:,iVarPnt) = BULKvar2D(:,iVarPnt)*area2D
   end do
#endif
#else
! BULK BUDGET (in budget mode, amount of substance of state variable is calculated directly in the budget check routine)
   if (BULKbudVar%type==TBNTpro2dVar.or.BULKbudVar%type==TBNTdia2dVar.or.BULKbudVar%type==TBNTdum2dVar) then ! 2D variable
      if (BULKbudVar%type==TBNTdum2dVar) then
         call get_nc_field(bulk_nc_unit, trim(BULKbudVar%auxVar%name), iReadRec, ierr, fieldData2D=BULKvar2D(:,1))
         BULKvar2D(:,1) = BULKvar2D(:,1) * BULKbudVar%auxVar%fac
      else
         call get_nc_field(bulk_nc_unit, trim(BULKbudVar%name), iReadRec, ierr, fieldData2D=BULKvar2D(:,1))
      end if
   else ! 3D variable
      if (BULKbudVar%type==TBNTdum3dVar) then
         call get_nc_field(bulk_nc_unit, trim(BULKbudVar%auxVar%name), iReadRec, ierr, fieldData3D=BULKvar3D(:,1))
         BULKvar3D(:,1) = BULKvar3D(:,1) * BULKbudVar%auxVar%fac
      else
         call get_nc_field(bulk_nc_unit, trim(BULKbudVar%name), iReadRec, ierr, fieldData3D=BULKvar3D(:,1))
      end if
   end if
#endif

   end subroutine get_bulk_vars
!=======================================================================
!
! !INTERFACE:
   subroutine get_bulk_fluxes(iReadRec, ierr)
#define SUBROUTINE_NAME 'get_bulk_fluxes'
!
! !DESCRIPTION:
! read bulk fluxes and calculate substance flux
!
! !USES:
   implicit none
!
! !OUTPUT PARAMETERS:
   integer, intent(inout) :: ierr
   integer, intent(in   ) :: iReadRec
!
! !LOCAL VARIABLES:
!
   real(dp), parameter    :: eps = 1.0e-06
   
   integer                :: flxType, auxFlxType, iiFlux
#ifdef TBNTonly_bulk_bud
   logical                :: isOutVar
#endif
   
!-----------------------------------------------------------------------
#include "call-trace.inc"
   
#ifndef TBNTonly_bulk_bud
! TBNT
   
   ! read fluxes
   do iFlux = 1,TBNTnFluxes3D+TBNTnLinkedFluxes3D ! 3D fluxes
      flxType = TBNTflux(iFlux)%type
      auxFlxType = TBNTflux(iFlux)%auxFlux%type
      ! get bulk flux
      if (flxType==TBNTder3DFlux) then ! 3D derived flux
         call get_nc_field(bulk_nc_unit, trim(TBNTflux(iFlux)%auxFlux%name),iReadRec, ierr, fieldData3D=BULKflux3D(:,iFlux))
         BULKflux3D(:,iFlux) = BULKflux3D(:,iFlux)*TBNTtimeFac * TBNTflux(iFlux)%auxFlux%fac
      else                             ! 3D prognostic fluxes
         call get_nc_field(bulk_nc_unit, trim(TBNTflux(iFlux)%name), iReadRec, ierr, fieldData3D=BULKflux3D(:,iFlux))
         BULKflux3D(:,iFlux) = BULKflux3D(:,iFlux)*TBNTtimeFac
      end if
#ifndef TBNTmass_fluxes
      ! calculate substance flux
      BULKflux3D(:,iFlux) = BULKflux3D(:,iFlux)*area3D
#endif
   end do
   do iFlux = 1,TBNTnFluxes2D+TBNTnLinkedFluxes2D ! 2D fluxes
      iiFlux = TBNTnFluxes3D+TBNTnLinkedFluxes3D + iFlux
      flxType = TBNTflux(iiFlux)%type
      auxFlxType = TBNTflux(iiFlux)%auxFlux%type
      if     (flxType==TBNTsed2DFlux) then ! 2D prognostic flux (sediment)
         call get_nc_field(bulk_nc_unit, trim(TBNTflux(iiFlux)%name), iReadRec, ierr, fieldData2D=BULKflux2D(:,iFlux))
         BULKflux2D(:,iFlux) = BULKflux2D(:,iFlux)*TBNTtimeFac
      elseif (flxType==TBNTder2DFlux) then ! 2D derived flux (sediment)
         call get_nc_field(bulk_nc_unit, trim(TBNTflux(iiFlux)%auxFlux%name),iReadRec, ierr, fieldData2D=BULKflux2D(:,iFlux))
         BULKflux2D(:,iFlux) = BULKflux2D(:,iFlux)*TBNTtimeFac * TBNTflux(iiFlux)%auxFlux%fac
      end if
#ifndef TBNTmass_fluxes
      ! calculate substance flux
      BULKflux2D(:,iFlux) = BULKflux2D(:,iFlux)*area2D
#endif
   end do
   
#else
! BULK BUDGET
   ! 3D fluxes
   do iFlux = 1,BULKnBudFluxes3D
      flxType = BULKbudFlux(iFlux)%type
      auxFlxType = BULKbudFlux(iFlux)%auxFlux%type
      if (flxType==TBNTder3DFlux) then ! 3D derived flux
         call get_nc_field(bulk_nc_unit,trim(BULKbudFlux(iFlux)%auxFlux%name),iReadRec,ierr,fieldData3D=BULKflux3D(:,iFlux))
         BULKflux3D(:,iFlux) = BULKflux3D(:,iFlux)*BULKbudFlux(iFlux)%auxFlux%fac ! calculation of derived flux from C-flux using redfield
      else                                 ! 3D prognostic fluxes
         call get_nc_field(bulk_nc_unit,trim(BULKbudFlux(iFlux)%name),iReadRec,ierr,fieldData3D=BULKflux3D(:,iFlux))
      end if
      ! define sign of flux (source (+) or sink (-)) in case that it's no exchange flux (no transport/mixing)
      if (flxType/=TBNTadvecFlux.and.auxFlxType/=TBNTadvecFlux.and. &
          flxType/=TBNTdiffuFlux.and.auxFlxType/=TBNTdiffuFlux) then
         isOutVar = .false.
         call string_check(BULKbudVar%name, BULKbudFlux(iFlux)%varOut%name, isOutVar)
         if (.not.isOutVar) BULKflux3D(:,iFlux) = -1.e0 * BULKflux3D(:,iFlux)
      end if
      !where (abs(BULKflux3D(:,iFlux))<eps) BULKflux3D(:,iFlux) = 0.0e0
#ifndef TBNTmass_fluxes
      ! calculate substance flux
      BULKflux3D(:,iFlux) = BULKflux3D(:,iFlux)*area3D
#endif
   end do
   ! 2D fluxes
   do iFlux = 1,BULKnBudFluxes2D
      iiFlux = BULKnBudFluxes3D + iFlux
      flxType = BULKbudFlux(iiFlux)%type
      auxFlxType = BULKbudFlux(iiFlux)%auxFlux%type
      if     (flxType==TBNTsed2DFlux) then ! 2D prognostic flux (sediment)
         call get_nc_field(bulk_nc_unit,trim(BULKbudFlux(iiFlux)%name),iReadRec,ierr,fieldData2D=BULKflux2D(:,iFlux))
      elseif (flxType==TBNTder2DFlux) then ! 2D derived flux (sediment)
         call get_nc_field(bulk_nc_unit,trim(BULKbudFlux(iiFlux)%auxFlux%name),iReadRec,ierr,fieldData2D=BULKflux2D(:,iFlux))
         BULKflux2D(:,iFlux) = BULKflux2D(:,iFlux)*BULKbudFlux(iiFlux)%auxFlux%fac ! calculation of derived flux from C-flux using redfield
      end if
      ! define sign of flux: source (+) or sink (-)
      isOutVar = .false.
      call string_check(BULKbudVar%name, BULKbudFlux(iiFlux)%varOut%name, isOutVar)
      if (.not.isOutVar) BULKflux2D(:,iFlux) = -1.e0 * BULKflux2D(:,iFlux)
      !where (abs(BULKflux2D(:,iFlux))<eps) BULKflux2D(:,iFlux) = 0.0e0
#ifndef TBNTmass_fluxes
      ! calculate substance flux
      BULKflux2D(:,iFlux) = BULKflux2D(:,iFlux)*area2D
#endif
   end do
#endif

   end subroutine get_bulk_fluxes
!=======================================================================
#ifndef TBNTonly_bulk_bud
!
! !INTERFACE:
   subroutine check_bulk_vars(positivity,ierr)
#define SUBROUTINE_NAME 'check_bulk_vars'
!
! !DESCRIPTION:
!  check bulk variable for positivity
!
! !USES:
   implicit none
!
! !OUTPUT PARAMETERS:
   integer, intent(inout) :: ierr
   logical, intent(in   ) :: positivity
!
! !LOCAL VARIABLES:
   integer                :: iiVar, i
   character(len=nameLen) :: srcName
!-----------------------------------------------------------------------
#include "call-trace.inc"

   iiVar = 0

   srcName = 'BULK'
   do iVar = 1,TBNTnVars3D ! 3D variables
      if (TBNTvar(iVar)%type==TBNTdum3dVar.and.TBNTvar(iVar)%auxVar%fac==0.e0) cycle
      do i = 1,nWetCells
         if (excludeCell3D(i)) cycle
         if (positivity) then
            call check_absVar(TBNTvar(iVar)%name,srcName,BULKvar3D(i,iVar),i,ierr)
         else
            call check_absVar(TBNTvar(iVar)%name,srcName,BULKvar3D(i,iVar),i,ierr, absBulk=BULKvar3D(i,iVar))
         end if
      end do
   end do
#ifndef TBNTnoVar2D
   do iVar = 1,TBNTnVars2D ! 2D variables
      iiVar = TBNTnVars3D + iVar
      if (TBNTvar(iiVar)%type==TBNTdum2dVar.and.TBNTvar(iiVar)%auxVar%fac==0.e0) cycle
      do i = 1,nBotCells
         if (excludeCell2D(i)) cycle
         if (positivity) then
            call check_absVar(TBNTvar(iiVar)%name,srcName,BULKvar2D(i,iVar),i,ierr)
         else
            call check_absVar(TBNTvar(iiVar)%name,srcName,BULKvar2D(i,iVar),i,ierr, absBulk=BULKvar2D(i,iVar))
         end if
      end do
   end do
#endif

   end subroutine check_bulk_vars
#endif
!=======================================================================
#ifndef TBNTonly_bulk_bud
!
! !INTERFACE:
   recursive subroutine tbnt_recursive_step(iReadRec, cumFac_, flxFac_, level_, ierr)
#define SUBROUTINE_NAME 'tbnt_recursive_step'
!
! !DESCRIPTION:
! recursive time step subdivision
!
! !USES:
   implicit none
!
! !OUTPUT PARAMETERS:
   integer , intent(in   ) :: iReadRec
   integer , intent(inout) :: level_, ierr
   real(dp), intent(in   ) :: flxFac_
   real(dp), intent(inout) :: cumFac_
!
! !LOCAL VARIABLES:
   integer, parameter      :: max_level = 50
   integer                 :: int_err
   real(dp)                :: flxFac_halved
!
!-----------------------------------------------------------------------
#include "call-trace.inc"

   if (level_ > max_level) then
      ierr = ierr + 1
      write(error_msg,'("Maximum level for time step subdivision exceeded.")') 
      call stop_tbnt(error_msg,ierr)
   end if
   
   if (level_==1) then      
      ! calculate fraction fluxes and variable fractions (TBNT core)
      call calc_tbnt_fractions

      ! bring fractions into system
      call get_tbnt_fractions(ierr)
   end if
   
   ! apply fluxes to variable fractions
   int_err = 0
   call update_tbnt_fractions(cumFac_, flxFac_, level_, int_err)
   
   if (int_err==0) then
      cumFac_ = cumFac_ + flxFac_
      
      ! calculate relative variable fractions & bulk variable
      call calc_tbnt_relative_fractions(ierr)
      
      if (level_>1.and.cumFac_<1.e0) then
         ! calculate fraction fluxes and variable fractions (TBNT core)
         call calc_tbnt_fractions

         ! bring fractions into system
         call get_tbnt_fractions(ierr)
      end if
   else
      flxFac_halved = 0.5 * flxFac_
      level_ = level_ + 1
      call tbnt_recursive_step(iReadRec, cumFac_, flxFac_halved, level_, ierr)
      call tbnt_recursive_step(iReadRec, cumFac_, flxFac_halved, level_, ierr)
      level_ = level_ - 1
   end if

   end subroutine tbnt_recursive_step
#endif
!=======================================================================
#ifndef TBNTonly_bulk_bud
!
! !INTERFACE:
   subroutine update_tbnt_fractions(cumFac_, flxFac_,level_,integ_err)
!
! !DESCRIPTION:
! update state variable fractions
!
! !USES:
   implicit none
!
! !OUTPUT PARAMETERS:
   real(dp), intent(in   ) :: cumFac_, flxFac_
   integer , intent(in   ) :: level_
   integer , intent(inout) :: integ_err
!
! !LOCAL VARIABLES:
   real(dp), parameter     :: eps = 1.e-10
   
   integer                 :: iiVar, i
   real(dp)                :: value, r_flxFac
   logical                 :: isVar3D
!
!-----------------------------------------------------------------------

   integ_err = 0
   isVar3D = .true.
   iiVar = 0
   
   r_flxFac = 1.e0/flxFac_
   
   ! check all fractions for positivity
   fractionLoop: do iFrac = 1,TBNTnFractions
      ! 3D variables
      do iVar = 1,TBNTnVars3D
         iVarPnt = TBNTvar3Dpnt(iVar,iFrac)
         do i = 1,nWetCells
            if (excludeCell3D(i)) cycle
            value = TBNTvar3D(i,iVarPnt) + varChange3D(i,iVarPnt)*flxFac_
            call check_absVar(TBNTvar(iVar)%name,TBNTsource(iFrac)%name,value,i,integ_err)
! FGrosse: Nov. 28, 2017: first "if line" replaced by if-block
!            if (integ_err/=0) exit fractionLoop
            if (integ_err/=0) then
               if (abs(value)>eps) exit fractionLoop
               ! correct variable change such that new variable fraction is ZERO
               integ_err = 0
               varChange3D(i,iVarPnt) = varChange3D(i,iVarPnt) - value*r_flxFac
            end if
         end do
      end do
      ! 2D variables
      do iVar = 1,TBNTnVars2D
         iiVar = iVar + TBNTnVars3D
         iVarPnt = TBNTvar2Dpnt(iVar,iFrac)
         do i = 1,nBotCells
            if (excludeCell2D(i)) cycle
            value = TBNTvar2D(i,iVarPnt) + varChange2D(i,iVarPnt)*flxFac_
            call check_absVar(TBNTvar(iiVar)%name,TBNTsource(iFrac)%name,value,i,integ_err)
! FGrosse: Nov. 28, 2017: first if-block replaced by second if-block
!            if (integ_err/=0) then
!               isVar3D = .false.
!               exit fractionLoop
!            end if
            if (integ_err/=0) then
               if (abs(value)>eps) then
                  isVar3D = .false.
                  exit fractionLoop
               end if
               ! correct variable change such that new variable fraction is ZERO
               integ_err = 0
               varChange2D(i,iVarPnt) = varChange2D(i,iVarPnt) - value*r_flxFac
            end if
         end do
      end do
   end do fractionLoop
   
   if (integ_err/=0) then
      if (isVar3D) iiVar = iVar
      call breakline(30)
      write(tbnt_log_unit, '("VARIABLE negative: ",a)') trim(TBNTvar(iiVar)%name)
      write(tbnt_log_unit, '("   frac => ",i14)') iFrac 
      write(tbnt_log_unit, '("   step => ",2i7)') iStep, iSubStep
      write(tbnt_log_unit, '("      i => ",i14)') i
#ifdef TBNTconvert_3Dto1D
      write(tbnt_log_unit, '("(i,j,k) => ",i4,2i5)') index1Dto3D(i,:)
#endif
      write(tbnt_log_unit, '("  level => ",i14)') level_
      write(tbnt_log_unit, '("    fac => ",E18.10)') flxFac_
      write(tbnt_log_unit, '(" cumFac => ",E18.10)') cumFac_
      write(tbnt_log_unit, '("  value => ",E18.10)') value
      if (isVar3D) then
         write(tbnt_log_unit, '("   BULK => ",E18.10)') BULKvar3D(i,iVar)
         write(tbnt_log_unit, '("    var => ",E18.10)') TBNTvar3D(i,TBNTvar3Dpnt(iVar,iFrac))
         write(tbnt_log_unit, '(" change => ",E18.10)') varChange3D(i,TBNTvar3Dpnt(iVar,iFrac))*flxFac_
#ifndef TBNTnoVar2D
      else
         write(tbnt_log_unit, '("    var => ",E18.10)') TBNTvar2D(i,TBNTvar2Dpnt(iVar,iFrac))
         write(tbnt_log_unit, '(" change => ",E18.10)') varChange2D(i,TBNTvar2Dpnt(iVar,iFrac))*flxFac_
#endif
      end if
      return
   end if
   
   ! if integ_err==0: update variables
   TBNTvar3D = TBNTvar3D + varChange3D*flxFac_
   TBNTvar2D = TBNTvar2D + varChange2D*flxFac_

   ! cumulate fraction fluxes
   if (TBNTabsFracOutStep>=0) then
      TBNTflux3D = TBNTflux3D + tempFlxFrac3D*flxFac_
      TBNTflux2D = TBNTflux2D + tempFlxFrac2D*flxFac_
   end if
   
   ! cumulate linked fraction fluxes
   if (TBNTlinkedFluxes) then
      TBNTlinkedFlux3D = TBNTlinkedFlux3D + lnkFlxFrac3D*flxFac_
      TBNTlinkedFlux2D = TBNTlinkedFlux2D + lnkFlxFrac2D*flxFac_
   end if

   end subroutine update_tbnt_fractions
#endif
!=======================================================================
#ifndef TBNTonly_bulk_bud
!
! !INTERFACE:
   subroutine get_tbnt_fractions(ierr)
#define SUBROUTINE_NAME 'get_tbnt_fractions'
!
! !DESCRIPTION:
! brings tracers into system
!
! !USES:
   implicit none
!
! !OUTPUT PARAMETERS:
   integer, intent(inout) :: ierr
!
! !LOCAL VARIABLES:
!
!-----------------------------------------------------------------------
#include "call-trace.inc"

   if (TBNTnRiverFractions>0.or.TBNTnDummyRivers>0) call get_tbnt_river_fracs(ierr)
   if (TBNTnOpenbFractions>0)                       call get_tbnt_openb_fracs(ierr)
   if (TBNTnAtmosFractions>0.or.TBNTnDummyAtmos >0) call get_tbnt_atmos_fracs(ierr)

   end subroutine get_tbnt_fractions
#endif
!=======================================================================
#ifndef TBNTonly_bulk_bud
!
! !INTERFACE:
   subroutine get_tbnt_river_fracs(ierr)
#define SUBROUTINE_NAME 'get_tbnt_river_fracs'
!
! !DESCRIPTION:
! bring river tracers into system (and dummy rivers)
!
! !USES:
   implicit none
!
! !OUTPUT PARAMETERS:
   integer, intent(inout) :: ierr
!
! !LOCAL VARIABLES:
!
   real(dp), parameter    :: eps = 1.e-15

   integer                :: iRvDisPnt, i
   integer                :: nInputs, ii
   integer, dimension(:), allocatable :: i_to
!-----------------------------------------------------------------------
#include "call-trace.inc"
   
   ! loop over variables and river fractions
   do iVar = 1,TBNTnVars3D
      if ((TBNTvar(iVar)%type==TBNTdum3dVar.and.TBNTvar(iVar)%auxVar%fac==0.e0).or.sum(iRvDisFlux(iVar,:))==0) cycle
      ! loop over fractions/sources
      do iFrac = 1,TBNTnFractions
         if (TBNTsource(iFrac)%type/=TBNTriverSource.and.TBNTsource(iFrac)%type/=TBNTdummySource) cycle
         iVarPnt = TBNTvar3Dpnt(iVar,iFrac)                         ! pointer to variable fraction
         do iFlux = 1,nRivFluxes
            iRvDisPnt = TBNTflux3Dpnt(iRvDisFlux(iVar,iFlux),iFrac) ! pointer to river discharge flux on variable fraction
            if (iRvDisPnt==0) cycle
            ! loop over subsources
            do iSource = 1,TBNTsource(iFrac)%nSubsources
               if (TBNTsource(iFrac)%subsource(iSource)%type/=TBNTriverSource) cycle
               ! river input cells
               nInputs = TBNTsource(iFrac)%subsource(iSource)%nInputs
               allocate ( i_to(nInputs) )
               i_to = TBNTsource(iFrac)%subsource(iSource)%i_to
               do i = 1,nInputs
                  ii = i_to(i)
                  if (excludeCell3D(ii)) cycle
                  tempFlxFrac3D(ii,iRvDisPnt) = BULKflux3D(ii,TBNTflux(iRvDisFlux(iVar,iFlux))%iTBNT)
                  if (tempFlxFrac3D(ii,iRvDisPnt)<eps) tempFlxFrac3D(ii,iRvDisPnt) = 0.e0
                  varChange3D(ii,iVarPnt) = varChange3D(ii,iVarPnt) + tempFlxFrac3D(ii,iRvDisPnt)
               end do
               deallocate ( i_to )
            end do
         end do
      end do
   end do
   
   end subroutine get_tbnt_river_fracs
#endif
!=======================================================================
#ifndef TBNTonly_bulk_bud
!
! !INTERFACE:
   subroutine get_tbnt_openb_fracs(ierr)
#define SUBROUTINE_NAME 'get_tbnt_openb_fracs'
!
! !DESCRIPTION:
! bring open boundary tracers into system
!
! !USES:
   implicit none
!
! !OUTPUT PARAMETERS:
   integer, intent(inout) :: ierr
!
! !LOCAL VARIABLES:
!
   real(dp), parameter    :: eps = 1.e-15
      
   integer                :: flxType, auxFlxType, flxComp
   real(dp)               :: flx, flxSign
   integer                :: i_fr, i_to, i_comp, i_flx
!-----------------------------------------------------------------------
#include "call-trace.inc"

   iVarPnt = 0
   iFlxPnt = 0
   iOff = TBNTnRiverFractions
   do iFrac = iOff+1,iOff+TBNTnOpenBFractions
      if (TBNTsource(iFrac)%type/=TBNTopenbSource) then
         ierr = ierr + 1
         write(error_msg,'(a)')'Selected source is no open boundary source.'
         call stop_tbnt(error_msg,ierr)
      end if
      ! only one subsource for open boundary sources => subsource(1)
      do iSource = 1,TBNTsource(iFrac)%subsource(1)%nInputs
         i_fr = TBNTsource(iFrac)%subsource(1)%i_from(iSource)
         i_to = TBNTsource(iFrac)%subsource(1)%i_to(iSource)
         i_comp = TBNTsource(iFrac)%subsource(1)%i_comp(iSource)
         flxSign = TBNTsource(iFrac)%subsource(1)%sign(iSource)
         if (flxSign*compSign(i_comp)<0.e0) then
            i_flx = i_fr
         else
            i_flx = i_to
         end if
         do iFlux = 1,TBNTnFluxes3D
            flxType = TBNTflux(iFlux)%type
            auxFlxType = TBNTflux(iFlux)%auxFlux%type
            flxComp = TBNTflux(iFlux)%i_comp
            if ((flxType/=TBNTadvecFlux.and.auxFlxType/=TBNTadvecFlux.and. &
                 flxType/=TBNTdiffuFlux.and.auxFlxType/=TBNTdiffuFlux).or.flxComp/=i_comp) cycle
            flx = flxSign*BULKflux3D(i_flx,iFlux)
            if (abs(flx)<eps.or.flx<0.e0) cycle
            iVar = TBNTflux(iFlux)%varOut%iTBNT
            iVarPnt = TBNTvar3Dpnt(iVar,iFrac)   ! pointer to variable fraction affected by open boundary input
            iFlxPnt = TBNTflux3Dpnt(iFlux,iFrac) ! pointer to open boundary flux on variable fraction
            tempFlxFrac3D(i_flx,iFlxPnt) = flx
            varChange3D(i_to,iVarPnt) = varChange3D(i_to,iVarPnt) + flx
         end do
      end do
   end do

   end subroutine get_tbnt_openb_fracs
#endif
!=======================================================================
#ifndef TBNTonly_bulk_bud
!
! !INTERFACE:
   subroutine get_tbnt_atmos_fracs(ierr)
#define SUBROUTINE_NAME 'get_tbnt_atmos_fracs'
!
! !DESCRIPTION:
! bring atmospheric tracers into system (and atmospheric dummy tracers)
!
! !USES:
   implicit none
!
! !OUTPUT PARAMETERS:
   integer, intent(inout) :: ierr
!
! !LOCAL VARIABLES:
!
   real(dp), parameter    :: eps = 1.e-15
   integer                :: i_to, i_fr, iSrc, iiFlux
!-----------------------------------------------------------------------
#include "call-trace.inc"
   
   iVarPnt = 0
   iFlxPnt = 0
   iOff = TBNTnRiverFractions+TBNTnOpenBFractions
   ! 3D fluxes
   do iFlux = 1,TBNTnFluxes3D
      if (TBNTflux(iFlux)%type/=TBNTatm3DFlux.and.TBNTflux(iFlux)%auxFlux%type/=TBNTatm3DFlux.and. &
          TBNTflux(iFlux)%type/=TBNTa2s3DFlux.and.TBNTflux(iFlux)%auxFlux%type/=TBNTa2s3DFlux) cycle
      iVar = TBNTflux(iFlux)%varOut%iTBNT
      do iFrac = iOff+1,TBNTnFractions
         if (TBNTsource(iFrac)%type/=TBNTatmosSource.and.TBNTsource(iFrac)%type/=TBNTdummySource) then
            ierr = ierr + 1
            write(error_msg,'(a)')'Selected source is no atmospheric deposition flux.'
            call stop_tbnt(error_msg,ierr)
         end if
         iVarPnt = TBNTvar3Dpnt(iVar,iFrac)   ! pointer to variable fraction affected by atmospheric deposition
         iFlxPnt = TBNTflux3Dpnt(iFlux,iFrac) ! pointer to atmospheric deposition flux on variable fraction
         do iSource = 1,TBNTsource(iFrac)%nSubsources
            if (TBNTsource(iFrac)%subsource(iSource)%type/=TBNTatmosSource) cycle
            ! set atmospheric deposition flux at affected grid points
            do iSrc = 1,TBNTsource(iFrac)%subsource(iSource)%nInputs
               i_to = TBNTsource(iFrac)%subsource(iSource)%i_to(iSrc)
               if (BULKflux3D(i_to,iFlux)<eps) cycle
               tempFlxFrac3D(i_to,iFlxPnt) = BULKflux3D(i_to,iFlux)
               varChange3D(i_to,iVarPnt) = varChange3D(i_to,iVarPnt) + BULKflux3D(i_to,iFlux)
            end do
         end do
      end do
   end do
   ! 2D fluxes
   do iFlux = 1,TBNTnFluxes2D
      iiFlux = TBNTnFluxes3D+TBNTnLinkedFluxes3D + iFlux
      if (TBNTflux(iiFlux)%type/=TBNTatm2DFlux.and.TBNTflux(iiFlux)%auxFlux%type/=TBNTatm2DFlux.and. &
          TBNTflux(iiFlux)%type/=TBNTa2s2DFlux.and.TBNTflux(iiFlux)%auxFlux%type/=TBNTa2s2DFlux) cycle
      iVar = TBNTflux(iiFlux)%varOut%iTBNT
      do iFrac = iOff+1,TBNTnFractions
         if (TBNTsource(iFrac)%type/=TBNTatmosSource.and.TBNTsource(iFrac)%type/=TBNTdummySource) then
            ierr = ierr + 1
            write(error_msg,'(a)')'Selected source is no atmospheric deposition flux.'
            call stop_tbnt(error_msg,ierr)
         end if
         iVarPnt = TBNTvar3Dpnt(iVar,iFrac)   ! pointer to variable fraction affected by atmospheric deposition
         iFlxPnt = TBNTflux2Dpnt(iFlux,iFrac) ! pointer to atmospheric deposition flux on variable fraction
         do iSource = 1,TBNTsource(iFrac)%nSubsources
            if (TBNTsource(iFrac)%subsource(iSource)%type/=TBNTatmosSource) cycle
            ! set atmospheric deposition flux at affected grid points
            do iSrc = 1,TBNTsource(iFrac)%subsource(iSource)%nInputs
               i_to = TBNTsource(iFrac)%subsource(iSource)%i_to(iSrc)
               do i_fr = 1,nSurfCells
                  if (surf2pelag(i_fr)==i_to) exit
               end do
               if (BULKflux2D(i_fr,iFlux)<eps) cycle
               tempFlxFrac2D(i_fr,iFlxPnt) = BULKflux2D(i_fr,iFlux)
               varChange3D(i_to,iVarPnt) = varChange3D(i_to,iVarPnt) + BULKflux2D(i_fr,iFlux)
            end do
         end do
      end do
   end do
   
   end subroutine get_tbnt_atmos_fracs
   
#endif
!=======================================================================
#ifndef TBNTonly_bulk_bud
!
! !INTERFACE:
   subroutine get_tbnt_vars
#define SUBROUTINE_NAME 'get_tbnt_vars'
!
! !DESCRIPTION:
! calculate initial distribution of variable fractions
!
! !USES:
   implicit none
!
! !LOCAL VARIABLES
!
   real(dp), parameter :: eps = 1.e-15
   
   integer :: iiVar
!
!-----------------------------------------------------------------------
#include "call-trace.inc"
   
   iiVar = 0

   do iFrac = 1,TBNTnFractions
      do iVar = 1,TBNTnVars3D ! 3D variables
         if (TBNTvar(iVar)%type==TBNTdum3dVar.and.TBNTvar(iVar)%auxVar%fac==0.e0) cycle ! skip dummy variables w/o auxiliary variable
         iVarPnt = TBNTvar3Dpnt(iVar,iFrac)
         TBNTvar3D(:,iVarPnt) = TBNTrelVar3D(:,iVarPnt) * BULKvar3D(:,iVar)
      end do
#ifndef TBNTnoVar2D
      do iVar = 1,TBNTnVars2D ! 2D variables
         iiVar = TBNTnVars3D + iVar
         if (TBNTvar(iiVar)%type==TBNTdum2dVar.and.TBNTvar(iiVar)%auxVar%fac==0.e0) cycle ! skip dummy variables w/o auxiliary variable
         iVarPnt = TBNTvar2Dpnt(iVar,iFrac)
         TBNTvar2D(:,iVarPnt) = TBNTrelVar2D(:,iVarPnt) * BULKvar2D(:,iVar)
      end do
#endif
   end do
   
   end subroutine get_tbnt_vars
#endif
!=======================================================================
#ifndef TBNTonly_bulk_bud
!
! !INTERFACE:
   subroutine calc_tbnt_relative_fractions(ierr)
#define SUBROUTINE_NAME 'calc_tbnt_relative_fractions'
!
! !DESCRIPTION:
! calculate relative variable fractions
!
! !USES:
   implicit none
!
! !OUTPUT PARAMETERS:
   integer, intent(inout) :: ierr
!
! !LOCAL VARIABLES:
   real(dp), parameter    :: eps = 1.e-12

   real(dp)               :: epsFrac
   integer                :: i
!-----------------------------------------------------------------------
#include "call-trace.inc"

   epsFrac = eps/TBNTnFractions
   
   ! calculate fraction sum (= BULKvar)
   BULKvar2D = 0.e0
   BULKvar3D = 0.e0
   do iFrac = 1,TBNTnFractions
      ! 3D variables 
      do iVar = 1,TBNTnVars3D
         iVarPnt = TBNTvar3Dpnt(iVar,iFrac)
         where (.not.excludeCell3D.and.abs(TBNTvar3D(:,iVarPnt))<epsFrac) TBNTvar3D(:,iVarPnt) = 0.e0
         where (.not.excludeCell3D) BULKvar3D(:,iVar) = BULKvar3D(:,iVar) + TBNTvar3D(:,iVarPnt)
      end do
      ! 2D variables
      do iVar = 1,TBNTnVars2D
         iVarPnt = TBNTvar2Dpnt(iVar,iFrac)
         where (.not.excludeCell2D.and.abs(TBNTvar2D(:,iVarPnt))<epsFrac) TBNTvar2D(:,iVarPnt) = 0.e0
         where (.not.excludeCell2D) BULKvar2D(:,iVar) = BULKvar2D(:,iVar) + TBNTvar2D(:,iVarPnt)
      end do
   end do
   
   ! calculate relative fractions
   do iFrac = 1,TBNTnFractions
      ! 3D variables
      do iVar = 1,TBNTnVars3D
         iVarPnt = TBNTvar3Dpnt(iVar,iFrac)
         do i = 1,nWetCells
            if (excludeCell3D(i)) cycle
            if (BULKvar3D(i,iVar)<eps) then
               if (TBNTiniFrac<=0) then
                  TBNTrelVar3D(i,iVarPnt) = 1.e0/real(TBNTnFractions,dp)
               else
                  if (iFrac==TBNTiniFrac) then
                     TBNTrelVar3D(i,iVarPnt) = 1.e0
                  else
                     TBNTrelVar3D(i,iVarPnt) = 0.e0
                  end if
               end if
            else
               TBNTrelVar3D(i,iVarPnt) = TBNTvar3D(i,iVarPnt)/BULKvar3D(i,iVar)
               call check_relVar(TBNTvar(iVar)%name,TBNTsource(iFrac)%name,TBNTrelVar3D(i,iVarPnt),i,ierr)
            end if
         end do
      end do
      ! 2D variables
      do iVar = 1,TBNTnVars2D
         iVarPnt = TBNTvar2Dpnt(iVar,iFrac)
         do i = 1,nBotCells
            if (excludeCell2D(i)) cycle
            if (BULKvar2D(i,iVar)<eps) then
               if (TBNTiniFrac<=0) then
                  TBNTrelVar2D(i,iVarPnt) = 1.e0/real(TBNTnFractions,dp)
               else
                  if (iFrac==TBNTiniFrac) then
                     TBNTrelVar2D(i,iVarPnt) = 1.e0
                  else
                     TBNTrelVar2D(i,iVarPnt) = 0.e0
                  end if
               end if
            else
               TBNTrelVar2D(i,iVarPnt) = TBNTvar2D(i,iVarPnt)/BULKvar2D(i,iVar)
               call check_relVar(TBNTvar(TBNTnVars3D+iVar)%name,TBNTsource(iFrac)%name, &
                                 TBNTrelVar2D(i,iVarPnt),i,ierr)
            end if
         end do
      end do
   end do
   
   end subroutine calc_tbnt_relative_fractions
#endif
!=======================================================================
#ifndef TBNTonly_bulk_bud
!
! !INTERFACE:
   subroutine calc_tbnt_fractions
#define SUBROUTINE_NAME 'calc_tbnt_fractions'
!
! !DESCRIPTION:
! NUMERICAL CORE of TBNT: calculates fraction fluxes and variable fractions
!
! BASIC IDEA (for fluxes within one grid cell or between pelagial and sediment):
! "relative variable fraction" = "absolute fraction" / "bulk variable" 
! "fraction flux" = "relative variable fraction" * "bulk flux"
! "new fraction input variable" = "old fraction input variable" - "fraction flux"
! "new fraction output variable" = "old fraction output variable" + "fraction flux"
!
! FOR EXCHANGE FLUXES: "input variable" and "output variable" are identical, but locations are different, i.e.:
! "relative variable fraction" = "absolute fraction in origin cell" / "bulk variable in origin cell" 
! "fraction flux" = "relative variable fraction" * "bulk flux between the two cells"
! "new fraction variable in origin cell" = "old fraction variable in origin cell" - "fraction flux"
! "new fraction variable in target cell" = "old fraction variable in target cell" + "fraction flux"
!
! !USES:
   implicit none
!
! !LOCAL VARIABLES
   real(dp), parameter :: eps = 1.e-12
   integer             :: flxType, auxFlxType, flxComp, iiFlux
   integer             :: iVari, iVariPnt, iVaro, iVaroPnt
   integer             :: i, i_fr, iNeigh
   real(dp)            :: flxSign, flx
   real(dp)            :: flx2D(nBotCells), flx3D(nWetCells)
   logical             :: isVari3D, isVaro3D, isExchFlx, isLinked
!
!-----------------------------------------------------------------------
#include "call-trace.inc"

   isExchFlx = .false.
   flxComp = 0
   
   isLinked = .false.

   ! empty fraction fluxes
   tempflxFrac3D = 0.e0
   tempFlxFrac2D = 0.e0
   if (TBNTlinkedFluxes) then
      lnkflxFrac3D = 0.e0
      lnkFlxFrac2D = 0.e0
   end if
      
   ! empty cumulated change applied to variable s
   varChange3D = 0.e0
   varChange2D = 0.e0

   i_fr = 0
   iNeigh = 0
   flxComp = 0
   
   ! 3D FLUXES
   do iFlux = 1,TBNTnFluxes3D+TBNTnLinkedFluxes3D
      if (iFlux>TBNTnFluxes3D) then
         isLinked = .true.
      else
         isLinked = .false.
      end if
      ! get flux type and dimensionality
      flxType = TBNTflux(iFlux)%type
      auxFlxType = TBNTflux(iFlux)%auxFlux%type
      ! cycle over external fluxes
      if (flxType==TBNTrvdisFlux.or.auxFlxType==TBNTrvdisFlux.or. &
          flxType==TBNTatm3DFlux.or.auxFlxType==TBNTatm3DFlux) cycle
      ! get dimensionality of involved variables
      iVari = TBNTflux(iFlux)%varIn%iTBNT                           ! pointer to BULK input variable
      iVaro = TBNTflux(iFlux)%varOut%iTBNT                          ! pointer to BULK output variable
      isVari3D = .true.                                             ! input variable is 3D (default)
      isVaro3D = .true.                                             ! output variable is 3D (default)
      if (flxType==TBNTsed3DFlux.or.auxFlxType==TBNTsed3DFlux.or. & ! 3D flux between pelagic and sediment OR
          flxType==TBNTa2s3DFlux.or.auxFlxType==TBNTa2s3DFlux) then ! 3D flux between pelagic and atmosphere
         if (TBNTflux(iFlux)%varIn%type==TBNTpro2dVar.or. &
             TBNTflux(iFlux)%varIn%type==TBNTdia2dVar.or. &
             TBNTflux(iFlux)%varIn%type==TBNTdum2dVar) then
            isVari3D = .false. ! input variable is 2D
         end if
         if (TBNTflux(iFlux)%varOut%type==TBNTpro2dVar.or. &
             TBNTflux(iFlux)%varOut%type==TBNTdia2dVar.or. &
             TBNTflux(iFlux)%varOut%type==TBNTdum2dVar) then
            isVaro3D = .false. ! output variable is 2D
         end if
      end if
      ! UPDATE VARIABLE FRACTIONS DEPENDING ON DIMENSIONALITY OF INVOLVED VARIABLES
      flx3D = BULKflux3D(:,iFlux) ! get bulk flux
      if (isVari3D.and.isVaro3D) then
         ! ========================================================
         ! 3D flux, 3D input variable, 3D output variable (pelagic)
         ! ========================================================
         ! check if flux is exchange flux
         isExchFlx = .false.
         if (flxType==TBNTadvecFlux.or.auxFlxType==TBNTadvecFlux.or. &
             flxType==TBNTdiffuFlux.or.auxFlxType==TBNTdiffuFlux) isExchFlx = .true.
         do iFrac = 1,TBNTnFractions
            if (isLinked) then
               iFlxPnt = TBNTlinkedFlux3Dpnt(iFlux-TBNTnFluxes3D,iFrac) ! 3D flux pointer
            else
               iFlxPnt = TBNTflux3Dpnt(iFlux,iFrac) ! 3D flux pointer
            end if
            iVariPnt = TBNTvar3Dpnt(iVari,iFrac) ! 3D input variable pointer
            iVaroPnt = TBNTvar3Dpnt(iVaro,iFrac) ! 3D output variable pointer
            if (.not.isExchFlx) then ! PELAGIC FLUXES WITHIN A SINGLE CELL
               do i = 1,nWetCells
                  if (excludeCell3D(i).or.abs(flx3D(i))<eps) cycle
                  if (isLinked) then
                     lnkFlxFrac3D(i,iFlxPnt) = TBNTrelVar3D(i,iVariPnt) * flx3D(i)
                     cycle
                  end if
                  ! update variable fraction
                  tempFlxFrac3D(i,iFlxPnt) = TBNTrelVar3D(i,iVariPnt) * flx3D(i)
                  varChange3D(i,iVariPnt) = varChange3D(i,iVariPnt) - tempFlxFrac3D(i,iFlxPnt)
                  varChange3D(i,iVaroPnt) = varChange3D(i,iVaroPnt) + tempFlxFrac3D(i,iFlxPnt)
               end do
            else ! PELAGIC FLUXES BETWEEN TWO CELLS (TRANSPORT/MIXING)
               flxComp = TBNTflux(iFlux)%i_comp
               do i = 1,nWetCells
                  iNeigh = neighbours(i,flxComp)
                  if (iNeigh<=0.or.abs(flx3D(i))<eps) cycle             ! skip cells without neighbour and zero fluxes
                  if (excludeCell3D(i).and.excludeCell3D(iNeigh)) cycle ! skip cells outside the domain
                  flxSign = compSign(flxComp)
                  if (flxSign*flx3D(i)<0.e0) then
                     i_fr = i
                  else
                     i_fr = iNeigh
                  end if
                  if (isLinked) then
                     lnkFlxFrac3D(i,iFlxPnt) = TBNTrelVar3D(i_fr,iVaroPnt) * flx3D(i)
                     cycle
                  end if
                  ! update variable fraction
                  tempFlxFrac3D(i,iFlxPnt) = TBNTrelVar3D(i_fr,iVaroPnt) * flx3D(i)
                  flx = flxSign*tempFlxFrac3D(i,iFlxPnt)
                  if (.not.excludeCell3D(i).and..not.excludeCell3D(iNeigh)) then ! update core area
                     varChange3D(i     ,iVaroPnt) = varChange3D(i     ,iVaroPnt) + flx
                     varChange3D(iNeigh,iVaroPnt) = varChange3D(iNeigh,iVaroPnt) - flx
                  elseif (.not.excludeCell3D(i).and.excludeCell3D(iNeigh).and.flx<0.e0) then ! outward flow across boundary in positive direction
                     varChange3D(i,iVaroPnt) = varChange3D(i,iVaroPnt) + flx
                  elseif (excludeCell3D(i).and..not.excludeCell3D(iNeigh).and.flx>0.e0) then ! outward flow across boundary in negative direction
                     varChange3D(iNeigh,iVaroPnt) = varChange3D(iNeigh,iVaroPnt) - flx
                  end if
               end do
            end if
         end do
      elseif (isVari3D.and..not.isVaro3D) then
         ! =====================================================================
         ! 3D flux, 3D input variable, 2D output variable (pelagial to sediment)
         ! =====================================================================
         ! update variable fractions
         do iFrac = 1,TBNTnFractions
            if (isLinked) then
               iFlxPnt = TBNTlinkedFlux3Dpnt(iFlux-TBNTnFluxes3D,iFrac) ! 3D flux pointer
            else
               iFlxPnt = TBNTflux3Dpnt(iFlux,iFrac) ! 3D flux pointer
            end if
            iVariPnt = TBNTvar3Dpnt(iVari,iFrac) ! 3D input variable pointer
            iVaroPnt = TBNTvar2Dpnt(iVaro,iFrac) ! 2D output variable pointer
            do i = 1,nBotCells
               if (excludeCell2D(i)) cycle
               iPut = bottom2pelag(i)
               if (abs(flx3D(iPut))<eps) cycle
               if (isLinked) then
                  lnkFlxFrac3D(iPut,iFlxPnt) = TBNTrelVar3D(iPut,iVariPnt) * flx3D(iPut)
                  cycle
               end if
               ! update variable fraction
               tempFlxFrac3D(iPut,iFlxPnt) = TBNTrelVar3D(iPut,iVariPnt) * flx3D(iPut)
               varChange3D(iPut,iVariPnt)  = varChange3D(iPut ,iVariPnt) - tempFlxFrac3D(iPut,iFlxPnt)
               varChange2D(i   ,iVaroPnt)  = varChange2D(i    ,iVaroPnt) + tempFlxFrac3D(iPut,iFlxPnt)
            end do
         end do
      elseif (.not.isVari3D.and.isVaro3D) then
         ! =====================================================================================
         ! 3D flux, 2D input variable, 3D output variable (sediment to pelagial OR air-sea flux)
         ! =====================================================================================
         do iFrac = 1,TBNTnFractions
            if (isLinked) then
               iFlxPnt = TBNTlinkedFlux3Dpnt(iFlux-TBNTnFluxes3D,iFrac) ! 3D flux pointer
            else
               iFlxPnt = TBNTflux3Dpnt(iFlux,iFrac) ! 3D flux pointer
            end if
            iVariPnt = TBNTvar2Dpnt(iVari,iFrac) ! 2D input variable pointer
            iVaroPnt = TBNTvar3Dpnt(iVaro,iFrac) ! 3D output variable pointer
            if (flxType==TBNTsed3DFlux.or.auxFlxType==TBNTsed3DFlux) then ! sediment flux
               do i = 1,nBotCells
                  if (excludeCell2D(i)) cycle
                  iPut = bottom2pelag(i)
                  if (abs(flx3D(iPut))<eps) cycle
                  if (isLinked) then
                     lnkFlxFrac3D(iPut,iFlxPnt) = TBNTrelVar2D(i,iVariPnt) * flx3D(iPut)
                     cycle
                  end if
                  ! update variable fraction
                  tempFlxFrac3D(iPut,iFlxPnt) = TBNTrelVar2D(i,iVariPnt) * flx3D(iPut)
                  varChange2D(i   ,iVariPnt)    = varChange2D(i   ,iVariPnt) - tempFlxFrac3D(iPut,iFlxPnt)
                  varChange3D(iPut,iVaroPnt)    = varChange3D(iPut,iVaroPnt) + tempFlxFrac3D(iPut,iFlxPnt)
               end do
            elseif (flxType==TBNTa2s3DFlux.or.auxFlxType==TBNTa2s3DFlux) then ! air-sea flux
               do i = 1,nSurfCells
                  iPut = surf2pelag(i)
                  if (excludeCell3D(iPut).or.abs(flx3D(iPut))<eps.or.flx3D(iPut)>=eps) cycle
                  ! 3rd condition: flux into sea => already brought into system as atmospheric source (get_tbnt_atmos_fracs)
                  if (isLinked) then
                     lnkFlxFrac3D(iPut,iFlxPnt) = TBNTrelVar3D(iPut,iVariPnt) * flx3D(iPut)
                     cycle
                  end if
                  ! update variable fraction
                  tempFlxFrac3D(iPut,iFlxPnt) = TBNTrelVar3D(iPut,iVaroPnt) * flx3D(iPut)
                  varChange2D(i   ,iVariPnt)  = varChange2D(i   ,iVariPnt) - tempFlxFrac3D(iPut,iFlxPnt)
                  varChange3D(iPut,iVaroPnt)  = varChange3D(iPut,iVaroPnt) + tempFlxFrac3D(iPut,iFlxPnt)
               end do
            end if
         end do
      end if
   end do

   ! 2D FLUXES (sediment to pelagial or vice versa)
   do iFlux = 1,TBNTnFluxes2D+TBNTnLinkedFluxes2D
      if (iFlux>TBNTnFluxes2D) then
         isLinked = .true.
      else
         isLinked = .false.
      end if
      iiFlux = TBNTnFluxes3D+TBNTnLinkedFluxes3D + iFlux
      ! get flux ype and dimensionality of involved variables
      flxType = TBNTflux(iiFlux)%type
      auxFlxType = TBNTflux(iiFlux)%auxFlux%type
      iVari = TBNTflux(iiFlux)%varIn%iTBNT           ! pointer to BULK input variable
      iVaro = TBNTflux(iiFlux)%varOut%iTBNT          ! pointer to BULK output variable
      isVari3D = .true.                              ! input variable is 3D (default)
      isVaro3D = .true.                              ! input variable is 3D (default)
      if (TBNTflux(iiFlux)%varIn%type==TBNTpro2dVar.or. &
          TBNTflux(iiFlux)%varIn%type==TBNTdia2dVar.or. &
          TBNTflux(iiFlux)%varIn%type==TBNTdum2dVar) isVari3D = .false.  ! input variable is 2D
      if (TBNTflux(iiFlux)%varOut%type==TBNTpro2dVar.or. &
          TBNTflux(iiFlux)%varOut%type==TBNTdia2dVar.or. &
          TBNTflux(iiFlux)%varOut%type==TBNTdum2dVar) isVaro3D = .false. ! output variable is 2D
      ! UPDATE BULK VARIABLES AND VARIABLE FRACTIONS DEPENDING ON DIMENSIONALITY OF INVOLVED VARIABLES
      flx2D = BULKflux2D(:,iFlux) ! get bulk flux
      if (isVari3D.and..not.isVaro3D) then
         ! =====================================================================
         ! 2D flux, 3D input variable, 2D output variable (pelagial to sediment)
         ! =====================================================================
         do iFrac = 1,TBNTnFractions
            if (isLinked) then
               iFlxPnt = TBNTlinkedFlux2Dpnt(iFlux-TBNTnFluxes2D,iFrac) ! 2D flux pointer
            else
               iFlxPnt = TBNTflux2Dpnt(iFlux,iFrac) ! 2D flux pointer
            end if
            iVariPnt = TBNTvar3Dpnt(iVari,iFrac) ! 3D input variable pointer
            iVaroPnt = TBNTvar2Dpnt(iVaro,iFrac) ! 2D output variable pointer
            do i = 1,nBotCells
               if (excludeCell2D(i).or.abs(flx2D(i))<eps) cycle
               iPut = bottom2pelag(i)
               if (isLinked) then
                  lnkFlxFrac2D(i,iFlxPnt) = TBNTrelVar3D(iPut,iVariPnt) * flx2D(i)
                  cycle   
               end if
               ! update variables fractions
               tempFlxFrac2D(i,iFlxPnt) = TBNTrelVar3D(iPut,iVariPnt) * flx2D(i)
               varChange3D(iPut,iVariPnt) = varChange3D(iPut,iVariPnt) - tempFlxFrac2D(i,iFlxPnt)
               varChange2D(i   ,iVaroPnt) = varChange2D(i   ,iVaroPnt) + tempFlxFrac2D(i,iFlxPnt)
            end do
         end do
      elseif (.not.isVari3D.and.isVaro3D) then
         ! =====================================================================
         ! 2D flux, 2D input variable, 3D output variable (sediment to pelagial)
         ! =====================================================================
         do iFrac = 1,TBNTnFractions
            if (isLinked) then
               iFlxPnt = TBNTlinkedFlux2Dpnt(iFlux-TBNTnFluxes2D,iFrac) ! 2D flux pointer
            else
               iFlxPnt = TBNTflux2Dpnt(iFlux,iFrac) ! 2D flux pointer
            end if
            iVariPnt = TBNTvar2Dpnt(iVari,iFrac) ! 2D input variable pointer
            iVaroPnt = TBNTvar3Dpnt(iVaro,iFrac) ! 3D output variable pointer
            do i = 1,nBotCells
               if (excludeCell2D(i).or.abs(flx2D(i))<eps) cycle
               iPut = bottom2pelag(i)
               if (isLinked) then
                  lnkFlxFrac2D(i,iFlxPnt) = TBNTrelVar2D(i,iVariPnt) * flx2D(i)
                  cycle   
               end if
               ! update variables fractions
               tempFlxFrac2D(i,iFlxPnt) = TBNTrelVar2D(i,iVariPnt) * flx2D(i)
               varChange2D(i   ,iVariPnt) = varChange2D(i   ,iVariPnt) - tempFlxFrac2D(i,iFlxPnt)
               varChange3D(iPut,iVaroPnt) = varChange3D(iPut,iVaroPnt) + tempFlxFrac2D(i,iFlxPnt)
            end do
         end do
      elseif  (.not.isVari3D.and..not.isVaro3D) then
         ! ==============================================================================================
         ! 2D flux, 2D input variable, 2D output variable (within sediment OR sediment to dummy variable)
         ! ==============================================================================================
         do iFrac = 1,TBNTnFractions
            if (isLinked) then
               iFlxPnt = TBNTlinkedFlux2Dpnt(iFlux-TBNTnFluxes2D,iFrac) ! 2D flux pointer
            else
               iFlxPnt = TBNTflux2Dpnt(iFlux,iFrac) ! 2D flux pointer
            end if
            iVariPnt = TBNTvar2Dpnt(iVari,iFrac) ! 2D input variable pointer
            iVaroPnt = TBNTvar2Dpnt(iVaro,iFrac) ! 2D output variable pointer
            do i = 1,nBotCells
               if (excludeCell2D(i).or.abs(flx2D(i))<eps) cycle
               if (isLinked) then
                  lnkFlxFrac2D(i,iFlxPnt) = TBNTrelVar2D(i,iVariPnt) * flx2D(i)
                  cycle   
               end if
               ! update variables fractions
               tempFlxFrac2D(i,iFlxPnt) = TBNTrelVar2D(i,iVariPnt) * flx2D(i)
               varChange2D(i,iVariPnt) = varChange2D(i,iVariPnt) - tempFlxFrac2D(i,iFlxPnt)
               varChange2D(i,iVaroPnt) = varChange2D(i,iVaroPnt) + tempFlxFrac2D(i,iFlxPnt)
            end do
         end do
      elseif  (isVari3D.and.isVaro3D) then
         ! ==============================================================================================
         ! 2D flux, 3D input variable, 3D output variable (processes via sediment or air-sea exchange)
         ! ==============================================================================================
         do iFrac = 1,TBNTnFractions
            if (isLinked) then
               iFlxPnt = TBNTlinkedFlux2Dpnt(iFlux-TBNTnFluxes2D,iFrac) ! 2D flux pointer
            else
               iFlxPnt = TBNTflux2Dpnt(iFlux,iFrac) ! 2D flux pointer
            end if
            iVariPnt = TBNTvar3Dpnt(iVari,iFrac) ! 3D input variable pointer
            iVaroPnt = TBNTvar3Dpnt(iVaro,iFrac) ! 3D output variable pointer
            if (flxType==TBNTsed2DFlux.or.auxFlxType==TBNTsed2DFlux) then ! sediment flux
               do i = 1,nBotCells
                  if (excludeCell2D(i).or.abs(flx2D(i))<eps) cycle
                  iPut = bottom2pelag(i)
                  if (isLinked) then
                     lnkFlxFrac2D(i,iFlxPnt) = TBNTrelVar3D(iPut,iVariPnt) * flx2D(i)
                     cycle   
                  end if
                  ! update variables fractions
                  tempFlxFrac2D(i,iFlxPnt) = TBNTrelVar3D(iPut,iVariPnt) * flx2D(i)
                  varChange3D(iPut,iVariPnt) = varChange3D(iPut,iVariPnt) - tempFlxFrac2D(i,iFlxPnt)
                  varChange3D(iPut,iVaroPnt) = varChange3D(iPut,iVaroPnt) + tempFlxFrac2D(i,iFlxPnt)
               end do
            else if (flxType==TBNTa2s2DFlux.or.auxFlxType==TBNTa2s2DFlux) then ! air-sea flux
               do i = 1,nSurfCells
                  iPut = surf2pelag(i)
                  if (excludeCell3D(iPut).or.abs(flx2D(iPut))<eps.or.flx2D(i)>=eps) cycle
                  ! 3rd condition: flux into sea => already brought into system as atmospheric source (get_tbnt_atmos_fracs)
                  if (isLinked) then
                     lnkFlxFrac3D(iPut,iFlxPnt) = TBNTrelVar3D(iPut,iVariPnt) * flx3D(iPut)
                     cycle
                  end if
                  ! update variable fraction
                  tempFlxFrac2D(iPut,iFlxPnt) = TBNTrelVar3D(iPut,iVaroPnt) * flx2D(i)
                  varChange3D(iPut,iVariPnt)  = varChange3D(iPut,iVariPnt) - tempFlxFrac2D(i,iFlxPnt)
                  varChange3D(iPut,iVaroPnt)  = varChange3D(iPut,iVaroPnt) + tempFlxFrac2D(i,iFlxPnt)
               end do
            end if
         end do
      end if
   end do
   
   end subroutine calc_tbnt_fractions
#endif
!=======================================================================
#ifdef TBNTonly_bulk_bud
!
! !INTERFACE:
   subroutine calc_bulk_bud(iReadRec, ierr)
#define SUBROUTINE_NAME 'calc_bulk_bud'
!
! !DESCRIPTION:
! calculate budget for selected position
!
! !USES:
   implicit none
!
! !OUTPUT PARAMETERS:
   integer, intent(inout) :: ierr
   integer, intent(in   ) :: iReadRec
!
! !LOCAL VARIABLES:
   real(dp), parameter    :: eps = 1.e-6
   integer                :: iiFlux, flxType, auxFlxType, flxComp
   real(dp)               :: varNewDaily
   logical                :: isVar3D
   integer                :: i
#ifdef TBNTbulk_bud_out
   integer                  :: ios, system
   character(len=fileLen)   :: filename, budget_dir
   character(len=formatLen) :: addStr
   character(len=formatLen) :: compStr
   logical                  :: found, lExist
#endif
   
!-----------------------------------------------------------------------
#include "call-trace.inc"

#ifdef TBNTbulk_bud_out
   if (iReadRec==TBNTiStart) then
      ! create output directory
      ios = 0
      lExist = .false.
      budget_dir = trim(output_dir)//trim(BULKbudVar%name)'_bulk_budget'
      inquire(file=trim(budget_dir),exist=lExist)
      if (.not.lExist) ios = system('mkdir -p '//trim(budget_dir))
      if (ios/=0) then
         ierr = ierr + ios
         write(error_msg,'(a)')'Error creating output directory: '//trim(budget_dir)
         call stop_tbnt(error_msg,ierr)
      end if
      write(addStr,'(i20)')iBud
      addStr = '_i'//adjustl(addStr)
      filename = trim(budget_dir)//trim(BULKbudVar%name)//'_budget'//trim(addStr)//'.dat'
      open(BULKbudUnit,file=trim(filename),form='formatted',status='replace')
   end if
#endif

   isVar3D = .true.
   if (BULKbudVar%type==TBNTpro2dVar.or.BULKbudVar%type==TBNTdia2dVar) isVar3D = .false.
   if (isVar3D) then ! 3D variables (pelagic variables)
      varNewDaily = BULKvar3D(iBud,1)*volOld(iBud)
      if (iReadRec==TBNTiStart) then
         BULKbudVarNew = varNewDaily
#ifdef TBNTbulk_bud_out
         write(BULKbudUnit,'(i5,2x,3(E25.17))') TBNTiStart-1,varNewDaily,varNewDaily,varNewDaily
#endif
      end if
      ! 2D fluxes
      if (bottomCell(iBud)) then
         do i = 1,nBotCells
            if (bottom2pelag(i)==iBud) exit
         end do
         do iFlux = 1,BULKnBudFluxes2D
            iiFlux = BULKnBudFluxes3D + iFlux
            flxType = BULKbudFlux(iiFlux)%type
            auxFlxType = BULKbudFlux(iiFlux)%auxFlux%type
            if (flxType==TBNTatm2DFlux.or.auxFlxType==TBNTatm2DFlux.or. &
                flxType==TBNTa2s3DFlux.or.auxFlxType==TBNTa2s2DFlux) cycle
            varNewDaily   = varNewDaily   + BULKflux2D(i,iFlux)
            BULKbudVarNew = BULKbudVarNew + BULKflux2D(i,iFlux)
            print*,BULKflux2D(i,iFlux), i, bottom2pelag(i), count(bottom2pelag==iBud)
#ifdef TBNTbulk_bud_out
            if (iReadrec==TBNTiStart) then
               filename = trim(budget_dir)//trim(BULKbudFlux(iiFlux)%name)//trim(addStr)//'.dat'
               open(BULKbudUnit+iiFlux,file=trim(filename), form='formatted', access='append',status='replace')
            end if
            write(BULKbudUnit+iiFlux,'(i5,E25.17)') iReadRec, BULKflux2D(i,iFlux)
            if (iReadRec==TBNTiStart+TBNTnSteps-1) close(BULKbudUnit+iiFlux)
#endif
         end do
      else if (surfaceCell(iBud)) then
         do i = 1,nSurfCells
            if (surf2pelag(i)==iBud) exit
         end do
         do iFlux = 1,BULKnBudFluxes2D
            iiFlux = BULKnBudFluxes3D + iFlux
            flxType = BULKbudFlux(iiFlux)%type
            auxFlxType = BULKbudFlux(iiFlux)%auxFlux%type
            if (flxType==TBNTsed2DFlux.or.auxFlxType==TBNTsed2DFlux) cycle
            varNewDaily   = varNewDaily   + BULKflux2D(i,iFlux)
            BULKbudVarNew = BULKbudVarNew + BULKflux2D(i,iFlux)
#ifdef TBNTbulk_bud_out
            if (iReadrec==TBNTiStart) then
               filename = trim(budget_dir)//trim(BULKbudFlux(iiFlux)%name)//trim(addStr)//'.dat'
               open(BULKbudUnit+iiFlux,file=trim(filename), form='formatted', access='append',status='replace')
            end if
            write(BULKbudUnit+iiFlux,'(i5,E25.17)') iReadRec, BULKflux2D(i,iFlux)
            if (iReadRec==TBNTiStart+TBNTnSteps-1) close(BULKbudUnit+iiFlux)
#endif
         end do
      end if
      
      ! 3D fluxes
      do iFlux = 1,BULKnBudFluxes3D
         flxType = BULKbudFlux(iFlux)%type
         auxFlxType = BULKbudFlux(iFlux)%auxFlux%type
         flxComp = BULKbudFlux(iFlux)%i_comp
         if (flxType/=TBNTadvecFlux.and.auxFlxType/=TBNTadvecFlux.and. &
             flxType/=TBNTdiffuFlux.and.auxFlxType/=TBNTdiffuFlux) then
            varNewDaily   = varNewDaily   + BULKflux3D(iBud,iFlux)
            BULKbudVarNew = BULKbudVarNew + BULKflux3D(iBud,iFlux)
#ifdef TBNTbulk_bud_out
            if (iReadrec==TBNTiStart) then
               filename = trim(budget_dir)//trim(BULKbudFlux(iFlux)%name)//trim(addStr)//'.dat'
               open(BULKbudUnit+iFlux,file=trim(filename), form='formatted', access='append',status='replace')
            end if
            write(BULKbudUnit+iFlux,'(i5,E25.17)') iReadRec, BULKflux3D(iBud,iFlux)
            if (iReadRec==TBNTiStart+TBNTnSteps-1) close(BULKbudUnit+iFlux)
#endif
         else
#ifdef TBNTbulk_bud_out
            if (iReadrec==TBNTiStart) then
               write(compStr,'(i20)')flxComp
               ! flux at edge positive direction
               filename = trim(budget_dir)//trim(BULKbudFlux(iFlux)%name)//  &
                          '_pos-comp'//trim(adjustl(compStr))//trim(addStr)//'.dat'
               open(BULKbudUnit+iFlux,file=trim(filename), form='formatted', access='append',status='replace')
               ! flux at edge negative direction
               filename = trim(budget_dir)//trim(BULKbudFlux(iFlux)%name)//  &
                          '_neg-comp'//trim(adjustl(compStr))//trim(addStr)//'.dat'
               open(BULKbudUnit+BULKnBudFluxes+iFlux,file=trim(filename), form='formatted', access='append',status='replace')
            end if
#endif
            if (neighbours(iBud,flxComp)>0) then ! flux at edge in positive index direction
               varNewDaily   = varNewDaily   + compSign(flxComp) * BULKflux3D(iBud,iFlux)
               BULKbudVarNew = BULKbudVarNew + compSign(flxComp) * BULKflux3D(iBud,iFlux)
#ifdef TBNTbulk_bud_out
               write(BULKbudUnit+iFlux,'(i5,E25.17)') iReadRec, compSign(flxComp)*BULKflux3D(iBud,iFlux)
            else
               write(BULKbudUnit+iFlux,'(i5,E25.17)') iReadRec, 0.e0
#endif
            end if
#ifdef TBNTbulk_bud_out
            found = .false.
#endif
            do i = 1,nWetCells
               if (neighbours(i,flxComp)==iBud) then ! flux at edge in negative index direction
                  varNewDaily   = varNewDaily   - compSign(flxComp) * BULKflux3D(i,iFlux)
                  BULKbudVarNew = BULKbudVarNew - compSign(flxComp) * BULKflux3D(i,iFlux)
#ifdef TBNTbulk_bud_out
                  write(BULKbudUnit+BULKnBudFluxes+iFlux,'(i5,E25.17)') iReadRec, -compSign(flxComp)*BULKflux3D(i,iFlux)
                  found = .true.
#endif
                  exit
               end if
            end do
#ifdef TBNTbulk_bud_out
            if (.not.found) write(BULKbudUnit+BULKnBudFluxes+iFlux,'(i5,E25.17)') iReadRec, 0.e0
            if (iReadRec==TBNTiStart+TBNTnSteps-1) then
               close(BULKbudUnit+iFlux)
               close(BULKbudUnit+BULKnBudFluxes+iFlux)
            end if
#endif
         end if
      end do
      call get_bulk_vars(iReadRec+1, ierr)
      BULKvar3D(iBud,1) = BULKvar3D(iBud,1)*volNew(iBud)
#ifdef TBNTbulk_bud_out
      write(BULKbudUnit,'(i5,2x,3(E25.17))')iReadRec, BULKvar3D(iBud,1), varNewDaily, BULKbudVarNew
#endif
      ! budget check
      !eps = epsFac*BULKvar3D(iBud,1)
      if (BULKvar3D(iBud,1)>eps) then
         if (abs(1.0e0-varNewDaily/BULKvar3D(iBud,1)) > eps) then
            ierr = ierr + 1
            write(error_msg,'(a,i4,2(a,E22.15))')'Daily budget failed at step: ', iStep, new_line('a')// &
                                                 '  >  reference = ', BULKvar3D(iBud,1), new_line('a')// &
                                                 '  > calculated = ', varNewDaily
            call stop_tbnt(error_msg,ierr)
         elseif (abs(1.0e0-BULKbudVarNew/BULKvar3D(iBud,1)) > eps) then
            ierr = ierr + 1
            write(error_msg,'(a,i4,2(a,E22.15))')'Cumulated budget failed at step: ', iStep, new_line('a')// &
                                                 '  >  reference = ', BULKvar3D(iBud,1),     new_line('a')// &
                                                 '  > calculated = ', BULKbudVarNew
            call stop_tbnt(error_msg,ierr)
         end if
      end if
#ifndef TBNTnoVar2D
   else ! 2D variables (sediment variables)
      varNewDaily = BULKvar2D(iBud,1)*area2D(iBud)
      if (iReadRec==TBNTiStart) then
         BULKbudVarNew = varNewDaily
#ifdef TBNTbulk_bud_out
         write(BULKbudUnit,'(i5,2x,3(E25.17))') TBNTiStart-1,varNewDaily,varNewDaily,varNewDaily
#endif
      end if
      ! 2D fluxes
      do iFlux = 1,BULKnBudFluxes2D
         iiFlux = BULKnBudFluxes3D + iFlux
         flxType = BULKbudFlux(iiFlux)%type
         auxFlxType = BULKbudFlux(iiFlux)%auxFlux%type
         if (flxType==TBNTatm2DFlux.or.auxFlxType==TBNTatm2DFlux.or. &
             flxType==TBNTa2s2DFlux.or.auxFlxType==TBNTa2s2DFlux) cycle
         varNewDaily = varNewDaily + BULKflux2D(iBud,iFlux)
         BULKbudVarNew = BULKbudVarNew + BULKflux2D(iBud,iFlux)
#ifdef TBNTbulk_bud_out
         if (iReadrec==TBNTiStart) then
            filename = trim(budget_dir)//trim(BULKbudFlux(iiFlux)%name)//trim(addStr)//'.dat'
            open(BULKbudUnit+iiFlux,file=trim(filename), form='formatted', access='append',status='replace')
         end if
         write(BULKbudUnit+iiFlux,'(i5,E25.17)') iReadRec, BULKflux2D(iBud,iFlux)
         if (iReadRec==TBNTiStart+TBNTnSteps-1) close(BULKbudUnit+iiFlux)
#endif
      end do
      ! 3D fluxes
      i = bottom2pelag(iBud)
      do iFlux = 1,BULKnBudFluxes3D
         flxType = BULKbudFlux(iFlux)%type
         auxFlxType = BULKbudFlux(iFlux)%auxFlux%type
         if (flxType/=TBNTsed3DFlux.and.auxFlxType/=TBNTsed3DFlux) cycle
         varNewDaily = varNewDaily + BULKflux3D(i,iFlux)
         BULKbudVarNew = BULKbudVarNew + BULKflux3D(i,iFlux)
#ifdef TBNTbulk_bud_out
         if (iReadrec==TBNTiStart) then
            filename = trim(budget_dir)//trim(BULKbudFlux(iFlux)%name)//trim(addStr)//'.dat'
            open(BULKbudUnit+iFlux,file=trim(filename), form='formatted', access='append',status='replace')
         end if
         write(BULKbudUnit+iFlux,'(i5,E25.17)') iReadRec, BULKflux3D(i,iFlux)
         if (iReadRec==TBNTiStart+TBNTnSteps-1) close(BULKbudUnit+iFlux)
#endif
      end do
      call get_bulk_vars(iReadRec+1, ierr)
      BULKvar2D(iBud,1) = BULKvar2D(iBud,1)*area2D(iBud)
#ifdef TBNTbulk_bud_out
      write(BULKbudUnit,'(i5,2x,3(E25.17))')iReadRec, BULKvar2D(iBud,1), varNewDaily, BULKbudVarNew
#endif
      ! budget check
      !eps = epsFac*BULKvar2D(iBud,1)
      if (BULKvar2D(iBud,1)>eps) then
         if (abs(1.0e0-varNewDaily/BULKvar2D(iBud,1)) > eps) then
            ierr = ierr + 1
            write(error_msg,'(a,i4,2(a,E22.15))')'Daily budget failed at step: ', iStep, new_line('a')// &
                                                 '  >  reference = ', BULKvar2D(iBud,1), new_line('a')// &
                                                 '  > calculated = ', varNewDaily
            call stop_tbnt(error_msg,ierr)
         elseif (abs(1.0e0-BULKbudVarNew/BULKvar2D(iBud,1)) > eps) then
            ierr = ierr + 1
            write(error_msg,'(a,i4,2(a,E22.15))')'Cumulated budget failed at step: ', iStep, new_line('a')// &
                                                 '  >  reference = ', BULKvar2D(iBud,1),     new_line('a')// &
                                                 '  > calculated = ', BULKbudVarNew
            call stop_tbnt(error_msg,ierr)
         end if
      end if
#endif
   end if
#ifdef TBNTbulk_bud_out
   if (iReadRec==TBNTiStart+TBNTnSteps-1) close(BULKbudUnit)
#endif

   end subroutine calc_bulk_bud
#endif
!=======================================================================
#ifndef TBNTonly_bulk_bud
!
! !INTERFACE:
   subroutine check_relVar(varName,srcName,relFrac, ic, ierr)
!
! !DESCRIPTION:
! checks if relative variable fraction is within range 0 <= relFrac <= 1
!
! !USES:

   implicit none
!
! !OUTPUT PARAMETERS:
!
   character(len=nameLen), intent(in   ) :: varName, srcName
   real(dp)              , intent(inout) :: relFrac
   integer               , intent(in   ) :: ic
   integer               , intent(inout) :: ierr
!
! !LOCAL VARIABLES:
!
   real(dp), parameter      :: epsRel = 1.e-20
   character(len=formatLen) :: icStr
!-----------------------------------------------------------------------
   
!   if (relFrac>0..and.relFrac<epsRel) then
!      relFrac = 0.
!      return
!   elseif (relFrac>1..and.relFrac-1.<epsRel) then
!      relFrac = 1.
!      return
!   end if
   if (relFrac>=0.e0.and.relFrac<=1.e0) return
   
   ierr = ierr + 1
   write(icStr,'(i20)')ic
   write(error_msg,'(a,a,a,E23.16,a,i4)') trim(varName)//' - '//trim(srcName)//': relative fraction out of range'// &
                                          new_line('a')//'  >       i = ', trim(adjustl(icStr)),                    &
                                          new_line('a')//'  >   value = ', relFrac,                                 &
                                          new_line('a')//'  > substep = ', iSubStep
   call stop_tbnt(error_msg,ierr)
   
   end subroutine check_relVar
#endif
!=======================================================================
#ifndef TBNTonly_bulk_bud
!
! !INTERFACE:
   subroutine check_absVar(varName, srcName, absFrac, ic, ierr, absBulk)
!
! !DESCRIPTION:
! checks if absolute variable fraction absFrac > 0             (if absBulk is not provided) or
!                                      0 <= absFrac <= absBulk (if absBulk is provided)
!
! !USES:

   implicit none
!
! !OUTPUT PARAMETERS:
!
   character(len=nameLen), intent(in   )           :: varName, srcName
   real(dp)              , intent(in   )           :: absFrac
   integer               , intent(in   )           :: ic
   integer               , intent(inout)           :: ierr
   real(dp)              , intent(in   ), optional :: absBulk
!
! !LOCAL VARIABLES:
!
   real(dp), parameter      :: epsFac = 1.e-08
   character(len=formatLen) :: icStr
!-----------------------------------------------------------------------

   if (.not.present(absBulk)) then
      if (absFrac<0.e0) ierr = 1
!      if (absFrac>=0.e0) return
!      ierr = ierr + 1
!      write(icStr,'(i20)')ic
!      write(error_msg,'(a,a,a,E23.16,a,i4)') trim(varName)//' - '//trim(srcName)//': variable became negative'// &
!                                             new_line('a')//'  >       i = ', trim(adjustl(icStr)),              &
!                                             new_line('a')//'  >   value = ', absFrac,                           &
!                                             new_line('a')//'  > substep = ', iSubStep
   else
      if ((absFrac>=0.e0.and.absFrac<=absBulk).or.(absFrac>absBulk.and.absFrac-absBulk<epsFac*absBulk)) return
      ierr = ierr + 1
      write(icStr,'(i20)')ic
      write(error_msg,'(a,a,2(a,E23.16),a,i4)') trim(varName)//' - '//trim(srcName)//': absolute fraction out of range'// &
                                                new_line('a')//'  >       i = ', trim(adjustl(icStr)),                    &
                                                new_line('a')//'  >    frac = ', absFrac,                                 &
                                                new_line('a')//'  >    bulk = ', absBulk,                                 &
                                                new_line('a')//'  > substep = ', iSubStep
      call stop_tbnt(error_msg,ierr)
   end if

   end subroutine check_absVar
#endif
!=======================================================================
#ifndef TBNTonly_bulk_bud
!
! !INTERFACE:
   subroutine do_tbnt_targets
#define SUBROUTINE_NAME 'do_tbnt_targets'
!
! !DESCRIPTION:
! do spatial aggregation to target areas
!
! !USES:
   implicit none
!
! !LOCAL VARIABLES:
!
   real(dp), parameter :: eps = 1.e-20

   integer :: iiVar, iTarget, iArea, varType
   logical :: mask3D(nWetCells)
#ifndef TBNTnoVar2D
   logical :: mask2D(nBotCells)
#endif
!
!-----------------------------------------------------------------------
#include "call-trace.inc"
   
   TBNTtargetVals = 0.e0
   TBNTtargetWeights = 0.e0

   ! spatial aggregation
   do iTarget = 1,TBNTnTargetVars
      do iVar = 1,TBNTtargetVars(iTarget)%nVars
         iiVar = TBNTtargetVars(iTarget)%vars(iVar)%iTBNT
         varType = TBNTtargetVars(iTarget)%vars(iVar)%type
         do iArea = 1,TBNTnTargetAreas
#ifndef TBNTnoVar2D
            if (varType==TBNTpro2dVar.or.varType==TBNTdia2dVar.or.varType==TBNTdum2dVar) then
               mask2D = TBNTtargetAreas(iArea)%isPart2D
               ! bulk mass
               TBNTtargetWeights(iArea,iTarget) = TBNTtargetWeights(iArea,iTarget) + sum(BULKvar2D(:,iiVar),mask2D)
               ! fraction masses
               do iFrac = 1,TBNTnFractions-1
                  iVarPnt = TBNTvar2Dpnt(iiVar,iFrac)
                  TBNTtargetVals(iArea, iTarget, iFrac) = TBNTtargetVals(iArea, iTarget, iFrac) + &
                                                          sum(TBNTvar2D(:,iVarPnt),mask2D)
               end do
            else
#endif
               mask3D = TBNTtargetAreas(iArea)%isPart3D
               ! bulk mass
               TBNTtargetWeights(iArea,iTarget) = TBNTtargetWeights(iArea,iTarget) + sum(BULKvar3D(:,iiVar),mask3D)
               ! fraction masses
               do iFrac = 1,TBNTnFractions-1
                  iVarPnt = TBNTvar3Dpnt(iiVar,iFrac)
                  TBNTtargetVals(iArea, iTarget, iFrac) = TBNTtargetVals(iArea, iTarget, iFrac) + &
                                                          sum(TBNTvar3D(:,iVarPnt),mask3D)
               end do
#ifndef TBNTnoVar2D
            end if
#endif
         end do
      end do
   end do
   
   ! calculate relative values
   do iFrac = 1,TBNTnFractions-1
      where (TBNTtargetWeights>eps) 
         TBNTtargetVals(:,:,iFrac) = TBNTtargetVals(:,:,iFrac)/TBNTtargetWeights
      else where
         TBNTtargetVals(:,:,iFrac) = 0.e0
      end where
   end do
   where (TBNTtargetVals<eps) TBNTtargetVals = 0.e0

   end subroutine do_tbnt_targets
#endif
!=======================================================================
#ifndef TBNTonly_bulk_bud
#ifndef TBNTonline_budget_check
!
! !INTERFACE:
   subroutine final_bulk_check(ierr)
#define SUBROUTINE_NAME 'final_bulk_check'
!
! !DESCRIPTION:
! do spatial aggregation to target areas
!
! !USES:
   implicit none
!
! !OUTPUT PARAMETERS:
   integer, intent(inout) :: ierr
!
! !LOCAL VARIABLES:
!
   integer                             :: i
   real(dp), dimension(:), allocatable :: BULKcalc, BULKread
!
!-----------------------------------------------------------------------
#include "call-trace.inc"

   i = 0

   allocate ( BULKcalc(TBNTnVars3D+TBNTnVars2D), BULKread(TBNTnVars3D+TBNTnVars2D) )
   BULKcalc = 0.e0
   BULKread = 0.e0

   ! calculate BULK masses within whole model domain at end of TBNT run
   do iVar = 1,TBNTnVars3D ! 3D variables
      BULKcalc(iVar) = sum(BULKvar3D(:,iVar),.not.excludeCell3D)
   end do
   do iVar = 1,TBNTnVars2D ! 2D variables
      BULKcalc(TBNTnVars3D+iVar) = sum(BULKvar2D(:,iVar),.not.excludeCell2D)
   end do

   ! get BULK masses at end of TBNT run from input NetCDF file
   call get_vols(TBNTiEnd, ierr)
   call get_bulk_vars(TBNTiEnd+1, ierr)
   do iVar = 1,TBNTnVars3D ! 3D variables
      BULKread(iVar) = sum(BULKvar3D(:,iVar)/volOldIn*volNewIn,.not.excludeCell3D)
   end do
   do iVar = 1,TBNTnVars2D ! 2D variables
      BULKread(TBNTnVars3D+iVar) = sum(BULKvar2D(:,iVar),.not.excLudeCell2D)
   end do
   
   ! write BULK variables to log file
   call breakLine(60)
   write(tbnt_log_unit,'("TOTAL MASS OF BULK VARIABLES AT END OF CALCULATION")')
   write(tbnt_log_unit,'(a60)') repeat('-',60)
   write(tbnt_log_unit,'(4a15)')'VAR NAME', 'CALC (C)', 'READ (R)', '(C-R)/MEAN'
   do iVar = 1,TBNTnVars3D+TBNTnVars2D
      write(tbnt_log_unit,'(a15,3E15.6)') trim(TBNTvar(iVar)%name), BULKcalc(iVar), BULKread(iVar), &
                                           2.e0*(BULKcalc(iVar)-BULKread(iVar))/(BULKcalc(iVar)+BULKread(iVar))
   end do
   call breakLine(60)
   
   deallocate ( BULKcalc, BULKread )
      
   end subroutine final_bulk_check   
#endif
#endif
!=======================================================================
end module mod_etrac_main
