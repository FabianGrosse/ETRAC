program ETRAC 
#define PROGRAM_NAME 'ETRAC'
!
! !USES:
   use mod_etrac_common
   use mod_etrac_init,  only: init_tbnt
#ifndef TBNTonly_bulk_bud
   use mod_etrac_main,  only: do_tbnt, get_bulk_vars, get_vols, do_tbnt_targets, get_tbnt_vars, final_bulk_check
#else
   use mod_etrac_main,  only: do_tbnt, get_bulk_vars, get_vols
#endif
   use mod_etrac_output
!
! !LOCAL VARIABLES
   integer :: iRec
!
!=======================================================================
! PROGRAM

   ierr = 0

   ! initialise TBNT
   call init_tbnt(ierr)

#ifndef TBNTonly_bulk_bud
   ! get bulk variable at start of simulation
   call get_vols(TBNTiStart, ierr)      ! get volume at beginning
   call get_bulk_vars(TBNTiStart, ierr) ! get bulk variables from BULK input file
   call get_tbnt_vars                   ! calculate initial distribution of fraction variables
   
   ! initialise output
   if (TBNTrelFracOutStep>=0.or.TBNTabsFracOutStep>=0.or.TBNTtargetOut.or.TBNTlinkedFluxes) call init_tbnt_output(ierr)
   
   write(tbnt_log_unit, '(a)') '==========================================='//new_line('a')// &
                               '============= TBNT SIMULATION ============='//new_line('a')// &
                               '==========================================='//new_line('a')
#else
   write(tbnt_log_unit, '(a)') '==========================================='//new_line('a')// &
                               '========== CALCULATE BULK BUDGET =========='//new_line('a')// &
                               '==========================================='//new_line('a')
#endif
   
   ! time loop based on time step of input NetCDF file
   do iStep = 1,TBNTnSteps
   
      iRec = TBNTiStart + iStep - 1
      
      ! do TBNT calculations
      call do_tbnt(iRec, ierr)

#ifndef TBNTonly_bulk_bud
      ! do spatial aggregation
      if (TBNTtargetOut.and.mod(iStep*TBNTstep,iMinPerDay)==0) call do_tbnt_targets

      ! write TBNT output
      if (TBNTrelFracOutStep>=0.or.TBNTabsFracOutStep>=0.or.TBNTtargetOut.or.TBNTlinkedFluxes) call do_tbnt_output(ierr)
#endif
   end do

#ifndef TBNTonly_bulk_bud
#ifndef TBNTonline_budget_check
   call final_bulk_check(ierr)
#endif
   
   ! finalize output
   if (TBNTrelFracOutStep>=0.or.TBNTabsFracOutStep>=0.or.TBNTtargetOut.or.TBNTlinkedFluxes) call finalize_tbnt_output(ierr,.true.)
#endif

   ! finalise TBNT
   call finalize_tbnt

   contains
!
!=======================================================================
!
! !INTERFACE:
   subroutine finalize_tbnt
#define SUBROUTINE_NAME 'finalize_tbnt'
!
! !DESCRIPTION:
! deallocate all allocated fields
!
! !USES:
!
   implicit none
!
! !LOCAL VARIABLES:
!
!-----------------------------------------------------------------------
#include "call-trace.inc"

   deallocate ( TBNTflux, TBNTvar, volOldIn, volOld, volNewIn, volNew )
   deallocate ( bottom2pelag, surf2pelag, compSign )
   deallocate ( BULKflux2D,  BULKflux3D )
#if !defined TBNTnoVar2D || !defined TBNTmass_fluxes
   if (allocated(area2D)) deallocate ( area2D )
#endif
#ifndef TBNTmass_fluxes
   deallocate ( area3D )
#endif

#ifndef TBNTnoVar2D
   if (allocated(BULKvar2D)) deallocate ( BULKvar2D )
#endif
   if (allocated(BULKvar3D)) deallocate ( BULKvar3D )
   
#ifdef TBNTonly_bulk_bud
   write(tbnt_log_unit,'(a)')new_line('a')// &
                             '==========================================='//new_line('a')// &
                             '==== finalized BULK BUDGET: normal end ===='//new_line('a')// &
                             '==========================================='
#else
   if (TBNTnRiverFractions>0) deallocate( iRvDisFlux )
   
   deallocate ( TBNTsource, TBNTflux2Dpnt, TBNTflux3Dpnt )
   deallocate ( TBNTvar3D, TBNTrelVar3D, TBNTvar3Dpnt )
#ifndef TBNTnoVar2D
   deallocate ( TBNTvar2D, TBNTrelVar2D, TBNTvar2Dpnt )
#endif

   if (allocated(TBNTflux2D)) deallocate ( TBNTflux2D )
   if (allocated(TBNTflux3D)) deallocate ( TBNTflux3D )
   deallocate ( tempFlxFrac2D, tempFlxFrac3D, varChange3D )
#ifndef TBNTnoVar2D
   deallocate ( varChange2D )
#endif

   if (TBNTlinkedFluxes) then
      deallocate ( lnkFlxFrac2D    , lnkFlxFrac3D )
      deallocate ( TBNTlinkedFlux2D, TBNTlinkedFlux3D )
   end if

   write(tbnt_log_unit,'(a)')new_line('a')// &
                             '==========================================='//new_line('a')// &
                             '======== finalized TBNT: normal end ======='//new_line('a')// &
                             '==========================================='
#endif
   close(tbnt_log_unit)

   end subroutine finalize_tbnt
!=======================================================================
end program ETRAC
