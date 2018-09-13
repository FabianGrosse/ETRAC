#ifdef PRINT_CALL_TRACE
! horizontal tab = char(9), vertical tab = char(11)

#  ifdef FUNCTION_NAME
   write (*,'(''###    function: '',a20,''  myID='',i2.2)') FUNCTION_NAME, myid
#  undef FUNCTION_NAME
#  elif defined SUBROUTINE_NAME
   if(ierr/=0) return
   write (*,'(''###  subroutine: '',a20,''  myID='',i2.2,''  ierr='',i1)') SUBROUTINE_NAME, myid, ierr
#  undef SUBROUTINE_NAME
#  elif defined PROGRAM_NAME
   write (*,'(''###     program: '',a20,''  myID='',i2.2)') PROGRAM_NAME, myid
#  undef PROGRAM_NAME
#  endif

#else

#  ifdef FUNCTION_NAME
#  undef FUNCTION_NAME
#  elif defined SUBROUTINE_NAME
   if(ierr/=0) return
#  undef SUBROUTINE_NAME
#  elif defined PROGRAM_NAME
#  undef PROGRAM_NAME
#  endif

#endif
