#################################################
#         makefile for ETRAC software          #
#################################################
#
# set name for executable
#------------------------------------------------------------------------
ifndef EXECUTABLE
   EXECUTABLE=ETRAC
endif
#
# set default compiler
#------------------------------------------------------------------------
ifndef FORTRAN_COMPILER
    FORTRAN_COMPILER=ifort
endif
#
# include compiler specific parts into the makefile
#------------------------------------------------------------------------
ifndef MK_CONFIG_FILE
  ifdef FORTRAN_COMPILER
    MAKE_INCLUDE=${PWD}/make-config/$(FORTRAN_COMPILER).config
  endif
else
  MAKE_INCLUDE=${PWD}/make-config/$(MK_CONFIG_FILE)
endif

$(info =========================================================)
ifeq ($(wildcard $(MAKE_INCLUDE)),)
  include ${PWD}/make-config/dummy.config
else
  $(info ### make include : $(MAKE_INCLUDE))
  include $(MAKE_INCLUDE)
endif
$(info =========================================================)
#
#
#------------------------------------------------------------------------
# preprocessor defines for program control
#------------------------------------------------------------------------
   DEFINES   += -DTBNTconvert_3Dto1D#        # use 3D input and convert to 1D vector => this is currently mandatory
   DEFINES   += -DTBNTmass_fluxes#           # fluxes in input file are mass fluxes, otherwise fluxes are vertically integrated, i.e., per unit area
#   DEFINES   += -DTBNTnoVar2D#               # model contains no 2D variables
#   DEFINES   += -DTBNTswitchNS#              # enable if indexing in y-direction in grid information files is opposite to indexing in NetCDF files (S-N)
#   DEFINES   += -DTBNTonline_budget_check#   # check bulk variable budget after each day
#   DEFINES   += -DTBNTonly_bulk_bud#         # only calculate balance for BULK variables and fluxes => no ETRAC calculation
#   DEFINES   += -DTBNTbulk_bud_out#          # write ascii output for BULK balance fluxes; only applicable if TBNTonly_bulk_bud is defined
#
#------------------------------------------------------------------------
#     adjust to IBM compiler syntax
#     -------------------------------------
ifeq ($(IBM_compiler),TRUE)
   DEFINES:=-WF,"$(DEFINES)"
endif
#$(info $(DEFINES))
#------------------------------------------------------------------------
#
###############################################
# SPECIFY COMPILATION OBJECTS                 #
###############################################

OBJS = \
       etrac_common.o   \
       etrac_init.o     \
       etrac_main.o     \
       etrac_output.o   \
       $(EXECUTABLE).o

###############################################
# BASIC COMPILE INSTRUCTIONS AND DEPENDENCIES #
###############################################
#.SUFFIXES:
#.SUFFIXES: .F90 .f90 .o

%.f90: %.F90
	$(CPP) $(CPPFLAGS) $(DEFINES) $< > $*.f90
	cp *.f90 ../$(BUILD_DIR) 

%.o : %.f90
	$(FC) $(FFLAGS) $(INCDIRS) $(LIBDIRS) $(EXTRA_LIBS) -c $<

$(EXECUTABLE):   $(OBJS)
	@echo "===== linking"
	$(FC) -o $@ $(LDFLAGS) $(INCDIRS) $(OBJS) $(LIBDIRS) $(EXTRA_LIBS)

clean:
	rm -f $(EXECUTABLE) *~ *.o *.f90 i.* *.mod *.lst

# dependencies
#--------------
etrac_init.o:              etrac_common.o
etrac_output.o:            etrac_common.o etrac_init.o
etrac_main.o:              etrac_common.o etrac_init.o
ETRAC.o:                   etrac_common.o etrac_init.o etrac_main.o etrac_output.o
