#makefile for make_TBNTcontrol4ROMS and make_nicemap4ROMS

#-----------------------------------
# set default compiler and other
#-----------------------------------

ifndef MK_CONFIG_FILE
  ifdef FORTRAN_COMPILER
    MAKE_INCLUDE=${PWD}/make-config/$(FORTRAN_COMPILER).config
  endif
else
  MAKE_INCLUDE=${PWD}/make-config/$(MK_CONFIG_FILE)
endif
ifeq ($(wildcard $(MAKE_INCLUDE)),)
  include ${PWD}/make-config/dummy.config
else
  $(info ### make include : $(MAKE_INCLUDE))
  include $(MAKE_INCLUDE)
endif


#-----------------------------------
# list of allowed suffixes (first cleared)
#-----------------------------------
.SUFFIXES: .f90 .o

#-----------------------------------
# general compilation rules
#-----------------------------------

%.o : %.f90
	$(FC) $(FFLAGS) $(INCDIRS) $(LIBDIRS) $(EXTRA_LIBS) -c $<

#-----------------------------------
# targets
#-----------------------------------
.PHONY: run clean

$(EXECUTABLE): $(EXECUTABLE).o
	$(FC) -o $(EXECUTABLE) $(EXECUTABLE).o $(INCDIRS) $(LIBDIRS) $(EXTRA_LIBS)

$(EXECUTABLE).f90:
	$(FC) $(FFLAGS) $(INCDIRS) $(LIBDIRS) $(EXTRA_LIBS) -c $(EXECUTABLE).f90

run: $(EXECUTABLE)
	./$(EXECUTABLE)

clean:
	rm -f $(EXECUTABLE) *.o *.mod

#-----------------------------------
# dependencies
#-----------------------------------
$(EXECUTABLE).o: $(EXECUTABLE).f90
$(EXECUTABLE): $(EXECUTABLE).o
