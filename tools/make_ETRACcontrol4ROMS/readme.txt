This tool is used to prepare the ASCII input files used by the ETRAC software. It consists of two different tools, which are described below.

For the compilation of either of the two just run "./compile.sh [name of the tool w/o file extension]".
Before compilation make sure that the correct compiler and compiler path is set in the compile.sh and in the makefile.

To run the compiled programs just execute "./[name of the tool]

1. make_nicemap4ROMS.f90
 - this program reads the ROMS grid file and creates an ASCII land-sea mask, which is used as template
   for the creation of the input files used by the 2nd script (make_ETRACcontrol4ROMS.f90)
 - its input is defined in the namelist file grid.nml

2. make_ETRACcontrol4ROMS.f90
 - this program reads the ROMS grid file, river forcing file, a HIS file, a nudging file (optional), and the ASCII nice-maps for ETRAC sources and creates standardized input files for the ETRAC application
 - its input is defined in the FORTRAN name list file setup.nml
 - the names of the grid, river forcing, HIS and nudging files are defined in the name lists
 - the names of the ASCII input files must be as follows:
    - etrac_openb_source_map_from_[DOMAIN].txt
    - etrac_openb_source_map_to_[DOMAIN].txt
    - etrac_exclude_map_[DOMAIN].txt
    - etrac_atmos_source_map_[DOMAIN].txt
    
   "[DOMAIN]" is a place holder for the "domain" name defined in the setup.nml
   
 - the output files are created in the "Output" subdirectory and must be copied to the user-defined ETRAC input folders
