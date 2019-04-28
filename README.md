# ETRAC
Element TRACing software for numerical models



Introduction

The ETRAC software is a software that allows for the tracing of elements from individual sources through all physical and biogeochemical processes represented by the physical-biogeochemical model, to which it is applied. The software was initially developed for the ECOHAM North Sea model (Große et al., 2017; doi: 10.3389/fmars.2017.00383), however, with the intention to make it applicable also to other models. It requires specific model output, which needs to be implemented to the model to be used if it is not implemented already. For ROMS, this output has been implemented based on ROMS 3.7. Detailed information on the theory and implementation of the element tracing, and the model output requirements can be found in the documentation.



Documentation

A detailed documentation can be found here: /docu/ETRAC_documentation.pdf (description of a test case still missing)


Source code (and compilation)

The ETRAC software is written in Fortran and requires Fortran 2000 compatibility of the compiler.
The source code and the scripts to compile the software can be found in /software.
To compile ETRAC, change into /software and run ./compile.sh
You will need a make configuration file for your architecture in /software/src/make-config/



Running ETRAC

To run ETRAC just modify the files the following files according to your setup (see comments in files):
 - run_ETRAC_CHAIN.sh
 - ETRAC_SINGLE_BASE.slurm
 - etrac_set_BASE.nml

Make sure that you have execution permission for run_ETRAC_CHAIN.sh.
Then run: nohup ./run_ETRAC_CHAIN.sh &



Tools for preparation of ETRAC setup files

ETRAC requires a number of specific input files that define the model structure (e.g. grid, state variable interactions) and the tracing setup. These files must comply with a defined format, and their preparation is model-dependent, as different models use different files, e.g. to describe their grids. The set of tools for preparing the input files for a ROMS-based application can be found in /tools/make_ETRACcontrol4ROMS. A short description of the tools is provided in the readme.txt files in the directory.
