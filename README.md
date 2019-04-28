# ETRAC
Element TRACing software for numerical models



Introduction

The ETRAC software allows for the tracing of elements from individual sources through all physical and biogeochemical processes represented by the physical-biogeochemical model, to which it is applied. The software was initially developed for the ECOHAM North Sea model (Große et al., 2017; doi: 10.3389/fmars.2017.00383), however, with the intention to make it applicable also to other models. It requires specific model output, which needs to be implemented to the model to be used if it is not implemented already. For ROMS, this output has been implemented based on ROMS 3.7. Detailed information on the theory and implementation of the element tracing, and the model output requirements can be found in the documentation.



Documentation

A detailed documentation can be found here: /docu/ETRAC_documentation.pdf (description of a test case still missing)


Source code (and compilation)

The ETRAC software is written in Fortran and requires Fortran 2000 compatibility of the compiler.
The source code and the scripts to compile the software can be found in /software.

To compile ETRAC, change into /software and run ./compile.sh

You will need a make configuration file for your architecture in /software/src/make-config/



Running ETRAC

To run ETRAC modify the following files according to your setup (see comments in files).

In the directory setup_files/model_setup
 - model_dummy_vars.txt     => list of variables that are used internally by the model and can be calculated from other available variables by simple multiplication with a constant factor, but are not stored in the model output
 - model_fluxes.txt         => a list of all model fluxes and the variables they link to each other
 - model_grid.txt           => the model grid described in terms of 1-D vectors of neighboring cells for each spatial dimension
 - model_iDep.txt           => an ASCII "map" of the model grid with the number of wet cells per horizontal grid cell
 - model_rivers.txt         => a list of all river input locations in the 1-D indexing scheme; river names used in this file are used to for definition of river source groups
 
The files model_grid.txt, model_idep.txt and model_rivers.txt can be generated with the tools provided in tools/make_ETRACcontrol4ROMS. All these files are mandatory for running ETRAC.
 
 
In the directory setup_files/etrac_setup
 - etrac_openb_source.txt   => a list of open boundary input sources defined as grid cell indices of source and target cells using the 1-D indexing scheme
 - etrac_atmos_source.txt   => a list of open boundary input sources defined as grid cell indices using the 1-D indexing scheme
 - linked_fluxes.txt        => a list of fluxes for which relative contributions can be calculated based on selected state variables (see Große et al. (2017))
 - target_areas.txt         => list of target/integration regions defined as 1-D index vectors; simple csv output will be created for these regions
 - target_variables.txt     => list of target variables to be stored in csv files for target areas
 
The files etrac_openb_source.txt and etrac_atmos_source.txt can be generated with the tools provided in tools/make_ETRACcontrol4ROMS
Only etrac_openb_source.txt is mandatory for running ETRAC. 
 
In the ETRAC main directory:
 - etrac_set_BASE.nml       => file in which your ETRAC setup is defined. It points to the setup files listed above. It contains placeholder strings
 - ETRAC_SINGLE_BASE.slurm  => the batch job script used to run ETRAC with SLURM; it contains placeholder strings
 - run_ETRAC_CHAIN.sh       => the main script to run ETRAC; it replaces the placeholders in the two previous files and deals with running a sequence of jobs etc.
 
After having prepared all your files, make sure that you have execution permission for run_ETRAC_CHAIN.sh.
Then run: nohup ./run_ETRAC_CHAIN.sh &



Tools for preparation of ETRAC setup files

ETRAC requires a number of specific input files that define the model structure (e.g. grid, state variable interactions) and the tracing setup. These files must comply with a defined format, and their preparation is model-dependent, as different models use different files, e.g. to describe their grids. The set of tools for preparing the input files for a ROMS-based application can be found in /tools/make_ETRACcontrol4ROMS. A short description of the tools is provided in the readme.txt files in the directory.
