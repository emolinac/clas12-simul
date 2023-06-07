# CLAS12 simulations
This branch is designed to obtain generated events with LEPTO at a 10.6GeV beam energy

## Prereqs
To use this software it is necesary to install https://github.com/utfsm-eg2-data-analysis/clas_software in your top directory.

## Necessary Steps Before Use
- In **JOB_run_lept_fullchain.sh** change the **out_dir** variable.
- 

## Key Codes
- **JOB_run_lepto_fullchain.sh** : bash script to be submitted in JLAB's farm to obtain the txt outputs used by lepto2dat. Modify the directories according to your needs
- **lepto2dat** : perl code modified from W. Brooks' original code. It transforms the output of lepto to a .dat file that can be easily read by *dat2tuple*.
    - *usage* : perl lepto2dat.pl z_vertex < lepto_original.out > lepto_out.dat
- **dat2tuple** : C++ code that takes a file formated by *lepto2dat* and outputs a root file with tuples for the electrons, hadrons, and raw.
    - *usage* :
       1. Execute *make*
       2. In bin folder: *./dat2tuple <input_file_name> <output_file_name>*
