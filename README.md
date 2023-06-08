# CLAS12 simulations
This branch is designed to obtain generated events with LEPTO at a 10.6GeV beam energy

## Prereqs
To use this software it is necesary to install https://github.com/utfsm-eg2-data-analysis/clas_software in your top directory (see line 100 of **JOB_run_lepto_fullchain.sh**)

## Necessary Steps Before Use
- In **JOB_run_lept_fullchain.sh** change the **out_dir** variable.
- In **JOB_run_lept_fullchain.sh**, specify the number of jobs in line 9; specifically, change the second number.

## Key Codes
- **send_jobs.sh**: bash script to submit jobs
- **JOB_run_lepto_fullchain.sh** : job script to be submitted to JLAB's farm.
- **lepto2dat** : perl code modified from W. Brooks' original code. It transforms the output of lepto to a .dat file that can be easily read by *dat2tuple*.
    - *usage* : perl lepto2dat.pl z_vertex < lepto_original.out > lepto_out.dat
- **dat2tuple** : C++ code that takes a file formated by *lepto2dat* and outputs a root file with tuples for the electrons, hadrons, and raw.
    - *usage* :
       1. Execute *make*
       2. In bin folder: *./dat2tuple <input_file_name> <output_file_name>*
