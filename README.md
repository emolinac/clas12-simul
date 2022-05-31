# CLAS12 simulations
Simple package set to obtain simulations in JLAB's farm.
## Thrown (LEPTO)
...Obtain the text files bla bla bla
- **lepto2dat** : perl code modified from W. Brooks' original code. It transforms the output of lepto to a .dat file that can be easily read by *dat2tuple*.
    - *usage* : perl lepto2dat.pl z_shift < lepto_original.out > lepto_out.dat
- **dat2tuple** : C++ code that takes a file formated by *lepto2dat* and outputs a root file with tuples for the electrons, hadrons, and raw.
    - *usage* :
       1. Execute *make*
       2. In bin folder: *./dat2tuple <input_file_name> <output_file_name>*
## Reconstructed (GEMC)
W.I.P.