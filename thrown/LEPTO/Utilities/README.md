## Brief description of added codes:
- leptoLUND.pl : Use it to format lepto's output to LUND format, which is used by GEMC
- leptodata.pl : Use it to format lepto's output to dat format. The idea is to then use MACRO_dat2root.cpp
- MACRO_dat2root.cpp : Use it to transform lepto's dat output into a ROOT's TTree

### Observation:
The codes to format lepto's output works with energy, momentum and PIDs. Calculations of hadronic or leptonic variables are not in the scope of these codes.