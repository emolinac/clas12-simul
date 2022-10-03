#!/bin/bash
# WARNING! : Do not move this script from this folder!

################################################################################################
########################           Checking variables          #################################
################################################################################################
directories_check(){
    # checking execution directories
    if [[ ! -d ${LEPTO_dir} || ! -d ${execution_dir} || ! -d ${lepto2dat_dir} || ! -d ${dat2tuple_dir} || ! -d ${rec_utils_dir} || ! -d ${out_dir_recon} || ! -d ${out_dir_lepto} ]]
    then
	echo "One of the necessary directories does not exist."
	exit 1
    fi
}
executables_check(){
    # checking executables existence
    if [[ ! -f ${dat2tuple_dir}/bin/dat2tuple || ! -f ${LEPTO_dir}/lepto.exe ]]
    then
	echo "One of the necessary executables does not exist."
	exit 1
    fi
}
errout_check(){
    # checking execution directories
    if [[ ! -d ${main_dir}/reconstructed-double-target-OFF/err || ! -d ${main_dir}/reconstructed-double-target-OFF/out ]]
    then
	echo "Making log out directories!"
	mkdir ${main_dir}/reconstructed-double-target-OFF/err
	mkdir ${main_dir}/reconstructed-double-target-OFF/out
    fi
}
################################################################################################
############################# Hermes-like script hehe ##########################################
################################################################################################
Njobs=1000
Njobsmax=25

################################################################################################
########################               Directories              ################################
################################################################################################
cd ..
main_dir=$(pwd)
LEPTO_dir=~/Lepto64Sim/bin ## CHECK THIS DIRECTORY!
execution_dir=/volatile/clas12/emolinac
lepto2dat_dir=${main_dir}/thrown/lepto2dat
dat2tuple_dir=${main_dir}/thrown/dat2tuple
rec_utils_dir=${main_dir}/reconstructed-double-target-OFF/utils

out_dir_lepto=/volatile/clas12/emolinac/lepto_files
out_dir_recon=/volatile/clas12/emolinac/hipo_files

################################################################################################
########################               Simul Specs             #################################
################################################################################################
# Use    : Sets the number of events per job (#electrons)
# Values : This is the sweet spot between quantity and performance
Nevents=500

# Use    : Sets the scaling of the magnetic fields
# Values : From -1 to 1
torus=-1
solenoid=-1

# Use    : Determine ID, vertex, and set the u/d ratio in LEPTO
# Values : D2
target=D2

# Use    : Determine the dt configuration present
# Values : lD2
target_variation=lD2

# Use    : Determine the cryotarget length
# Values : 2, 3, 5 (just the number! do not write points or things like that)
lD2_length=3

################################################################################################
########################                SHOWTIME               #################################
################################################################################################
directories_check
executables_check
errout_check

cd ${main_dir}/reconstructed-double-target-OFF
sbatch --array=1-${Njobs}%${Njobsmax} run_full_reconstruction_fmt_cryoresize_fullD2vertex.sh \
${LEPTO_dir} ${execution_dir} ${lepto2dat_dir} ${dat2tuple_dir} ${rec_utils_dir} ${out_dir_lepto} ${out_dir_recon} \
${Nevents} ${torus} ${solenoid} ${target} ${target_variation} ${lD2_length}
