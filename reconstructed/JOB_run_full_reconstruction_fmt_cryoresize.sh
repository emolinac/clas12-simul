#!/bin/bash

#SBATCH --account=clas12
#SBATCH --partition=production
#SBATCH --job-name=gemc-rec
#SBATCH --output=./out/%x.%j.array%a.out
#SBATCH --error=./err/%x.%j.array%a.err
#SBATCH --time=01:30:00
#SBATCH --mem=2G
#SBATCH --array=1-100

#--output=./out/%x.%j.array%a.out
#--error=./err/%x.%j.array%a.err


echo "This is JOB ${SLURM_ARRAY_JOB_ID} task ${SLURM_ARRAY_TASK_ID}"
echo "Its name is ${SLURM_JOB_NAME} and its ID is ${SLURM_JOB_ID}"

###########################################################################
###########################      FUNCTIONS      ###########################
###########################################################################

AZ_assignation(){
    # Function to assignate A and Z according to the target specified
    if [[ "$1" == "D" || "$1" == "D2" ]]
    then
	echo "Using Deuterium"
	A=2
	Z=1
    elif [[ "$1" == "C" ]]
    then
	echo "Using Carbon"
	A=12
	Z=6
    elif [[ "$1" == "Al" ]]
    then
	echo "Using Aluminum"
	A=27
	Z=13
    elif [[ "$1" == "Cu" ]]
    then
	echo "Using Copper"
	A=64
	Z=29
    elif [[ "$1" == "Sn" ]]
    then
	echo "Using Tin"
	A=119
	Z=50
    elif [[ "$1" == "Pb" ]]
    then
	echo "Using Lead"
	A=207
	Z=82
    else
	echo "No target input!"
	exit 1
    fi
    echo "A=${A} and Z=${Z}"    
}

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

###########################################################################
###########################     DIRECTORIES     ###########################
###########################################################################

cd ..
main_dir=$(pwd)
echo "main dire is ${main_dir}"
LEPTO_dir=~/Lepto64Sim/bin ## CHECK THIS DIRECTORY!
execution_dir=/volatile/clas12/emolinac
lepto2dat_dir=${main_dir}/thrown/lepto2dat
dat2tuple_dir=${main_dir}/thrown/dat2tuple
rec_utils_dir=${main_dir}/reconstructed/utils

out_dir_lepto=/work/clas12/rg-e/emolinac/lepto_11gev_fullchain
out_dir_recon=/work/clas12/rg-e/emolinac/hipo_11gev_fullchain

###########################################################################
###########################      VARIABLES      ###########################
###########################################################################

Nevents=750
target=D
torus=1
solenoid=-1
z_shift=0.
lD2_length=2.

id=${target}_${lD2_length}cmlD2_${SLURM_ARRAY_JOB_ID}${SLURM_ARRAY_TASK_ID}
temp_dir=${execution_dir}/${id}

###########################################################################
###########################       PREAMBLE      ###########################
###########################################################################

echo "Checking directories"
# Directory and files check
directories_check
executables_check

# Create folder in volatile to not interfere with other lepto executions
mkdir ${temp_dir}
cd ${temp_dir}

###########################################################################
###########################       LEPTO         ###########################
###########################################################################

echo "Running LEPTO"
# Prereqs setting
if [ -z "${CERN}" ]
then
    source ~/software/env_scripts/set_all.sh
fi

lepto_out=lepto_out_${id}

# Copy lepto executable to temp folder
cp ${LEPTO_dir}/lepto.exe ${temp_dir}/lepto_${id}.exe
echo "Copying LEPTO to temp dir"

# Assign the targets A and Z numbers
AZ_assignation ${target}
echo "${Nevents} ${A} ${Z}" > lepto_input.txt
echo "Target assignation done!"
# EXECUTE LEPTO
lepto_${id}.exe < lepto_input.txt > ${lepto_out}.txt
echo "LEPTO execution done"
# Transform lepto's output to dat files
cp ${lepto2dat_dir}/lepto2dat.pl ${temp_dir}/
perl lepto2dat.pl ${z_shift} < ${lepto_out}.txt > ${lepto_out}.dat
echo "lepto2dat done"
# Transform's dat files into ROOT NTuples
echo "dat2tuple start"
cp ${dat2tuple_dir}/bin/dat2tuple ${temp_dir}/
./dat2tuple ${lepto_out}.dat ${lepto_out}_ntuple.root
echo "Finished LEPTO"

###########################################################################
###########################       GEMC          ###########################
###########################################################################

# Prereqs setting
echo "Running GEMC"
if [ -z "${GEMC_DATA_DIR}" ]
then
    source /group/clas12/packages/setup.sh
    module load clas12
fi

gemc_out=gemc_out_${id}_${target}_s${solenoid}_t${torus}
gcard_name=clas12_fmt_cryoresize

# Transform lepto's output to LUND format
LUND_lepto_out=LUND${lepto_out}
cp ${rec_utils_dir}/leptoLUND.pl ${temp_dir}/
perl leptoLUND.pl ${z_shift} < ${lepto_out}.txt > ${LUND_lepto_out}.dat

# Copy the gcard you'll use into the temp folder and set the torus value
cp ${rec_utils_dir}/${gcard_name}.gcard ${temp_dir}/
cp /group/clas12/gemc/4.4.2/experiments/clas12/micromegas/micromegas__bank.txt ${temp_dir}/

# Copy the cryotarget model into the execution folder
cp -r ${rec_utils_dir}/targets/${lD2_length}cmlD2/ ${temp_dir}/

sed -i "s/TORUS_VALUE/${torus}/g" ${gcard_name}.gcard
sed -i "s/SOLENOID_VALUE/${solenoid}/g" ${gcard_name}.gcard
sed -i "s/CRYOTARGET_VARIATION/${lD2_length}cmlD2/g" ${gcard_name}.gcard

# EXECUTE GEMC
gemc ${gcard_name}.gcard -INPUT_GEN_FILE="LUND, ${LUND_lepto_out}.dat" -OUTPUT="evio, ${gemc_out}.ev" -USE_GUI="0"
echo "GEMC execution succesful!"

# Transform to HIPO
evio2hipo -t ${torus} -s ${solenoid} -r 11 -o ${gemc_out}.hipo -i ${gemc_out}.ev
rm ${gemc_out}.ev
echo "Evio 2 HIPO transformation successful"

# EXECUTE RECONSTRUCTION
cp ${rec_utils_dir}/clas12.yaml ${temp_dir}/
recon-util -y clas12.yaml -i ${gemc_out}.hipo -o ${gemc_out}.rec.hipo
echo "Reconstruction successful"
rm ${gemc_out}.hipo

# Move output to its folder
#mv ${lepto_out}.txt ${lepto_out}.dat ${lepto_out}_ntuple.root ${LUND_lepto_out}.dat ${out_dir_lepto}/
mv ${lepto_out}_ntuple.root ${out_dir_lepto}/
mv ${gemc_out}.rec.hipo ${out_dir_recon}/

# Remove folder
rm -rf ${temp_dir}
echo "Done!"
