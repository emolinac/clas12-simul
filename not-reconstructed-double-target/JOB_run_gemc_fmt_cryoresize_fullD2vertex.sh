#!/bin/bash

#SBATCH --account=clas12
#SBATCH --partition=production
#SBATCH --job-name=gemc
#SBATCH --output=./out/%x.%j.array%a.out
#SBATCH --error=./err/%x.%j.array%a.err
#SBATCH --time=01:30:00
#SBATCH --mem=2G
#SBATCH --array=1-10

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
    if [[ ! -d ${LEPTO_dir} || ! -d ${execution_dir} || ! -d ${lepto2dat_dir} || ! -d ${dat2tuple_dir} || ! -d ${not-rec_utils_dir} || ! -d ${out_dir_gemc} || ! -d ${out_dir_lepto} ]]
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
not-rec_utils_dir=${main_dir}/not-reconstructed-double-target/utils

out_dir_lepto=/work/clas12/rg-e/emolinac/lepto_11gev_fullchain
out_dir_gemc=/work/clas12/rg-e/emolinac/hipo_11gev_notrecon

###########################################################################
###########################      VARIABLES      ###########################
###########################################################################
Nevents=750
torus=1
solenoid=-1
# Use    : determine ID, vertex, and set the u/d ratio in LEPTO
# Values : D2, C, Al, Cu, Sn, Pb
target=D2
# Use    : determine the dt configuration present
# Values : lD2, eg2-X, eg2-X-lD2, where X = {C,Al,Cu,Sn,Pb}
target_variation=lD2
lD2_length=2 # write just the number!
cryotarget_variation=${lD2_length}cmlD2

id=${target}_${cryotarget_variation}_${SLURM_ARRAY_JOB_ID}${SLURM_ARRAY_TASK_ID}
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

# Setting the vertex
cp ${not-rec_utils_dir}/*.py .
rdm=$(python random_gen.py)
z_vertex=$(python vertex.py ${lD2_length} ${rdm} ${target})

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
perl lepto2dat.pl ${z_vertex} < ${lepto_out}.txt > ${lepto_out}.dat
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

gemc_out=gemc_out_Zeq8_${id}_${target_variation}_s${solenoid}_t${torus}
gcard_name=clas12_fmt_cryoresize

# Transform lepto's output to LUND format
# Yes, It is necessary to specify the same z_vertex again
LUND_lepto_out=LUND${lepto_out}
cp ${not-rec_utils_dir}/leptoLUND.pl ${temp_dir}/
perl leptoLUND.pl ${z_vertex} < ${lepto_out}.txt > ${LUND_lepto_out}.dat

# Copy the gcard you'll use into the temp folder and set the torus value
cp ${not-rec_utils_dir}/${gcard_name}.gcard ${temp_dir}/
cp /group/clas12/gemc/4.4.2/experiments/clas12/micromegas/micromegas__bank.txt ${temp_dir}/
cp /group/clas12/gemc/4.4.2/experiments/clas12/dc/dc__bank.txt ${temp_dir}/
cp /group/clas12/gemc/4.4.2/experiments/clas12/bst/bst__bank.txt ${temp_dir}/

# Copy the cryotarget model into the execution folder and 
cp -r ${not-rec_utils_dir}/targets/${cryotarget_variation}/ ${temp_dir}/
cd ${not-rec_utils_dir}/targets/${cryotarget_variation}/
./targets.pl config.dat

# Change some variables in the gcard
cd ${temp_dir}

sed -i "s/TORUS_VALUE/${torus}/g" ${gcard_name}.gcard
sed -i "s/SOLENOID_VALUE/${solenoid}/g" ${gcard_name}.gcard
sed -i "s/CRYOTARGET_VARIATION/${cryotarget_variation}/g" ${gcard_name}.gcard
sed -i "s/TARGET_VARIATION/${target_variation}/g" ${gcard_name}.gcard

# EXECUTE GEMC
gemc ${gcard_name}.gcard -INTEGRATEDRAW="bst,dc,micromegas" -INPUT_GEN_FILE="LUND, ${LUND_lepto_out}.dat" -OUTPUT="evio, ${gemc_out}.ev" -USE_GUI="0" -FASTMCMODE="10"
echo "GEMC execution succesful!"

# Transform to HIPO
evio2hipo -t ${torus} -s ${solenoid} -r 11 -o ${gemc_out}.hipo -i ${gemc_out}.ev
#rm ${gemc_out}.ev
echo "Evio 2 HIPO transformation successful"
# Move output to its folder
mv ${gemc_out}.hipo ${out_dir_gemc}/
mv ${gemc_out}.ev ${out_dir_gemc}/
# Remove folder
rm -rf ${temp_dir}
echo "Done!"
