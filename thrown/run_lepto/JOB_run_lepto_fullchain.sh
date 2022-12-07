#!/bin/bash

#SBATCH --account=clas12
#SBATCH --partition=production
#SBATCH --job-name=lepto
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null
#SBATCH --time=00:05:00
#SBATCH --array=1-10000


# --output=./out/%x.%j.array%a.out
# --error=./err/%x.%j.array%a.err

echo "This is JOB ${SLURM_ARRAY_JOB_ID} task ${SLURM_ARRAY_TASK_ID}"

### FUNCTIONS
AZ_assignation(){
    # Function to assignate A and Z according to the target specified
    if [[ "$1" == "D" || "$1" == "D2" ]]
    then
	echo "Using Deuterium"
	A=2
	Z=1
	echo "A=${A} and Z=${Z}"
    elif [[ "$1" == "C" ]]
    then
	echo "Using Carbon"
	A=12
	Z=6
	echo "A=${A} and Z=${Z}"
    elif [[ "$1" == "Fe" ]]
    then
	echo "Using Iron"
	A=56
	Z=26
	echo "A=${A} and Z=${Z}"
    elif [[ "$1" == "Pb" ]]
    then
	echo "Using Lead"
	A=207
	Z=82
	echo "A=${A} and Z=${Z}"
    else
	echo "No target input!"
	exit 1
    fi    
}

directory_files_check(){
    # checking execution directories
    if [[ ! -d ${LEPTO_dir} || ! -d ${execution_dir} || ! -d ${out_dir} ]]
    then
	echo "One of the necessary directories does not exist."
	exit 1
    fi
    # checking lepto executable existence
    if [[ ! -f ${LEPTO_dir}/lepto.exe ]]
    then
	echo "The executable does not exist"
	exit 1
    fi
}

executable_file_check(){
    # checking dat2tuple executable existence
    if [[ ! -f ${dat2tuple_dir}/bin/dat2tuple ]]
    then
	echo "The dat2tuple executable does not exist"
	exit 1
    fi
}

## DIRECTORIES
main_dir=$(pwd)
LEPTO_dir=~/Lepto64Sim/bin ## CHECK THIS DIRECTORY!
execution_dir=/volatile/clas12/emolinac
lepto2dat_dir=${main_dir}/thrown/lepto2dat
dat2tuple_dir=${main_dir}/thrown/dat2tuple
rec_utils_dir=${main_dir}/reconstructed-double-target/utils

out_dir=/volatile/clas12/emolinac/only-lepto

## VARIABLES
Nevents=1000
target=D
id=${target}_${SLURM_ARRAY_JOB_ID}${SLURM_ARRAY_TASK_ID}
temp_dir=${execution_dir}/${id}
lepto_out=lepto_out_${id}
beam_energy=11
z_vertex=0

# Directory and files check
directory_files_check

# Set required variables for LEPTO execution
if [ -z "${CERN}" ]
then
    source ~/software/env_scripts/set_all.sh
fi

# Create folder in volatile to not interfere with other lepto executions
mkdir ${temp_dir}
cd ${temp_dir}

# Copy the executable to a different folder 
cp ${LEPTO_dir}/lepto.exe ${temp_dir}/lepto_${id}.exe

# Execution 
AZ_assignation ${target}
echo "${Nevents} ${A} ${Z}" > lepto_input.txt
lepto_${id}.exe < lepto_input.txt > ${lepto_out}.txt

cp ${lepto2dat_dir}/lepto2dat.pl ${temp_dir}/
perl lepto2dat.pl ${z_vertex} < ${lepto_out}.txt > ${lepto_out}.dat

executable_file_check
cp ${dat2tuple_dir}/bin/dat2tuple ${temp_dir}/
./dat2tuple ${lepto_out}.dat ${lepto_out}_ntuple.root

# Obtain LUND formated output
LUND_lepto_out=LUND${lepto_out}
cp ${rec_utils_dir}/leptoLUND.pl ${temp_dir}/
perl leptoLUND.pl ${z_vertex} ${beam_energy} < ${lepto_out}.txt > ${LUND_lepto_out}.dat

# Move output to its folder
#mv ${LUND_lepto_out}.dat ${lepto_out}_ntuple.root ${out_dir}/
#mv ${lepto_out}.dat ${lepto_out}_ntuple.root ${out_dir}/
mv ${lepto_out}_ntuple.root ${out_dir}/

# Remove folder
rm -rf ${temp_dir}
echo "Done!"
