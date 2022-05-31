#!/bin/bash

#SBATCH --account=clas12
#SBATCH --partition=production
#SBATCH --job-name=lepto
#SBATCH --output=./out/%x.%j.array%a.out
#SBATCH --error=./err/%x.%j.array%a.err
#SBATCH --time=00:05:00
#SBATCH --array=1-200

echo "This is JOB ${SLURM_ARRAY_JOB_ID} task ${SLURM_ARRAY_TASK_ID}"

### FUNCTIONS
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
    elif [[ "$1" == "Fe" ]]
    then
	echo "Using Iron"
	A=56
	Z=26
    elif [[ "$1" == "Pb" ]]
    then
	echo "Using Lead"
	A=207
	Z=82
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
LEPTO_dir=~/clas12_simul/Lepto64Sim/bin ## CHECK THIS DIRECTORY!
execution_dir=/volatile/clas12/emolinac
lepto2dat_dir=~/clas12_simul/thrown/lepto2dat
dat2tuple_dir=~/clas12_simul/thrown/dat2tuple

out_dir=/work/clas12/rg-e/emolinac/lepto

## VARIABLES
Nevents=1000
target=D
id=${target}_${SLURM_ARRAY_JOB_ID}${SLURM_ARRAY_TASK_ID}
temp_dir=${execution_dir}/${id}

# Directory and files check
directory_files_check

# Set required variables for LEPTO execution
source ~/software/env_scripts/set_all.sh

# Create folder in volatile to not interfere with other lepto executions
mkdir ${temp_dir}
cd ${temp_dir}

# Copy the executable to a different folder 
cp ${LEPTO_dir}/lepto.exe ${temp_dir}/lepto_${id}.exe

# Execution 
lepto_out=lepto_out_${id}
z_shift=0.

AZ_assignation ${target}
echo "${Nevents} ${A} ${Z}" > lepto_input.txt
lepto_${id}.exe < lepto_input.txt > ${lepto_out}.txt

cp ${lepto2dat_dir}/lepto2dat.pl ${temp_dir}/
perl lepto2dat.pl ${z_shift} < ${lepto_out}.txt > ${lepto_out}.dat

executable_file_check
cp ${dat2tuple_dir}/bin/dat2tuple ${temp_dir}/
./dat2tuple ${lepto_out}.dat ${lepto_out}_ntuple.root

# Move output to its folder
mv ${lepto_out}.dat ${lepto_out}_ntuple.root ${out_dir}/

# Remove folder
rm -rf ${temp_dir}
echo "Done!"
