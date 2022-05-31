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

## DIRECTORIES
LEPTO_dir=~/clas12_simul/Lepto64Sim/bin
execution_dir=/volatile/clas12/emolinac
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

# Create folder in volatile to not interferewith other lepto executions
mkdir ${temp_dir}
cd ${temp_dir}

# Copy the executable to a different folder 
cp ${LEPTO_dir}/lepto.exe ${temp_dir}/lepto_${id}.exe

# Execution 
AZ_assignation ${target}
echo "${Nevents} ${A} ${Z}" > lepto_input.txt
lepto_${id}.exe < lepto_input.txt > lepto_out_${id}.txt

# Move output to its folder
mv lepto_out_${id}.txt ${out_dir}/

# Remove folder
rm -rf ${temp_dir}
echo "Done!"
