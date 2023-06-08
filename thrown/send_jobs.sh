#!/bin/bash

# Necessary directories
main_dir=$(pwd)
LEPTO_dir=${main_dir}/LEPTO
dat2tuple_dir=${main_dir}/dat2tuple

# Set necessary environmental variables to compile the codes
if [ -z "${CERN}" ]
then
    source ~/software/env_scripts/set_all.sh
fi

# Compile LEPTO from scratch
cd ${LEPTO_dir}
make clean
make

# Compile dat2tuple from scratch
cd ${dat2tuple_dir}
make clean
make

# Send jobs
cd ${main_dir}/run_lepto
sbatch JOB_run_lepto_fullchain.sh

