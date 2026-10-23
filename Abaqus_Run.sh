#!/bin/bash
#SBATCH --partition=general
#SBATCH -N 1
#SBATCH --ntasks-per-node=126
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --output=Abaqus_Run.out
#SBATCH --open-mode=truncate
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=david.o.cohen@uconn.edu

module purge
module load abaqus/2022
module load intel/oneapi/2022.3
WORK_DIR=/home/tan22006/Tony_Nguyen/ABAQUS/2_8_2024
HOME_DIR=/home/tan22006/Tony_Nguyen/ABAQUS/2_8_2024

# mkdir $WORK_DIR
# cp Heat_Analysis.inp $WORK_DIR
# cp welding-CuNi-heat-constant.for $WORK_DIR
# cp Stress_Analysis.inp $WORK_DIR
# cp welding-CuNi-stress-constant.for $WORK_DIR

cd $WORK_DIR
#abq2022 job=Heat_Analysis user=welding-CuNi-heat-constant.for cpus=$SLURM_NTASKS interactive
abq2022 job=Stress_Analysis user=welding-CuNi-stress-constant.for cpus=$SLURM_NTASKS interactive
