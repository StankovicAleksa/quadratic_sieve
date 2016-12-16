#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=4096
#SBATCH --time=0:10:00
#SBATCH --job-name=nbody4

echo "$N ${SLURM_JOB_ID}" >> jobs.txt

#/usr/bin/time -f "Time=%e"

module purge
module load intelmpi/5.1.3  

srun ./NBODY.x $1 $2
srun ./NBODY.x $1 $2
srun ./NBODY.x $1 $2
srun ./NBODY.x $1 $2
srun ./NBODY.x $1 $2
srun ./NBODY.x $1 $2
srun ./NBODY.x $1 $2
srun ./NBODY.x $1 $2
srun ./NBODY.x $1 $2
srun ./NBODY.x $1 $2

# according to man time:
# -f specifise the format





