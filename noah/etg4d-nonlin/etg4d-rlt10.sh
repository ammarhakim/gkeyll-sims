#!/bin/bash -l
#SBATCH -q regular
#SBATCH -N 2
#SBATCH -t 01:00:00
#SBATCH -J my_job
#SBATCH -o my_job.o%j
#SBATCH -L SCRATCH,project

#Edison has 24 cores per compute node
srun -n 48 ~/gkylsoft/gkyl/bin/gkyl etg4d-rlt10.lua
