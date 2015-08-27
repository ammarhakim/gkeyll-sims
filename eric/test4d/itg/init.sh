#!/bin/bash 
#PBS -N init-4d
#-- mail on execution("b"), termination ("e"), or interruption ("a")
#PBS -m ae
#PBS -M eshi@pppl.gov 
#PBS -l nodes=4:ppn=8
#PBS -l mem=96000mb
#PBS -l walltime=12:00:00
#PBS -r n
#PBS -V 
#PBS -j oe 
NPROCS=`wc -l < $PBS_NODEFILE`

CMD="/p/gke/eshi/gkeyllall/par-opt/gkeyll/gkeyll -i init.lua -pc_type lu -pc_factor_mat_solver_package superlu_dist"
cd $PBS_O_WORKDIR 
mpiexec -np $NPROCS $CMD 
exit
