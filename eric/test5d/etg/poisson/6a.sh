#!/bin/bash 
#PBS -N par-5d-6a
#-- mail on execution("b"), termination ("e"), or interruption ("a")
#PBS -m ae
#PBS -M eshi@pppl.gov 
#PBS -l nodes=2:ppn=32
#PBS -l mem=220000mb
#PBS -l walltime=48:00:00
#PBS -r n
#PBS -V 
#PBS -j oe 
#PBS -q dtest
NPROCS=`wc -l < $PBS_NODEFILE`

CMD="/p/gke/eshi/gkeyllall/par-opt/gkeyll/gkeyll -i 6a.lua -pc_type lu -pc_factor_mat_solver_package superlu_dist" 
cd $PBS_O_WORKDIR 
mpiexec -np $NPROCS $CMD 
exit
