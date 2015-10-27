#!/bin/bash 
#PBS -N itg-4d-4a
#-- mail on execution("b"), termination ("e"), or interruption ("a")
#PBS -m ae
#PBS -M eshi@pppl.gov 
#PBS -l nodes=1:ppn=32
#PBS -l mem=62000mb
#PBS -l walltime=24:00:00
#PBS -r n
#PBS -V 
#PBS -j oe 
#PBS -q dtest
NPROCS=`wc -l < $PBS_NODEFILE`

CMD="/p/gke/eshi/gkeyllall/par-opt/gkeyll/gkeyll -i 4a.lua -pc_type lu -pc_factor_mat_solver_package superlu_dist -r 15"
cd $PBS_O_WORKDIR 
mpiexec -np $NPROCS $CMD 
exit
