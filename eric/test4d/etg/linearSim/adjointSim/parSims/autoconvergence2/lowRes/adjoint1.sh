#!/bin/bash 
#PBS -N auto-adj-1
#-- mail on execution("b"), termination ("e"), or interruption ("a")
#PBS -m ae
#PBS -M eshi@pppl.gov 
#PBS -l nodes=1:ppn=2
#PBS -l mem=48000mb
#PBS -l walltime=24:00:00
#PBS -r n
#PBS -V 
#PBS -j oe 
NPROCS=`wc -l < $PBS_NODEFILE`

CMD="/p/gke/eshi/gkeyllall/par-opt/gkeyll/gkeyll -i adjoint1.lua -pc_type lu" 
cd $PBS_O_WORKDIR 
mpiexec -np $NPROCS $CMD 
exit
