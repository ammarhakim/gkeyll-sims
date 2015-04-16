#!/bin/bash 
#PBS -N par-5d
#-- mail on execution("b"), termination ("e"), or interruption ("a")
#PBS -m ae
#PBS -M eshi@pppl.gov 
#PBS -l nodes=8:ppn=8
#PBS -l mem=16gb 
#PBS -l walltime=96:00:00
#PBS -r n
#PBS -V 
#PBS -j oe 
 
NPROCS=`wc -l < $PBS_NODEFILE`

CMD="/p/gke/eshi/gkeyllall/par-opt/gkeyll/gkeyll -i test0.lua -pc_type lu" 
cd $PBS_O_WORKDIR 
mpiexec -np $NPROCS $CMD 
exit
