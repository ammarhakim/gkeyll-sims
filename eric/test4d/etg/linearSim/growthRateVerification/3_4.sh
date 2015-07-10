#!/bin/bash 
#PBS -N par-3-4
#-- mail on execution("b"), termination ("e"), or interruption ("a")
#PBS -m ae
#PBS -M eshi@pppl.gov 
#PBS -l nodes=1:ppn=8
#PBS -l mem=48000mb
#PBS -l walltime=12:00:00
#PBS -r n
#PBS -V 
#PBS -j oe 
NPROCS=`wc -l < $PBS_NODEFILE`

CMD="/p/gke/eshi/gkeyllall/par-opt/gkeyll/gkeyll -i 3_4.lua -pc_type lu" 
cd $PBS_O_WORKDIR 
mpiexec -np $NPROCS $CMD 
exit
