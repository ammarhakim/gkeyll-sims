#!/bin/bash 
#PBS -N adj-2
#-- mail on execution("b"), termination ("e"), or interruption ("a")
#PBS -m ae
#PBS -M eshi@pppl.gov 
#PBS -l nodes=1:ppn=1 
#PBS -l mem=8000mb 
#PBS -l walltime=12:00:00
#PBS -r n
#PBS -V 
#PBS -j oe 
 
CMD="/p/gke/eshi/gkeyllall/ser-opt/gkeyll/gkeyllser -i adjoint2.lua -pc_type lu" 
cd $PBS_O_WORKDIR 
$CMD 
exit
