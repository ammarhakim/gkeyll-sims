#!/bin/bash 

#PBS -N sol-rec-0
#-- mail on execution("b"), termination ("e"), or interruption ("a")
#PBS -m ae
#PBS -M eshi@pppl.gov 
#PBS -l nodes=1:ppn=1 
#PBS -l mem=2000mb 
#PBS -l walltime=96:00:00
#PBS -r n
#PBS -V 
#PBS -j oe 
 
CMD="/p/gke/eshi/gkeyllall/ser-opt/gkeyll/gkeyllser -i 0.lua -pc_type lu -r 40" 
cd $PBS_O_WORKDIR 
$CMD 
exit
