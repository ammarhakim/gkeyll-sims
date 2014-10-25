#!/bin/bash 

#PBS -N rs4-24
#-- mail on execution("b"), termination ("e"), or interruption ("a")
#PBS -m ae
#PBS -M eshi@pppl.gov 
#PBS -l nodes=1:ppn=1 
#PBS -l mem=4gb 
#PBS -l walltime=32:00:00
#PBS -r n
#PBS -V 
#PBS -j oe 
 
CMD="/p/gke/eshi/gkeyllall/ser-opt/gkeyll/gkeyllser -i 24.lua -pc_type lu" 
cd $PBS_O_WORKDIR 
$CMD 
exit
