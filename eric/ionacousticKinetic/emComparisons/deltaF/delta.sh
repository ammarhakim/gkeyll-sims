#!/bin/bash 

#PBS -N gke-kinetic-test
#-- mail on execution("b"), termination ("e"), or interruption ("a")
#PBS -m ab 
#PBS -M eshi@pppl.gov 
#PBS -l nodes=1:ppn=1 
#PBS -l mem=4gb 
#PBS -l walltime=16:00:00
#PBS -r n
#PBS -V 
#PBS -j oe 
 
CMD="/p/gke/eshi/gkeyllall/ser-opt/gkeyll/gkeyllser -i delta.lua -pc_type lu" 
cd $PBS_O_WORKDIR 
$CMD 
exit