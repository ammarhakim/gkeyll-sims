#!/bin/bash 

#PBS -N mtest
#-- mail on execution("b"), termination ("e"), or interruption ("a")
#PBS -m ae
#PBS -M eshi@pppl.gov 
#PBS -l nodes=2:ppn=4 
#PBS -l mem=1000mb
#PBS -l walltime=1:00:00
#PBS -r n
#PBS -V 
#PBS -j oe 
CMD="matlab -nodisplay -nosplash < hello.m" 
cd $PBS_O_WORKDIR 
$CMD 
exit
