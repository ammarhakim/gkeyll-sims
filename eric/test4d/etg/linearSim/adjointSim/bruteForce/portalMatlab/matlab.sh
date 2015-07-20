#!/bin/bash 

#PBS -N mtest
#-- mail on execution("b"), termination ("e"), or interruption ("a")
#PBS -m ae
#PBS -M eshi@pppl.gov 
#PBS -l nodes=1:ppn=8 
#PBS -l mem=1000mb
#PBS -l walltime=1:00:00
#PBS -r n
#PBS -V 
#PBS -j oe 
cd $PBS_O_WORKDIR 
matlab -nodisplay -nosplash < portalProcessMatrices.m > run.loog
echo ""
echo "Done at "`date`
