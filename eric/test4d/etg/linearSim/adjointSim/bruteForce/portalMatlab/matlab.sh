#!/bin/bash
#
#-- mail on execution("b"), termination ("e"), or interruption ("a")
#PBS -m ae
#PBS -M eshi@pppl.gov 
#PBS -l nodes=1:ppn=8 
#PBS -l mem=16000mb 
#PBS -l walltime=24:00:00
#PBS -r n
#PBS -V 
#PBS -j oe 

cd $PBS_O_WORKDIR
matlab -nosplash -nodisplay <<EOF
matlabpool close force local
matlabpool open local 8   
% call the function
portalProcessMatrices
matlabpool close
exit  
EOF
   
echo ""
echo "Done at " `date`
