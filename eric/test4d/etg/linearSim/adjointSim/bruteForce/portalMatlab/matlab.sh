#!/bin/bash
#
#PBS -l walltime=8:00:00,nodes=1:ppn=8,mem=16gb
#PBS -j oe
#
#PBS -m ae
#PBS -M eshi@pppl.gov
#

cd $PBS_O_WORKDIR
matlab -nosplash -nodisplay <<EOF
matlabpool open local 8   
% call the function
portalProcessMatrices
matlabpool close
   
EOF
   
echo ""
echo "Done at " `date`
