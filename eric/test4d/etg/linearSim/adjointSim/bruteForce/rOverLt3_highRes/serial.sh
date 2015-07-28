#!/bin/bash
#
#PBS -l walltime=24:00:00,nodes=1:ppn=1,mem=16000mb
#PBS -j oe
#
#PBS -m ae
#PBS -M eshi@pppl.gov
#

cd $PBS_O_WORKDIR
matlab -nosplash -nodisplay <<EOF
% call the function
processMatrices_3
EOF
   
echo ""
echo "Done at " `date`
