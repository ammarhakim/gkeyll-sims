#!/bin/bash
#
#PBS -l walltime=48:00:00,nodes=1:ppn=1,mem=16000mb
#PBS -j oe
#PBS -N matlab-svd
#
#PBS -m ae
#PBS -M eshi@pppl.gov
#

cd $PBS_O_WORKDIR
matlab -nosplash -nodisplay <<EOF
% call the function
processSVD
EOF
   
echo ""
echo "Done at " `date`
