#!/bin/bash

#SBATCH -J nstx-no-colls
#SBATCH -A TG-PHY160007       # Project
#SBATCH -p skx-dev     # Queue (partition) name -- normal, development, etc.
#SBATCH -N 4
#SBATCH -n 192 
#SBATCH -t 2:00:00
#SBATCH --mail-user=tnbernard@utexas.edu
#SBATCH --mail-type=ALL

mpirun ~/gkylsoft/gkyl/bin/gkyl nstx-Lz12-no-collisions.lua
