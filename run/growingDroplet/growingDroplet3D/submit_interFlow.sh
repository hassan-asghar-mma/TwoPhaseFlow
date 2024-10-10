#!/bin/bash -l

#allocate 1 node for 24 hour on LB2 partition
#SBATCH -A special00005 
#SBATCH -n 4 
#SBATCH --mem-per-cpu=3000 
#SBATCH -t 10:30:00     
#SBATCH --error=error.err
#SBATCH -J gD2D 

#interFlow > growingDroplet3D.o
srun -n 4 interFlow -parallel  > growingDroplet3D.o
./combineCaseData.sh
#reconstructPar -fileHandler collated
