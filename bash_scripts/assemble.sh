#!/bin/bash
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=36:00:00
#SBATCH --mem=100G
#SBATCH --job-name=rnaSpades
#SBATCH --output=%x.out
#SBATCH --error=%x.err

# activate conda environment spades 
# this script loops through each organism number, and uses the respective .yml file
# containing sequence file names, locations, and sequence types (Illumina vs PacBio or forward vs reverse read)

samples="04 06 08 13"

for s in `echo $samples`; do
	echo ${s}
 
	rnaspades.py --dataset input${s}.yaml -o /path/to/Assembly/${s}/
done;

