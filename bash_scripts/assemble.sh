#!/bin/bash
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=36:00:00
#SBATCH --mem=100G
#SBATCH --job-name=rnaSpades
#SBATCH --output=%x.out
#SBATCH --error=%x.err

#using conda environment spades 

samples="04 06 08 13"

for s in `echo $samples`; do
	echo ${s}
 
	rnaspades.py --dataset input${r}${s}.yaml -o /path/to/Assembly/${s}/rnaSpades/
done;

