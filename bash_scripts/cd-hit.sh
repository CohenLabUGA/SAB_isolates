#!/bin/bash
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --time=6:00:00
#SBATCH --mem=10G
#SBATCH --output=%x.out
#SBATCH --error=%x.err

## load the CD-HIT module
module load CD-HIT

samples='04 06 08 13'

for s in $samples; do

  pathin="/path/to/Assembly/${s}/rnaSpades"
  pathout="/path/to/Assembly/${s}/cd-hit"

  echo "collaping transcripts in ${s} with over 94% identity"
#we collaps at 95% identity match
  cd-hit-est -i $pathin/transcripts.fasta -o $pathout/${s}transcript95.fasta -c 0.95 -n 9
  
done;



