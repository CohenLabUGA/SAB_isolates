#!/bin/bash
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --time=6:00:00
#SBATCH --mem=10G
#SBATCH --output=%x.out
#SBATCH --error=%x.err

module load CD-HIT

samples='04 06 08 13"
for s in $samples; do

  path="/work/nclab/lucy/SAB/Assembly/${s}/rnaSpades"
  pathout="/work/nclab/lucy/SAB/Assembly/${s}/cd-hit"

  echo "collaping transcripts in ${s} with over 94% identity"
#we collaps at 100% identity match
  cd-hit-est -i $path/transcripts.fasta -o $pathout/${s}transcript95.fasta -c 0.95 -n 9
  
done;



