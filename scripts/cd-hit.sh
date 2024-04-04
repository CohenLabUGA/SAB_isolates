#!/bin/bash
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --time=6:00:00
#SBATCH --mem=10G
#SBATCH --output=%x.out
#SBATCH --error=%x.err

module load CD-HIT

path="/work/nclab/lucy/SAB/Assembly/$1/rnaSpades"
pathout="/work/nclab/lucy/SAB/Assembly/$1/cd-hit"

#we collaps at 100% identity match
cd-hit-est -i $path/transcripts.fasta -o $pathout/$1transcript95.fasta -c 0.95 -n 9



