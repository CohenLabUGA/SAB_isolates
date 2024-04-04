#!/bin/bash
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --job-name=test
#SBATCH --output=%x.out
#SBATCH --error=%x.err
#SBATCH --mem=1G 
#SBATCH --time=00:10:00

indir="/work/nclab/lucy/SAB/Assembly/ribodetector/nonrna/"

files=`ls $indir | grep -oP $indir.*A.*R1'

echo $files

names= fopen("names_match.txt","w")
fputs($files,names)
fclose(names)
