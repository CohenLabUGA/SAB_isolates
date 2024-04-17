#!/bin/bash
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --mem=10G
#SBATCH --job-name=trimmomatic
#SBATCH --output=%x.out
#SBATCH --error=%x.err

#load module trimmomatic
module load Trimmomatic

#name indir and outdir
indir="/work/nclab/data/SAB_isolates/Illumina/raw"
outdir="/work/nclab/data/SAB_isolates/Illumina/trimmed"

#need to add loop here for naming reverse and forward read seqs


trimmomatic PE <forward read sequence> <reverse read> <output forward read> <output reverse read> $adapters 

# $adapters is a file that lists the adapters which we will trim from our reads

