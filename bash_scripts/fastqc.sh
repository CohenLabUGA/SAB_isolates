#!/bin/bash
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --time=6:00:00
#SBATCH --mem=10G
#SBATCH --job-name=fastqc
#SBATCH --output=%x.out
#SBATCH --error=%x.err

#load modules FastQC and MultiQC before running

indir=/path/to/raw/reads
outdir=/path/to/output/fastqc


fastqc $indir/*.fq -o $outdir

multiqc $outdir/ 
