#!/bin/bash
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --time=6:00:00
#SBATCH --mem=10G
#SBATCH --job-name=fastqc
#SBATCH --output=%x.out
#SBATCH --error=%x.err

indir=/work/nclab/lucy/SAB/Assembly/ribodetector/nonrna
outdir=/work/nclab/lucy/SAB/Assembly/ribodetector/fastqc


#load module fast qc

indir=/path/to/raw/reads
outdir=/path/to/output/fastqc


fastqc $indir/*.fq -o $outdir
