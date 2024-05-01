#!/bin/bash
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --time=36:00:00
#SBATCH --mem=100G
#SBATCH --job-name=rnaQuast
#SBATCH --output=%x.out
#SBATCH --error=%x.err

###  activate env rnaQuast for rnaQUAST

samples="04 06 08 13"

for s in `echo $samples`; do

indir=/path/to/${s}/assembly
outdir=/path/to/${s}/quast

rnaQUAST.py --transcripts $indir/transcripts.fasta --output_dir $outdir/rnaQUAST_${s} --busco /path/to/.conda/env/rnaQuast/bin/busco_downloads/eukaryota --gene_mark -d

done;
