#!/bin/bash
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --time=36:00:00
#SBATCH --mem=100G
#SBATCH --job-name=rnaQuast
#SBATCH --output=%x.out
#SBATCH --error=%x.err

###  activate env-quast1 for rnaQUAST

samples="04 06 08 13"

for s in `echo $samples`; do

indir=/work/nclab/lucy/SAB/Assembly/${s}/rnaSpades/post_ribo
outdir=/work/nclab/lucy/SAB/Assembly/${s}/quast/post_ribo

rnaQUAST.py --transcripts $indir/transcripts.fasta --output_dir $outdir/rnaQUAST_${s} --busco /home/leq40065/.conda/env/env-quast1/bin/busco_downloads/eukaryota --gene_mark -d
#gmst.pl Euka_transcripts.fasta

done;
