#!/bin/bash
#SBATCH --partition=batch
#SBATCH --time=24:00:00
#SBATCH --job-name=transdecoder
#SBATCH --ntasks=1
#SBATCH --mem=30G
#SBATCH --output=transdecoder.out
#SBATCH --error=transdecoder.err
#SBATCH --verbose

## activate env transdecoder
#translate transcriptomes before running through eggnogmapper

samples="04 08 06 13"

for s in $samples; do
	indir=/path/to/${s}/assembly
	outdir=/path/to/transdecoder/output
	
	echo "extracting long open reading frames from ${s}"
	TransDecoder.LongOrfs -t $indir/transcripts.fasta --gene_trans_map --output_dir $outdir

	echo "predicitng proteins in ${s}"
	TransDecoder.Predict -t $indir/transcripts.fasta --output_dir $outdir
done;
