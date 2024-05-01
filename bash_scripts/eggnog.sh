#!/bin/bash
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --time=36:00:00
#SBATCH --mem=100G
#SBATCH --job-name=eggnog
#SBATCH --output=%x.out
#SBATCH --error=%x.err

## using conda environment eggnog
#map the 'pep' file from transdecoder which is in 'protein space' 

samples="04 06 08 13"

for s in $samples; do

	echo ${s}
	indir=/work/nclab/lucy/SAB/Annotation/transdecoder/${s}
	outdir=/work/nclab/lucy/SAB/Annotation/eggnog/gff/


	emapper.py -i $indir/transcripts.fasta.transdecoder.pep -o ${s} --output_dir $outdir

done;
