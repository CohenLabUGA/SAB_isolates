#!/bin/bash
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --time=36:00:00
#SBATCH --mem=100G
#SBATCH --job-name=eggnog
#SBATCH --output=%x.out
#SBATCH --error=%x.err

## using conda environmnet eggnog
#map the 'pep' file from transdecoder which is in 'protein space' 

samples="04 06 08 13"

for s in $samples; do

	echo ${s}
	indir=/work/nclab/lucy/SAB/Annotation/transdecoder/${s}
	outdir=/work/nclab/lucy/SAB/Annotation/eggnog/gff/

## original run:
#	emapper.py -i $indir/transcripts.fasta.transdecoder.pep -o ${s} --output_dir $outdir

## run to get a gff file of annotations:
	emapper.py -i $indir/transcripts.fasta.transdecoder.pep -o ${s} --output_dir $outdir --decorate_gff yes --decorate_gff_ID_field GeneID --override

## run which translates as it annotates:
#	emapper.py -m diamond --itype CDS --translate -i $indir/${s}transcripts.fast -o ${s} --output_dir $outdir
#	emapper.py -m hmmer -d Eukaryotic -i $indir/${s}transcripts.fasta


done;
