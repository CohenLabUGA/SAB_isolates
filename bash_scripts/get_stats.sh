#!/bin/bash
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --time=6:00:00
#SBATCH --mem=10G
#SBATCH --job-name=getstat
#SBATCH --output=getstat.out
#SBATCH --error=getstat.err

## Using assembly stats package to get the number and % of transcriptome which was translated into predicted proteins by transdecoder. 
## assembly-stats downloaded via conda

outdir='/work/nclab/lucy/SAB/Annotation/transdecoder/stats/'

samples="04 06 08 13"

for s in $samples; do

echo ${s}

org_assembly=/work/nclab/lucy/SAB/Assembly/${s}/rnaSpades/post_ribo/transcripts.fasta
transdecoder=/work/nclab/lucy/SAB/Annotation/transdecoder/${s}/transcripts.fasta.transdecoder.pep

assembly-stats -t $org_assembly $transdecoder > $outdir/${s}stats

done;
