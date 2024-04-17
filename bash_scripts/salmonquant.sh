#!/bin/bash
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=36:00:00
#SBATCH --mem=100G
#SBATCH --job-name=quant
#SBATCH --output=%x.out
#SBATCH --error=%x.err

transcripts=/work/nclab/lucy/SAB/Assembly/ribodetector/nonrna

source trimmed_names.sh

for s in `echo $samples`; do
        echo ${s}
        transcriptome=/work/nclab/lucy/SAB/Assembly/$1/rnaSpades/post_ribo/transcripts.fasta
        outdir=/work/nclab/lucy/SAB/Assembly/$1/salmon

	R1=`ls -1 $transcripts | grep -o ${s}.*R1.fq`
	R2=`ls -1 $transcripts | grep -o ${s}.*R2.fq`
	base=${s}

	echo ${R1}
	echo ${R2}

salmon quant -i $outdir/$1_index -l A -1 $transcripts/${R1} -2 $transcripts/${R2} -p 16 --validate Mappings -o $outdir/${base}_quant

done;

