#!/bin/bash
#SBATCH --partition=batch
#SBATCH --time=24:00:00
#SBATCH --job-name=transdecoder
#SBATCH --ntasks=1
#SBATCH --mem=30G
#SBATCH --output=transdecoder.out
#SBATCH --error=transdecoder.err
#SBATCH --verbose


#translate transcriptomes before running through eggnogmapper

samples="04 08 06 13"

for s in $samples; do
	indir=/work/nclab/lucy/SAB/Assembly/${s}/rnaSpades/post_ribo
	outdir=/work/nclab/lucy/SAB/Annotation/transdecoder/${s}/

pushd ${outdir}

	echo "extracting long open reading frames from ${s}"
	TransDecoder.LongOrfs -t $indir/transcripts.fasta

	echo "predicitng proteins in ${s}"
	TransDecoder.Predict -t $indir/transcripts.fasta

	popd

done;
#viewing ORF predictions?..
#java -jar $GENOMEVIEW/genomeview.jar transcripts.fasta transcripts.fasta.transdecoder.bed
