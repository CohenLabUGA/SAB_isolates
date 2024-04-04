#!/bin/bash
#SBATCH --partition=batch
#SBATCH --time=24:00:00
#SBATCH --job-name=transdecoder
#SBATCH --ntasks=1
#SBATCH --mem=30G
#SBATCH --output=cox1Translate.out
#SBATCH --error=cox1Translate.err
#SBATCH --verbose


#translate nucleotide cox1 matches, output in cox1 folder under Assembly/

samples="04 08 06 13"

for s in $samples; do
        indir=/work/nclab/lucy/SAB/Assembly/${s}/cox1

pushd ${indir}

        echo "extracting long open reading frames from ${s} cox1 gene"
        TransDecoder.LongOrfs -t ${s}cox1-match.fasta

        echo "predicitng proteins in ${s} cox1 gene"
        TransDecoder.Predict -t ${s}cox1-match.fasta

        popd

done;

