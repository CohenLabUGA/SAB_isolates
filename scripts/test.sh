#!/bin/bash
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=36:00:00
#SBATCH --mem=10G
#SBATCH --job-name=test
#SBATCH --output=%x.out
#SBATCH --error=%x.err

indir="/work/nclab/lucy/SAB/Annotation/transdecoder/04/transcripts.fasta.transdecoder_dir"
    echo ${indir}
    sed 's/\_length.*i0//g' $indir/shorter.pep.fasta >> /work/nclab/lucy/SAB/shortest.pep.fasta
