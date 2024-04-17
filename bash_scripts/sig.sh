#!/bin/bash
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --time=36:00:00
#SBATCH --mem=20G
#SBATCH --job-name=sig.sourmash
#SBATCH --output=%x.out
#SBATCH --error=%x.err
#---------------------------------------
#Lucy Quirk UGA 2023-04-19

#activate environment called sourmash
#make a signature for each assembly file
#use annother script to compare all of these signatures to one annother and see similarity

#----------------------------------------

#samples="04 06 08 13"
#reps="A B C"

## location of the database which I will use to query against my samples
#database=/work/nclab/lucy/SAB/src/sourmash_euk_rna_db     

#for s in $samples; do

## define directory from which you are pulling sample and where the output will be put
## these locations are defined at each iteration of the loop to follow the folder 
## structure based on my 4 organisms

#	echo ${s}

#	for r in $reps; do
#		echo ${r}
#		
#		indir=/work/nclab/lucy/SAB/Assembly/${s}/rnaSpades/replicates/${r}
#
#		outdir=/work/nclab/lucy/SAB/Assembly/${s}/sourmash
#
#		sourmash sketch dna -p k=31,scaled=10000 $indir/transcripts.fasta -o $outdir/sig/${s}${r}.sig
#
## match the sig of the transcriptome to a database of eukaryotic microbes downloaded from sourmash website
#
#		sourmash gather -k 31 --scaled 10000 -o $outdir/query/${s}${r}query.csv $outdir/sig/${s}${r}.sig $database/*
#
#	done;
#done;


# -------------------------------------------------------
## creating a sig for transcriptome of each individual sample
samples='04 08 06 13'
for i in $samples; do
echo${i}
indir='/work/nclab/lucy/SAB/Assembly/'${i}'/rnaSpades/post_ribo'
outdir="/work/nclab/lucy/SAB/Assembly/sourmash/full_transcriptome/sig"

sourmash sketch dna -p k=31,scaled=10000 $indir/transcripts.fasta -o $outdir/${i}.sig

done;
