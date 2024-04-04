#!/bin/bash
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --time=36:00:00
#SBATCH --mem=20G
#SBATCH --job-name=compare.plot
#SBATCH --output=compare.out
#SBATCH --error=compare.err

## using conda env sourmash:
## match the sig of the transcriptome to a database of eukaryotic microbes. 

# -------------------------------------------------------------------------
## comparing replicates within organisms ## 

samples="04 08 06 13"

for s in $samples; do

	indir=/work/nclab/lucy/SAB/Assembly/${s}/sourmash/sig
	outdir=/work/nclab/lucy/SAB/Assembly/${s}/sourmash

	sourmash compare $indir/${s}{A,B,C}.sig -o $outdir/${s}compare -k 31 --csv $outdir/${s}compare.csv
	
	sourmash plot --pdf --labels $outdir/${s}compare 
	echo ${s}
done;

# ---------------------------------------------------------------------
## comparing replicates between taxa groups

taxa="diatoms coccolithophores"
diatoms="04 08"
coccolithophores="06 13"

for t in $taxa; do

echo ${t}[1]

done;






## separate use case ---------------------------
#diatoms=`ls -1 $indir/* | grep -oP '(04|08).*'`
#cocco=`ls -1 $indir/* | grep -oP '(06|^13).*'`

#echo "diatoms files are:
#$diatoms"
#echo "coccoliths are:
#$cocco"

#sourmash compare $indir/$diatoms -o $outdir/diatoms_similarity -k 31
#sourmash compare $indir/$coccoliths -o $outdir/coccoliths_similarity -k 31

#sourmash plot --pdf --labels $outdir/all_transcripts_similarity.mat
