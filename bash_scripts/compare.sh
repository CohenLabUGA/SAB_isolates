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
#
#samples="04 08 06 13"
#
#for s in $samples; do
#
#	indir=/work/nclab/lucy/SAB/Assembly/${s}/sourmash/sig
#	outdir=/work/nclab/lucy/SAB/Assembly/${s}/sourmash
#
#	sourmash compare $indir/${s}{A,B,C}.sig -o $outdir/${s}compare -k 31 --csv $outdir/${s}compare.csv
#	
#	sourmash plot --pdf --labels $outdir/${s}compare
#	echo ${s}
#done;

# ---------------------------------------------------------------------
## comparing replicates between taxa groups

#sourmash compare /work/nclab/lucy/SAB/Assembly/{04,08}/sourmash/sig/*{A,B,C}.sig -o /work/nclab/lucy/SAB/Assembly/sourmash/diatoms_compare -k 31 --csv /work/nclab/lucy/SAB/Assembly/sourmash/diatoms_compare.csv

#sourmash plot --pdf --labels /work/nclab/lucy/SAB/Assembly/sourmash/diatoms_compare

## coccolithophores:

#sourmash compare /work/nclab/lucy/SAB/Assembly/{06,13}/sourmash/sig/*{A,B,C}.sig -o /work/nclab/lucy/SAB/Assembly/sourmash/coccolithophores_compare -k 31 --csv /work/nclab/lucy/SAB/Assembly/sourmash/coccolithophores_compare.csv

#sourmash plot --pdf --labels /work/nclab/lucy/SAB/Assembly/sourmash/coccolithophores_compare

# --------------------------------------------------------------------
## all taxa together:

#sourmash compare /work/nclab/lucy/SAB/Assembly/{04,08,06,13}/sourmash/sig/*{A,B,C}.sig -o /work/nclab/lucy/SAB/Assembly/sourmash/all_taxa_compare -k 31 --csv /work/nclab/lucy/SAB/Assembly/sourmash/all_taxa_compare.csv

#sourmash plot --pdf --labels /work/nclab/lucy/SAB/Assembly/sourmash/all_taxa_compare

# -------------------------------------------------------------------
## using transcriptome constructed for each individual sample comparing across taxa groups and all organisms

## separate use case ---------------------------
#diatoms=`ls -1 $indir/* | grep -oP '(04|08).*'`
#cocco=`ls -1 $indir/* | grep -oP '(06|^13).*'`

indir="/work/nclab/lucy/SAB/Assembly/sourmash/full_transcriptome/sig"
outdir="/work/nclab/lucy/SAB/Assembly/sourmash/full_transcriptome/compare"

#comparing the diatoms:

#sourmash compare $indir/{04,08}*.sig -o $outdir/diatoms_by_sample_compare -k 31 --csv $outdir/diatoms_by_sample_compare.csv

#comparing the coccolithophores:

#sourmash compare $indir/{06,13}*.sig -o $outdir/coccoliths_by_sample_compare -k 31 --csv $outdir/coccoliths_by_sample_compare.csv

#comparing all samples:

#sourmash compare $indir/*.sig -o $outdir/all_samples_compare -k 31 --csv $outdir/all_samples_compare.csv

#sourmash plot --pdf --labels $outdir/all_transcripts_similarity.mat
### compare complete transcriptomes

sourmash compare $indir/*.sig -o $outdir/transcriptome_compare -k 31 --csv $outdir/transcriptome_compare.csv
