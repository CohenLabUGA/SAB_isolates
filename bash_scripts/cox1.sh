#!/bin/bash
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --time=6:00:00
#SBATCH --mem=10G
#SBATCH --job-name=createDB
#SBATCH --output=cox1blast.out
#SBATCH --error=cox1blast.err

# activate ncbi_database conda env

#samples is a space delinated list of organism names used as the folder naming scheme
samples="04 06 08 13"

for s in `echo $samples`; do
    echo ${s}
## create indirectory and outdirectory for each iteration (s) and name them as variables to be called throughout the script 
    indir="/work/nclab/lucy/SAB/Annotation/transdecoder/${s}/transcripts.fasta.transdecoder_dir"
    outdir="/work/nclab/lucy/SAB/Assembly/${s}/cox1/prot"
    echo ${indir}
    echo ${outdir}

## if you have not already, make a folder for the database to go
#mkdir $outdir/${s}db
###

### 1. location and name of protein file from ncbi database
if [${s} = '04'] || [${s} = '08']; then
	cox_gene="/path/to/cylindrotheca_cox1_gene.faa"
 else
	cox_gene="/path/to/coccolith_cox1_gene/${s}_cox1.fasta"
 fi

### 2. make a blast database from the transcriptome of each sample (s)
# awk -F " " '{if ($0 ~ /^>/) {print $1} else {print $0}}' $indir/longest_orfs.pep > $indir/shorter.pep.fasta
	 sed 's/\_length.*g[[:digit:]]*//g' $indir/shorter.pep.fasta > $outdir/shortest.pep.fasta

echo head $outdir/shortest.pep.fasta

        makeblastdb -in $outdir/shortest.pep.fasta -out $outdir/${s}db/ -dbtype prot

### 3. query the database (transcriptome) with each isip gene uploaded
## we will use tblastp because the query is a protein and the db is a
## translated nucleotide

	echo 'blasting transcriptome for cox1'
       blastp -query $cox_gene -db ${outdir}/${s}db -out ${outdir}/${s}blastp.txt -evalue 1e-30
        #this finds unique orfs listed in match and writes the orf id's into a file call hits.txt
        cat $outdir/${s}blastp.txt | grep -oP 'NODE\S*'| sort -u > ${outdir}/${s}cox1_hits.txt
        echo 'orf_ids:'
	cat $outdir/${s}cox1_hits.txt

        blastdbcmd -db ${outdir}/${s}db -dbtype prot -entry_batch $outdir/${s}cox1_hits.txt -outfmt %f -out $outdir/${s}cox1.fasta

       grep -wf $outdir/${s}cox1_hits.txt -A1 $outdir/shortest.pep.fasta > $outdir/${s}cox1_hits.fasta #| grep -v '^>' > $outdir/${s}cox1_hits.fasta


done;
