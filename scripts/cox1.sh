#!/bin/bash
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --time=6:00:00
#SBATCH --mem=10G
#SBATCH --job-name=createDB
#SBATCH --output=cox1blast.out
#SBATCH --error=cox1blast.err

#make sure to activate ncbi_database env
#samples is a space delinated list of organism names used as the folder naming scheme
samples="04 08"

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

#location and name of protein file from ncbi database
cox_gene="/work/nclab/lucy/SAB/src/cylindrotheca_cox1_gene/ncbi_dataset/data/protein.faa"
#cox_gene="/work/nclab/lucy/SAB/src/coccolith_cox1_gene/${s}_cox1.fasta"

## make a blast database from the transcriptome of each sample (s)
# awk -F " " '{if ($0 ~ /^>/) {print $1} else {print $0}}' $indir/longest_orfs.pep > $indir/shorter.pep.fasta
	 sed 's/\_length.*g[[:digit:]]*//g' $indir/shorter.pep.fasta > $outdir/shortest.pep.fasta
#when using the .pep files, I have to remove -parse_seqids
head $outdir/shortest.pep.fasta

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

#        blastdbcmd -db ${outdir}/${s}db -dbtype prot -entry_batch $outdir/${s}cox1_hits.txt -outfmt %f -out $outdir/${s}cox1.fasta

       grep -wf $outdir/${s}cox1_hits.txt -A1 $outdir/shortest.pep.fasta > $outdir/${s}cox1_hits.fasta #| grep -v '^>' > $outdir/${s}cox1_hits.fasta

#	makeblastdb -in $indir/transcripts.fasta.transdecoder.pep -parse_seqids -title ${s} -out $outdir/${s}db/${s} -dbtype prot

####
## query the database (transcriptome) with the cox1 gene uploaded
##we will use tblastn because the query is a protein and the db is a nucleotide which needs to be translated
#	tblastn -query $cox_gene -db ${outdir}/${s}db/${s} -out ${outdir}/${s}cox1_tblastn.txt -evalue 1e-30 
	##using transdecorder translated transcriptome:
#	blastp -query $cox_gene -db ${outdir}/${s}db/${s} -out ${outdir}/${s}cox1_blastp.txt -evalue 1e-30 
####
## extract the transcripts in the database which the protein matched to 
## use grep to greedy match Node id's listed in the output {s}cox1_tblastn.txt file and append that list of matches into a txt file called hits.txt we will use the hits.txt to extract the matched transcripts from the database using the blastdbcmd command

#	cat $outdir/${s}cox1_blastp.txt | grep -oP 'NODE\S*' > $outdir/hits.txt
#	cat $outdir/hits.txt

#	blastdbcmd -db ${outdir}/${s}db/${s} -dbtype prot -entry_batch $outdir/hits.txt -outfmt %f -out $outdir/${s}cox1-match.fasta

# -entry_batch: file containing sequence names/id's that matched the protien query
# -outfmt: output format (here as fasta files through %f)
done;
