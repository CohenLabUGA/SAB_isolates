#!/bin/bash
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --time=6:00:00
#SBATCH --mem=10G
#SBATCH --job-name=isip_db_search
#SBATCH --output=isip.out
#SBATCH --error=isip.err

#### activate ncbi_database env 
## samples is a space delinated list of organism names used as the folder naming scheme.
## We will loop through each organism, and within each protein to blast against 
## the organism's translated transcriptome

samples="04 08 06 13"
proteins="isip1a isip1b isip1 isip2a isip2b isip3"
proteins_dir="/path/to/isip_prot"

for s in `echo $samples`; do
    echo ${s}

### 1. name in and out dir as variables to be called throughout the script 

    indir="/path/to/${s}/transdecoder/ouput"

    outdir="/path/to/isip/output"

 	echo ${indir}
        echo ${outdir}

### 2. make a blast database from the transcriptome of each sample (s) 
## fasta headings are too long for blast, this will shorten the file header

	awk -F " " '{if ($0 ~ /^>/) {print $1} else {print $0}}' $indir/longest_orfs.pep > $indir/shorter.pep.fasta

	makeblastdb -in $indir/shorter.pep.fasta -title ${s} -out $outdir/pep_db/${s}db/${s} -dbtype prot

### 3. query the database (transcriptome) with each isip gene uploaded 
## we will use tblastp because the query is a protein and the db is a 
## translated nucleotide

	for p in `echo $proteins`; do
		echo 'blasting transcriptome for' ${p}

		blastp -query $proteins_dir/${p}.fasta -db ${outdir}/pep_db/${s}db/${s} -out ${outdir}/pep_hit/${s}${p}.txt -evalue 1e-30 


### 4. extract the transcripts in the database which the protein matched 
## to use grep to greedy match Node id's listed in the output {s}{p}.txt 
## file and append that list of matches into a txt file called hits.txt 
## we will use the hits.txt to extract the matched transcripts from the 
## database using the blastdbcmd command. if a node is matched to multiple 
## proteins in the file, remove any duplicates

			#this finds unique orfs listed in match and writes the orf id's into a file call hits.txt
		orf_id=($(cat $outdir/pep_hit/${s}${p}.txt | grep -oP 'NODE\S*'| sort -u ))
		echo 'orf_ids:' $orf_id

		echo $orf_id > $outdir/pep_hit/${s}${p}_hits.txt
	done
done;

