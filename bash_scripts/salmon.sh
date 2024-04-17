#!/bin/bash
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=36:00:00
#SBATCH --mem=100G
#SBATCH --job-name=Salmon.2
#SBATCH --output=%x.out
#SBATCH --error=%x.err

# first we will index the data and then we will map the transcripts to the index to quantify contigs or transcripts
# index in created from the .cds file from transdecoder. the .cds file works in nucleotide space so that we map with the transcripts (nucleotides)
## define variables ##

source trimmed_names.sh #used to call transcript files
transcripts=/work/nclab/lucy/SAB/Assembly/ribodetector/nonrna
taxa="04 06 08 13"
patern='(?:[[:alnum:]]+_)+'

for s in `echo $taxa`; do
	echo ${s}
	#locate transdecoder file '.cds' for each ${s} or taxa
	transcriptome=/work/nclab/lucy/SAB/Annotation/transdecoder/${s}/transcripts.fasta.transdecoder.cds
	outdir=/work/nclab/lucy/SAB/Assembly/${s}/salmon/run_2
	echo $transcriptome
	#index each transcriptome
	salmon index -t $transcriptome -i $outdir/${s}_index

	file_list=$(grep -Po ${s}$patern <<< $samples)
	echo "file list:  ${file_list}"

	for f in ${file_list}; do

		R1=`ls -1 $transcripts | grep -o ${f}.*R1.fq`
		R2=`ls -1 $transcripts | grep -o ${f}.*R2.fq`
		base=${f}
		echo "analysing files:
		${R1}

		${R2}
	base file name is:
		${f}"
	#quantify reads using indexed transcriptome and transcripts with rRNA removed
		salmon quant -i $outdir/${s}_index -l A -1 $transcripts/${R1} -2 $transcripts/${R2} -p 16 --validate Mappings -o $outdir/${base}_quant

	done;
done;
