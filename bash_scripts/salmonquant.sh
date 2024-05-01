#!/bin/bash
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=36:00:00
#SBATCH --mem=100G
#SBATCH --job-name=quant
#SBATCH --output=%x.out
#SBATCH --error=%x.err

transcripts=/path/to/ribodetector/output

# file listing file names to loop through
source trimmed_names.sh
isolates="04 06 08 13"
for s in `echo $samples`; do
for i in $isolates; do
        echo ${s}
        transcriptome=/path/to/${i}/rnaSpades/assemvly
        outdir=/path/to/${i}/salmon

	R1=`ls -1 $transcripts | grep -o ${s}.*R1.fq`
	R2=`ls -1 $transcripts | grep -o ${s}.*R2.fq`
	base=${s}

	echo ${R1}
	echo ${R2}

salmon quant -i $outdir/$1_index -l A -1 $transcripts/${R1} -2 $transcripts/${R2} -p 16 --validate Mappings -o $outdir/${base}_quant

done;
done;
