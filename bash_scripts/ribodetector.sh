#!/bin/bash
#SBATCH --partition=batch
#SBATCH --time=24:00:00
#SBATCH --job-name=ribodetector
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=100G		#memory per node
#SBATCH --threads-per-core=1
#SBATCH --output=ribodetector.out
#SBATCH --error=ribodetector.err
#SBATCH --verbose

export OMP_NUM_THREADS=20

# trimmed_names.sh is a file listing illumina sequence file names shortened 
source trimmed_names.sh

indir=/path/to/trimmed/Illumina/reads
outdir=/path/to/ribodetector/output

for s in `echo $samples`; do
	echo ${s}
	R1=`ls -1 $indir | grep -o ${s}R1.*P[^\.fastq.gz]*`
	R2=`ls -1 $indir | grep -o ${s}R2.*P[^\.fastq.gz]*`
echo $R1
echo $R2

#rrna=/work/nclab/lucy/SAB/Assembly/ribodetector/rna


ribodetector_cpu -t 20 -l 151 -i $indir/${R1}.fastq.gz $indir/${R2}.fastq.gz -e rrna --chunk_size 128 -o $outdir/${s}.R1.fq $outdir/${s}.R2.fq

done;
