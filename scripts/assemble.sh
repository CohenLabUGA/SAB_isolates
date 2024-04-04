#!/bin/bash
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=36:00:00
#SBATCH --mem=100G
#SBATCH --job-name=rnaSpades
#SBATCH --output=%x.out
#SBATCH --error=%x.err

#using conda environment spades or, module load SPAdes

samples="04 06 08 13"
reps="A B C"

for s in `echo $samples`; do
	echo ${s}

	for r in `echo $reps`; do
	rnaspades.py --dataset input${r}${s}.yaml -o /work/nclab/lucy/SAB/Assembly/${s}/rnaSpades/replicates/${r}/
done;
done;

#creating a transcriptome for each sample indivdually
# for running sourmash individually on each sample


source sample_names.sh
indir="/work/nclab/lucy/SAB/Assembly/ribodetector/nonrna"

for s in `echo $samples`; do
	echo ${s}
	R1=`ls -1 $indir/* | grep -Po ${s}.*R1.*`
	R2=`ls -1 $indir/* | grep -Po ${s}.*R2.*`

	rnaspades.py --pe1-1 $indir/${R1} --pe1-2 $indir/${R2} -o /work/nclab/lucy/SAB/Assembly/ribodetector/samples_spades/${s}

done;
