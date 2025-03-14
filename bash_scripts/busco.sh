#!/bin/bash
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=36:00:00
#SBATCH --mem=30G
#SBATCH --job-name=busco
#SBATCH --output=%x.out
#SBATCH --error=%x.err

## activate busco environment called busco

#busco Eukaryote files were downloaded:
busco_files=/path/to/.conda/envs/busco/bin/busco_downloads

samples="04 06 08 13"

for s in `echo $samples`; do
	echo ${s}
	indir=/path/to/Assembly/${s}/rnaSpades/post_ribo
	outdir=/path/to/Assembly/${s}/busco/

	busco -i $indir/transcripts.fasta -l eukaryota_odb10 -o ${s}busco_post_ribo --out_path $outdir -m transcriptome --offline --download_path $busco_files -f
done;
