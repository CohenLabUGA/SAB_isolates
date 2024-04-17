#!/bin/bash
#SBATCH --partition=batch
#SBATCH --job-name=psi-04-hit-est
#SBATCH --ntasks=1
#SBATCH --time=36:00:00
#SBATCH --mem=30G
#SBATCH --output=psi-hit-04.out

module load CD-HIT
module load BLAST+

path=$"/work/nclab/data/SAB_isolates/Assembly/rnaSpades/04"

psi-cd-hit.pl -i $path/transcripts.fasta -o $path/transcripts.psi.fasta -c 95 
./psi-cd-hit.pl -i $path/transcripts95.fasta -o $path/transcripts95.psi.fasta -c 0.9 -G 1 -g 1 -prog blastn -circle 1 -exec local -para 8 -blp 4



