#!/bin/bash
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --time=6:00:00
#SBATCH --mem=10G
#SBATCH --job-name=storage
#SBATCH --output=%x.out
#SBATCH --error=%x.err

du /work/nclab/
