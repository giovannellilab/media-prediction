#!/bin/bash

#SBATCH --job-name="mmsearch"
#SBATCH --time=160:00:00
#SBATCH --cpus-per-task=80
#SBATCH --mem=200G
#SBATCH --partition=parallel

mmseqs search --threads $SLURM_CPUS_PER_TASK queryDB /path/to/UniRef90 resultDB tmp
