#!/bin/bash

#SBATCH -p long
#SBATCH --job-name=HANK
#SBATCH --cpus-per-task=36
#SBATCH --mem=200G
#SBATCH --output=Output.txt
#SBATCH --error=Errors.txt
#SBATCH --mail-type=END

matlab -r Main
