#!/bin/bash

#SBATCH --job-name=HANK
#SBATCH --cpus-per-task=33
#SBATCH --mem=30G
#SBATCH --output=../mycode/TRANS/Output.txt
#SBATCH --error=../mycode/TRANS/Errors.txt
#SBATCH --mail-type=END

./Main
