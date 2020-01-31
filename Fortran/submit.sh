#!/bin/bash

#SBATCH --job-name=HANK
#SBATCH --cpus-per-task=36
#SBATCH --mem=20G
#SBATCH --mail-type=END

./Main
