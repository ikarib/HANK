#!/bin/bash
#SBATCH --time=00:010:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=40G
##SBATCH --gres=gpu:pascal:1

./Main
