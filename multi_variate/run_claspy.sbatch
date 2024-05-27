#!/bin/bash
#SBATCH --partition=THIN
#SBATCH --job-name=claspy
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBARCH --cpus-per-task=4
#SBATCH --mem 100g
#SBATCH --time=10:00:00

conda activate claspy

cd /orfeo/LTS/LADE/LT_storage/lvaleriani/CNA/segmentation/claspy/multi_variate

python main.py -mode mult
python main.py -mode sum
python main.py -mode max 
