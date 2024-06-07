#!/bin/bash
#SBATCH --partition=THIN
#SBATCH --job-name=clasp
#SBATCH --nodes=1
###SBATCH --ntasks=1
#SBATCH --array=1-10
#SBARCH --cpus-per-task=4
#SBATCH --mem 50g
#SBATCH --time=10:00:00
#SBATCH --output=clasp_%a.out

s=${SLURM_ARRAY_TASK_ID}
sim=$(($s + 10))
sim_n="sim_${sim}"
echo $sim_n

export PATH="/orfeo/scratch/area/lvaleriani/myconda/bin:$PATH"
source activate claspy
cd /orfeo/LTS/LADE/LT_storage/lvaleriani/CNA/segmentation/claspy/multi_variate

for mode in sum mult max
do
	python main.py -mode $mode -sim $sim_n
done
#python read_results.py 
