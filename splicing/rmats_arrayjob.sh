#!/bin/bash
#
#SBATCH --array=3-365
#SBATCH --cpus-per-task 10
#SBATCH --mem=16GB


. "/home/tbrittoborges/miniconda3/etc/profile.d/conda.sh"
conda activate rmats

FILES=(mappings/*.bam)
NUMFILES=${#FILES[@]}
f=${FILES[$SLURM_ARRAY_TASK_ID]}
echo "$f"
out=$(basename $f .bam)
echo $f >> rmats/prep/$out.txt

python ~/rmats-turbo/rmats.py --b1 rmats/prep/$out.txt --gtf /beegfs/biodb/genomes/homo_sapiens/GRCh38_96/GRCh38.96.gtf --variable-read-length --readLength 100 --nthread 10 --novelSS --libType fr-secondstrand --allow-clipping --od rmats/ --tmp rmats/tmp/$out --task prep 

