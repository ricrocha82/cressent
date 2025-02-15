#!/bin/bash

#SBATCH --account=PAS1117
#SBATCH --job-name=align_test_job.sh
#SBATCH --time=10:00:00
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=40
#SBATCH --error=/fs/project/PAS1117/ricardo/ssDNA_tool/jobs/%x_%j.err
#SBATCH --output=/fs/project/PAS1117/ricardo/ssDNA_tool/jobs/%x_%j.out

source ~/miniconda3/bin/activate
conda activate mafft

NUM_THREADS=$SLURM_NTASKS

python /fs/project/PAS1117/ricardo/ssDNA_tool/ssDNA_annotator/modules/align.py  --threads $NUM_THREADS \
                    --input_fasta /fs/project/PAS1117/ricardo/ssDNA_tool/ssDNA_annotator/test_data/sub_reps.fa \
                    --db_family Circoviridae Adamaviridae Microviridae \
                    --db_path /fs/scratch/Sullivan_Lab/Ricardo/ssDNA_db/concat/reps \
                    -d /fs/project/PAS1117/ricardo/ssDNA_tool/test_data/output_4/