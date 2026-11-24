#!/bin/bash
#SBATCH --job-name=phylog_prep
#SBATCH --output=slurm-log/%x_%j.out
#SBATCH --error=slurm-log/%x_%j.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00

set -euo pipefail

PHYL_DIR="$HOME/PhyLoad"
INPUT="$PHYL_DIR/input/ref-region.tsv"
OUTPUT="$PHYL_DIR/output"
MAIL="no@thankyou.com"

mkdir -p "$OUTPUT"

# Set install source for conda
source $HOME/miniforge3/etc/profile.d/conda.sh
# module load ...
# conda init
# Activate conda environment :
conda activate colletotrichum-ncbi


echo "Job $SLURM_JOB_ID démarré sur $(hostname) à $(date)"

python3 "$PHYL_DIR/phyload.py" \
    --email "$MAIL" \
    --input "$INPUT" \
    --outdir "$OUTPUT"


echo "Finished at $(date)"